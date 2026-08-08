#!/usr/bin/env python3
"""Qualify serial, two-rank MPI, and two-thread OpenMP zone semantics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from typing import Mapping, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "test" / "regression"))

from xnet_regression import (  # noqa: E402
    ALPHA_SPECIES,
    FinalState,
    ParsingFailure,
    parse_diagnostic,
)


EXPECTED_ZONES = tuple(range(1, 11))
FIXTURE_DIRECTORY = Path(__file__).with_name("parallel_zones")
NETWORK_INPUTS = ("sunet", "netsu", "netweak", "netwinv")
PROCESS_ARTIFACTS = (
    "xnet.command.json",
    "xnet.stdout.txt",
    "xnet.stderr.txt",
    "xnet.status.txt",
)


class QualificationFailure(RuntimeError):
    """The qualification setup, execution, inventory, or comparison failed."""


@dataclass(frozen=True)
class ProcessResult:
    command: tuple[str, ...]
    work_directory: Path
    return_code: int
    stdout: str
    stderr: str


def _resolve_executable(path: Path, label: str) -> Path:
    executable = path.expanduser().resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise QualificationFailure(f"{label} is not an executable file: {executable}")
    return executable


def _validate_fixture() -> None:
    required = [FIXTURE_DIRECTORY / "control"]
    for zone in EXPECTED_ZONES:
        required.extend(
            (
                FIXTURE_DIRECTORY / f"abundance_{zone:02d}",
                FIXTURE_DIRECTORY / f"thermo_{zone:02d}",
            )
        )
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise QualificationFailure(
            "qualification fixture is incomplete: " + ", ".join(map(str, missing))
        )

    for prefix in ("abundance", "thermo"):
        payloads = {
            (FIXTURE_DIRECTORY / f"{prefix}_{zone:02d}").read_bytes()
            for zone in EXPECTED_ZONES
        }
        if len(payloads) != len(EXPECTED_ZONES):
            raise QualificationFailure(
                f"qualification {prefix} inputs are not all distinguishable"
            )


def prepare_work_directory(work_directory: Path) -> Path:
    """Stage one isolated, writable execution directory."""

    _validate_fixture()
    work_directory = work_directory.resolve()
    if work_directory.exists():
        if not work_directory.is_dir() or any(work_directory.iterdir()):
            raise QualificationFailure(
                f"work directory must be absent or empty: {work_directory}"
            )
    else:
        work_directory.mkdir(parents=True)

    shutil.copy2(FIXTURE_DIRECTORY / "control", work_directory / "control")
    data_directory = work_directory / "Data_alpha"
    data_directory.mkdir()
    source_data = REPOSITORY_ROOT / "test" / "Data_alpha"
    for filename in NETWORK_INPUTS:
        source = source_data / filename
        if not source.is_file():
            raise QualificationFailure(f"required network input is missing: {source}")
        (data_directory / filename).symlink_to(source.resolve())

    helm_table = REPOSITORY_ROOT / "tools" / "starkiller-helmholtz" / "helm_table.dat"
    if not helm_table.is_file():
        raise QualificationFailure(f"required Helmholtz table is missing: {helm_table}")
    (work_directory / helm_table.name).symlink_to(helm_table.resolve())

    input_directory = work_directory / "inputs"
    input_directory.mkdir()
    for zone in EXPECTED_ZONES:
        for prefix in ("abundance", "thermo"):
            source = FIXTURE_DIRECTORY / f"{prefix}_{zone:02d}"
            (input_directory / source.name).symlink_to(source.resolve())
    return work_directory


def _write_process_artifacts(
    work_directory: Path,
    command: Sequence[str],
    stdout: str,
    stderr: str,
    status: str,
) -> None:
    (work_directory / PROCESS_ARTIFACTS[0]).write_text(
        json.dumps(list(command), indent=2) + "\n", encoding="utf-8"
    )
    (work_directory / PROCESS_ARTIFACTS[1]).write_text(stdout, encoding="utf-8")
    (work_directory / PROCESS_ARTIFACTS[2]).write_text(stderr, encoding="utf-8")
    (work_directory / PROCESS_ARTIFACTS[3]).write_text(status + "\n", encoding="utf-8")


def run_process(
    command: Sequence[Path | str],
    work_directory: Path,
    *,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
    expect_success: bool = True,
) -> ProcessResult:
    """Run a process group with a bounded timeout and preserve its artifacts."""

    if timeout_seconds <= 0.0 or not math.isfinite(timeout_seconds):
        raise QualificationFailure("timeout must be a positive finite value")
    normalized_command = tuple(str(item) for item in command)
    try:
        process = subprocess.Popen(
            normalized_command,
            cwd=work_directory,
            env=None if environment is None else dict(environment),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
    except OSError as error:
        _write_process_artifacts(
            work_directory, normalized_command, "", str(error), "launch-error"
        )
        raise QualificationFailure(
            f"could not launch {' '.join(normalized_command)}: {error}; "
            f"artifacts: {work_directory}"
        ) from error

    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as error:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=2.0)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
        _write_process_artifacts(
            work_directory, normalized_command, stdout, stderr, "timeout"
        )
        raise QualificationFailure(
            f"process timed out after {timeout_seconds:g} seconds; "
            f"artifacts: {work_directory}"
        ) from error

    _write_process_artifacts(
        work_directory,
        normalized_command,
        stdout,
        stderr,
        f"return_code={process.returncode}",
    )
    if expect_success and process.returncode != 0:
        raise QualificationFailure(
            f"process returned {process.returncode}; artifacts: {work_directory}"
        )
    if not expect_success and process.returncode == 0:
        raise QualificationFailure(
            f"failure probe returned zero; artifacts: {work_directory}"
        )
    return ProcessResult(
        normalized_command, work_directory, process.returncode, stdout, stderr
    )


def validate_output_inventory(work_directory: Path, worker_count: int) -> tuple[Path, ...]:
    """Require exactly ten ASCII/binary zone outputs and one diagnostic per worker."""

    expected_ev = {f"ev_parallel_zones_{zone:02d}" for zone in EXPECTED_ZONES}
    expected_ts = {f"ts_parallel_zones_{zone:02d}" for zone in EXPECTED_ZONES}
    actual_ev = {path.name for path in work_directory.glob("ev_parallel_zones_*")}
    actual_ts = {path.name for path in work_directory.glob("ts_parallel_zones_*")}
    if actual_ev != expected_ev:
        raise QualificationFailure(
            f"ASCII output inventory mismatch: {sorted(actual_ev)} != {sorted(expected_ev)}"
        )
    if actual_ts != expected_ts:
        raise QualificationFailure(
            f"binary output inventory mismatch: {sorted(actual_ts)} != {sorted(expected_ts)}"
        )
    for filename in (*actual_ev, *actual_ts):
        if (work_directory / filename).stat().st_size == 0:
            raise QualificationFailure(f"output is empty: {filename}")

    diagnostics = tuple(sorted(work_directory.glob("net_diag*")))
    if len(diagnostics) != worker_count or any(path.stat().st_size == 0 for path in diagnostics):
        raise QualificationFailure(
            f"expected {worker_count} nonempty worker diagnostics, found "
            f"{[path.name for path in diagnostics]}"
        )
    return diagnostics


def _discover_diagnostic_groups(text: str, path: Path) -> tuple[tuple[int, ...], ...]:
    groups: list[tuple[int, ...]] = []
    current: list[int] = []
    for line in text.splitlines():
        if line.startswith("End"):
            fields = line.split()
            if len(fields) < 2 or not fields[1].isdigit():
                raise QualificationFailure(f"malformed End record in {path}: {line}")
            current.append(int(fields[1]))
        elif line.startswith("Counters:"):
            if not current:
                raise QualificationFailure(f"Counters record without final states in {path}")
            groups.append(tuple(current))
            current = []
    if current or not groups:
        raise QualificationFailure(f"incomplete or empty diagnostic grouping in {path}")
    return tuple(groups)


def normalize_worker_states(
    worker_states: Sequence[Sequence[FinalState]],
    expected_zones: Sequence[int] = EXPECTED_ZONES,
) -> tuple[FinalState, ...]:
    """Merge rank/thread results by global zone without assuming worker ownership."""

    expected_zones = tuple(expected_zones)
    by_zone: dict[int, FinalState] = {}
    duplicates: set[int] = set()
    for states in worker_states:
        for state in states:
            if state.zone in by_zone:
                duplicates.add(state.zone)
            else:
                by_zone[state.zone] = state
    if duplicates:
        raise QualificationFailure(f"duplicate final zones: {sorted(duplicates)}")
    actual = set(by_zone)
    expected = set(expected_zones)
    if actual != expected:
        raise QualificationFailure(
            f"final-zone inventory mismatch; missing={sorted(expected - actual)}, "
            f"unexpected={sorted(actual - expected)}"
        )
    return tuple(by_zone[zone] for zone in expected_zones)


def parse_worker_diagnostics(paths: Sequence[Path]) -> tuple[FinalState, ...]:
    worker_states: list[tuple[FinalState, ...]] = []
    for path in paths:
        try:
            text = path.read_text(encoding="utf-8")
            groups = _discover_diagnostic_groups(text, path)
            zones = tuple(zone for group in groups for zone in group)
            worker_states.append(parse_diagnostic(text, zones, ALPHA_SPECIES, groups))
        except (OSError, UnicodeError, ParsingFailure) as error:
            raise QualificationFailure(f"could not parse worker diagnostic {path}: {error}") from error
    return normalize_worker_states(worker_states)


def _endpoint_signature(state: FinalState) -> tuple[object, ...]:
    return (
        state.target_time,
        state.time,
        state.temperature_gk,
        state.density,
        state.electron_fraction,
        tuple(state.mass_fractions.items()),
    )


def compare_endpoint_states(
    actual: Sequence[FinalState], reference: Sequence[FinalState], label: str
) -> None:
    """Apply the batch regression's exact printed-endpoint equivalence policy."""

    actual_by_zone = {state.zone: state for state in actual}
    reference_by_zone = {state.zone: state for state in reference}
    if set(actual_by_zone) != set(reference_by_zone):
        raise QualificationFailure(f"{label} endpoint zone inventory differs from serial")
    differences = [
        zone
        for zone in reference_by_zone
        if _endpoint_signature(actual_by_zone[zone])
        != _endpoint_signature(reference_by_zone[zone])
    ]
    if differences:
        raise QualificationFailure(
            f"{label} endpoint values differ from serial in zones {differences}"
        )


def _close(actual: float, expected: float, relative: float) -> bool:
    return math.isclose(actual, expected, rel_tol=relative, abs_tol=5.0e-99)


def validate_ascii_association(
    work_directory: Path, states: Sequence[FinalState]
) -> None:
    """Match each filename's final ASCII row to that global zone's diagnostic state."""

    for state in states:
        path = work_directory / f"ev_parallel_zones_{state.zone:02d}"
        lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(lines) < 2:
            raise QualificationFailure(f"ASCII history has no final row: {path}")
        fields = lines[-1].split()
        expected_fields = 7 + len(ALPHA_SPECIES) + 2
        if len(fields) != expected_fields:
            raise QualificationFailure(f"malformed final ASCII row in {path}: {lines[-1]}")
        try:
            step = int(fields[0])
            time = float(fields[1].replace("D", "E"))
            temperature = float(fields[2].replace("D", "E"))
            density = float(fields[3].replace("D", "E"))
            mass_fractions = tuple(
                float(token.replace("D", "E"))
                for token in fields[7 : 7 + len(ALPHA_SPECIES)]
            )
        except ValueError as error:
            raise QualificationFailure(f"non-numeric final ASCII row in {path}") from error
        # The diagnostic End step is the batch-level maximum, while the ASCII
        # history step is the per-zone time-step counter.
        if step != state.counters.ts or not _close(time, state.time, 5.1e-9):
            raise QualificationFailure(f"time/step association mismatch for zone {state.zone}")
        if not _close(temperature, state.temperature_gk, 5.1e-4) or not _close(
            density, state.density, 5.1e-4
        ):
            raise QualificationFailure(f"thermodynamic association mismatch for zone {state.zone}")
        for species, actual, expected in zip(
            ALPHA_SPECIES, mass_fractions, state.mass_fractions.values(), strict=True
        ):
            if not _close(actual, expected, 5.1e-3):
                raise QualificationFailure(
                    f"composition association mismatch for zone {state.zone} species {species}"
                )


def _run_success(
    label: str,
    command: Sequence[Path | str],
    work_directory: Path,
    worker_count: int,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
) -> tuple[FinalState, ...]:
    prepare_work_directory(work_directory)
    run_process(
        command,
        work_directory,
        timeout_seconds=timeout_seconds,
        environment=environment,
    )
    diagnostics = validate_output_inventory(work_directory, worker_count)
    states = parse_worker_diagnostics(diagnostics)
    validate_ascii_association(work_directory, states)
    if len({_endpoint_signature(state) for state in states}) != len(EXPECTED_ZONES):
        raise QualificationFailure(f"{label} did not retain distinguishable per-zone results")
    return states


def run_qualification(arguments: argparse.Namespace) -> Path:
    serial = _resolve_executable(arguments.serial_executable, "serial executable")
    mpi = _resolve_executable(arguments.mpi_executable, "MPI executable")
    openmp = _resolve_executable(arguments.openmp_executable, "OpenMP executable")
    mpi_launcher = _resolve_executable(arguments.mpi_launcher, "MPI launcher")

    work_root = arguments.work_root.expanduser().resolve()
    if work_root.exists():
        if not work_root.is_dir() or any(work_root.iterdir()):
            raise QualificationFailure(f"work root must be absent or empty: {work_root}")
    else:
        work_root.mkdir(parents=True)

    serial_states = _run_success(
        "serial", (serial,), work_root / "serial", 1, arguments.timeout
    )
    mpi_command = (
        mpi_launcher,
        *arguments.mpi_launcher_argument,
        "-n",
        "2",
        mpi,
    )
    mpi_states = _run_success(
        "MPI", mpi_command, work_root / "mpi", 2, arguments.timeout
    )
    compare_endpoint_states(mpi_states, serial_states, "MPI")

    failure_directory = prepare_work_directory(work_root / "mpi-nonroot-failure")
    # With three rank-strided batches, two ranks assign zones 5-8 to rank 1.
    (failure_directory / "inputs" / "thermo_05").unlink()
    run_process(
        mpi_command,
        failure_directory,
        timeout_seconds=arguments.timeout,
        expect_success=False,
    )

    openmp_environment = os.environ.copy()
    openmp_environment.update({"OMP_NUM_THREADS": "2", "OMP_DYNAMIC": "FALSE"})
    openmp_states: list[tuple[FinalState, ...]] = []
    for repetition in range(1, 4):
        states = _run_success(
            f"OpenMP repetition {repetition}",
            (openmp,),
            work_root / f"openmp-{repetition}",
            2,
            arguments.timeout,
            openmp_environment,
        )
        compare_endpoint_states(states, serial_states, f"OpenMP repetition {repetition}")
        openmp_states.append(states)
    for repetition, states in enumerate(openmp_states[1:], start=2):
        compare_endpoint_states(
            states, openmp_states[0], f"OpenMP repetition {repetition} repeatability"
        )

    summary = work_root / "qualification_summary.json"
    summary.write_text(
        json.dumps(
            {
                "fixture": "ten distinguishable zones, nzbatchmx=4",
                "comparison": "exact normalized printed endpoints by global zone",
                "serial": {"runs": 1, "workers": 1},
                "mpi": {"runs": 1, "ranks": 2, "nonroot_failure_probe": "nonzero"},
                "openmp": {"runs": 3, "threads": 2, "OMP_DYNAMIC": "FALSE"},
                "inactive_final_batch_lanes": 2,
                "status": "passed",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return summary


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial-executable", type=Path, required=True)
    parser.add_argument("--mpi-executable", type=Path, required=True)
    parser.add_argument("--openmp-executable", type=Path, required=True)
    parser.add_argument(
        "--mpi-launcher", type=Path, default=Path(shutil.which("mpiexec") or "mpiexec")
    )
    parser.add_argument(
        "--mpi-launcher-argument",
        action="append",
        default=[],
        help="repeatable launcher option placed before -n (for example --oversubscribe)",
    )
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=60.0)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    try:
        summary = run_qualification(parse_arguments(argv))
    except (QualificationFailure, OSError) as error:
        print(f"parallel-zone qualification failed: {error}", file=sys.stderr)
        return 1
    print(f"parallel-zone qualification passed; summary: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
