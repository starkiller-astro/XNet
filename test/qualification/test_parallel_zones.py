"""Focused tests for the parallel-zone qualification runner."""

from dataclasses import replace
import os
from pathlib import Path
import signal
import sys
import time

import pytest

from parallel_zones import (
    ALPHA_SPECIES,
    AsciiEndpoint,
    EXPECTED_ZONES,
    QualificationFailure,
    compare_ascii_endpoints,
    compare_endpoint_states,
    normalize_worker_states,
    run_process,
    validate_ascii_association,
    validate_mpi_topology,
    validate_openmp_topology,
    validate_output_inventory,
)
from xnet_regression import FinalState, SolverCounters


def _state(zone: int, value: float | None = None) -> FinalState:
    value = float(zone) if value is None else value
    return FinalState(
        zone=zone,
        step=zone + 10,
        target_time=value * 1.0e-3,
        time=value * 1.0e-3,
        temperature_gk=1.0 + value,
        density=1.0e7 * value,
        electron_fraction=0.5,
        mass_fractions={species: value * 1.0e-3 for species in ALPHA_SPECIES},
        counters=SolverCounters(1, 2, 3, 4, 5),
    )


def _complete_inventory(directory: Path, worker_count: int = 1) -> None:
    for zone in EXPECTED_ZONES:
        (directory / f"ev_parallel_zones_{zone:02d}").write_text("ev\n", encoding="utf-8")
        (directory / f"ts_parallel_zones_{zone:02d}").write_bytes(b"ts")
    for worker in range(worker_count):
        (directory / f"net_diag{worker:02d}").write_text("diag\n", encoding="utf-8")


def _write_ascii_state(
    directory: Path,
    state: FinalState,
    filename_zone: int,
    *,
    neutrino_loss_rate: str = "0.0",
    time_value: str | None = None,
) -> None:
    values = " ".join(f"{value:.8E}" for value in state.mass_fractions.values())
    row = (
        f"{state.counters.ts} {time_value or f'{state.time:.8E}'} "
        f"{state.temperature_gk:.3E} "
        f"{state.density:.3E} 0.0 {neutrino_loss_rate} 1.0E-6 {values} 1 1\n"
    )
    (directory / f"ev_parallel_zones_{filename_zone:02d}").write_text(
        "header\n" + row, encoding="utf-8"
    )


def test_normalization_rejects_duplicate_and_missing_zones() -> None:
    with pytest.raises(QualificationFailure, match="duplicate final zones"):
        normalize_worker_states(((_state(1),), (_state(1),)), expected_zones=(1, 2))
    with pytest.raises(QualificationFailure, match=r"missing=\[2\]"):
        normalize_worker_states(((_state(1),),), expected_zones=(1, 2))


def test_exact_endpoint_comparison_detects_state_leakage() -> None:
    reference = (_state(1), _state(2))
    leaked = (reference[0], replace(reference[1], mass_fractions=reference[0].mass_fractions))
    with pytest.raises(QualificationFailure, match=r"zones \[2\]"):
        compare_endpoint_states(leaked, reference, "mutated")


@pytest.mark.parametrize("kind", ("missing", "off-by-one"))
def test_inventory_rejects_missing_or_off_by_one_zone(tmp_path: Path, kind: str) -> None:
    _complete_inventory(tmp_path)
    if kind == "missing":
        (tmp_path / "ts_parallel_zones_10").unlink()
    else:
        (tmp_path / "ev_parallel_zones_11").write_text("ev\n", encoding="utf-8")
    with pytest.raises(QualificationFailure, match="output inventory mismatch"):
        validate_output_inventory(tmp_path)


def test_ascii_filename_association_rejects_swapped_zone_content(tmp_path: Path) -> None:
    states = (_state(1), _state(2))
    _write_ascii_state(tmp_path, states[1], filename_zone=1)
    _write_ascii_state(tmp_path, states[0], filename_zone=2)
    with pytest.raises(QualificationFailure, match="association mismatch for zone 1"):
        validate_ascii_association(tmp_path, states)


def test_ascii_parser_accepts_fortran_omitted_exponent_letter(tmp_path: Path) -> None:
    state = _state(1)
    _write_ascii_state(
        tmp_path,
        state,
        filename_zone=1,
        neutrino_loss_rate="6.95-310",
    )
    endpoint = validate_ascii_association(tmp_path, (state,))[0]
    assert endpoint.neutrino_loss_rate == float("6.95e-310")


def test_ascii_association_accepts_difference_from_documented_time_formats(
    tmp_path: Path,
) -> None:
    state = replace(_state(1), time=1.2345678)
    _write_ascii_state(
        tmp_path,
        state,
        filename_zone=1,
        time_value="1.23456784E+00",
    )
    endpoint = validate_ascii_association(tmp_path, (state,))[0]
    assert endpoint.zone == 1


def test_expected_failure_requires_nonzero_status(tmp_path: Path) -> None:
    result = run_process(
        (sys.executable, "-c", "raise SystemExit(3)"),
        tmp_path,
        timeout_seconds=2.0,
        expect_success=False,
    )
    assert result.return_code == 3
    with pytest.raises(QualificationFailure, match="failure probe returned zero"):
        run_process(
            (sys.executable, "-c", "raise SystemExit(0)"),
            tmp_path,
            timeout_seconds=2.0,
            expect_success=False,
        )


def test_ascii_endpoint_comparison_detects_energy_leakage() -> None:
    reference = (
        AsciiEndpoint(1, 1.0, 2.0, 3.0),
        AsciiEndpoint(2, 4.0, 5.0, 6.0),
    )
    leaked = (reference[0], replace(reference[1], energy_generation_rate=1.0))
    with pytest.raises(QualificationFailure, match=r"zones \[2\]"):
        compare_ascii_endpoints(leaked, reference, "mutated")


def test_parallel_topology_rejects_serial_worker_headers(tmp_path: Path) -> None:
    serial_diagnostic = tmp_path / "net_diag"
    serial_diagnostic.write_text(" MyId    0    1\n", encoding="utf-8")
    with pytest.raises(QualificationFailure, match="MPI topology mismatch"):
        validate_mpi_topology((serial_diagnostic,), 2)
    with pytest.raises(QualificationFailure, match="OpenMP topology mismatch"):
        validate_openmp_topology((serial_diagnostic,), 2)


def test_parallel_topology_accepts_two_reported_workers(tmp_path: Path) -> None:
    mpi_diagnostics = tuple(tmp_path / f"mpi-{rank}" for rank in range(2))
    for rank, path in enumerate(mpi_diagnostics):
        path.write_text(f" MyId{rank:5d}{2:5d}\n", encoding="utf-8")
    validate_mpi_topology(mpi_diagnostics, 2)

    openmp_diagnostics = tuple(tmp_path / f"openmp-{thread}" for thread in range(1, 3))
    for thread, path in enumerate(openmp_diagnostics, start=1):
        path.write_text(
            f" MyId{0:5d}{1:5d}\nThread {thread:4d} of {2:4d}\n",
            encoding="utf-8",
        )
    validate_openmp_topology(openmp_diagnostics, 2)


def test_timeout_terminates_the_process_group(tmp_path: Path) -> None:
    child_script = tmp_path / "child.py"
    child_script.write_text(
        """\
import os
from pathlib import Path
import signal
import sys
import time

marker = Path(sys.argv[1])
pid_file = Path(sys.argv[2])

def terminate(signum, frame):
    marker.write_text("terminated\\n", encoding="utf-8")
    raise SystemExit(0)

signal.signal(signal.SIGTERM, terminate)
pid_file.write_text(str(os.getpid()), encoding="utf-8")
while True:
    time.sleep(1)
""",
        encoding="utf-8",
    )
    marker = tmp_path / "child-terminated.txt"
    pid_file = tmp_path / "child.pid"
    leader = (
        "import subprocess, sys, time; "
        f"subprocess.Popen([sys.executable, {str(child_script)!r}, "
        f"{str(marker)!r}, {str(pid_file)!r}], "
        "stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); "
        "time.sleep(10)"
    )

    child_pid: int | None = None
    try:
        with pytest.raises(QualificationFailure, match="timed out"):
            run_process(
                (sys.executable, "-c", leader),
                tmp_path,
                timeout_seconds=0.5,
            )
        deadline = time.monotonic() + 2.0
        while not marker.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        assert marker.read_text(encoding="utf-8") == "terminated\n"
    finally:
        if pid_file.exists():
            child_pid = int(pid_file.read_text(encoding="utf-8"))
        if child_pid is not None:
            try:
                os.kill(child_pid, 0)
            except ProcessLookupError:
                pass
            else:
                os.kill(child_pid, signal.SIGKILL)
    assert (tmp_path / "xnet.status.txt").read_text(encoding="utf-8") == "timeout\n"
