#!/usr/bin/env python3
"""Run and validate the manual Frontier GPU correctness qualification."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import PurePosixPath
from typing import Mapping, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT / "test" / "qualification"))
sys.path.insert(0, str(REPOSITORY_ROOT / "test" / "regression"))

from parallel_zones import (  # noqa: E402
    AsciiEndpoint,
    EXPECTED_ZONES,
    QualificationFailure as ParallelQualificationFailure,
    run_configuration,
    validate_ascii_association,
)
from xnet_regression import (  # noqa: E402
    FinalState,
    RegressionFailure,
    heat_sn160_case,
    parse_diagnostic,
    prepare_work_directory as prepare_regression_work_directory,
    run_xnet,
)


MANIFEST_SCHEMA = "xnet-frontier-qualification-v2"
POLICY_SCHEMA = "xnet-frontier-comparison-v1"
GPU_LINALG_RESIDUAL_LIMIT = 1.0e-12
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
SOURCE_SHA_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SLURM_TIME_PATTERN = re.compile(r"^(?:\d+-)?\d{1,2}:\d{2}:\d{2}$")
LINALG_BATCH_PATTERN = re.compile(
    r"^XNET_GPU_LINALG batch\s+(\d+)\s+info\s+(-?\d+)\s+"
    r"relative_residual\s+([^\s]+)$"
)
REQUIRED_MODULE_MARKERS = (
    "PrgEnv-cray",
    "rocm",
    "craype-accel-amd-gfx90a",
    "hipfort",
    "cray-python",
)
FAILURE_CATEGORIES = (
    "source",
    "environment",
    "submission",
    "queue",
    "allocation",
    "facility",
    "build",
    "test",
    "comparison",
)
CPU_BUILD_VARIABLES = {
    "BUILD_NAME": "frontier-cpu",
    "CMODE": "OPT",
    "PE_ENV": "CRAY",
    "MPI_MODE": "OFF",
    "OPENMP_MODE": "OFF",
    "GPU_MODE": "OFF",
    "OPENACC_MODE": "OFF",
    "OPENMP_OL_MODE": "OFF",
    "EOS": "STARKILLER",
    "MATRIX_SOLVER": "dense",
    "LAPACK_VER": "LIBSCI",
}
GPU_BUILD_VARIABLES = {
    "BUILD_NAME": "frontier-gpu",
    "CMODE": "OPT",
    "PE_ENV": "CRAY",
    "MPI_MODE": "OFF",
    "OPENMP_MODE": "OFF",
    "GPU_MODE": "ON",
    "GPU_BACKEND": "HIP",
    "OPENACC_MODE": "OFF",
    "OPENMP_OL_MODE": "ON",
    "GPU_LAPACK_VER": "ROCM",
    "EOS": "STARKILLER",
    "MATRIX_SOLVER": "dense",
    "LAPACK_VER": "LIBSCI",
}
RESOLVED_BUILD_VARIABLES = (
    "FC",
    "LDR",
    "FFLAGS",
    "LDFLAGS",
    "CRAY_OMP_PREPROCESS",
    "XNET_CPP",
)
EXPECTED_EXECUTABLES = {
    "cpu": ("xnet",),
    "gpu": ("xnet", "frontier_gpu_linalg_probe"),
}
HEAT_SN160_ASCII_SPECIES = (
    "p",
    "he4",
    "c12",
    "o16",
    "ne20",
    "mg24",
    "si28",
    "ca40",
    "ti44",
    "cr48",
    "fe52",
    "fe54",
    "ni56",
    "zn60",
)


class FrontierFailure(RuntimeError):
    """A classified source, environment, build, execution, or comparison failure."""

    def __init__(self, category: str, phase: str, message: str):
        super().__init__(message)
        self.category = category
        self.phase = phase


@dataclass(frozen=True)
class Bounds:
    atol: float
    rtol: float


@dataclass(frozen=True)
class NumericalPolicy:
    status: str
    target_time_exact: bool
    scalar_fields: Mapping[str, Bounds]
    ascii_fields: Mapping[str, Bounds]
    reported_ascii_fields: tuple[str, ...]
    material_threshold: float
    anchors: tuple[str, ...]
    selected: Bounds
    complete_vector_l1_limit: float
    complete_vector_linf_limit: float
    normalization_atol: float


def _finite_nonnegative(value: object, context: str) -> float:
    if not isinstance(value, (int, float)):
        raise FrontierFailure("source", "comparison-policy", f"{context} is not numeric")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise FrontierFailure(
            "source", "comparison-policy", f"{context} must be finite and nonnegative"
        )
    return normalized


def _load_bounds(document: object, context: str) -> Bounds:
    if not isinstance(document, dict) or set(document) != {"atol", "rtol"}:
        raise FrontierFailure(
            "source", "comparison-policy", f"{context} must contain only atol and rtol"
        )
    return Bounds(
        _finite_nonnegative(document["atol"], f"{context}.atol"),
        _finite_nonnegative(document["rtol"], f"{context}.rtol"),
    )


def load_policy(path: Path) -> NumericalPolicy:
    """Load the reviewed CPU/GPU numerical limits."""

    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise FrontierFailure(
            "source", "comparison-policy", f"could not read policy {path}: {error}"
        ) from error
    required = {
        "schema",
        "status",
        "target_time_exact",
        "scalar_fields",
        "ascii_fields",
        "reported_ascii_fields",
        "mass_fractions",
    }
    if not isinstance(document, dict) or set(document) != required:
        raise FrontierFailure(
            "source", "comparison-policy", "comparison policy has unexpected fields"
        )
    if document["schema"] != POLICY_SCHEMA:
        raise FrontierFailure("source", "comparison-policy", "comparison policy schema differs")
    scalar_names = {
        "achieved_time",
        "temperature_gk",
        "density",
        "electron_fraction",
    }
    ascii_names = {"neutrino_loss_rate", "timestep"}
    reported_ascii_names = ["energy_generation_rate"]
    scalar_document = document["scalar_fields"]
    ascii_document = document["ascii_fields"]
    if not isinstance(scalar_document, dict) or set(scalar_document) != scalar_names:
        raise FrontierFailure("source", "comparison-policy", "scalar field policy is incomplete")
    if not isinstance(ascii_document, dict) or set(ascii_document) != ascii_names:
        raise FrontierFailure("source", "comparison-policy", "ASCII field policy is incomplete")
    reported_ascii_fields = document["reported_ascii_fields"]
    if (
        not isinstance(reported_ascii_fields, list)
        or reported_ascii_fields != reported_ascii_names
    ):
        raise FrontierFailure(
            "source", "comparison-policy", "reported ASCII field policy is incomplete"
        )
    mass = document["mass_fractions"]
    mass_names = {
        "material_threshold",
        "anchors",
        "selected",
        "complete_vector_l1_limit",
        "complete_vector_linf_limit",
        "normalization_atol",
    }
    if not isinstance(mass, dict) or set(mass) != mass_names:
        raise FrontierFailure("source", "comparison-policy", "mass-fraction policy is incomplete")
    anchors = mass["anchors"]
    if (
        not isinstance(anchors, list)
        or not anchors
        or any(not isinstance(item, str) or not item for item in anchors)
        or len(set(anchors)) != len(anchors)
    ):
        raise FrontierFailure("source", "comparison-policy", "anchors are invalid")
    if not isinstance(document["status"], str) or not document["status"].strip():
        raise FrontierFailure("source", "comparison-policy", "policy status is empty")
    if document["target_time_exact"] is not True:
        raise FrontierFailure("source", "comparison-policy", "target time must remain exact")
    return NumericalPolicy(
        status=document["status"],
        target_time_exact=True,
        scalar_fields={
            name: _load_bounds(scalar_document[name], f"scalar_fields.{name}")
            for name in sorted(scalar_names)
        },
        ascii_fields={
            name: _load_bounds(ascii_document[name], f"ascii_fields.{name}")
            for name in sorted(ascii_names)
        },
        reported_ascii_fields=tuple(reported_ascii_fields),
        material_threshold=_finite_nonnegative(
            mass["material_threshold"], "mass_fractions.material_threshold"
        ),
        anchors=tuple(anchors),
        selected=_load_bounds(mass["selected"], "mass_fractions.selected"),
        complete_vector_l1_limit=_finite_nonnegative(
            mass["complete_vector_l1_limit"],
            "mass_fractions.complete_vector_l1_limit",
        ),
        complete_vector_linf_limit=_finite_nonnegative(
            mass["complete_vector_linf_limit"],
            "mass_fractions.complete_vector_linf_limit",
        ),
        normalization_atol=_finite_nonnegative(
            mass["normalization_atol"], "mass_fractions.normalization_atol"
        ),
    )


def _difference(actual: float, reference: float, bounds: Bounds) -> tuple[bool, float, float]:
    difference = abs(actual - reference)
    allowed = bounds.atol + bounds.rtol * abs(reference)
    return difference <= allowed, difference, allowed


def compare_endpoint_states(
    actual: Sequence[FinalState],
    reference: Sequence[FinalState],
    policy: NumericalPolicy,
    label: str,
) -> dict[str, object]:
    """Compare normalized per-zone endpoints and return observed-difference evidence."""

    actual_by_zone = {state.zone: state for state in actual}
    reference_by_zone = {state.zone: state for state in reference}
    if len(actual_by_zone) != len(actual) or len(reference_by_zone) != len(reference):
        raise FrontierFailure("test", label, "duplicate endpoint zone")
    if set(actual_by_zone) != set(reference_by_zone):
        raise FrontierFailure(
            "test", label, "CPU/GPU endpoint zone inventory differs"
        )

    failures: list[str] = []
    observations: list[dict[str, object]] = []
    maximum_selected_ratio = 0.0
    for zone in sorted(reference_by_zone):
        candidate = actual_by_zone[zone]
        baseline = reference_by_zone[zone]
        if tuple(candidate.mass_fractions) != tuple(baseline.mass_fractions):
            failures.append(f"zone {zone} species identity/order differs")
            continue
        if candidate.target_time != baseline.target_time:
            failures.append(f"zone {zone} target_time is not exact")

        scalar_differences: dict[str, dict[str, float]] = {}
        for policy_name, attribute in (
            ("achieved_time", "time"),
            ("temperature_gk", "temperature_gk"),
            ("density", "density"),
            ("electron_fraction", "electron_fraction"),
        ):
            observed = getattr(candidate, attribute)
            expected = getattr(baseline, attribute)
            passed, difference, allowed = _difference(
                observed, expected, policy.scalar_fields[policy_name]
            )
            scalar_differences[policy_name] = {
                "observed": observed,
                "expected": expected,
                "difference": difference,
                "allowed": allowed,
            }
            if not passed:
                failures.append(
                    f"zone {zone} {policy_name} difference {difference:.3e} "
                    f"exceeds {allowed:.3e}"
                )

        selected = tuple(
            species
            for species, value in baseline.mass_fractions.items()
            if species in policy.anchors or value >= policy.material_threshold
        )
        if not selected:
            failures.append(f"zone {zone} selects no material species")
        vector_differences = {
            species: abs(candidate.mass_fractions[species] - expected)
            for species, expected in baseline.mass_fractions.items()
        }
        linf_species = max(vector_differences, key=vector_differences.__getitem__)
        linf = vector_differences[linf_species]
        l1 = math.fsum(vector_differences.values())
        if l1 > policy.complete_vector_l1_limit:
            failures.append(
                f"zone {zone} composition L1 {l1:.3e} exceeds "
                f"{policy.complete_vector_l1_limit:.3e}"
            )
        if linf > policy.complete_vector_linf_limit:
            failures.append(
                f"zone {zone} composition Linf {linf:.3e} at {linf_species} exceeds "
                f"{policy.complete_vector_linf_limit:.3e}"
            )
        selected_differences: dict[str, dict[str, float]] = {}
        for species in selected:
            observed = candidate.mass_fractions[species]
            expected = baseline.mass_fractions[species]
            passed, difference, allowed = _difference(
                observed,
                expected,
                policy.selected,
            )
            selected_differences[species] = {
                "observed": observed,
                "expected": expected,
                "difference": difference,
                "allowed": allowed,
            }
            ratio = difference / allowed if allowed > 0.0 else math.inf
            maximum_selected_ratio = max(maximum_selected_ratio, ratio)
            if not passed:
                failures.append(
                    f"zone {zone} {species} difference {difference:.3e} exceeds "
                    f"{allowed:.3e}"
                )
        for state_name, state in (("CPU", baseline), ("GPU", candidate)):
            values = tuple(state.mass_fractions.values())
            if any(not math.isfinite(value) or value < 0.0 for value in values):
                failures.append(f"zone {zone} {state_name} composition is non-finite or negative")
            normalization_error = abs(math.fsum(values) - 1.0)
            if normalization_error > policy.normalization_atol:
                failures.append(
                    f"zone {zone} {state_name} normalization error "
                    f"{normalization_error:.3e} exceeds {policy.normalization_atol:.3e}"
                )
        observations.append(
            {
                "zone": zone,
                "scalar_differences": scalar_differences,
                "selected_species": list(selected),
                "selected_differences": selected_differences,
                "observed_mass_fractions": dict(candidate.mass_fractions),
                "reference_mass_fractions": dict(baseline.mass_fractions),
                "composition_l1": l1,
                "composition_linf": linf,
                "composition_linf_species": linf_species,
            }
        )
    if failures:
        raise FrontierFailure(
            "comparison", label, "CPU/GPU endpoint comparison failed:\n  " + "\n  ".join(failures)
        )
    return {
        "status": "passed",
        "policy_status": policy.status,
        "maximum_selected_fraction_of_allowed": maximum_selected_ratio,
        "zones": observations,
    }


def compare_ascii_endpoints(
    actual: Sequence[AsciiEndpoint],
    reference: Sequence[AsciiEndpoint],
    policy: NumericalPolicy,
    label: str,
) -> dict[str, object]:
    actual_by_zone = {item.zone: item for item in actual}
    reference_by_zone = {item.zone: item for item in reference}
    if set(actual_by_zone) != set(reference_by_zone):
        raise FrontierFailure("test", label, "ASCII endpoint zone inventory differs")
    failures: list[str] = []
    maximum_ratio = 0.0
    observations: list[dict[str, object]] = []
    for zone in sorted(reference_by_zone):
        field_differences: dict[str, dict[str, object]] = {}
        for field in policy.ascii_fields:
            observed = getattr(actual_by_zone[zone], field)
            expected = getattr(reference_by_zone[zone], field)
            passed, difference, allowed = _difference(
                observed, expected, policy.ascii_fields[field]
            )
            field_differences[field] = {
                "comparison": "bounded",
                "observed": observed,
                "expected": expected,
                "difference": difference,
                "allowed": allowed,
            }
            maximum_ratio = max(
                maximum_ratio, difference / allowed if allowed > 0.0 else math.inf
            )
            if not passed:
                failures.append(
                    f"zone {zone} {field}: observed {observed:.8e}, "
                    f"expected {expected:.8e}, difference {difference:.3e} "
                    f"exceeds {allowed:.3e}"
                )
        for field in policy.reported_ascii_fields:
            observed = getattr(actual_by_zone[zone], field)
            expected = getattr(reference_by_zone[zone], field)
            field_differences[field] = {
                "comparison": "reported_only",
                "observed": observed,
                "expected": expected,
                "difference": abs(observed - expected),
            }
        observations.append({"zone": zone, "field_differences": field_differences})
    if failures:
        raise FrontierFailure(
            "comparison", label, "CPU/GPU ASCII comparison failed:\n  " + "\n  ".join(failures)
        )
    return {
        "status": "passed",
        "maximum_fraction_of_allowed": maximum_ratio,
        "reported_only_fields": list(policy.reported_ascii_fields),
        "zones": observations,
    }


def parse_linalg_probe(
    text: str, residual_limit: float = GPU_LINALG_RESIDUAL_LIMIT
) -> dict[str, object]:
    """Require real device execution, two successful factors, and small residuals."""

    device_count: int | None = None
    offloaded: bool | None = None
    data_present: bool | None = None
    status: str | None = None
    batches: list[dict[str, object]] = []
    for line in text.splitlines():
        fields = line.split()
        if line.startswith("XNET_GPU_LINALG device_count ") and len(fields) == 3:
            device_count = int(fields[2])
        elif line.startswith("XNET_GPU_LINALG offloaded ") and len(fields) == 3:
            offloaded = fields[2].upper() == "T"
        elif line.startswith("XNET_GPU_LINALG data_present ") and len(fields) == 3:
            data_present = fields[2].upper() == "T"
        elif line.startswith("XNET_GPU_LINALG status ") and len(fields) >= 3:
            status = fields[2]
        else:
            match = LINALG_BATCH_PATTERN.match(line)
            if match:
                residual = float(match.group(3).replace("D", "E"))
                batches.append(
                    {
                        "batch": int(match.group(1)),
                        "info": int(match.group(2)),
                        "relative_residual": residual,
                    }
                )
    if device_count is None or device_count < 1:
        raise FrontierFailure("allocation", "gpu-linalg", "probe did not find a GPU")
    if offloaded is not True:
        raise FrontierFailure("allocation", "gpu-linalg", "OpenMP target stayed on the host")
    if data_present is not True:
        raise FrontierFailure("test", "gpu-linalg", "mapped solve data is not present")
    if status != "passed" or [item["batch"] for item in batches] != [1, 2]:
        raise FrontierFailure("test", "gpu-linalg", "probe report is incomplete or failed")
    for item in batches:
        residual = item["relative_residual"]
        if (
            item["info"] != 0
            or not isinstance(residual, float)
            or not math.isfinite(residual)
            or residual > residual_limit
        ):
            raise FrontierFailure(
                "test", "gpu-linalg", f"factor/solve failure in batch {item['batch']}"
            )
    return {
        "status": "passed",
        "device_count": device_count,
        "offloaded": True,
        "data_present": True,
        "residual_limit": residual_limit,
        "batches": batches,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_relative_path(value: str, context: str) -> PurePosixPath:
    path = PurePosixPath(value)
    raw_parts = value.split("/")
    if (
        not value
        or "\\" in value
        or path.is_absolute()
        or any(part in ("", ".", "..") for part in raw_parts)
    ):
        raise FrontierFailure("source", "manifest", f"unsafe {context} path")
    return path


def _source_tree_records_from_archive(archive: Path) -> list[tuple[str, str, int, str]]:
    records: list[tuple[str, str, int, str]] = []
    try:
        with tarfile.open(archive, "r:") as stream:
            for member in stream.getmembers():
                if member.isdir():
                    continue
                relative = _safe_relative_path(member.name, "source archive")
                mode = member.mode & 0o111
                if member.isfile():
                    extracted = stream.extractfile(member)
                    if extracted is None:
                        raise FrontierFailure(
                            "source", "source-archive", f"could not read {relative}"
                        )
                    digest = hashlib.sha256(extracted.read()).hexdigest()
                    records.append((relative.as_posix(), "file", mode, digest))
                elif member.issym():
                    records.append((relative.as_posix(), "symlink", mode, member.linkname))
                else:
                    raise FrontierFailure(
                        "source",
                        "source-archive",
                        f"unsupported archive member type for {relative}",
                    )
    except (OSError, tarfile.TarError) as error:
        raise FrontierFailure(
            "source", "source-archive", f"could not inspect source archive: {error}"
        ) from error
    return sorted(records)


def _source_tree_records_from_directory(
    source_root: Path,
) -> list[tuple[str, str, int, str]]:
    records: list[tuple[str, str, int, str]] = []
    for path in sorted(source_root.rglob("*")):
        relative = path.relative_to(source_root).as_posix()
        mode = path.lstat().st_mode & 0o111
        if path.is_symlink():
            records.append((relative, "symlink", mode, os.readlink(path)))
        elif path.is_file():
            records.append((relative, "file", mode, _sha256(path)))
    return sorted(records)


def _source_tree_sha256(records: Sequence[tuple[str, str, int, str]]) -> str:
    digest = hashlib.sha256()
    for path, kind, mode, payload in records:
        digest.update(f"{path}\0{kind}\0{mode:o}\0{payload}\n".encode("utf-8"))
    return digest.hexdigest()


def verify_source_binding(
    source_root: Path,
    artifact_root: Path,
    expected_source_sha: str,
    expected_archive_sha256: str,
) -> dict[str, object]:
    """Bind the extracted build tree to the staged git archive before compilation."""

    archive = artifact_root / "source.tar"
    if not archive.is_file() or _sha256(archive) != expected_archive_sha256:
        raise FrontierFailure(
            "source", "source-archive", "staged source archive SHA-256 differs"
        )
    try:
        with archive.open("rb") as stream:
            commit = subprocess.run(
                ["git", "get-tar-commit-id"],
                stdin=stream,
                capture_output=True,
                text=True,
                check=False,
            )
    except OSError as error:
        raise FrontierFailure(
            "source", "source-archive", f"could not inspect archive commit: {error}"
        ) from error
    archive_commit_sha = commit.stdout.strip()
    if commit.returncode != 0 or archive_commit_sha != expected_source_sha:
        raise FrontierFailure(
            "source", "source-archive", "source archive commit differs from requested SHA"
        )

    archive_records = _source_tree_records_from_archive(archive)
    extracted_records = _source_tree_records_from_directory(source_root)
    if archive_records != extracted_records:
        raise FrontierFailure(
            "source", "source-tree", "extracted source tree differs from verified archive"
        )
    evidence = {
        "archive_sha256": expected_archive_sha256,
        "archive_commit_sha": archive_commit_sha,
        "tree_sha256": _source_tree_sha256(archive_records),
        "verified_before_build": True,
    }
    (artifact_root / "source-verification.json").write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="utf-8"
    )
    return evidence


def inventory_regular_files(
    directory: Path,
    *,
    exclude: Sequence[Path] = (),
    exclude_trees: Sequence[Path] = (),
) -> list[dict[str, object]]:
    """Inventory every regular non-symlink file below a directory."""

    excluded = {path.as_posix() for path in exclude}
    excluded_trees = tuple(path.as_posix() for path in exclude_trees)
    inventory: list[dict[str, object]] = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(directory).as_posix()
        if relative in excluded or any(
            relative == tree or relative.startswith(tree + "/")
            for tree in excluded_trees
        ):
            continue
        inventory.append(
            {"path": relative, "size": path.stat().st_size, "sha256": _sha256(path)}
        )
    return inventory


def _record_command(
    command: Sequence[str],
    directory: Path,
    artifact_directory: Path,
    label: str,
    *,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
) -> tuple[subprocess.CompletedProcess[str], float]:
    artifact_directory.mkdir(parents=True, exist_ok=True)
    (artifact_directory / f"{label}.command.json").write_text(
        json.dumps(list(command), indent=2) + "\n", encoding="utf-8"
    )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            command,
            cwd=directory,
            env=None if environment is None else dict(environment),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        elapsed = time.monotonic() - started
        status = f"return_code={completed.returncode}"
    except subprocess.TimeoutExpired as error:
        elapsed = time.monotonic() - started
        completed = subprocess.CompletedProcess(
            command,
            124,
            "" if error.stdout is None else str(error.stdout),
            "" if error.stderr is None else str(error.stderr),
        )
        status = "timeout"
    (artifact_directory / f"{label}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (artifact_directory / f"{label}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    (artifact_directory / f"{label}.status.json").write_text(
        json.dumps({"status": status, "runtime_seconds": elapsed}, indent=2) + "\n",
        encoding="utf-8",
    )
    return completed, elapsed


def _make_command(
    source_root: Path,
    variables: Mapping[str, str],
    jobs: int,
    *targets: str,
) -> list[str]:
    return [
        "make",
        "-C",
        str(source_root),
        "-f",
        "test/qualification/frontier/Makefile",
        f"-j{jobs}",
        *(f"{name}={value}" for name, value in variables.items()),
        *targets,
    ]


def _neutralize_paths(value: str, source_root: Path, artifact_root: Path) -> str:
    replacements = {
        str(source_root.resolve()): "<source>",
        str(artifact_root.resolve()): "<artifact>",
    }
    for variable, marker in (
        ("ROCM_PATH", "<rocm>"),
        ("OLCF_HIPFORT_ROOT", "<hipfort>"),
    ):
        if os.environ.get(variable):
            replacements[os.environ[variable]] = marker
    normalized = value
    for path, marker in sorted(replacements.items(), key=lambda item: -len(item[0])):
        normalized = normalized.replace(path, marker)
    return normalized


def _parse_resolved_variables(
    text: str,
    names: Sequence[str],
    source_root: Path,
    artifact_root: Path,
) -> dict[str, str]:
    resolved: dict[str, str] = {}
    expected = set(names)
    for line in text.splitlines():
        if " = " not in line:
            continue
        name, value = line.split(" = ", 1)
        if name in expected:
            resolved[name] = _neutralize_paths(value, source_root, artifact_root)
    if set(resolved) != expected:
        missing = sorted(expected - set(resolved))
        raise FrontierFailure(
            "build", "build-evidence", "missing resolved variables: " + ", ".join(missing)
        )
    return resolved


def _linked_library_names(text: str) -> list[str]:
    libraries: list[str] = []
    for line in text.splitlines():
        fields = line.strip().split()
        if not fields:
            continue
        candidate = fields[0]
        if candidate.startswith("/"):
            candidate = Path(candidate).name
        if candidate not in libraries:
            libraries.append(candidate)
    return libraries


def _build_configuration(
    source_root: Path,
    artifact_root: Path,
    label: str,
    variables: Mapping[str, str],
    jobs: int,
    targets: Sequence[str],
) -> dict[str, object]:
    build_directory = source_root / "build" / variables["BUILD_NAME"]
    build_evidence_directory = artifact_root / "build" / label
    clean, _ = _record_command(
        _make_command(source_root, variables, 1, "clean"),
        source_root,
        build_evidence_directory,
        "clean",
        timeout_seconds=300.0,
    )
    if clean.returncode != 0:
        raise FrontierFailure("build", f"{label}-clean", "make clean failed")
    completed, runtime = _record_command(
        _make_command(source_root, variables, jobs, *targets),
        source_root,
        build_evidence_directory,
        "build",
        timeout_seconds=1800.0,
    )
    if completed.returncode != 0:
        raise FrontierFailure("build", f"{label}-build", "configuration build failed")
    resolved_names = (*variables, *RESOLVED_BUILD_VARIABLES)
    print_targets = [f"print-{name}" for name in resolved_names]
    resolved, _ = _record_command(
        _make_command(source_root, variables, 1, *print_targets),
        source_root,
        build_evidence_directory,
        "resolved-variables",
        timeout_seconds=120.0,
    )
    if resolved.returncode != 0:
        raise FrontierFailure("build", f"{label}-manifest", "could not resolve build variables")
    resolved_values = _parse_resolved_variables(
        resolved.stdout, resolved_names, source_root, artifact_root
    )

    binary_directory = artifact_root / "bin"
    binary_directory.mkdir(exist_ok=True)
    executables: dict[str, dict[str, object]] = {}
    for target in targets:
        source = build_directory / "bin" / target
        if not source.is_file() or not os.access(source, os.X_OK):
            raise FrontierFailure("build", f"{label}-build", f"missing executable {target}")
        destination_name = f"{target}-{label}"
        destination = binary_directory / destination_name
        shutil.copy2(source, destination)
        link, _ = _record_command(
            ["ldd", str(destination)],
            source_root,
            build_evidence_directory,
            f"{target}-link",
            timeout_seconds=120.0,
        )
        if link.returncode != 0:
            raise FrontierFailure("build", f"{label}-link", f"ldd failed for {target}")
        linked_libraries = _linked_library_names(link.stdout)
        if not linked_libraries:
            raise FrontierFailure(
                "build", f"{label}-link", f"ldd reported no libraries for {target}"
            )
        executables[target] = {
            "artifact": f"bin/{destination_name}",
            "size": destination.stat().st_size,
            "sha256": _sha256(destination),
            "link_evidence": f"build/{label}/{target}-link.stdout.txt",
            "linked_libraries": linked_libraries,
        }
    return {
        "status": "passed",
        "variables": dict(variables),
        "resolved_variables": f"build/{label}/resolved-variables.stdout.txt",
        "resolved_variable_values": resolved_values,
        "timeout_seconds": 1800.0,
        "runtime_seconds": runtime,
        "executables": executables,
    }


def _capture_environment(artifact_root: Path) -> dict[str, object]:
    environment_directory = artifact_root / "environment"
    environment_directory.mkdir(parents=True, exist_ok=True)
    modules = tuple(filter(None, os.environ.get("LOADEDMODULES", "").split(":")))
    missing = [
        marker for marker in REQUIRED_MODULE_MARKERS if not any(marker in item for item in modules)
    ]
    if missing:
        raise FrontierFailure(
            "environment", "modules", "required modules are missing: " + ", ".join(missing)
        )
    required_environment = ("ROCM_PATH", "OLCF_HIPFORT_ROOT")
    missing_environment = [name for name in required_environment if not os.environ.get(name)]
    if missing_environment:
        raise FrontierFailure(
            "environment",
            "modules",
            "required module variables are missing: " + ", ".join(missing_environment),
        )
    commands = {
        "compiler": ["ftn", "--version"],
        "preprocessor": ["cpp", "--version"],
        "rocm": ["rocminfo"],
        "gpu": ["rocm-smi", "--showproductname"],
    }
    command_summaries: dict[str, dict[str, object]] = {}
    for label, command in commands.items():
        completed, runtime = _record_command(
            command,
            REPOSITORY_ROOT,
            environment_directory,
            label,
            timeout_seconds=120.0,
        )
        if completed.returncode != 0:
            category = "allocation" if label in ("rocm", "gpu") else "environment"
            raise FrontierFailure(category, f"environment-{label}", f"{label} probe failed")
        command_summaries[label] = {
            "artifact": f"environment/{label}.stdout.txt",
            "sha256": _sha256(environment_directory / f"{label}.stdout.txt"),
            "runtime_seconds": runtime,
            "summary": next(
                (line.strip() for line in completed.stdout.splitlines() if line.strip()),
                "",
            ),
        }
    gpu_text = (environment_directory / "gpu.stdout.txt").read_text(encoding="utf-8")
    gpu_model = next(
        (line.strip() for line in gpu_text.splitlines() if "MI250" in line.upper()),
        "",
    )
    if not gpu_model:
        raise FrontierFailure("allocation", "environment-gpu", "MI250X GPU was not reported")
    return {
        "modules": list(modules),
        "compiler": command_summaries["compiler"],
        "preprocessor": command_summaries["preprocessor"],
        "rocm": command_summaries["rocm"],
        "gpu": command_summaries["gpu"],
        "gpu_model": gpu_model,
        "rocm_version": Path(os.environ["ROCM_PATH"]).name,
        "hipfort_module": next(item for item in modules if "hipfort" in item),
    }


def _input_inventory(source_root: Path) -> list[dict[str, object]]:
    partial = source_root / "test" / "qualification" / "parallel_zones"
    paths = [partial / "control"]
    for zone in EXPECTED_ZONES:
        paths.extend((partial / f"abundance_{zone:02d}", partial / f"thermo_{zone:02d}"))
    paths.extend(
        source_root / "test" / "Data_alpha" / name
        for name in ("sunet", "netsu", "netweak", "netwinv")
    )
    heat = heat_sn160_case(source_root)
    paths.extend((heat.control, heat.helm_table, *heat.trajectories))
    paths.extend(heat.network_data / name for name in heat.network_inputs)
    unique = sorted(set(path.resolve() for path in paths))
    missing = [path for path in unique if not path.is_file()]
    if missing:
        raise FrontierFailure("source", "input-inventory", f"missing input {missing[0]}")
    return [
        {
            "path": path.relative_to(source_root.resolve()).as_posix(),
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in unique
    ]


def _run_linalg_probe(executable: Path, artifact_root: Path) -> dict[str, object]:
    directory = artifact_root / "runs" / "gpu-linalg"
    timeout_seconds = 60.0
    environment = os.environ.copy()
    environment.update({"OMP_NUM_THREADS": "1", "OMP_TARGET_OFFLOAD": "MANDATORY"})
    completed, runtime = _record_command(
        [str(executable)],
        directory,
        directory,
        "probe",
        timeout_seconds=timeout_seconds,
        environment=environment,
    )
    if completed.returncode != 0:
        raise FrontierFailure(
            "test", "gpu-linalg", f"GPU linear-algebra probe returned {completed.returncode}"
        )
    result = parse_linalg_probe(completed.stdout)
    result["timeout_seconds"] = timeout_seconds
    result["runtime_seconds"] = runtime
    result["output_inventory"] = inventory_regular_files(directory)
    return result


def _run_partial_batch(
    cpu_executable: Path,
    gpu_executable: Path,
    artifact_root: Path,
    policy: NumericalPolicy,
) -> dict[str, object]:
    root = artifact_root / "runs" / "partial-batch"
    timeout_seconds = 180.0
    try:
        cpu_started = time.monotonic()
        cpu = run_configuration(
            "Frontier CPU",
            (cpu_executable,),
            root / "cpu",
            timeout_seconds,
        )
        cpu_runtime = time.monotonic() - cpu_started
        gpu_environment = os.environ.copy()
        gpu_environment.update({"OMP_NUM_THREADS": "1", "OMP_TARGET_OFFLOAD": "MANDATORY"})
        gpu_started = time.monotonic()
        gpu = run_configuration(
            "Frontier GPU",
            (gpu_executable,),
            root / "gpu",
            timeout_seconds,
            gpu_environment,
        )
        gpu_runtime = time.monotonic() - gpu_started
    except ParallelQualificationFailure as error:
        raise FrontierFailure("test", "partial-batch", str(error)) from error
    endpoint = compare_endpoint_states(gpu.states, cpu.states, policy, "partial-batch")
    ascii_result = compare_ascii_endpoints(
        gpu.ascii_endpoints, cpu.ascii_endpoints, policy, "partial-batch-ascii"
    )
    return {
        "status": "passed",
        "fixture": "ten distinguishable zones, nzbatchmx=4",
        "zones": list(EXPECTED_ZONES),
        "inactive_final_batch_lanes": 2,
        "timeout_seconds": timeout_seconds,
        "cpu_runtime_seconds": cpu_runtime,
        "gpu_runtime_seconds": gpu_runtime,
        "endpoint_comparison": endpoint,
        "ascii_comparison": ascii_result,
        "cpu_output_inventory": inventory_regular_files(
            root / "cpu", exclude=(Path("control"),)
        ),
        "gpu_output_inventory": inventory_regular_files(
            root / "gpu", exclude=(Path("control"),)
        ),
    }


def _run_heat_configuration(
    executable: Path,
    source_root: Path,
    directory: Path,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
) -> tuple[tuple[FinalState, ...], tuple[AsciiEndpoint, ...], float]:
    case = heat_sn160_case(source_root)
    try:
        prepared = prepare_regression_work_directory(case, directory)
        started = time.monotonic()
        run_xnet(
            executable,
            case,
            prepared,
            timeout_seconds=timeout_seconds,
            environment=environment,
        )
        runtime = time.monotonic() - started
        diagnostic_path = prepared / "net_diag01"
        diagnostic = diagnostic_path.read_text(encoding="utf-8")
        states = parse_diagnostic(
            diagnostic,
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
        ascii_endpoints = validate_ascii_association(
            prepared,
            states,
            filename_root="ev_heat_sn160_",
            output_species=HEAT_SN160_ASCII_SPECIES,
            zone_width=1,
        )
        return states, ascii_endpoints, runtime
    except (OSError, UnicodeError, RegressionFailure, ParallelQualificationFailure) as error:
        raise FrontierFailure("test", "heat-sn160", str(error)) from error


def _run_heat_sn160(
    cpu_executable: Path,
    gpu_executable: Path,
    source_root: Path,
    artifact_root: Path,
    policy: NumericalPolicy,
) -> dict[str, object]:
    root = artifact_root / "runs" / "heat-sn160"
    timeout_seconds = 600.0
    cpu, cpu_ascii, cpu_runtime = _run_heat_configuration(
        cpu_executable, source_root, root / "cpu", timeout_seconds
    )
    gpu_environment = os.environ.copy()
    gpu_environment.update({"OMP_NUM_THREADS": "1", "OMP_TARGET_OFFLOAD": "MANDATORY"})
    gpu, gpu_ascii, gpu_runtime = _run_heat_configuration(
        gpu_executable,
        source_root,
        root / "gpu",
        timeout_seconds,
        gpu_environment,
    )
    endpoint = compare_endpoint_states(gpu, cpu, policy, "heat-sn160")
    ascii_result = compare_ascii_endpoints(
        gpu_ascii, cpu_ascii, policy, "heat-sn160-ascii"
    )
    nonzero_neutrino_loss_zones = [
        endpoint.zone for endpoint in cpu_ascii if endpoint.neutrino_loss_rate != 0.0
    ]
    if not nonzero_neutrino_loss_zones:
        raise FrontierFailure(
            "test",
            "heat-sn160-ascii",
            "fixture did not exercise a nonzero neutrino-loss diagnostic",
        )
    return {
        "status": "passed",
        "case": "heat_sn160",
        "zones": [state.zone for state in cpu],
        "timeout_seconds": timeout_seconds,
        "cpu_runtime_seconds": cpu_runtime,
        "gpu_runtime_seconds": gpu_runtime,
        "endpoint_comparison": endpoint,
        "ascii_comparison": ascii_result,
        "nonzero_neutrino_loss_zones": nonzero_neutrino_loss_zones,
        "cpu_output_inventory": inventory_regular_files(
            root / "cpu", exclude=(Path("control"),)
        ),
        "gpu_output_inventory": inventory_regular_files(
            root / "gpu", exclude=(Path("control"),)
        ),
    }


def _slurm_manifest(time_limit: str) -> dict[str, object]:
    return {
        "job_id": os.environ.get("SLURM_JOB_ID", "unknown"),
        "account_supplied": True,
        "partition": os.environ.get("SLURM_JOB_PARTITION", "unknown"),
        "qos_supplied": bool(os.environ.get("SLURM_JOB_QOS")),
        "reservation_supplied": bool(os.environ.get("SLURM_JOB_RESERVATION")),
        "nodes": int(os.environ.get("SLURM_JOB_NUM_NODES", "1")),
        "tasks": int(os.environ.get("SLURM_NTASKS", "1")),
        "cpus_per_task": int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
        "gpus_per_task": 1,
        "time_limit": time_limit,
    }


def run_qualification(arguments: argparse.Namespace) -> Path:
    source_root = arguments.source_root.resolve()
    artifact_root = arguments.artifact_root.resolve()
    artifact_root.mkdir(parents=True, exist_ok=True)
    report_path = artifact_root / "qualification_manifest.json"
    started = datetime.now(timezone.utc)
    report: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "status": "failed",
        "failure": None,
        "source": {
            "sha": arguments.source_sha,
            "worktree_clean": True,
            "archive_sha256": arguments.archive_sha256,
            "archive_commit_sha": arguments.source_sha,
            "tree_sha256": None,
            "verified_before_build": False,
        },
        "environment": {},
        "slurm": _slurm_manifest(arguments.time_limit),
        "builds": {},
        "inputs": [],
        "checks": {},
        "artifact_inventory": [],
        "started_at_utc": started.isoformat(),
        "finished_at_utc": None,
        "runtime_seconds": None,
    }
    monotonic_start = time.monotonic()
    try:
        if not SOURCE_SHA_PATTERN.fullmatch(arguments.source_sha):
            raise FrontierFailure("source", "source-sha", "source SHA is not a full commit")
        if not SHA256_PATTERN.fullmatch(arguments.archive_sha256):
            raise FrontierFailure("source", "source-archive", "archive SHA-256 is invalid")
        source_evidence = verify_source_binding(
            source_root,
            artifact_root,
            arguments.source_sha,
            arguments.archive_sha256,
        )
        report["source"] = {
            "sha": arguments.source_sha,
            "worktree_clean": True,
            **source_evidence,
        }
        policy = load_policy(
            source_root
            / "test"
            / "qualification"
            / "frontier"
            / "comparison_policy.json"
        )
        report["environment"] = _capture_environment(artifact_root)
        report["inputs"] = _input_inventory(source_root)
        builds = report["builds"]
        assert isinstance(builds, dict)
        builds["cpu"] = _build_configuration(
            source_root,
            artifact_root,
            "cpu",
            CPU_BUILD_VARIABLES,
            arguments.build_jobs,
            ("xnet",),
        )
        builds["gpu"] = _build_configuration(
            source_root,
            artifact_root,
            "gpu",
            GPU_BUILD_VARIABLES,
            arguments.build_jobs,
            ("xnet", "frontier_gpu_linalg_probe"),
        )
        cpu_executable = artifact_root / "bin" / "xnet-cpu"
        gpu_executable = artifact_root / "bin" / "xnet-gpu"
        probe_executable = artifact_root / "bin" / "frontier_gpu_linalg_probe-gpu"
        checks = report["checks"]
        assert isinstance(checks, dict)
        checks["gpu_linalg"] = _run_linalg_probe(probe_executable, artifact_root)
        checks["partial_batch"] = _run_partial_batch(
            cpu_executable, gpu_executable, artifact_root, policy
        )
        checks["heat_sn160"] = _run_heat_sn160(
            cpu_executable, gpu_executable, source_root, artifact_root, policy
        )
        report["status"] = "passed"
    except FrontierFailure as error:
        report["failure"] = {
            "category": error.category,
            "phase": error.phase,
            "message": str(error),
        }
    except (OSError, ValueError) as error:
        report["failure"] = {
            "category": "test",
            "phase": "qualification-runner",
            "message": str(error),
        }
    finally:
        report["artifact_inventory"] = inventory_regular_files(
            artifact_root,
            exclude=(Path("qualification_manifest.json"), Path("source.tar")),
            exclude_trees=(Path("source"),),
        )
        report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        report["runtime_seconds"] = time.monotonic() - monotonic_start
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report_path


def write_failure_manifest(arguments: argparse.Namespace) -> Path:
    """Write classified evidence when Slurm cannot start the qualification step."""

    artifact_root = arguments.artifact_root.resolve()
    artifact_root.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    report = {
        "schema": MANIFEST_SCHEMA,
        "status": "failed",
        "failure": {
            "category": arguments.category,
            "phase": arguments.phase,
            "message": arguments.message,
        },
        "source": {
            "sha": arguments.source_sha,
            "worktree_clean": True,
            "archive_sha256": arguments.archive_sha256,
            "archive_commit_sha": arguments.source_sha,
            "tree_sha256": None,
            "verified_before_build": False,
        },
        "environment": {},
        "slurm": _slurm_manifest(arguments.time_limit),
        "builds": {},
        "inputs": [],
        "checks": {},
        "artifact_inventory": inventory_regular_files(
            artifact_root,
            exclude=(Path("qualification_manifest.json"), Path("source.tar")),
            exclude_trees=(Path("source"),),
        ),
        "started_at_utc": now,
        "finished_at_utc": now,
        "runtime_seconds": 0.0,
    }
    report_path = artifact_root / "qualification_manifest.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report_path


def _manifest_mapping(
    value: object, fields: set[str], context: str, category: str = "test"
) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != fields:
        raise FrontierFailure(category, "manifest", f"{context} fields are incomplete")
    return value


def _manifest_number(
    value: object, context: str, *, positive: bool = False
) -> float:
    if type(value) not in (int, float):
        raise FrontierFailure("test", "manifest", f"{context} is not numeric")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0 or (positive and normalized <= 0.0):
        raise FrontierFailure("test", "manifest", f"{context} is invalid")
    return normalized


def _manifest_close(actual: float, expected: float) -> bool:
    return math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-300)


def _validate_bounded_values(
    value: object,
    bounds: Bounds,
    context: str,
) -> float:
    evidence = _manifest_mapping(
        value,
        {"observed", "expected", "difference", "allowed"},
        context,
        "comparison",
    )
    observed = _manifest_number(evidence["observed"], f"{context} observed")
    expected = _manifest_number(evidence["expected"], f"{context} expected")
    difference = _manifest_number(evidence["difference"], f"{context} difference")
    allowed = _manifest_number(evidence["allowed"], f"{context} allowed")
    recalculated_difference = abs(observed - expected)
    recalculated_allowed = bounds.atol + bounds.rtol * abs(expected)
    if (
        not _manifest_close(difference, recalculated_difference)
        or not _manifest_close(allowed, recalculated_allowed)
        or difference > recalculated_allowed
    ):
        raise FrontierFailure("comparison", "manifest", f"{context} failed")
    if recalculated_allowed == 0.0:
        return 0.0
    return difference / recalculated_allowed


def _exact_integer_sequence(value: object, expected: Sequence[int]) -> bool:
    return (
        isinstance(value, list)
        and len(value) == len(expected)
        and all(type(item) is int and item == reference for item, reference in zip(value, expected, strict=False))
    )


def _validate_inventory(
    value: object, context: str, *, allow_empty: bool = False
) -> dict[str, tuple[int, str]]:
    if not isinstance(value, list) or (not value and not allow_empty):
        raise FrontierFailure("test", "manifest", f"{context} is empty")
    artifacts: dict[str, tuple[int, str]] = {}
    for item in value:
        artifact = _manifest_mapping(
            item, {"path", "size", "sha256"}, f"{context} item"
        )
        path = artifact["path"]
        if not isinstance(path, str):
            raise FrontierFailure("test", "manifest", f"invalid {context} path")
        _safe_relative_path(path, context)
        if path in artifacts:
            raise FrontierFailure("test", "manifest", f"duplicate {context} path")
        if type(artifact["size"]) is not int or artifact["size"] < 0:
            raise FrontierFailure("test", "manifest", f"invalid {context} size")
        if not SHA256_PATTERN.fullmatch(str(artifact["sha256"])):
            raise FrontierFailure("test", "manifest", f"invalid {context} hash")
        artifacts[path] = (artifact["size"], str(artifact["sha256"]))
    return artifacts


def _validate_endpoint_evidence(
    value: object,
    expected_zones: Sequence[int],
    policy: NumericalPolicy,
    context: str,
) -> None:
    evidence = _manifest_mapping(
        value,
        {
            "status",
            "policy_status",
            "maximum_selected_fraction_of_allowed",
            "zones",
        },
        context,
    )
    maximum = _manifest_number(
        evidence["maximum_selected_fraction_of_allowed"], f"{context} maximum"
    )
    if evidence["status"] != "passed" or evidence["policy_status"] != policy.status or maximum > 1.0:
        raise FrontierFailure("comparison", "manifest", f"{context} did not pass policy")
    zones = evidence["zones"]
    if (
        not isinstance(zones, list)
        or len(zones) != len(expected_zones)
        or any(
            not isinstance(item, dict)
            or type(item.get("zone")) is not int
            or item["zone"] != expected
            for item, expected in zip(zones, expected_zones, strict=False)
        )
    ):
        raise FrontierFailure("comparison", "manifest", f"{context} zones are incomplete")
    recalculated_maximum = 0.0
    for zone in zones:
        observation = _manifest_mapping(
            zone,
            {
                "zone",
                "scalar_differences",
                "selected_species",
                "selected_differences",
                "observed_mass_fractions",
                "reference_mass_fractions",
                "composition_l1",
                "composition_linf",
                "composition_linf_species",
            },
            f"{context} zone",
        )
        scalar = _manifest_mapping(
            observation["scalar_differences"],
            set(policy.scalar_fields),
            f"{context} scalar differences",
        )
        for name, difference in scalar.items():
            _validate_bounded_values(
                difference, policy.scalar_fields[name], f"{context} {name}"
            )
        observed_fractions = observation["observed_mass_fractions"]
        reference_fractions = observation["reference_mass_fractions"]
        if (
            not isinstance(observed_fractions, dict)
            or not isinstance(reference_fractions, dict)
            or not reference_fractions
            or set(observed_fractions) != set(reference_fractions)
            or any(not isinstance(name, str) or not name for name in reference_fractions)
        ):
            raise FrontierFailure(
                "comparison", "manifest", f"{context} composition inventory is invalid"
            )
        normalized_observed = {
            name: _manifest_number(value, f"{context} observed {name}")
            for name, value in observed_fractions.items()
        }
        normalized_reference = {
            name: _manifest_number(value, f"{context} reference {name}")
            for name, value in reference_fractions.items()
        }
        expected_selected = [
            name
            for name, value in normalized_reference.items()
            if name in policy.anchors or value >= policy.material_threshold
        ]
        selected = observation["selected_species"]
        if (
            not isinstance(selected, list)
            or selected != expected_selected
        ):
            raise FrontierFailure("comparison", "manifest", f"{context} species are invalid")
        selected_differences = _manifest_mapping(
            observation["selected_differences"],
            set(selected),
            f"{context} selected differences",
            "comparison",
        )
        for name in selected:
            selected_evidence = selected_differences[name]
            values = _manifest_mapping(
                selected_evidence,
                {"observed", "expected", "difference", "allowed"},
                f"{context} selected {name}",
                "comparison",
            )
            observed = _manifest_number(
                values["observed"], f"{context} selected {name} observed"
            )
            expected = _manifest_number(
                values["expected"], f"{context} selected {name} expected"
            )
            if (
                not _manifest_close(observed, normalized_observed[name])
                or not _manifest_close(expected, normalized_reference[name])
            ):
                raise FrontierFailure(
                    "comparison", "manifest", f"{context} selected {name} is inconsistent"
                )
            recalculated_maximum = max(
                recalculated_maximum,
                _validate_bounded_values(
                    selected_evidence, policy.selected, f"{context} selected {name}"
                ),
            )
        vector_differences = {
            name: abs(normalized_observed[name] - expected)
            for name, expected in normalized_reference.items()
        }
        recalculated_l1 = math.fsum(vector_differences.values())
        recalculated_linf_species = max(
            vector_differences, key=vector_differences.__getitem__
        )
        recalculated_linf = vector_differences[recalculated_linf_species]
        l1 = _manifest_number(observation["composition_l1"], f"{context} L1")
        linf = _manifest_number(observation["composition_linf"], f"{context} Linf")
        if (
            not _manifest_close(l1, recalculated_l1)
            or not _manifest_close(linf, recalculated_linf)
            or observation["composition_linf_species"] != recalculated_linf_species
            or l1 > policy.complete_vector_l1_limit
            or linf > policy.complete_vector_linf_limit
        ):
            raise FrontierFailure("comparison", "manifest", f"{context} composition failed")
        for label, fractions in (
            ("observed", normalized_observed),
            ("reference", normalized_reference),
        ):
            normalization_error = abs(math.fsum(fractions.values()) - 1.0)
            if normalization_error > policy.normalization_atol:
                raise FrontierFailure(
                    "comparison", "manifest", f"{context} {label} normalization failed"
                )
    if not _manifest_close(maximum, recalculated_maximum):
        raise FrontierFailure("comparison", "manifest", f"{context} maximum is inconsistent")


def _validate_ascii_evidence(
    value: object,
    expected_zones: Sequence[int],
    policy: NumericalPolicy,
    context: str,
) -> None:
    evidence = _manifest_mapping(
        value,
        {"status", "maximum_fraction_of_allowed", "reported_only_fields", "zones"},
        context,
    )
    maximum = _manifest_number(evidence["maximum_fraction_of_allowed"], f"{context} maximum")
    if (
        evidence["status"] != "passed"
        or maximum > 1.0
        or evidence["reported_only_fields"] != list(policy.reported_ascii_fields)
    ):
        raise FrontierFailure("comparison", "manifest", f"{context} did not pass policy")
    zones = evidence["zones"]
    if (
        not isinstance(zones, list)
        or len(zones) != len(expected_zones)
        or any(
            not isinstance(item, dict)
            or type(item.get("zone")) is not int
            or item["zone"] != expected
            for item, expected in zip(zones, expected_zones, strict=False)
        )
    ):
        raise FrontierFailure("comparison", "manifest", f"{context} zones are incomplete")
    expected_fields = set(policy.ascii_fields) | set(policy.reported_ascii_fields)
    recalculated_maximum = 0.0
    for zone in zones:
        observation = _manifest_mapping(zone, {"zone", "field_differences"}, f"{context} zone")
        fields = _manifest_mapping(
            observation["field_differences"], expected_fields, f"{context} fields"
        )
        for name in policy.ascii_fields:
            bounded = _manifest_mapping(
                fields[name],
                {"comparison", "observed", "expected", "difference", "allowed"},
                f"{context} {name}",
            )
            if bounded["comparison"] != "bounded":
                raise FrontierFailure("comparison", "manifest", f"{context} {name} failed")
            recalculated_maximum = max(
                recalculated_maximum,
                _validate_bounded_values(
                    {key: bounded[key] for key in ("observed", "expected", "difference", "allowed")},
                    policy.ascii_fields[name],
                    f"{context} {name}",
                ),
            )
        for name in policy.reported_ascii_fields:
            reported = _manifest_mapping(
                fields[name],
                {"comparison", "observed", "expected", "difference"},
                f"{context} {name}",
            )
            observed_raw = reported["observed"]
            expected_raw = reported["expected"]
            if type(observed_raw) not in (int, float) or type(expected_raw) not in (int, float):
                raise FrontierFailure("comparison", "manifest", f"{context} {name} is not numeric")
            observed = float(observed_raw)
            expected = float(expected_raw)
            difference = _manifest_number(reported["difference"], f"{context} {name} difference")
            if (
                reported["comparison"] != "reported_only"
                or not math.isfinite(observed)
                or not math.isfinite(expected)
                or not _manifest_close(difference, abs(observed - expected))
            ):
                raise FrontierFailure("comparison", "manifest", f"{context} {name} is invalid")
    if not _manifest_close(maximum, recalculated_maximum):
        raise FrontierFailure("comparison", "manifest", f"{context} maximum is inconsistent")


def validate_manifest(document: object, *, require_pass: bool = True) -> None:
    """Validate every required result from a successful Frontier run."""

    required = {
        "schema", "status", "failure", "source", "environment", "slurm",
        "builds", "inputs", "checks", "artifact_inventory", "started_at_utc",
        "finished_at_utc", "runtime_seconds",
    }
    manifest = _manifest_mapping(document, required, "manifest", "source")
    if manifest["schema"] != MANIFEST_SCHEMA:
        raise FrontierFailure("source", "manifest", "manifest schema differs")
    if manifest["status"] not in ("passed", "failed"):
        raise FrontierFailure("test", "manifest", "qualification status is invalid")
    if require_pass and manifest["status"] != "passed":
        raise FrontierFailure("test", "manifest", "qualification status is not passed")

    source = _manifest_mapping(
        manifest["source"],
        {
            "sha", "worktree_clean", "archive_sha256", "archive_commit_sha",
            "tree_sha256", "verified_before_build",
        },
        "source evidence",
        "source",
    )
    if (
        not SOURCE_SHA_PATTERN.fullmatch(str(source["sha"]))
        or source["worktree_clean"] is not True
        or not SHA256_PATTERN.fullmatch(str(source["archive_sha256"]))
        or not SOURCE_SHA_PATTERN.fullmatch(str(source["archive_commit_sha"]))
        or not isinstance(source["verified_before_build"], bool)
        or (source["tree_sha256"] is not None and not SHA256_PATTERN.fullmatch(str(source["tree_sha256"])))
    ):
        raise FrontierFailure("source", "manifest", "source evidence is incomplete")

    _manifest_number(manifest["runtime_seconds"], "runtime seconds")
    for field in ("started_at_utc", "finished_at_utc"):
        if not isinstance(manifest[field], str):
            raise FrontierFailure("test", "manifest", f"{field} is invalid")
        try:
            timestamp = datetime.fromisoformat(manifest[field])
        except ValueError as error:
            raise FrontierFailure("test", "manifest", f"{field} is invalid") from error
        if timestamp.tzinfo is None:
            raise FrontierFailure("test", "manifest", f"{field} lacks a timezone")

    if manifest["status"] == "failed":
        failure = _manifest_mapping(
            manifest["failure"], {"category", "phase", "message"}, "failure classification"
        )
        if (
            failure["category"] not in FAILURE_CATEGORIES
            or not all(isinstance(failure[name], str) and failure[name] for name in failure)
        ):
            raise FrontierFailure("test", "manifest", "failure classification is incomplete")
        return
    if manifest["failure"] is not None:
        raise FrontierFailure("test", "manifest", "passed manifest contains a failure")
    if (
        source["verified_before_build"] is not True
        or source["archive_commit_sha"] != source["sha"]
        or not SHA256_PATTERN.fullmatch(str(source["tree_sha256"]))
    ):
        raise FrontierFailure("source", "manifest", "source binding was not verified")

    artifact_claims: list[tuple[str, int | None, str, str]] = []
    environment = _manifest_mapping(
        manifest["environment"],
        {"modules", "compiler", "preprocessor", "rocm", "gpu", "gpu_model", "rocm_version", "hipfort_module"},
        "environment evidence",
        "environment",
    )
    modules = environment["modules"]
    if not isinstance(modules, list) or any(not isinstance(item, str) or not item for item in modules):
        raise FrontierFailure("environment", "manifest", "module evidence is incomplete")
    if any(not any(marker in item for item in modules) for marker in REQUIRED_MODULE_MARKERS):
        raise FrontierFailure("environment", "manifest", "module evidence is incomplete")
    for label in ("compiler", "preprocessor", "rocm", "gpu"):
        summary = _manifest_mapping(
            environment[label], {"artifact", "sha256", "runtime_seconds", "summary"}, f"{label} evidence", "environment"
        )
        if not isinstance(summary["artifact"], str) or not isinstance(summary["summary"], str) or not summary["summary"].strip():
            raise FrontierFailure("environment", "manifest", f"{label} summary is incomplete")
        _safe_relative_path(summary["artifact"], f"{label} evidence")
        if not SHA256_PATTERN.fullmatch(str(summary["sha256"])):
            raise FrontierFailure("environment", "manifest", f"{label} hash is incomplete")
        artifact_claims.append(
            (
                summary["artifact"],
                None,
                str(summary["sha256"]),
                f"{label} evidence",
            )
        )
        _manifest_number(summary["runtime_seconds"], f"{label} runtime")
    if not all(isinstance(environment[name], str) and environment[name] for name in ("gpu_model", "rocm_version", "hipfort_module")):
        raise FrontierFailure("environment", "manifest", "environment summary is incomplete")

    slurm = _manifest_mapping(
        manifest["slurm"],
        {"job_id", "account_supplied", "partition", "qos_supplied", "reservation_supplied", "nodes", "tasks", "cpus_per_task", "gpus_per_task", "time_limit"},
        "Slurm evidence",
        "allocation",
    )
    if (
        not isinstance(slurm["job_id"], str) or not slurm["job_id"].isdigit()
        or slurm["account_supplied"] is not True
        or not isinstance(slurm["partition"], str) or not slurm["partition"] or slurm["partition"] == "unknown"
        or not isinstance(slurm["qos_supplied"], bool)
        or not isinstance(slurm["reservation_supplied"], bool)
        or type(slurm["nodes"]) is not int or slurm["nodes"] != 1
        or type(slurm["tasks"]) is not int or slurm["tasks"] != 1
        or type(slurm["cpus_per_task"]) is not int or slurm["cpus_per_task"] < 1
        or type(slurm["gpus_per_task"]) is not int or slurm["gpus_per_task"] != 1
        or not isinstance(slurm["time_limit"], str) or not SLURM_TIME_PATTERN.fullmatch(slurm["time_limit"])
    ):
        raise FrontierFailure("allocation", "manifest", "Slurm evidence is incomplete")

    builds = _manifest_mapping(manifest["builds"], {"cpu", "gpu"}, "build evidence", "build")
    required_artifacts = {"source-verification.json"}
    for label, expected in (("cpu", CPU_BUILD_VARIABLES), ("gpu", GPU_BUILD_VARIABLES)):
        build = _manifest_mapping(
            builds[label],
            {"status", "variables", "resolved_variables", "resolved_variable_values", "timeout_seconds", "runtime_seconds", "executables"},
            f"{label} build evidence",
            "build",
        )
        if build["status"] != "passed" or build["variables"] != expected:
            raise FrontierFailure("build", "manifest", f"{label} build variables differ")
        resolved_path = build["resolved_variables"]
        if not isinstance(resolved_path, str):
            raise FrontierFailure("build", "manifest", f"{label} resolved-variable path is invalid")
        _safe_relative_path(resolved_path, f"{label} resolved variables")
        required_artifacts.add(resolved_path)
        resolved = build["resolved_variable_values"]
        expected_resolved = set(expected) | set(RESOLVED_BUILD_VARIABLES)
        if not isinstance(resolved, dict) or set(resolved) != expected_resolved or any(resolved[name] != value for name, value in expected.items()):
            raise FrontierFailure("build", "manifest", f"{label} resolved variables differ")
        if not all(isinstance(resolved[name], str) for name in RESOLVED_BUILD_VARIABLES):
            raise FrontierFailure("build", "manifest", f"{label} compiler evidence is invalid")
        if label == "gpu" and resolved["CRAY_OMP_PREPROCESS"] != "yes":
            raise FrontierFailure("build", "manifest", "Cray GPU preprocessing is not selected")
        _manifest_number(build["timeout_seconds"], f"{label} build timeout", positive=True)
        _manifest_number(build["runtime_seconds"], f"{label} build runtime")
        executables = build["executables"]
        if not isinstance(executables, dict) or tuple(executables) != EXPECTED_EXECUTABLES[label]:
            raise FrontierFailure("build", "manifest", f"{label} executable set differs")
        for target, value in executables.items():
            executable = _manifest_mapping(
                value, {"artifact", "size", "sha256", "link_evidence", "linked_libraries"}, f"{label} {target} executable", "build"
            )
            if not isinstance(executable["artifact"], str) or not isinstance(executable["link_evidence"], str):
                raise FrontierFailure("build", "manifest", f"{label} executable paths are invalid")
            _safe_relative_path(executable["artifact"], f"{label} executable")
            _safe_relative_path(executable["link_evidence"], f"{label} link evidence")
            required_artifacts.update((executable["artifact"], executable["link_evidence"]))
            if type(executable["size"]) is not int or executable["size"] <= 0 or not SHA256_PATTERN.fullmatch(str(executable["sha256"])):
                raise FrontierFailure("build", "manifest", f"{label} executable evidence is invalid")
            artifact_claims.append(
                (
                    executable["artifact"],
                    executable["size"],
                    str(executable["sha256"]),
                    f"{label} {target} executable",
                )
            )
            libraries = executable["linked_libraries"]
            if not isinstance(libraries, list) or not libraries or any(not isinstance(name, str) or not name for name in libraries) or len(set(libraries)) != len(libraries):
                raise FrontierFailure("build", "manifest", f"{label} link summary is incomplete")

    checks = _manifest_mapping(manifest["checks"], {"gpu_linalg", "partial_batch", "heat_sn160"}, "check evidence")
    linalg = _manifest_mapping(
        checks["gpu_linalg"],
        {"status", "device_count", "offloaded", "data_present", "residual_limit", "batches", "timeout_seconds", "runtime_seconds", "output_inventory"},
        "GPU linear-algebra evidence",
    )
    residual_limit = _manifest_number(linalg["residual_limit"], "residual limit", positive=True)
    if residual_limit != GPU_LINALG_RESIDUAL_LIMIT:
        raise FrontierFailure(
            "test",
            "manifest",
            f"GPU residual limit must be {GPU_LINALG_RESIDUAL_LIMIT:.1e}",
        )
    if linalg["status"] != "passed" or type(linalg["device_count"]) is not int or linalg["device_count"] < 1 or linalg["offloaded"] is not True or linalg["data_present"] is not True:
        raise FrontierFailure("test", "manifest", "GPU linear-algebra evidence failed")
    batches = linalg["batches"]
    if (
        not isinstance(batches, list)
        or len(batches) != 2
        or any(
            not isinstance(item, dict)
            or type(item.get("batch")) is not int
            or item["batch"] != expected
            for item, expected in zip(batches, (1, 2), strict=False)
        )
    ):
        raise FrontierFailure("test", "manifest", "GPU batch evidence is incomplete")
    for batch in batches:
        result = _manifest_mapping(batch, {"batch", "info", "relative_residual"}, "GPU batch result")
        residual = _manifest_number(result["relative_residual"], "GPU relative residual")
        if type(result["info"]) is not int or result["info"] != 0 or residual > residual_limit:
            raise FrontierFailure("test", "manifest", "GPU factor/solve evidence failed")
    _manifest_number(linalg["timeout_seconds"], "GPU probe timeout", positive=True)
    _manifest_number(linalg["runtime_seconds"], "GPU probe runtime")
    probe_paths = _validate_inventory(linalg["output_inventory"], "GPU probe inventory")
    required_artifacts.update(f"runs/gpu-linalg/{path}" for path in probe_paths)
    artifact_claims.extend(
        (f"runs/gpu-linalg/{path}", size, sha256, "GPU probe inventory")
        for path, (size, sha256) in probe_paths.items()
    )

    policy = load_policy(REPOSITORY_ROOT / "test" / "qualification" / "frontier" / "comparison_policy.json")
    partial = _manifest_mapping(
        checks["partial_batch"],
        {"status", "fixture", "zones", "inactive_final_batch_lanes", "timeout_seconds", "cpu_runtime_seconds", "gpu_runtime_seconds", "endpoint_comparison", "ascii_comparison", "cpu_output_inventory", "gpu_output_inventory"},
        "partial-batch evidence",
    )
    if (
        partial["status"] != "passed"
        or not _exact_integer_sequence(partial["zones"], EXPECTED_ZONES)
        or type(partial["inactive_final_batch_lanes"]) is not int
        or partial["inactive_final_batch_lanes"] != 2
        or not isinstance(partial["fixture"], str)
        or not partial["fixture"]
    ):
        raise FrontierFailure("test", "manifest", "partial-batch evidence failed")
    _validate_endpoint_evidence(partial["endpoint_comparison"], EXPECTED_ZONES, policy, "partial-batch endpoint")
    _validate_ascii_evidence(partial["ascii_comparison"], EXPECTED_ZONES, policy, "partial-batch ASCII")
    for name in ("timeout_seconds", "cpu_runtime_seconds", "gpu_runtime_seconds"):
        _manifest_number(partial[name], f"partial-batch {name}", positive=name == "timeout_seconds")
    partial_cpu_paths = _validate_inventory(
        partial["cpu_output_inventory"], "partial-batch CPU inventory"
    )
    partial_gpu_paths = _validate_inventory(
        partial["gpu_output_inventory"], "partial-batch GPU inventory"
    )
    required_artifacts.update(
        f"runs/partial-batch/cpu/{path}" for path in partial_cpu_paths
    )
    required_artifacts.update(
        f"runs/partial-batch/gpu/{path}" for path in partial_gpu_paths
    )
    artifact_claims.extend(
        (f"runs/partial-batch/cpu/{path}", size, sha256, "partial-batch CPU inventory")
        for path, (size, sha256) in partial_cpu_paths.items()
    )
    artifact_claims.extend(
        (f"runs/partial-batch/gpu/{path}", size, sha256, "partial-batch GPU inventory")
        for path, (size, sha256) in partial_gpu_paths.items()
    )

    heat_zones = list(range(1, 7))
    heat = _manifest_mapping(
        checks["heat_sn160"],
        {"status", "case", "zones", "timeout_seconds", "cpu_runtime_seconds", "gpu_runtime_seconds", "endpoint_comparison", "ascii_comparison", "nonzero_neutrino_loss_zones", "cpu_output_inventory", "gpu_output_inventory"},
        "heat_sn160 evidence",
    )
    nonzero_zones = heat["nonzero_neutrino_loss_zones"]
    if (
        heat["status"] != "passed"
        or heat["case"] != "heat_sn160"
        or not _exact_integer_sequence(heat["zones"], heat_zones)
        or not isinstance(nonzero_zones, list)
        or not nonzero_zones
        or any(type(zone) is not int for zone in nonzero_zones)
        or len(set(nonzero_zones)) != len(nonzero_zones)
        or not set(nonzero_zones).issubset(heat_zones)
    ):
        raise FrontierFailure("test", "manifest", "heat_sn160 evidence failed")
    _validate_endpoint_evidence(heat["endpoint_comparison"], heat_zones, policy, "heat_sn160 endpoint")
    _validate_ascii_evidence(heat["ascii_comparison"], heat_zones, policy, "heat_sn160 ASCII")
    for name in ("timeout_seconds", "cpu_runtime_seconds", "gpu_runtime_seconds"):
        _manifest_number(heat[name], f"heat_sn160 {name}", positive=name == "timeout_seconds")
    heat_ascii_zones = {
        item["zone"]: item["field_differences"]["neutrino_loss_rate"]["expected"]
        for item in heat["ascii_comparison"]["zones"]
    }
    if any(heat_ascii_zones[zone] == 0.0 for zone in nonzero_zones):
        raise FrontierFailure(
            "test", "manifest", "nonzero neutrino-loss evidence is inconsistent"
        )
    heat_cpu_paths = _validate_inventory(
        heat["cpu_output_inventory"], "heat_sn160 CPU inventory"
    )
    heat_gpu_paths = _validate_inventory(
        heat["gpu_output_inventory"], "heat_sn160 GPU inventory"
    )
    required_artifacts.update(
        f"runs/heat-sn160/cpu/{path}" for path in heat_cpu_paths
    )
    required_artifacts.update(
        f"runs/heat-sn160/gpu/{path}" for path in heat_gpu_paths
    )
    artifact_claims.extend(
        (f"runs/heat-sn160/cpu/{path}", size, sha256, "heat_sn160 CPU inventory")
        for path, (size, sha256) in heat_cpu_paths.items()
    )
    artifact_claims.extend(
        (f"runs/heat-sn160/gpu/{path}", size, sha256, "heat_sn160 GPU inventory")
        for path, (size, sha256) in heat_gpu_paths.items()
    )

    if len(manifest["inputs"]) != 38:
        raise FrontierFailure("source", "manifest", "input inventory is incomplete")
    _validate_inventory(manifest["inputs"], "input inventory")
    artifacts = _validate_inventory(manifest["artifact_inventory"], "artifact inventory")
    for label in ("compiler", "preprocessor", "rocm", "gpu"):
        required_artifacts.add(environment[label]["artifact"])
    if not required_artifacts.issubset(artifacts):
        raise FrontierFailure("test", "manifest", "required artifacts are not inventoried")
    for path, claimed_size, claimed_sha256, context in artifact_claims:
        actual_size, actual_sha256 = artifacts[path]
        if (
            (claimed_size is not None and claimed_size != actual_size)
            or claimed_sha256 != actual_sha256
        ):
            raise FrontierFailure(
                "test", "manifest", f"{context} disagrees with artifact inventory"
            )


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="build and run inside the Frontier allocation")
    run.add_argument("--source-root", type=Path, required=True)
    run.add_argument("--artifact-root", type=Path, required=True)
    run.add_argument("--source-sha", required=True)
    run.add_argument("--archive-sha256", required=True)
    run.add_argument("--build-jobs", type=int, default=8)
    run.add_argument("--time-limit", required=True)
    validate = subparsers.add_parser("validate", help="validate a retained manifest")
    validate.add_argument("manifest", type=Path)
    validate.add_argument("--allow-failure", action="store_true")
    failure = subparsers.add_parser("failure", help="record a pre-run allocation failure")
    failure.add_argument("--artifact-root", type=Path, required=True)
    failure.add_argument("--source-sha", required=True)
    failure.add_argument("--archive-sha256", required=True)
    failure.add_argument("--time-limit", required=True)
    failure.add_argument("--category", required=True, choices=FAILURE_CATEGORIES)
    failure.add_argument("--phase", required=True)
    failure.add_argument("--message", required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        if arguments.command == "run":
            report_path = run_qualification(arguments)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            validate_manifest(report, require_pass=False)
            if report["status"] != "passed":
                failure = report["failure"]
                print(
                    f"Frontier qualification failed [{failure['category']}/"
                    f"{failure['phase']}]: {failure['message']}; manifest: {report_path}",
                    file=sys.stderr,
                )
                return 1
            validate_manifest(report)
            print(f"Frontier qualification passed; manifest: {report_path}")
            return 0
        if arguments.command == "failure":
            report_path = write_failure_manifest(arguments)
            report = json.loads(report_path.read_text(encoding="utf-8"))
            validate_manifest(report, require_pass=False)
            print(f"Frontier qualification failure manifest: {report_path}")
            return 0
        report = json.loads(arguments.manifest.read_text(encoding="utf-8"))
        validate_manifest(report, require_pass=not arguments.allow_failure)
    except (FrontierFailure, OSError, UnicodeError, json.JSONDecodeError) as error:
        print(f"Frontier qualification validation failed: {error}", file=sys.stderr)
        return 1
    print("Frontier qualification manifest is valid")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
