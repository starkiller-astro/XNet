"""Focused tests for the parallel-zone qualification runner."""

from dataclasses import replace
from pathlib import Path
import sys

import pytest

from parallel_zones import (
    ALPHA_SPECIES,
    EXPECTED_ZONES,
    QualificationFailure,
    compare_endpoint_states,
    normalize_worker_states,
    run_process,
    validate_ascii_association,
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


def _write_ascii_state(directory: Path, state: FinalState, filename_zone: int) -> None:
    values = " ".join(f"{value:.8E}" for value in state.mass_fractions.values())
    row = (
        f"{state.counters.ts} {state.time:.8E} {state.temperature_gk:.3E} "
        f"{state.density:.3E} 0.0 0.0 1.0E-6 {values} 1 1\n"
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
        validate_output_inventory(tmp_path, 1)


def test_ascii_filename_association_rejects_swapped_zone_content(tmp_path: Path) -> None:
    states = (_state(1), _state(2))
    _write_ascii_state(tmp_path, states[1], filename_zone=1)
    _write_ascii_state(tmp_path, states[0], filename_zone=2)
    with pytest.raises(QualificationFailure, match="association mismatch for zone 1"):
        validate_ascii_association(tmp_path, states)


def test_expected_failure_requires_nonzero_status(tmp_path: Path) -> None:
    result = run_process(
        (sys.executable, "-c", "raise SystemExit(3)"),
        tmp_path,
        timeout_seconds=2.0,
        expect_success=False,
    )
    assert result.return_code == 3


def test_timeout_terminates_the_process_group(tmp_path: Path) -> None:
    with pytest.raises(QualificationFailure, match="timed out"):
        run_process(
            (sys.executable, "-c", "import time; time.sleep(10)"),
            tmp_path,
            timeout_seconds=0.1,
        )
    assert (tmp_path / "xnet.status.txt").read_text(encoding="utf-8") == "timeout\n"
