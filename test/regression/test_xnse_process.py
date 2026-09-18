"""Bounded process contracts for the standalone NSE executable."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CASE_DIRECTORY = REPOSITORY_ROOT / "test" / "regression" / "cases" / "xnse_sn160"
NETWORK_INPUTS = ("sunet", "netsu", "netweak", "netwinv")
ELEMENT_Z = {
    symbol: charge
    for charge, symbol in enumerate(
        (
            "",
            "h",
            "he",
            "li",
            "be",
            "b",
            "c",
            "n",
            "o",
            "f",
            "ne",
            "na",
            "mg",
            "al",
            "si",
            "p",
            "s",
            "cl",
            "ar",
            "k",
            "ca",
            "sc",
            "ti",
            "v",
            "cr",
            "mn",
            "fe",
            "co",
            "ni",
            "cu",
            "zn",
            "ga",
            "ge",
        )
    )
}
STATE_ROW = re.compile(
    r"^\s*NSE solved\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+",
    re.MULTILINE,
)
COUNTER_ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$", re.MULTILINE)
TIMER_HEADING = re.compile(r"^Timers:", re.MULTILINE)
TIMER_ROW = re.compile(
    r"^\s*([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s+"
    r"([+-]?\d+\.\d+E[+-]\d+)\s*$",
    re.MULTILINE,
)
ABUNDANCE_ROW = re.compile(
    r"^\s*([a-z][a-z0-9]*)\s+([+-]?\d+\.\d+E[+-]\d+)\s*$", re.MULTILINE
)


def _mass_and_charge(species: str) -> tuple[int, int]:
    if species == "n":
        return 1, 0
    if species == "p":
        return 1, 1
    if species == "d":
        return 2, 1
    match = re.fullmatch(r"([a-z]+)(\d+)", species)
    assert match is not None
    return int(match.group(2)), ELEMENT_Z[match.group(1)]


def _prepare_work_directory(tmp_path: Path) -> Path:
    work_directory = tmp_path / "xnse"
    work_directory.mkdir()
    shutil.copy2(CASE_DIRECTORY / "control", work_directory / "control")
    network_directory = work_directory / "Data_SN160"
    network_directory.mkdir()
    source_directory = REPOSITORY_ROOT / "test" / "Data_SN160"
    for filename in NETWORK_INPUTS:
        shutil.copy2(source_directory / filename, network_directory / filename)
    shutil.copy2(
        REPOSITORY_ROOT / "tools" / "starkiller-helmholtz" / "helm_table.dat",
        work_directory / "helm_table.dat",
    )
    return work_directory


def _run_xnse(
    executable: Path, work_directory: Path, input_text: str, timeout: float = 15.0
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            [str(executable.resolve())],
            cwd=work_directory,
            input=input_text,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(f"xnse exceeded the {timeout:g}-second process timeout", pytrace=False)


def test_xnse_multirow_output_association(
    xnse_executable: Path, tmp_path: Path
) -> None:
    work_directory = _prepare_work_directory(tmp_path)
    input_text = (CASE_DIRECTORY / "input").read_text(encoding="utf-8")
    result = _run_xnse(xnse_executable, work_directory, input_text)

    assert result.returncode == 0, result.stdout + result.stderr
    diagnostic = (work_directory / "nse_diag01").read_text(encoding="utf-8")
    expected_states = [
        tuple(float(value) for value in line.split()[:3])
        for line in input_text.splitlines()
    ]
    state_matches = tuple(STATE_ROW.finditer(diagnostic))
    actual_states = [
        tuple(float(value) for value in match.groups()) for match in state_matches
    ]
    assert actual_states == expected_states
    expected_species = tuple(
        line.strip()
        for line in (REPOSITORY_ROOT / "test" / "Data_SN160" / "sunet")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    )
    counters = []
    timers = []
    for index, match in enumerate(state_matches):
        state_end = (
            state_matches[index + 1].start()
            if index + 1 < len(state_matches)
            else len(diagnostic)
        )
        counter_heading = diagnostic.index("NSE Counters:", match.end(), state_end)
        abundances = tuple(
            (name, float(value))
            for name, value in ABUNDANCE_ROW.findall(
                diagnostic, match.end(), counter_heading
            )
        )
        assert tuple(name for name, _ in abundances) == expected_species
        assert all(math.isfinite(value) and value >= 0.0 for _, value in abundances)
        assert math.isclose(
            math.fsum(value for _, value in abundances),
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-6,
        )
        reconstructed_ye = math.fsum(
            charge / mass * value
            for name, value in abundances
            for mass, charge in (_mass_and_charge(name),)
        )
        assert math.isclose(
            reconstructed_ye,
            expected_states[index][2],
            rel_tol=0.0,
            abs_tol=1.0e-6,
        )
        counter_match = COUNTER_ROW.search(diagnostic, counter_heading, state_end)
        assert counter_match is not None
        counters.append(tuple(int(value) for value in counter_match.groups()))
        timer_headings = tuple(
            TIMER_HEADING.finditer(diagnostic, counter_match.end(), state_end)
        )
        assert len(timer_headings) == 1
        timer_match = TIMER_ROW.search(
            diagnostic, timer_headings[0].end(), state_end
        )
        assert timer_match is not None
        timers.append(tuple(float(value) for value in timer_match.groups()))
    assert [row[0] for row in counters] == [1, 2, 3]
    assert all(all(value >= 0 for value in row[1:]) and row[3] > 0 for row in counters)
    assert all(
        all(math.isfinite(value) and value >= 0.0 for value in row) and row[0] > 0.0
        for row in timers
    )
    assert diagnostic.count("NSE Counters:") == 3


@pytest.mark.parametrize(
    ("input_text", "expected_message"),
    [
        ("1.0E+07 9.0\n", "Not enough inputs"),
        ("1.0E+07 9.0 1.01\n", "NSE ERROR: NR Failed"),
    ],
)
def test_xnse_rejects_bad_state_rows(
    xnse_executable: Path,
    tmp_path: Path,
    input_text: str,
    expected_message: str,
) -> None:
    work_directory = _prepare_work_directory(tmp_path)
    result = _run_xnse(xnse_executable, work_directory, input_text)

    assert result.returncode != 0
    assert expected_message in result.stdout + result.stderr
    diagnostic_path = work_directory / "nse_diag01"
    if diagnostic_path.exists():
        assert not STATE_ROW.search(diagnostic_path.read_text(encoding="utf-8"))
