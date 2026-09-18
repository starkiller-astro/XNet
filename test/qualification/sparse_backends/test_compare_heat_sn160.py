"""Focused checks for the optional sparse-backend qualification runner."""

from pathlib import Path
import subprocess
import sys


SCRIPT = Path(__file__).with_name("compare_heat_sn160.py")


def test_identical_executable_content_is_rejected(tmp_path: Path) -> None:
    dense = tmp_path / "xnetd"
    sparse = tmp_path / "xnetm"
    content = b"#!/bin/sh\nexit 0\n"
    dense.write_bytes(content)
    sparse.write_bytes(content)

    result = subprocess.run(
        (
            sys.executable,
            str(SCRIPT),
            "--provider=ma48",
            f"--dense-executable={dense}",
            f"--sparse-executable={sparse}",
            f"--work-directory={tmp_path / 'work'}",
        ),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert "identical SHA-256 hashes" in result.stderr
    assert not (tmp_path / "work").exists()
