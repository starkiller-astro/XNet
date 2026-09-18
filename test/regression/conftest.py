"""pytest configuration for external XNet regression cases."""

from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    group = parser.getgroup("xnet regression")
    group.addoption(
        "--xnet-executable",
        metavar="PATH",
        help="explicit path to the XNet executable under test (required by regression cases)",
    )
    group.addoption(
        "--xnet-timeout",
        type=float,
        default=30.0,
        metavar="SECONDS",
        help="per-case XNet process timeout (default: 30)",
    )
    group.addoption(
        "--xnse-executable",
        metavar="PATH",
        help="explicit path to the standalone xnse executable under test",
    )


@pytest.fixture
def xnet_executable(pytestconfig: pytest.Config) -> Path:
    value = pytestconfig.getoption("xnet_executable")
    if value is None:
        pytest.fail(
            "setup failure: --xnet-executable=PATH is required; "
            "no historical executable name is selected implicitly",
            pytrace=False,
        )
    return Path(value)


@pytest.fixture
def xnet_timeout(pytestconfig: pytest.Config) -> float:
    return float(pytestconfig.getoption("xnet_timeout"))


@pytest.fixture
def xnse_executable(pytestconfig: pytest.Config) -> Path:
    value = pytestconfig.getoption("xnse_executable")
    if value is None:
        pytest.fail(
            "setup failure: --xnse-executable=PATH is required by xnse process tests",
            pytrace=False,
        )
    return Path(value)
