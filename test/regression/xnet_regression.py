"""Execution and comparison helpers for XNet regression tests."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import date
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Mapping, Sequence


ALPHA_SPECIES = (
    "he4",
    "c12",
    "o16",
    "ne20",
    "mg24",
    "si28",
    "s32",
    "ar36",
    "ca40",
    "ti44",
    "cr48",
    "fe52",
    "ni56",
    "zn60",
)

TORCH47_SPECIES = (
    "n",
    "p",
    "d",
    "t",
    "he3",
    "he4",
    "li7",
    "be7",
    "b8",
    "c12",
    "c13",
    "n13",
    "n14",
    "n15",
    "o14",
    "o15",
    "o16",
    "o17",
    "o18",
    "f17",
    "f18",
    "f19",
    "ne18",
    "ne19",
    "ne20",
    "na23",
    "mg23",
    "mg24",
    "al27",
    "si27",
    "si28",
    "p30",
    "p31",
    "s31",
    "s32",
    "cl35",
    "ar36",
    "k39",
    "ca40",
    "sc43",
    "ti44",
    "v47",
    "cr48",
    "mn51",
    "fe52",
    "co55",
    "ni56",
)

SN160_SPECIES = (
    "n",
    "p",
    "d",
    "he3",
    "he4",
    "li6",
    "li7",
    "be7",
    "be9",
    "b8",
    "b10",
    "b11",
    "c12",
    "c13",
    "c14",
    "n13",
    "n14",
    "n15",
    "o14",
    "o15",
    "o16",
    "o17",
    "o18",
    "f17",
    "f18",
    "f19",
    "ne18",
    "ne19",
    "ne20",
    "ne21",
    "ne22",
    "na21",
    "na22",
    "na23",
    "mg23",
    "mg24",
    "mg25",
    "mg26",
    "al25",
    "al26",
    "al27",
    "si28",
    "si29",
    "si30",
    "si31",
    "si32",
    "p29",
    "p30",
    "p31",
    "p32",
    "p33",
    "s32",
    "s33",
    "s34",
    "s35",
    "s36",
    "cl33",
    "cl34",
    "cl35",
    "cl36",
    "cl37",
    "ar36",
    "ar37",
    "ar38",
    "ar39",
    "ar40",
    "k37",
    "k38",
    "k39",
    "k40",
    "k41",
    "ca40",
    "ca41",
    "ca42",
    "ca43",
    "ca44",
    "ca45",
    "ca46",
    "ca47",
    "ca48",
    "sc43",
    "sc44",
    "sc45",
    "sc46",
    "sc47",
    "sc48",
    "sc49",
    "ti44",
    "ti45",
    "ti46",
    "ti47",
    "ti48",
    "ti49",
    "ti50",
    "ti51",
    "v46",
    "v47",
    "v48",
    "v49",
    "v50",
    "v51",
    "v52",
    "cr48",
    "cr49",
    "cr50",
    "cr51",
    "cr52",
    "cr53",
    "cr54",
    "mn50",
    "mn51",
    "mn52",
    "mn53",
    "mn54",
    "mn55",
    "fe52",
    "fe53",
    "fe54",
    "fe55",
    "fe56",
    "fe57",
    "fe58",
    "co53",
    "co54",
    "co55",
    "co56",
    "co57",
    "co58",
    "co59",
    "ni56",
    "ni57",
    "ni58",
    "ni59",
    "ni60",
    "ni61",
    "ni62",
    "ni63",
    "ni64",
    "cu57",
    "cu58",
    "cu59",
    "cu60",
    "cu61",
    "cu62",
    "cu63",
    "cu64",
    "cu65",
    "zn59",
    "zn60",
    "zn61",
    "zn62",
    "zn63",
    "zn64",
    "zn65",
    "zn66",
    "ga62",
    "ga63",
    "ga64",
    "ge63",
    "ge64",
)

# These established products remain comparison anchors for continuity even
# when one falls below the general material endpoint threshold.
SILICON_BURNING_COMPARISON_SPECIES = (
    "si28",
    "s32",
    "ar36",
    "ca40",
    "ti44",
    "cr48",
    "fe52",
    "ni56",
)
MATERIAL_MASS_FRACTION_THRESHOLD = 1.0e-4
COMPARISON_SCHEMA = "xnet-comparison-v1"


def comparison_species_for_zone(
    expected_species: Sequence[str],
    mass_fractions: Mapping[str, float],
) -> tuple[str, ...]:
    """Select retained anchors and material abundances in complete-vector order."""

    return tuple(
        species
        for species in expected_species
        if species in SILICON_BURNING_COMPARISON_SPECIES
        or mass_fractions[species] >= MATERIAL_MASS_FRACTION_THRESHOLD
    )


FLOAT_TOKEN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
END_RECORD = re.compile(
    rf"^End\s+(\d+)\s+(\d+)\s+({FLOAT_TOKEN})\s+({FLOAT_TOKEN})\s+"
    rf"({FLOAT_TOKEN})\s+({FLOAT_TOKEN})\s+({FLOAT_TOKEN})\s*$"
)
ABUNDANCE_PAIR = re.compile(rf"([A-Za-z][A-Za-z0-9]*)\s+({FLOAT_TOKEN})")
TIMER_HEADING = re.compile(r"^Timers Summary:\s*$")
TIMER_ROW = re.compile(rf"^\s+[A-Za-z][A-Za-z0-9_/-]*\s+{FLOAT_TOKEN}\s*$")
COUNTER_HEADING = re.compile(
    r"^Counters:\s+Zone\s+TS\s+NR\s+Jacobian\s+Deriv\s+CrossSect\s*$"
)
COUNTER_ROW = re.compile(r"^\s*\d+(?:\s+\d+){5}\s*$")
NONFINITE_TOKEN = re.compile(
    r"(?i)(?<![A-Za-z0-9_])[+-]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])"
)


class RegressionFailure(RuntimeError):
    """Base class for a classified regression failure."""


class SetupFailure(RegressionFailure):
    """The executable, case definition, inputs, or work directory is invalid."""


class ExecutionFailure(RegressionFailure):
    """XNet did not complete and create all required output."""


class ParsingFailure(RegressionFailure):
    """A required diagnostic record is absent or malformed."""


class ComparisonFailure(RegressionFailure):
    """Parsed output disagrees with the characterization reference."""


@dataclass(frozen=True)
class StagedInput:
    """One tracked file linked into a relative path in an isolated run."""

    source: Path
    destination: Path


@dataclass(frozen=True)
class RegressionCase:
    name: str
    control: Path
    network_data: Path
    trajectories: tuple[Path, ...]
    helm_table: Path
    reference: Path
    expected_zones: tuple[int, ...]
    expected_species: tuple[str, ...]
    network_inputs: tuple[str, ...]
    reference_schema: str | None = None
    comparison_schema: str | None = None
    expected_legacy_id: int | None = None
    expected_solver: str | None = None
    expected_effective_controls: Mapping[str, float] | None = None
    expected_characterization_metadata: Mapping[str, object] | None = None
    staged_inputs: tuple[StagedInput, ...] = ()
    diagnostic_groups: tuple[tuple[int, ...], ...] | None = None
    equivalent_zone_groups: tuple[tuple[int, ...], ...] = ()
    required_diagnostic_markers: tuple[str, ...] = ()

    @property
    def required_outputs(self) -> tuple[str, ...]:
        zone_width = len(str(max(self.expected_zones)))
        zone_outputs = tuple(
            filename
            for zone in self.expected_zones
            for filename in (
                f"ev_{self.name}_{zone:0{zone_width}d}",
                f"ts_{self.name}_{zone:0{zone_width}d}",
            )
        )
        return ("net_diag01", *zone_outputs)

    @property
    def expected_diagnostic_groups(self) -> tuple[tuple[int, ...], ...]:
        """Return the declared End/counter/timer grouping for this case."""

        if self.diagnostic_groups is not None:
            return self.diagnostic_groups
        return tuple((zone,) for zone in self.expected_zones)


@dataclass(frozen=True)
class ProcessResult:
    executable: Path
    work_directory: Path
    return_code: int
    stdout: str
    stderr: str


@dataclass(frozen=True)
class FinalState:
    zone: int
    step: int
    target_time: float
    time: float
    temperature_gk: float
    density: float
    electron_fraction: float
    mass_fractions: Mapping[str, float]
    counters: SolverCounters


@dataclass(frozen=True)
class SolverCounters:
    """Source-labeled values from one XNet Counters record."""

    ts: int
    nr: int
    jacobian: int
    derivative: int
    cross_section: int


@dataclass(frozen=True)
class Tolerance:
    value: float
    atol: float
    rtol: float
    exact: bool = False


@dataclass(frozen=True)
class ToleranceBounds:
    """Absolute and relative bounds without a duplicated expected value."""

    atol: float
    rtol: float
    exact: bool = False


@dataclass(frozen=True)
class CompositionNormLimits:
    """Optional pass/fail limits for a complete composition vector."""

    l1: float | None = None
    linf: float | None = None


@dataclass(frozen=True)
class CharacterizationReference:
    """Complete expected final state plus selected pass/fail tolerances."""

    case_name: str
    expected_zones: tuple[int, ...]
    final_steps: Mapping[int, int]
    fields: Mapping[int, Mapping[str, Tolerance]]
    # The complete vector supplies expected values and norm diagnostics.
    mass_fractions: Mapping[int, Mapping[str, float]]
    # This subset selects species for field-aware pass/fail comparison.
    mass_fraction_tolerances: Mapping[int, Mapping[str, ToleranceBounds]]
    mass_fraction_sum_atols: Mapping[int, float]
    # The sum of the emitted complete vector is compared to its canonical sum.
    mass_fraction_printed_sum_tolerances: Mapping[int, Tolerance] | None = None
    composition_norm_limits: Mapping[int, CompositionNormLimits] | None = None
    solver_counters: Mapping[int, SolverCounters] | None = None
    reference_schema: str | None = None
    comparison_schema: str | None = None
    metadata_inputs: Mapping[str, object] | None = None
    input_sha256: Mapping[str, str] | None = None
    legacy_provenance: Mapping[str, object] | None = None
    characterization_metadata: Mapping[str, object] | None = None


@dataclass(frozen=True)
class CompositionNorms:
    zone: int
    l1: float
    l2: float
    linf: float
    linf_species: str | None


def tnsn_alpha_case(repository_root: Path) -> RegressionCase:
    case_directory = repository_root / "test" / "regression" / "cases" / "tnsn_alpha"
    return RegressionCase(
        name="tnsn_alpha",
        control=case_directory / "control",
        network_data=repository_root / "test" / "Data_alpha",
        trajectories=(
            repository_root / "test" / "Test_Problems" / "th_sn1aflame",
        ),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=tuple(range(1, 11)),
        expected_species=ALPHA_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv", "ab_co"),
        comparison_schema=COMPARISON_SCHEMA,
    )


def heat_alpha_case(repository_root: Path) -> RegressionCase:
    case_directory = (
        repository_root / "test" / "regression" / "cases" / "heat_alpha"
    )
    return RegressionCase(
        name="heat_alpha",
        control=case_directory / "control",
        network_data=repository_root / "test" / "Data_alpha",
        trajectories=tuple(
            repository_root / "test" / "Test_Problems" / f"th_co_burn_{zone}"
            for zone in range(1, 7)
        ),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=tuple(range(1, 7)),
        expected_species=ALPHA_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv", "ab_co"),
        comparison_schema=COMPARISON_SCHEMA,
    )


def tnsn_torch47_case(repository_root: Path) -> RegressionCase:
    case_directory = (
        repository_root / "test" / "regression" / "cases" / "tnsn_torch47"
    )
    return RegressionCase(
        name="tnsn_torch47",
        control=case_directory / "control",
        network_data=repository_root / "test" / "Data_torch47",
        trajectories=(
            repository_root / "test" / "Test_Problems" / "th_sn1aflame",
        ),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=(1,),
        expected_species=TORCH47_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv", "ab_co"),
        comparison_schema=COMPARISON_SCHEMA,
    )


def heat_sn160_case(repository_root: Path) -> RegressionCase:
    case_directory = (
        repository_root / "test" / "regression" / "cases" / "heat_sn160"
    )
    return RegressionCase(
        name="heat_sn160",
        control=case_directory / "control",
        network_data=repository_root / "test" / "Data_SN160",
        trajectories=tuple(
            repository_root / "test" / "Test_Problems" / f"th_co_burn_{zone}"
            for zone in range(1, 7)
        ),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=tuple(range(1, 7)),
        expected_species=SN160_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv", "ab_co"),
        comparison_schema=COMPARISON_SCHEMA,
    )


def batch_alpha_case(repository_root: Path) -> RegressionCase:
    """Legacy ID 61: four serial batches of four alpha-network zones."""

    case_directory = (
        repository_root / "test" / "regression" / "cases" / "batch_alpha"
    )
    network_data = repository_root / "test" / "Data_alpha"
    return RegressionCase(
        name="batch_alpha",
        control=case_directory / "control",
        network_data=network_data,
        trajectories=(),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=tuple(range(1, 17)),
        expected_species=ALPHA_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv"),
        reference_schema="xnet-characterization-v1",
        comparison_schema=COMPARISON_SCHEMA,
        expected_legacy_id=61,
        expected_solver="Backward Euler (isolv = 1)",
        expected_effective_controls={},
        expected_characterization_metadata={
            "baseline_kind": "characterization",
            "baseline_status": (
                "characterization-only; not independently validated scientific truth"
            ),
            "generated_from_revision": "47befd90fff3c05828bf2f1b528cd74d7aa6a36c",
            "generated_on": "2026-08-07",
            "platform": "macOS 26.6 arm64",
            "compiler": "GNU Fortran (Homebrew GCC 16.1.0) 16.1.0",
            "python": "Python 3.13.0",
            "pytest": "pytest 9.1.1",
            "build": {
                "executable": "source/xnet",
                "CMODE": "OPT",
                "PE_ENV": "GNU",
                "MPI_MODE": "OFF",
                "OPENMP_MODE": "OFF",
                "GPU_MODE": "OFF",
                "EOS": "STARKILLER",
                "MATRIX_SOLVER": "dense",
                "LAPACK_VER": "NETLIB",
            },
            "legacy_provenance": {
                "legacy_id": 61,
                "assembly": [
                    "test/test_settings_batch",
                    "test/Test_Problems/setup_batch_alpha",
                ],
                "maintained_solver": "Backward Euler (isolv = 1)",
                "normalized_changes": [
                    "remove Test_Results/ from ASCII and binary output roots",
                    "preserve nested abundance and trajectory prefixes in isolated staging",
                ],
                "effective_runtime_controls": {},
            },
        },
        staged_inputs=tuple(
            StagedInput(
                network_data / "ab_batch" / f"ab_batch_{zone:02d}",
                Path("Data_alpha") / "ab_batch" / f"ab_batch_{zone:02d}",
            )
            for zone in range(1, 17)
        )
        + tuple(
            StagedInput(
                repository_root
                / "test"
                / "Test_Problems"
                / "th_batch"
                / f"th_batch_{zone:02d}",
                Path("Test_Problems") / "th_batch" / f"th_batch_{zone:02d}",
            )
            for zone in range(1, 17)
        ),
        diagnostic_groups=((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)),
        equivalent_zone_groups=((1, 5, 9, 13), (2, 6, 10, 14), (3, 7, 11, 15), (4, 8, 12, 16)),
    )


def bdf_sn160_case(repository_root: Path) -> RegressionCase:
    case_directory = (
        repository_root / "test" / "regression" / "cases" / "bdf_sn160"
    )

    return RegressionCase(
        name="bdf_sn160",
        control=case_directory / "control",
        network_data=repository_root / "test" / "Data_SN160",
        trajectories=tuple(
            repository_root / "test" / "Test_Problems" / f"th_co_burn_{zone}"
            for zone in range(1, 7)
        ),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=tuple(range(1, 7)),
        expected_species=SN160_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv", "ab_co"),
        reference_schema="xnet-characterization-v1",
        comparison_schema=COMPARISON_SCHEMA,
        expected_legacy_id=54,
        expected_solver="Backward Differentiation Formula (isolv = 3)",
        expected_effective_controls={
            "maximum_abundance_change": 1.0e10,
            "maximum_temperature_change": 1.0e10,
        },
        expected_characterization_metadata={
            "baseline_kind": "characterization",
            "baseline_status": (
                "characterization-only; not independently validated scientific truth"
            ),
            "generated_from_revision": (
                "a8b64764a6d614f406da6c897e6b051fb3e1972d"
            ),
            "generated_on": "2026-08-05",
            "platform": "macOS 26.6 arm64",
            "compiler": "GNU Fortran (Homebrew GCC 16.1.0) 16.1.0",
            "python": "Python 3.13.0",
            "pytest": "pytest 9.1.1",
            "build": {
                "executable": "source/xnet",
                "CMODE": "OPT",
                "PE_ENV": "GNU",
                "MPI_MODE": "OFF",
                "OPENMP_MODE": "OFF",
                "GPU_MODE": "OFF",
                "EOS": "STARKILLER",
                "MATRIX_SOLVER": "dense",
                "LAPACK_VER": "NETLIB",
            },
            "legacy_provenance": {
                "legacy_id": 54,
                "assembly": [
                    "test/test_settings_bdf",
                    "test/Test_Problems/setup_bdf_sn160",
                ],
                "maintained_solver": (
                    "Backward Differentiation Formula (isolv = 3)"
                ),
                "normalized_changes": [
                    "remove Test_Results/ from ASCII and binary output roots",
                    "remove Test_Problems/ from six trajectory paths",
                    "change zone block size from legacy 4 to regression 1",
                ],
                "effective_runtime_controls": {
                    "maximum_abundance_change": 1.0e10,
                    "maximum_temperature_change": 1.0e10,
                },
            },
        },
    )


def nse_sn160_case(repository_root: Path) -> RegressionCase:
    """One ordinary trajectory initialized through the production NSE branch."""

    case_directory = (
        repository_root / "test" / "regression" / "cases" / "nse_sn160"
    )
    network_data = repository_root / "test" / "Data_SN160"
    return RegressionCase(
        name="nse_sn160",
        control=case_directory / "control",
        network_data=network_data,
        trajectories=(),
        helm_table=(
            repository_root / "tools" / "starkiller-helmholtz" / "helm_table.dat"
        ),
        reference=case_directory / "reference" / "final_state.json",
        expected_zones=(1,),
        expected_species=SN160_SPECIES,
        network_inputs=("sunet", "netsu", "netweak", "netwinv"),
        reference_schema="xnet-characterization-v1",
        comparison_schema=COMPARISON_SCHEMA,
        expected_legacy_id=41,
        expected_solver="Backward Euler (isolv = 1)",
        expected_effective_controls={},
        expected_characterization_metadata={
            "baseline_kind": "characterization",
            "baseline_status": (
                "characterization-only; not independently validated scientific truth"
            ),
            "generated_from_revision": "e5c01dded7ccae53f04122599eb0b2face773c2c",
            "generated_on": "2026-08-07",
            "platform": "macOS 26.6 arm64",
            "compiler": "GNU Fortran (Homebrew GCC 16.1.0) 16.1.0",
            "python": "Python 3.13.0",
            "pytest": "pytest 9.1.1",
            "build": {
                "executable": "source/xnet",
                "CMODE": "OPT",
                "PE_ENV": "GNU",
                "MPI_MODE": "OFF",
                "OPENMP_MODE": "OFF",
                "GPU_MODE": "OFF",
                "EOS": "STARKILLER",
                "MATRIX_SOLVER": "dense",
                "LAPACK_VER": "NETLIB",
            },
            "legacy_provenance": {
                "legacy_id": 41,
                "assembly": [
                    "test/test_settings_nse",
                    "test/Test_Problems/setup_nse_ccsn_nup",
                    "issue #40 SN160 characterization definition",
                ],
                "maintained_solver": "Backward Euler (isolv = 1)",
                "normalized_changes": [
                    "replace Data_nup with the maintained SN160 network",
                    "use ab_ye49 so the file state differs from trajectory Ye=0.55",
                    "remove Test_Results/ and Test_Problems/ runtime path prefixes",
                ],
                "effective_runtime_controls": {},
            },
        },
        staged_inputs=(
            StagedInput(
                network_data / "ab_ye49", Path("Data_SN160") / "ab_ye49"
            ),
            StagedInput(
                repository_root
                / "test"
                / "Test_Problems"
                / "th_frohlich2006_nse",
                Path("th_frohlich2006_nse"),
            ),
        ),
        required_diagnostic_markers=("Initial abundances from NSE",),
    )


def _read_sunet_species(path: Path) -> tuple[str, ...]:
    try:
        species = tuple(
            line.strip().lower()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    except (OSError, UnicodeError) as error:
        raise SetupFailure(
            f"could not read network species input {path}: {error}"
        ) from error
    if not species or len(set(species)) != len(species):
        raise SetupFailure(
            f"network species input is empty or contains duplicates: {path}"
        )
    return species


def validate_executable(executable: Path) -> Path:
    executable = executable.expanduser().resolve()
    if not executable.exists():
        raise SetupFailure(f"XNet executable does not exist: {executable}")
    if not executable.is_file():
        raise SetupFailure(f"XNet executable is not a regular file: {executable}")
    if not os.access(executable, os.X_OK):
        raise SetupFailure(f"XNet executable is not executable: {executable}")
    return executable


def _validate_staged_destination(destination: Path, context: str) -> Path:
    """Require a non-empty, traversal-free path below the work directory."""

    if (
        destination.is_absolute()
        or not destination.parts
        or any(part in ("", ".", "..") for part in destination.parts)
    ):
        raise SetupFailure(
            f"invalid staged destination for {context}: {destination}"
        )
    return destination


def _case_staged_inputs(case: RegressionCase) -> tuple[StagedInput, ...]:
    """Return all default and explicitly declared isolated-workdir links."""

    inputs = [
        StagedInput(case.network_data / filename, Path(case.network_data.name) / filename)
        for filename in case.network_inputs
    ]
    inputs.extend(
        StagedInput(trajectory, Path(trajectory.name)) for trajectory in case.trajectories
    )
    inputs.extend(case.staged_inputs)
    inputs.append(StagedInput(case.helm_table, Path(case.helm_table.name)))
    return tuple(inputs)


def _validate_diagnostic_groups(case: RegressionCase) -> None:
    groups = case.expected_diagnostic_groups
    if not groups or any(not group for group in groups):
        raise SetupFailure(f"invalid diagnostic group definition for {case.name}")
    flattened = tuple(zone for group in groups for zone in group)
    if flattened != case.expected_zones or len(set(flattened)) != len(flattened):
        raise SetupFailure(
            f"diagnostic groups must contain each expected zone once in order for {case.name}"
        )
    for group in case.equivalent_zone_groups:
        if len(group) < 2 or any(zone not in case.expected_zones for zone in group):
            raise SetupFailure(f"invalid equivalent-zone group definition for {case.name}")


def _reserved_work_directory_paths(case: RegressionCase) -> set[Path]:
    """Paths owned by the runner or XNet rather than case-declared staging."""

    return {
        Path("control"),
        *map(Path, case.required_outputs),
        Path("xnet.stdout.txt"),
        Path("xnet.stderr.txt"),
        Path("xnet.status.txt"),
        Path("composition_error_norms.json"),
    }


def _validate_case_inputs(case: RegressionCase) -> None:
    required = {
        "complete control input": case.control,
        "network data directory": case.network_data,
        "Helmholtz EOS table": case.helm_table,
        "characterization reference": case.reference,
    }
    required.update(
        {
            f"thermodynamic trajectory {index}": path
            for index, path in enumerate(case.trajectories, start=1)
        }
    )
    required.update(
        {
            f"declared staged input {index}": item.source
            for index, item in enumerate(case.staged_inputs, start=1)
        }
    )
    missing = [
        f"{description}: {path}"
        for description, path in required.items()
        if not path.exists()
    ]
    if missing:
        raise SetupFailure("required case input is missing:\n  " + "\n  ".join(missing))
    if not case.network_data.is_dir():
        raise SetupFailure(f"network data path is not a directory: {case.network_data}")
    if not case.expected_zones or len(set(case.expected_zones)) != len(case.expected_zones):
        raise SetupFailure(f"invalid expected zone definition for {case.name}")
    if not case.expected_species or len(set(case.expected_species)) != len(
        case.expected_species
    ):
        raise SetupFailure(f"invalid expected species definition for {case.name}")
    if (
        len(set(case.required_diagnostic_markers))
        != len(case.required_diagnostic_markers)
        or any(not marker.strip() for marker in case.required_diagnostic_markers)
    ):
        raise SetupFailure(
            f"invalid required diagnostic marker definition for {case.name}"
        )
    trajectory_names = tuple(path.name for path in case.trajectories)
    if case.trajectories and len(set(trajectory_names)) != len(trajectory_names):
        raise SetupFailure(
            f"invalid thermodynamic trajectory definition for {case.name}"
        )
    if (
        not case.network_inputs
        or len(set(case.network_inputs)) != len(case.network_inputs)
        or any(Path(filename).name != filename for filename in case.network_inputs)
    ):
        raise SetupFailure(f"invalid network input definition for {case.name}")
    missing_network_inputs = [
        case.network_data / filename
        for filename in case.network_inputs
        if not (case.network_data / filename).is_file()
    ]
    if missing_network_inputs:
        raise SetupFailure(
            "required network source input is missing:\n  "
            + "\n  ".join(str(path) for path in missing_network_inputs)
        )
    sunet_species = _read_sunet_species(case.network_data / "sunet")
    if sunet_species != case.expected_species:
        raise SetupFailure(
            f"network species input does not match the case definition for {case.name}: "
            f"{sunet_species} != {case.expected_species}"
        )
    _validate_diagnostic_groups(case)
    staged_destinations: set[Path] = set()
    reserved_paths = _reserved_work_directory_paths(case)
    for item in _case_staged_inputs(case):
        destination = _validate_staged_destination(item.destination, str(item.source))
        if any(
            reserved == destination
            or reserved in destination.parents
            or destination in reserved.parents
            for reserved in reserved_paths
        ):
            raise SetupFailure(f"staged destination is reserved: {destination}")
        if destination in staged_destinations:
            raise SetupFailure(f"duplicate staged destination: {destination}")
        if any(
            existing in destination.parents or destination in existing.parents
            for existing in staged_destinations
        ):
            raise SetupFailure(
                f"overlapping staged destinations: {destination}"
            )
        staged_destinations.add(destination)


def prepare_work_directory(case: RegressionCase, work_directory: Path) -> Path:
    """Create one empty isolated XNet working directory."""

    _validate_case_inputs(case)
    work_directory = work_directory.resolve()
    try:
        if work_directory.exists():
            if not work_directory.is_dir():
                raise SetupFailure(f"work path is not a directory: {work_directory}")
            if any(work_directory.iterdir()):
                raise SetupFailure(
                    f"work directory is not empty; refusing stale artifacts: {work_directory}"
                )
        else:
            work_directory.mkdir(parents=True)

        shutil.copy2(case.control, work_directory / "control")
        for item in _case_staged_inputs(case):
            destination = work_directory / item.destination
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() or destination.is_symlink():
                raise SetupFailure(f"refusing to replace staged input: {destination}")
            destination.symlink_to(item.source.resolve())
    except RegressionFailure:
        raise
    except OSError as error:
        raise SetupFailure(
            f"could not prepare work directory {work_directory}: {error}"
        ) from error
    return work_directory


def _captured_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def _write_process_artifacts(
    work_directory: Path, stdout: str, stderr: str, status: str
) -> None:
    try:
        (work_directory / "xnet.stdout.txt").write_text(stdout, encoding="utf-8")
        (work_directory / "xnet.stderr.txt").write_text(stderr, encoding="utf-8")
        (work_directory / "xnet.status.txt").write_text(status + "\n", encoding="utf-8")
    except OSError as error:
        raise ExecutionFailure(
            f"could not preserve process artifacts in {work_directory}: {error}"
        ) from error


def run_xnet(
    executable: Path,
    case: RegressionCase,
    work_directory: Path,
    *,
    timeout_seconds: float,
    environment: Mapping[str, str] | None = None,
) -> ProcessResult:
    """Run XNet once and require fresh output from the isolated directory."""

    executable = validate_executable(executable)
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise SetupFailure(f"timeout must be a positive finite value: {timeout_seconds}")

    try:
        completed = subprocess.run(
            [str(executable)],
            cwd=work_directory,
            env=None if environment is None else dict(environment),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        stdout = _captured_text(error.stdout)
        stderr = _captured_text(error.stderr)
        _write_process_artifacts(work_directory, stdout, stderr, "timeout")
        raise ExecutionFailure(
            f"XNet timed out after {timeout_seconds:g} seconds; "
            f"artifacts: {work_directory}"
        ) from error
    except OSError as error:
        _write_process_artifacts(work_directory, "", str(error), "launch-error")
        raise ExecutionFailure(
            f"could not launch XNet executable {executable}: {error}; "
            f"artifacts: {work_directory}"
        ) from error

    _write_process_artifacts(
        work_directory,
        completed.stdout,
        completed.stderr,
        f"return_code={completed.returncode}",
    )
    if completed.returncode < 0:
        raise ExecutionFailure(
            f"XNet was terminated by signal {-completed.returncode}; "
            f"artifacts: {work_directory}"
        )
    if completed.returncode != 0:
        raise ExecutionFailure(
            f"XNet returned nonzero status {completed.returncode}; "
            f"artifacts: {work_directory}"
        )

    missing = [
        filename
        for filename in case.required_outputs
        if not (work_directory / filename).is_file()
        or (work_directory / filename).stat().st_size == 0
    ]
    if missing:
        raise ExecutionFailure(
            "XNet returned zero but required fresh output is missing or empty: "
            f"{', '.join(missing)}; artifacts: {work_directory}"
        )

    return ProcessResult(
        executable=executable,
        work_directory=work_directory,
        return_code=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )


def strip_timer_sections(text: str) -> tuple[str, int]:
    """Remove only delimited timer headings and timer name/value rows."""

    lines = text.splitlines()
    retained: list[str] = []
    timer_count = 0
    index = 0
    while index < len(lines):
        if not TIMER_HEADING.match(lines[index]):
            retained.append(lines[index])
            index += 1
            continue

        timer_count += 1
        index += 1
        timer_rows = 0
        while index < len(lines) and TIMER_ROW.match(lines[index]):
            timer_rows += 1
            index += 1
        if timer_rows == 0:
            raise ParsingFailure(
                f"timer section {timer_count} contains no recognized timer rows"
            )
    return "\n".join(retained), timer_count


def _as_finite_float(token: str, context: str) -> float:
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as error:
        raise ParsingFailure(f"malformed numerical token for {context}: {token}") from error
    if not math.isfinite(value):
        raise ParsingFailure(f"non-finite numerical value for {context}: {token}")
    return value


def parse_diagnostic(
    text: str,
    expected_zones: Sequence[int],
    expected_species: Sequence[str],
    expected_groups: Sequence[Sequence[int]] | None = None,
) -> tuple[FinalState, ...]:
    """Parse declared End/composition/counter/timer groups in diagnostic order."""

    if NONFINITE_TOKEN.search(text):
        token = NONFINITE_TOKEN.search(text)
        assert token is not None
        raise ParsingFailure(f"diagnostic contains non-finite token: {token.group(0)}")

    expected_zones = tuple(expected_zones)
    expected_species = tuple(expected_species)
    if not expected_species or len(set(expected_species)) != len(expected_species):
        raise ParsingFailure("expected species definition is empty or contains duplicates")
    groups = (
        tuple(tuple(group) for group in expected_groups)
        if expected_groups is not None
        else tuple((zone,) for zone in expected_zones)
    )
    if not groups or any(not group for group in groups):
        raise ParsingFailure("expected diagnostic groups are empty or malformed")
    if tuple(zone for group in groups for zone in group) != expected_zones:
        raise ParsingFailure(
            "expected diagnostic groups do not contain ordered expected zones"
        )

    lines = text.splitlines()
    states: dict[int, FinalState] = {}
    index = 0
    for group_number, group in enumerate(groups, start=1):
        while index < len(lines) and not lines[index].startswith("End"):
            if (
                COUNTER_HEADING.match(lines[index])
                or TIMER_HEADING.match(lines[index])
                or COUNTER_ROW.match(lines[index])
                or (group_number > 1 and TIMER_ROW.match(lines[index]))
            ):
                raise ParsingFailure(
                    f"unexpected diagnostic structure before group {group_number}: "
                    f"{lines[index]}"
                )
            index += 1
        group_states: list[FinalState] = []
        for expected_zone in group:
            if index >= len(lines):
                raise ParsingFailure(
                    f"expected ordered final records for zone {expected_zone}; found none"
                )
            if COUNTER_HEADING.match(lines[index]):
                raise ParsingFailure(
                    f"missing End record for zone {expected_zone} before Counters"
                )
            match = END_RECORD.match(lines[index])
            if match is None:
                raise ParsingFailure(f"malformed End record: {lines[index]}")
            zone = int(match.group(1))
            if zone != expected_zone:
                raise ParsingFailure(
                    f"group {group_number} expected End zone {expected_zone}, found {zone}"
                )
            if zone in states:
                raise ParsingFailure(f"duplicate End record for zone {zone}")
            values = [
                _as_finite_float(token, f"zone {zone} End record")
                for token in match.groups()[2:]
            ]
            index += 1
            mass_fractions: dict[str, float] = {}
            while len(mass_fractions) < len(expected_species):
                if index >= len(lines) or lines[index].startswith(("End", "Counters:")):
                    raise ParsingFailure(
                        f"incomplete abundance record for zone {zone}: "
                        f"found {len(mass_fractions)} of {len(expected_species)} species"
                    )
                pairs = ABUNDANCE_PAIR.findall(lines[index])
                if not pairs or ABUNDANCE_PAIR.sub("", lines[index]).strip():
                    raise ParsingFailure(
                        f"malformed abundance row for zone {zone}: {lines[index]}"
                    )
                for species, token in pairs:
                    species = species.lower()
                    if species in mass_fractions:
                        raise ParsingFailure(
                            f"duplicate abundance for species {species} in zone {zone}"
                        )
                    mass_fractions[species] = _as_finite_float(
                        token, f"zone {zone} species {species}"
                    )
                index += 1
            if tuple(mass_fractions) != expected_species:
                raise ParsingFailure(
                    f"unexpected species structure for zone {zone}: {', '.join(mass_fractions)}"
                )
            group_states.append(
                FinalState(zone, int(match.group(2)), values[0], values[1], values[2], values[3], values[4], mass_fractions, SolverCounters(0, 0, 0, 0, 0))
            )
        if index >= len(lines) or not COUNTER_HEADING.match(lines[index]):
            raise ParsingFailure(f"missing Counters record after group {group_number} final states")
        index += 1
        for state in group_states:
            if index >= len(lines):
                raise ParsingFailure(f"missing counter values for zone {state.zone}")
            counters = lines[index].split()
            if len(counters) != 6 or not all(token.isdigit() for token in counters):
                raise ParsingFailure(f"malformed counter values for zone {state.zone}: {lines[index]}")
            counter_zone, *counter_values = (int(token) for token in counters)
            if counter_zone != state.zone:
                raise ParsingFailure(
                    f"zone {state.zone} counter record does not match its End record: {lines[index]}"
                )
            states[state.zone] = FinalState(
                state.zone, state.step, state.target_time, state.time,
                state.temperature_gk, state.density, state.electron_fraction,
                state.mass_fractions, SolverCounters(*counter_values)
            )
            index += 1
        if index >= len(lines) or not TIMER_HEADING.match(lines[index]):
            raise ParsingFailure(f"missing timer section for group {group_number}")
        index += 1
        timer_rows = 0
        timer_names: set[str] = set()
        while index < len(lines) and TIMER_ROW.match(lines[index]):
            timer_name = lines[index].split()[0]
            if timer_name in timer_names:
                raise ParsingFailure(
                    f"duplicate timer row {timer_name} in group {group_number}"
                )
            timer_names.add(timer_name)
            timer_rows += 1
            index += 1
        if timer_rows == 0:
            raise ParsingFailure(f"timer section {group_number} contains no recognized timer rows")

    for line in lines[index:]:
        if line.startswith("End") or COUNTER_HEADING.match(line) or TIMER_HEADING.match(line):
            raise ParsingFailure(f"unexpected extra diagnostic structure: {line}")

    actual_zones = tuple(states)
    if actual_zones != expected_zones:
        raise ParsingFailure(
            f"expected ordered final records for zones {expected_zones}, found {actual_zones}"
        )
    return tuple(states[zone] for zone in expected_zones)


def _expand_zone_items(
    item: object, expected_zones: Sequence[int], context: str
) -> Mapping[int, object]:
    if not isinstance(item, dict):
        return {zone: item for zone in expected_zones}

    values: dict[int, object] = {}
    for key, value in item.items():
        if not isinstance(key, str) or not key.isdigit():
            raise SetupFailure(
                f"reference entry {context} has a non-zone key: {key!r}"
            )
        zone = int(key)
        if str(zone) != key or zone in values:
            raise SetupFailure(
                f"reference entry {context} has an invalid zone key: {key!r}"
            )
        values[zone] = value
    expected = set(expected_zones)
    if set(values) != expected:
        raise SetupFailure(
            f"reference entry {context} must define zones {sorted(expected)}; "
            f"found {sorted(values)}"
        )
    return values


def _load_zone_floats(
    item: object, expected_zones: Sequence[int], context: str
) -> Mapping[int, float]:
    values = _expand_zone_items(item, expected_zones, context)
    if not all(
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        for value in values.values()
    ):
        raise SetupFailure(f"reference entry {context} contains a non-finite value")
    return {zone: float(value) for zone, value in values.items()}


def _load_zone_integers(
    item: object, expected_zones: Sequence[int], context: str
) -> Mapping[int, int]:
    values = _expand_zone_items(item, expected_zones, context)
    if not all(
        isinstance(value, int) and not isinstance(value, bool) and value >= 0
        for value in values.values()
    ):
        raise SetupFailure(
            f"reference entry {context} must contain nonnegative integers"
        )
    return {zone: int(value) for zone, value in values.items()}


def _load_tolerance(
    item: object,
    expected_zones: Sequence[int],
    context: str,
    *,
    require_explicit_exact: bool,
) -> Mapping[int, Tolerance]:
    required = {"value", "atol", "rtol"}
    if require_explicit_exact:
        required.add("exact")
    if not isinstance(item, dict) or set(item) != required:
        raise SetupFailure(
            f"reference entry {context} must contain exactly "
            + ", ".join(sorted(required))
        )
    exact = item.get("exact", False)
    if not isinstance(exact, bool):
        raise SetupFailure(f"reference entry {context}.exact must be a boolean")
    values = _load_zone_floats(item["value"], expected_zones, f"{context}.value")
    atols = _load_zone_floats(item["atol"], expected_zones, f"{context}.atol")
    rtols = _load_zone_floats(item["rtol"], expected_zones, f"{context}.rtol")
    if any(value < 0 for value in (*atols.values(), *rtols.values())):
        raise SetupFailure(f"reference entry {context} contains a negative tolerance")
    if exact and any(value != 0.0 for value in (*atols.values(), *rtols.values())):
        raise SetupFailure(
            f"reference entry {context} declares exact comparison with nonzero tolerance"
        )
    return {
        zone: Tolerance(values[zone], atols[zone], rtols[zone], exact)
        for zone in expected_zones
    }


def _load_tolerance_bounds(
    item: object,
    expected_zones: Sequence[int],
    context: str,
    *,
    require_explicit_exact: bool,
) -> Mapping[int, ToleranceBounds]:
    required = {"atol", "rtol"}
    if require_explicit_exact:
        required.add("exact")
    if not isinstance(item, dict) or set(item) != required:
        raise SetupFailure(
            f"reference entry {context} must contain exactly "
            + ", ".join(sorted(required))
        )
    exact = item.get("exact", False)
    if not isinstance(exact, bool):
        raise SetupFailure(f"reference entry {context}.exact must be a boolean")
    atols = _load_zone_floats(item["atol"], expected_zones, f"{context}.atol")
    rtols = _load_zone_floats(item["rtol"], expected_zones, f"{context}.rtol")
    if any(value < 0 for value in (*atols.values(), *rtols.values())):
        raise SetupFailure(f"reference entry {context} contains a negative tolerance")
    if exact and any(value != 0.0 for value in (*atols.values(), *rtols.values())):
        raise SetupFailure(
            f"reference entry {context} declares exact comparison with nonzero tolerance"
        )
    return {
        zone: ToleranceBounds(atols[zone], rtols[zone], exact)
        for zone in expected_zones
    }


def _load_composition_norm_limits(
    item: object, expected_zones: Sequence[int]
) -> Mapping[int, CompositionNormLimits]:
    if not isinstance(item, dict) or not item or set(item).difference({"l1", "linf"}):
        raise SetupFailure(
            "reference entry composition_norm_limits must contain one or both of l1 and linf"
        )
    limits = {
        name: _load_zone_floats(value, expected_zones, f"composition_norm_limits.{name}")
        for name, value in item.items()
    }
    if any(value < 0.0 for values in limits.values() for value in values.values()):
        raise SetupFailure(
            "reference entry composition_norm_limits contains a negative tolerance"
        )
    return {
        zone: CompositionNormLimits(
            l1=limits.get("l1", {}).get(zone),
            linf=limits.get("linf", {}).get(zone),
        )
        for zone in expected_zones
    }


def _load_composition_selection(
    item: object, expected_zones: Sequence[int]
) -> Mapping[int, tuple[str, ...]]:
    context = "mass_fraction_selection"
    if not isinstance(item, dict):
        raise SetupFailure(f"reference entry {context} must be a zone mapping")
    zone_items: dict[int, object] = {}
    for key, value in item.items():
        if not isinstance(key, str) or not key.isdigit():
            raise SetupFailure(
                f"reference entry {context} has a non-zone key: {key!r}"
            )
        zone = int(key)
        if str(zone) != key:
            raise SetupFailure(
                f"reference entry {context} has an invalid zone key: {key!r}"
            )
        zone_items[zone] = value
    expected = set(expected_zones)
    actual = set(zone_items)
    missing_zones = sorted(expected.difference(actual))
    unknown_zones = sorted(actual.difference(expected))
    if missing_zones or unknown_zones:
        raise SetupFailure(
            f"reference entry {context} must define exactly zones {sorted(expected)}; "
            f"missing {missing_zones}; unknown {unknown_zones}"
        )
    selections: dict[int, tuple[str, ...]] = {}
    for zone, value in zone_items.items():
        if not isinstance(value, list) or not value:
            raise SetupFailure(
                f"reference entry {context}.{zone} must be a nonempty species list"
            )
        if not all(isinstance(species, str) and species for species in value):
            raise SetupFailure(
                f"reference entry {context}.{zone} contains an invalid species"
            )
        duplicates = sorted(
            species for species in set(value) if value.count(species) > 1
        )
        if duplicates:
            raise SetupFailure(
                f"reference entry {context}.{zone} contains duplicate species: "
                + ", ".join(duplicates)
            )
        selections[zone] = tuple(value)
    return selections


def _unique_json_object(pairs: Sequence[tuple[str, object]]) -> dict[str, object]:
    document: dict[str, object] = {}
    for key, value in pairs:
        if key in document:
            raise ValueError(f"duplicate JSON object key: {key}")
        document[key] = value
    return document


def _load_reference_metadata(
    document: Mapping[str, object], reference_schema: str
) -> tuple[
    Mapping[str, object],
    Mapping[str, str],
    Mapping[str, object],
    Mapping[str, object],
]:
    if reference_schema != "xnet-characterization-v1":
        raise SetupFailure(f"unsupported reference_schema: {reference_schema!r}")

    required_strings = (
        "baseline_kind",
        "baseline_status",
        "generated_from_revision",
        "generated_on",
        "platform",
        "compiler",
        "python",
        "pytest",
    )
    for name in required_strings:
        value = document.get(name)
        if not isinstance(value, str) or not value:
            raise SetupFailure(
                f"reference metadata {name} must be a nonempty string"
            )
    if document["baseline_kind"] != "characterization":
        raise SetupFailure("reference baseline_kind must be 'characterization'")
    if "characterization-only" not in str(document["baseline_status"]):
        raise SetupFailure(
            "reference baseline_status must identify characterization-only status"
        )
    if (
        re.fullmatch(
            r"[0-9a-f]{40}", str(document["generated_from_revision"])
        )
        is None
    ):
        raise SetupFailure(
            "reference generated_from_revision must be a 40-character SHA"
        )
    try:
        date.fromisoformat(str(document["generated_on"]))
    except ValueError as error:
        raise SetupFailure(
            "reference generated_on must be an ISO calendar date"
        ) from error

    build = document.get("build")
    required_build = {
        "executable",
        "CMODE",
        "PE_ENV",
        "MPI_MODE",
        "OPENMP_MODE",
        "GPU_MODE",
        "EOS",
        "MATRIX_SOLVER",
        "LAPACK_VER",
    }
    if (
        not isinstance(build, dict)
        or set(build) != required_build
        or not all(isinstance(value, str) and value for value in build.values())
    ):
        raise SetupFailure(
            "reference build metadata must contain the complete resolved selection"
        )

    inputs = document.get("inputs")
    standard_inputs = {
        "control",
        "network_data",
        "network_source_inputs",
        "initial_abundance",
        "trajectories",
        "eos_table",
    }
    staged_inputs = {
        "control",
        "network_data",
        "network_source_inputs",
        "staged_inputs",
        "eos_table",
    }
    if not isinstance(inputs, dict) or set(inputs) not in (standard_inputs, staged_inputs):
        raise SetupFailure(
            "reference inputs metadata must contain the complete input provenance"
        )
    scalar_input_names = ["control", "network_data", "eos_table"]
    if set(inputs) == standard_inputs:
        scalar_input_names.append("initial_abundance")
    if not all(
        isinstance(inputs[name], str) and inputs[name]
        for name in scalar_input_names
    ):
        raise SetupFailure("reference inputs metadata contains an invalid path")
    for name in ("network_source_inputs",):
        values = inputs[name]
        if (
            not isinstance(values, list)
            or not values
            or not all(isinstance(value, str) and value for value in values)
        ):
            raise SetupFailure(
                f"reference inputs metadata {name} must be a nonempty string list"
            )
    if set(inputs) == standard_inputs:
        values = inputs["trajectories"]
        if (
            not isinstance(values, list)
            or not values
            or not all(isinstance(value, str) and value for value in values)
        ):
            raise SetupFailure(
                "reference inputs metadata trajectories must be a nonempty string list"
            )
    else:
        values = inputs["staged_inputs"]
        if not isinstance(values, list) or not values:
            raise SetupFailure("reference staged_inputs metadata must be a nonempty list")
        destinations: set[str] = set()
        for item in values:
            if (
                not isinstance(item, dict)
                or set(item) != {"source", "destination"}
                or not isinstance(item["source"], str)
                or not item["source"]
                or not isinstance(item["destination"], str)
            ):
                raise SetupFailure("reference staged_inputs metadata entry is invalid")
            destination = _validate_staged_destination(
                Path(item["destination"]), "reference metadata"
            ).as_posix()
            if destination in destinations:
                raise SetupFailure("reference staged_inputs metadata duplicates a destination")
            destinations.add(destination)

    legacy = document.get("legacy_provenance")
    required_legacy = {
        "legacy_id",
        "assembly",
        "maintained_solver",
        "normalized_changes",
        "effective_runtime_controls",
    }
    if not isinstance(legacy, dict) or set(legacy) != required_legacy:
        raise SetupFailure(
            "reference legacy_provenance must contain the complete control provenance"
        )
    if (
        not isinstance(legacy["legacy_id"], int)
        or isinstance(legacy["legacy_id"], bool)
        or legacy["legacy_id"] < 0
        or not isinstance(legacy["maintained_solver"], str)
        or not legacy["maintained_solver"]
    ):
        raise SetupFailure("reference legacy_provenance identity is invalid")
    for name in ("assembly", "normalized_changes"):
        values = legacy[name]
        if (
            not isinstance(values, list)
            or not values
            or not all(isinstance(value, str) and value for value in values)
        ):
            raise SetupFailure(
                f"reference legacy_provenance {name} must be a nonempty string list"
            )
    effective_controls = legacy["effective_runtime_controls"]
    if (
        not isinstance(effective_controls, dict)
        or any(not isinstance(name, str) or not name for name in effective_controls)
        or not all(
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(value)
            for value in effective_controls.values()
        )
    ):
        raise SetupFailure(
            "reference effective_runtime_controls metadata is invalid"
        )

    input_sha256 = document.get("input_sha256")
    if (
        not isinstance(input_sha256, dict)
        or not input_sha256
        or not all(
            isinstance(path, str)
            and path
            and isinstance(digest, str)
            and re.fullmatch(r"[0-9a-f]{64}", digest) is not None
            for path, digest in input_sha256.items()
        )
    ):
        raise SetupFailure("reference input_sha256 metadata is invalid")
    characterization_metadata = {
        name: document[name] for name in required_strings
    }
    characterization_metadata["build"] = build
    characterization_metadata["legacy_provenance"] = legacy
    return inputs, input_sha256, legacy, characterization_metadata


def load_reference(path: Path) -> CharacterizationReference:
    try:
        document = json.loads(
            path.read_text(encoding="utf-8"), object_pairs_hook=_unique_json_object
        )
    except OSError as error:
        raise SetupFailure(f"could not read characterization reference {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise SetupFailure(f"malformed characterization reference {path}: {error}") from error
    except ValueError as error:
        raise SetupFailure(f"invalid characterization reference {path}: {error}") from error

    if not isinstance(document, dict):
        raise SetupFailure(
            f"invalid characterization reference structure in {path}: "
            "top level must be an object"
        )
    reference_schema = document.get("reference_schema")
    if reference_schema is None:
        metadata_inputs = None
        input_sha256 = None
        characterization_metadata = None
    else:
        if not isinstance(reference_schema, str) or not reference_schema:
            raise SetupFailure("reference_schema must be a nonempty string")
        (
            metadata_inputs,
            input_sha256,
            legacy_provenance,
            characterization_metadata,
        ) = (
            _load_reference_metadata(document, reference_schema)
        )
    if reference_schema is None:
        legacy_provenance = None
    comparison_schema = document.get("comparison_schema")
    if comparison_schema is not None:
        if comparison_schema != COMPARISON_SCHEMA:
            raise SetupFailure(
                f"unsupported comparison_schema: {comparison_schema!r}"
            )
    require_explicit_exact = comparison_schema == COMPARISON_SCHEMA

    try:
        case_name = document["case"]
        if not isinstance(case_name, str) or not case_name:
            raise ValueError("case must be a nonempty string")
        zone_items = document["expected_zones"]
        if not isinstance(zone_items, list) or not all(
            isinstance(zone, int) and not isinstance(zone, bool)
            for zone in zone_items
        ):
            raise ValueError("expected_zones must be a list of integers")
        expected_zones = tuple(zone_items)
        if not expected_zones or len(set(expected_zones)) != len(expected_zones):
            raise ValueError("expected_zones is empty or contains duplicates")
        final_steps = _load_zone_integers(
            document["final_step"], expected_zones, "final_step"
        )
        fields_by_name = {
            name: _load_tolerance(
                item,
                expected_zones,
                f"fields.{name}",
                require_explicit_exact=require_explicit_exact,
            )
            for name, item in document["fields"].items()
        }
        mass_fractions_by_species = {
            name: _load_zone_floats(
                value, expected_zones, f"mass_fractions.{name}"
            )
            for name, value in document["mass_fractions"].items()
        }
        negative_mass_fractions = [
            f"{species} (zone {zone})"
            for species, zone_values in mass_fractions_by_species.items()
            for zone, value in zone_values.items()
            if value < 0.0
        ]
        if negative_mass_fractions:
            raise ValueError(
                "mass_fractions contains negative values: "
                + ", ".join(negative_mass_fractions)
            )
        mass_fraction_selection = _load_composition_selection(
            document["mass_fraction_selection"], expected_zones
        )
        tolerance_items = document["mass_fraction_tolerances"]
        if not isinstance(tolerance_items, dict):
            raise ValueError("mass_fraction_tolerances must be an object")
        all_selected_tolerances = (
            _load_tolerance_bounds(
                tolerance_items["all_selected"],
                expected_zones,
                "mass_fraction_tolerances.all_selected",
                require_explicit_exact=require_explicit_exact,
            )
            if set(tolerance_items) == {"all_selected"}
            else None
        )
        mass_fraction_tolerances_by_species = (
            {}
            if all_selected_tolerances is not None
            else {
                name: _load_tolerance_bounds(
                    item,
                    expected_zones,
                    f"mass_fraction_tolerances.{name}",
                    require_explicit_exact=require_explicit_exact,
                )
                for name, item in tolerance_items.items()
            }
        )
        norm_limit_item = document.get("composition_norm_limits")
        composition_norm_limits = (
            None
            if norm_limit_item is None
            else _load_composition_norm_limits(norm_limit_item, expected_zones)
        )
        mass_fraction_sum_atols = _load_zone_floats(
            document["mass_fraction_sum_atol"],
            expected_zones,
            "mass_fraction_sum_atol",
        )
        printed_sum_item = document.get("mass_fraction_printed_sum")
        if printed_sum_item is None:
            if require_explicit_exact:
                raise ValueError(
                    "mass_fraction_printed_sum is required by the comparison schema"
                )
            mass_fraction_printed_sum_tolerances = None
        else:
            mass_fraction_printed_sum_tolerances = _load_tolerance(
                printed_sum_item,
                expected_zones,
                "mass_fraction_printed_sum",
                require_explicit_exact=require_explicit_exact,
            )
        counter_items = document.get("solver_counters")
        if counter_items is None:
            if reference_schema is not None:
                raise ValueError(
                    "solver_counters is required by the reference schema"
                )
            solver_counters = None
        else:
            if not isinstance(counter_items, dict) or set(counter_items) != {
                "ts",
                "nr",
                "jacobian",
                "derivative",
                "cross_section",
            }:
                raise ValueError(
                    "solver_counters must contain exactly ts, nr, jacobian, "
                    "derivative, and cross_section"
                )
            counters_by_name = {
                name: _load_zone_integers(
                    values, expected_zones, f"solver_counters.{name}"
                )
                for name, values in counter_items.items()
            }
            solver_counters = {
                zone: SolverCounters(
                    ts=counters_by_name["ts"][zone],
                    nr=counters_by_name["nr"][zone],
                    jacobian=counters_by_name["jacobian"][zone],
                    derivative=counters_by_name["derivative"][zone],
                    cross_section=counters_by_name["cross_section"][zone],
                )
                for zone in expected_zones
            }
    except (AttributeError, KeyError, TypeError, ValueError) as error:
        raise SetupFailure(f"invalid characterization reference structure in {path}: {error}") from error

    fields = {
        zone: {name: policies[zone] for name, policies in fields_by_name.items()}
        for zone in expected_zones
    }
    mass_fractions = {
        zone: {
            species: values[zone]
            for species, values in mass_fractions_by_species.items()
        }
        for zone in expected_zones
    }
    required_fields = {
        "target_time",
        "temperature_gk",
        "density",
        "electron_fraction",
    }
    fields_with_achieved_time = required_fields | {"achieved_time"}
    if set(fields_by_name) not in (required_fields, fields_with_achieved_time):
        raise SetupFailure(
            "reference fields must contain the established scalar fields, "
            "optionally including achieved_time; "
            f"found {sorted(fields_by_name)}"
        )
    if reference_schema is not None and "achieved_time" not in fields_by_name:
        raise SetupFailure(
            "reference field achieved_time is required by the reference schema"
        )
    if not mass_fractions_by_species or not all(
        isinstance(species, str) and species
        for species in mass_fractions_by_species
    ):
        raise SetupFailure("reference composition is empty or invalid")
    if composition_norm_limits is not None and not mass_fractions_by_species:
        raise SetupFailure(
            "composition norm limits require a complete reference vector"
        )
    if not mass_fraction_tolerances_by_species and all_selected_tolerances is None:
        raise SetupFailure("reference mass-fraction tolerance selection is empty")
    selected_species = {
        species
        for zone_selection in mass_fraction_selection.values()
        for species in zone_selection
    }
    unknown_selected_species = selected_species.difference(mass_fractions_by_species)
    if unknown_selected_species:
        raise SetupFailure(
            "reference mass-fraction selection has no matching composition value: "
            + ", ".join(sorted(unknown_selected_species))
        )
    if all_selected_tolerances is None:
        unknown_tolerance_species = set(mass_fraction_tolerances_by_species).difference(
            mass_fractions_by_species
        )
        if unknown_tolerance_species:
            raise SetupFailure(
                "reference mass-fraction tolerances have no matching composition value: "
                + ", ".join(sorted(unknown_tolerance_species))
            )
        missing_tolerances = selected_species.difference(
            mass_fraction_tolerances_by_species
        )
        if missing_tolerances:
            raise SetupFailure(
                "reference mass-fraction selection has no matching tolerance: "
                + ", ".join(sorted(missing_tolerances))
            )
        unused_tolerances = set(mass_fraction_tolerances_by_species).difference(
            selected_species
        )
        if unused_tolerances:
            raise SetupFailure(
                "reference mass-fraction tolerances are not selected in any zone: "
                + ", ".join(sorted(unused_tolerances))
            )
    mass_fraction_tolerances = {
        zone: {
            species: (
                all_selected_tolerances[zone]
                if all_selected_tolerances is not None
                else mass_fraction_tolerances_by_species[species][zone]
            )
            for species in mass_fraction_selection[zone]
        }
        for zone in expected_zones
    }
    if any(value < 0 for value in mass_fraction_sum_atols.values()):
        raise SetupFailure("reference mass_fraction_sum_atol must be finite and nonnegative")
    if mass_fraction_printed_sum_tolerances is not None:
        for zone in expected_zones:
            canonical_sum = math.fsum(mass_fractions[zone].values())
            policy_value = mass_fraction_printed_sum_tolerances[zone].value
            if policy_value != canonical_sum:
                raise SetupFailure(
                    "reference mass_fraction_printed_sum value does not match "
                    f"the canonical composition sum for zone {zone}: "
                    f"{policy_value:.17g} != {canonical_sum:.17g}"
                )
    return CharacterizationReference(
        case_name=case_name,
        expected_zones=expected_zones,
        final_steps=final_steps,
        fields=fields,
        mass_fractions=mass_fractions,
        mass_fraction_tolerances=mass_fraction_tolerances,
        mass_fraction_sum_atols=mass_fraction_sum_atols,
        mass_fraction_printed_sum_tolerances=mass_fraction_printed_sum_tolerances,
        composition_norm_limits=composition_norm_limits,
        solver_counters=solver_counters,
        reference_schema=reference_schema,
        comparison_schema=comparison_schema,
        metadata_inputs=metadata_inputs,
        input_sha256=input_sha256,
        legacy_provenance=legacy_provenance,
        characterization_metadata=characterization_metadata,
    )


def validate_reference_for_case(
    case: RegressionCase, reference: CharacterizationReference
) -> None:
    """Require one reference to match its case and per-zone selection policy."""

    if reference.case_name != case.name:
        raise SetupFailure(
            "reference case does not match the case definition: "
            f"{reference.case_name!r} != {case.name!r}"
        )
    if reference.reference_schema != case.reference_schema:
        raise SetupFailure(
            "reference schema does not match the case definition: "
            f"{reference.reference_schema!r} != {case.reference_schema!r}"
        )
    if reference.comparison_schema != case.comparison_schema:
        raise SetupFailure(
            "reference comparison schema does not match the case definition: "
            f"{reference.comparison_schema!r} != {case.comparison_schema!r}"
        )
    if reference.expected_zones != case.expected_zones:
        raise SetupFailure(
            "reference expected_zones does not match the case definition: "
            f"{reference.expected_zones} != {case.expected_zones}"
        )
    for zone in case.expected_zones:
        if tuple(reference.mass_fractions[zone]) != case.expected_species:
            raise SetupFailure(
                "composition reference does not match the case species for "
                f"zone {zone}: {tuple(reference.mass_fractions[zone])} "
                f"!= {case.expected_species}"
            )
    if case.reference_schema is not None:
        if (
            reference.metadata_inputs is None
            or reference.input_sha256 is None
            or reference.legacy_provenance is None
            or reference.characterization_metadata is None
        ):
            raise SetupFailure(
                "reference schema requires complete provenance and input hashes"
            )
        if (
            reference.legacy_provenance["legacy_id"]
            != case.expected_legacy_id
            or reference.legacy_provenance["maintained_solver"]
            != case.expected_solver
            or reference.legacy_provenance["effective_runtime_controls"]
            != case.expected_effective_controls
        ):
            raise SetupFailure(
                "reference solver provenance does not match the case definition"
            )
        if (
            reference.characterization_metadata
            != case.expected_characterization_metadata
        ):
            raise SetupFailure(
                "reference characterization metadata does not match the case definition"
            )
        provenance_paths = (
            case.control,
            *(case.network_data / filename for filename in case.network_inputs),
            *case.trajectories,
            *(item.source for item in case.staged_inputs),
            case.helm_table,
        )
        repository_root = Path(
            os.path.commonpath(str(path.resolve()) for path in provenance_paths)
        )

        def input_label(path: Path) -> str:
            return path.resolve().relative_to(repository_root).as_posix()

        path_labels = {
            input_label(path): path for path in provenance_paths
        }
        expected_inputs: dict[str, object] = {
            "control": input_label(case.control),
            "network_data": input_label(case.network_data),
            "network_source_inputs": list(case.network_inputs),
            "eos_table": input_label(case.helm_table),
        }
        if case.staged_inputs:
            expected_inputs["staged_inputs"] = [
                {
                    "source": input_label(item.source),
                    "destination": item.destination.as_posix(),
                }
                for item in case.staged_inputs
            ]
        else:
            expected_inputs["initial_abundance"] = input_label(
                case.network_data / "ab_co"
            )
            expected_inputs["trajectories"] = [
                input_label(path) for path in case.trajectories
            ]
        if dict(reference.metadata_inputs) != expected_inputs:
            raise SetupFailure(
                "reference input provenance does not match the case definition"
            )
        if set(reference.input_sha256) != set(path_labels):
            raise SetupFailure(
                "reference input hashes do not match the complete case input set"
            )
        for label, source in path_labels.items():
            try:
                actual_digest = hashlib.sha256(source.read_bytes()).hexdigest()
            except OSError as error:
                raise SetupFailure(
                    f"could not hash reference input {source}: {error}"
                ) from error
            expected_digest = reference.input_sha256[label]
            if actual_digest != expected_digest:
                raise SetupFailure(
                    f"reference input hash does not match {label}: "
                    f"{actual_digest} != {expected_digest}"
                )
    for zone in case.expected_zones:
        comparison_species = comparison_species_for_zone(
            case.expected_species, reference.mass_fractions[zone]
        )
        selected_species = tuple(reference.mass_fraction_tolerances[zone])
        if selected_species != comparison_species:
            raise SetupFailure(
                "reference pass/fail species do not match the per-zone "
                f"composition policy for zone {zone}: "
                f"{selected_species} != {comparison_species}"
            )


def _difference(actual: float, reference: Tolerance) -> tuple[bool, float, float]:
    absolute_difference = abs(actual - reference.value)
    if reference.exact:
        allowed = 0.0
    else:
        allowed = reference.atol + reference.rtol * abs(reference.value)
        if allowed > 0.0:
            # Decimal policy values and parsed endpoints are represented as
            # binary floats.  Cover only the rounding in forming their
            # subtraction and the allowance; a zero allowance remains exact.
            allowed += (
                math.ulp(actual)
                + math.ulp(reference.value)
                + math.ulp(allowed)
            )
    return absolute_difference <= allowed, absolute_difference, allowed


def calculate_composition_norms(
    states: Sequence[FinalState], reference: CharacterizationReference
) -> tuple[CompositionNorms, ...]:
    """Calculate diagnostic-only norms over the complete composition vector."""

    diagnostics: list[CompositionNorms] = []
    for state in states:
        expected_mass_fractions = reference.mass_fractions[state.zone]
        differences = {
            species: abs(state.mass_fractions[species] - expected)
            for species, expected in expected_mass_fractions.items()
        }
        maximum_species = max(differences, key=differences.__getitem__)
        linf = differences[maximum_species]
        diagnostics.append(
            CompositionNorms(
                zone=state.zone,
                l1=math.fsum(differences.values()),
                l2=math.sqrt(math.fsum(value * value for value in differences.values())),
                linf=linf,
                linf_species=maximum_species if linf > 0.0 else None,
            )
        )
    return tuple(diagnostics)


def reanchor_reference_to_states(
    reference: CharacterizationReference,
    states: Sequence[FinalState],
    *,
    case_name: str,
) -> CharacterizationReference:
    """Reuse a comparison policy with values generated by a same-source run."""

    states_by_zone = {state.zone: state for state in states}
    zones = tuple(state.zone for state in states)
    if zones != reference.expected_zones or len(states_by_zone) != len(states):
        raise SetupFailure(
            f"comparison baseline zones {zones} do not match "
            f"{reference.expected_zones}"
        )

    fields: dict[int, dict[str, Tolerance]] = {}
    mass_fractions: dict[int, dict[str, float]] = {}
    for zone in reference.expected_zones:
        state = states_by_zone[zone]
        expected_species = set(reference.mass_fractions[zone])
        if set(state.mass_fractions) != expected_species:
            raise SetupFailure(
                f"comparison baseline species for zone {zone} do not match the policy"
            )
        fields[zone] = {}
        for name, policy in reference.fields[zone].items():
            value = state.time if name == "achieved_time" else getattr(state, name)
            fields[zone][name] = replace(policy, value=value)
        mass_fractions[zone] = dict(state.mass_fractions)

    printed_sum_tolerances = reference.mass_fraction_printed_sum_tolerances
    if printed_sum_tolerances is not None:
        printed_sum_tolerances = {
            zone: replace(
                printed_sum_tolerances[zone],
                value=math.fsum(states_by_zone[zone].mass_fractions.values()),
            )
            for zone in reference.expected_zones
        }

    return replace(
        reference,
        case_name=case_name,
        fields=fields,
        mass_fractions=mass_fractions,
        final_steps={zone: states_by_zone[zone].step for zone in reference.expected_zones},
        mass_fraction_printed_sum_tolerances=printed_sum_tolerances,
        solver_counters={
            zone: states_by_zone[zone].counters for zone in reference.expected_zones
        },
    )


def scale_comparison_tolerances(
    reference: CharacterizationReference, factor: float
) -> CharacterizationReference:
    """Scale non-exact cross-implementation bounds by one explicit factor."""

    if not math.isfinite(factor) or factor < 1.0:
        raise SetupFailure(f"comparison tolerance factor must be finite and >= 1: {factor}")

    def scale_tolerance(policy: Tolerance) -> Tolerance:
        if policy.exact:
            return policy
        return replace(policy, atol=factor * policy.atol, rtol=factor * policy.rtol)

    fields = {
        zone: {name: scale_tolerance(policy) for name, policy in policies.items()}
        for zone, policies in reference.fields.items()
    }
    mass_fraction_tolerances = {
        zone: {
            species: (
                bounds
                if bounds.exact
                else replace(
                    bounds,
                    atol=factor * bounds.atol,
                    rtol=factor * bounds.rtol,
                )
            )
            for species, bounds in policies.items()
        }
        for zone, policies in reference.mass_fraction_tolerances.items()
    }
    printed_sum_tolerances = reference.mass_fraction_printed_sum_tolerances
    if printed_sum_tolerances is not None:
        printed_sum_tolerances = {
            zone: scale_tolerance(policy)
            for zone, policy in printed_sum_tolerances.items()
        }
    composition_norm_limits = reference.composition_norm_limits
    if composition_norm_limits is not None:
        composition_norm_limits = {
            zone: CompositionNormLimits(
                l1=None if limits.l1 is None else factor * limits.l1,
                linf=None if limits.linf is None else factor * limits.linf,
            )
            for zone, limits in reference.composition_norm_limits.items()
        }

    return replace(
        reference,
        fields=fields,
        mass_fraction_tolerances=mass_fraction_tolerances,
        mass_fraction_printed_sum_tolerances=printed_sum_tolerances,
        composition_norm_limits=composition_norm_limits,
    )


def _write_composition_diagnostics(
    work_directory: Path,
    diagnostics: Sequence[CompositionNorms],
    reference: CharacterizationReference,
) -> None:
    limits = reference.composition_norm_limits
    document = {
        "status": (
            "L2 is diagnostic-only; L1 and L-infinity are pass/fail gates "
            "only when this reference supplies limits"
        ),
        "vector": "absolute mass-fraction errors for every species in the case",
        "zones": [
            {
                "zone": item.zone,
                "l1": item.l1,
                "l2": item.l2,
                "linf": item.linf,
                "linf_species": item.linf_species,
                "l1_limit": None if limits is None else limits[item.zone].l1,
                "linf_limit": None if limits is None else limits[item.zone].linf,
            }
            for item in diagnostics
        ],
    }
    path = work_directory / "composition_error_norms.json"
    try:
        path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    except OSError as error:
        raise ExecutionFailure(
            f"could not preserve composition diagnostics in {path}: {error}"
        ) from error


def compare_final_states(
    states: Sequence[FinalState], reference: CharacterizationReference
) -> None:
    failures: list[str] = []
    zones = tuple(state.zone for state in states)
    if zones != reference.expected_zones:
        failures.append(
            f"case {reference.case_name} zone records {zones} != expected "
            f"{reference.expected_zones}"
        )

    for state in states:
        if state.zone not in reference.fields:
            failures.append(
                f"case {reference.case_name} zone {state.zone} has no characterization reference"
            )
            continue

        field_policies = reference.fields[state.zone]
        target_policy = field_policies["target_time"]
        if "achieved_time" not in field_policies:
            completion_policy = Tolerance(
                state.target_time,
                target_policy.atol,
                target_policy.rtol,
                target_policy.exact,
            )
            completed, difference, allowed = _difference(
                state.time, completion_policy
            )
            if not completed:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} achieved_time: "
                    f"actual={state.time:.9e}, reference={completion_policy.value:.9e}, "
                    f"|difference|={difference:.3e}, allowed={allowed:.3e}"
                )

        for field, policy in field_policies.items():
            actual = (
                state.time if field == "achieved_time" else getattr(state, field)
            )
            passed, difference, allowed = _difference(actual, policy)
            if not passed:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} {field}: actual={actual:.9e}, "
                    f"reference={policy.value:.9e}, |difference|={difference:.3e}, "
                    f"allowed={allowed:.3e} (atol={policy.atol:.3e}, "
                    f"rtol={policy.rtol:.3e})"
                )

        expected_mass_fractions = reference.mass_fractions[state.zone]
        for species, bounds in reference.mass_fraction_tolerances[state.zone].items():
            actual = state.mass_fractions[species]
            policy = Tolerance(
                expected_mass_fractions[species],
                bounds.atol,
                bounds.rtol,
                bounds.exact,
            )
            passed, difference, allowed = _difference(actual, policy)
            if not passed:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} {species} mass fraction: actual={actual:.9e}, "
                    f"reference={policy.value:.9e}, |difference|={difference:.3e}, "
                    f"allowed={allowed:.3e} (atol={policy.atol:.3e}, "
                    f"rtol={policy.rtol:.3e})"
                )

        norms = calculate_composition_norms((state,), reference)[0]
        limits = reference.composition_norm_limits
        if limits is not None:
            zone_limits = limits[state.zone]
            if zone_limits.l1 is not None and norms.l1 > zone_limits.l1:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} L1: "
                    f"actual norm={norms.l1:.3e}, allowed norm={zone_limits.l1:.3e}"
                )
            if zone_limits.linf is not None and norms.linf > zone_limits.linf:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} L-infinity: "
                    f"actual norm={norms.linf:.3e}, allowed norm={zone_limits.linf:.3e}, "
                    f"species={norms.linf_species}"
                )

        negative_species = [
            species for species, value in state.mass_fractions.items() if value < 0
        ]
        if negative_species:
            failures.append(
                f"zone {state.zone} has negative mass fractions: {', '.join(negative_species)}"
            )
        printed_sum = math.fsum(state.mass_fractions.values())
        printed_sum_tolerances = reference.mass_fraction_printed_sum_tolerances
        if printed_sum_tolerances is not None:
            policy = printed_sum_tolerances[state.zone]
            passed, difference, allowed = _difference(printed_sum, policy)
            if not passed:
                failures.append(
                    f"case {reference.case_name} zone {state.zone} printed mass-fraction sum: "
                    f"actual={printed_sum:.12e}, reference={policy.value:.12e}, "
                    f"|difference|={difference:.3e}, allowed={allowed:.3e} "
                    f"(atol={policy.atol:.3e}, rtol={policy.rtol:.3e})"
                )
        mass_fraction_sum = math.fsum(state.mass_fractions.values())
        sum_atol = reference.mass_fraction_sum_atols[state.zone]
        normalization_passed, _, normalization_allowed = _difference(
            mass_fraction_sum,
            Tolerance(1.0, sum_atol, 0.0, exact=sum_atol == 0.0),
        )
        if not normalization_passed:
            failures.append(
                f"case {reference.case_name} zone {state.zone} recomputed mass-fraction normalization={mass_fraction_sum:.12e}; "
                f"allowed |sum - 1| <= {normalization_allowed:.3e}"
            )

    if failures:
        raise ComparisonFailure("numerical characterization failed:\n  " + "\n  ".join(failures))


def compare_equivalent_zone_groups(
    states: Sequence[FinalState], groups: Sequence[Sequence[int]]
) -> None:
    """Require declared repeated-input endpoints to be identical when characterized so."""

    states_by_zone = {state.zone: state for state in states}
    failures: list[str] = []
    for group in groups:
        first, *rest = group
        baseline = states_by_zone[first]
        for zone in rest:
            candidate = states_by_zone[zone]
            if (
                candidate.target_time != baseline.target_time
                or candidate.time != baseline.time
                or candidate.temperature_gk != baseline.temperature_gk
                or candidate.density != baseline.density
                or candidate.electron_fraction != baseline.electron_fraction
                or candidate.mass_fractions != baseline.mass_fractions
            ):
                failures.append(
                    f"equivalent zones {tuple(group)} have non-identical endpoint values"
                )
                break
    if failures:
        raise ComparisonFailure("repeated-zone equivalence failed:\n  " + "\n  ".join(failures))


def run_and_compare(
    executable: Path,
    case: RegressionCase,
    work_directory: Path,
    *,
    timeout_seconds: float,
) -> ProcessResult:
    reference = load_reference(case.reference)
    validate_reference_for_case(case, reference)
    prepared = prepare_work_directory(case, work_directory)
    result = run_xnet(
        executable,
        case,
        prepared,
        timeout_seconds=timeout_seconds,
    )
    diagnostic_path = prepared / "net_diag01"
    try:
        diagnostic = diagnostic_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise ParsingFailure(
            f"could not read diagnostic {diagnostic_path}: {error}; artifacts: {prepared}"
        ) from error
    try:
        missing_markers = tuple(
            marker
            for marker in case.required_diagnostic_markers
            if marker not in diagnostic
        )
        if missing_markers:
            raise ParsingFailure(
                "required diagnostic marker is missing: "
                + ", ".join(repr(marker) for marker in missing_markers)
            )
        states = parse_diagnostic(
            diagnostic,
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
        diagnostics = calculate_composition_norms(states, reference)
        _write_composition_diagnostics(prepared, diagnostics, reference)
        compare_final_states(states, reference)
        compare_equivalent_zone_groups(states, case.equivalent_zone_groups)
    except (ParsingFailure, ComparisonFailure) as error:
        raise type(error)(f"{error}; artifacts: {prepared}") from None
    return result
