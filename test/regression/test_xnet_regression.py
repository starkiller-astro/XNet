"""Focused tests for the XNet regression runner and comparator."""

from dataclasses import replace
import base64
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest

from xnet_regression import (
    ALPHA_SPECIES,
    batch_alpha_case,
    bdf_sn160_case,
    CharacterizationReference,
    CompositionNormLimits,
    CompositionNorms,
    ComparisonFailure,
    ExecutionFailure,
    FinalState,
    ParsingFailure,
    RegressionCase,
    SN160_SPECIES,
    SolverCounters,
    StagedInput,
    SetupFailure,
    Tolerance,
    ToleranceBounds,
    TORCH47_SPECIES,
    calculate_composition_norms,
    comparison_species_for_zone,
    compare_final_states,
    heat_alpha_case,
    heat_sn160_case,
    load_reference,
    nse_sn160_case,
    parse_diagnostic,
    prepare_work_directory,
    reanchor_reference_to_states,
    run_and_compare,
    run_xnet,
    scale_comparison_tolerances,
    strip_timer_sections,
    tnsn_alpha_case,
    tnsn_torch47_case,
    validate_reference_for_case,
    validate_executable,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _fabricated_final_diagnostic(*, timer_total: str = "1.000E-02") -> str:
    """Return invented parser input with no physical or tnsn_alpha meaning."""

    return f"""MyId    0    1
End     1    42 2.0000000E+00 2.0000000E+00 2.0000000E+00 4.0000000E+06 5.0000000E-01
  he4 1.0000000E-01   c12 2.0000000E-01   o16 3.0000000E-01  ne20 0.0000000E+00
 mg24 0.0000000E+00  si28 1.5000000E-01   s32 2.5000000E-01  ar36 0.0000000E+00
 ca40 0.0000000E+00  ti44 0.0000000E+00  cr48 0.0000000E+00  fe52 0.0000000E+00
 ni56 0.0000000E+00  zn60 0.0000000E+00
Counters:  Zone        TS        NR  Jacobian     Deriv CrossSect
              1        42        42        42        43        43
Timers Summary:
        Total      {timer_total}
        Solver     2.000E-03
"""


def _fabricated_diagnostic_with_species(species: tuple[str, ...]) -> str:
    lines = _fabricated_final_diagnostic().splitlines()
    abundance_start = next(
        index for index, line in enumerate(lines) if line.startswith("End")
    ) + 1
    abundance_end = next(
        index for index, line in enumerate(lines) if line.startswith("Counters:")
    )
    abundance_rows = [
        " ".join(
            f"{name:>5} {1.0 / len(species):.7E}"
            for name in species[index : index + 4]
        )
        for index in range(0, len(species), 4)
    ]
    return "\n".join(
        lines[:abundance_start] + abundance_rows + lines[abundance_end:]
    )


def _matching_unit_reference() -> CharacterizationReference:
    """Return expectations independently matching the fabricated parser input."""

    fields = {
        name: Tolerance(value, atol=1e-6, rtol=1e-6)
        for name, value in {
            "target_time": 2.0,
            "temperature_gk": 2.0,
            "density": 4.0e6,
            "electron_fraction": 0.5,
        }.items()
    }
    mass_fractions = {
        species: {
            "he4": 0.1,
            "c12": 0.2,
            "o16": 0.3,
            "si28": 0.15,
            "s32": 0.25,
        }.get(species, 0.0)
        for species in (
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
    }
    return CharacterizationReference(
        case_name="fake",
        expected_zones=(1,),
        final_steps={1: 42},
        fields={1: fields},
        mass_fractions={1: mass_fractions},
        mass_fraction_tolerances={
            1: {
                "si28": ToleranceBounds(atol=1e-8, rtol=1e-8),
                "s32": ToleranceBounds(atol=1e-8, rtol=1e-8),
            }
        },
        mass_fraction_sum_atols={1: 1e-8},
    )


def _state_from_reference(
    reference: CharacterizationReference, zone: int
) -> FinalState:
    fields = reference.fields[zone]
    return FinalState(
        zone=zone,
        step=reference.final_steps[zone],
        target_time=fields["target_time"].value,
        time=fields["target_time"].value,
        temperature_gk=fields["temperature_gk"].value,
        density=fields["density"].value,
        electron_fraction=fields["electron_fraction"].value,
        mass_fractions=reference.mass_fractions[zone],
        counters=SolverCounters(reference.final_steps[zone], 0, 0, 0, 0),
    )


def _fake_case(tmp_path: Path) -> RegressionCase:
    return RegressionCase(
        name="fake",
        control=tmp_path / "unused-control",
        network_data=tmp_path / "unused-data",
        trajectories=(tmp_path / "unused-trajectory",),
        helm_table=tmp_path / "unused-table",
        reference=tmp_path / "unused-reference",
        expected_zones=(1,),
        expected_species=("he4",),
        network_inputs=("sunet",),
    )


def _make_executable(tmp_path: Path, body: str) -> Path:
    executable = tmp_path / "fake_xnet.py"
    executable.write_text(f"#!{sys.executable}\n{body}\n", encoding="utf-8")
    executable.chmod(0o755)
    return executable


def _diagnostic_for_states(
    states: tuple[FinalState, ...],
    groups: tuple[tuple[int, ...], ...] | None = None,
) -> str:
    """Build a structurally valid XNet diagnostic for registered-path tests."""

    lines: list[str] = []
    states_by_zone = {state.zone: state for state in states}
    if groups is None:
        groups = tuple((state.zone,) for state in states)
    for group in groups:
        for zone in group:
            state = states_by_zone[zone]
            lines.append(
                "End "
                f"{state.zone} {state.step} {state.target_time:.16E} "
                f"{state.time:.16E} {state.temperature_gk:.16E} "
                f"{state.density:.16E} {state.electron_fraction:.16E}"
            )
            for species, value in state.mass_fractions.items():
                lines.append(f" {species} {value:.16E}")
        lines.extend(
            (
                "Counters:  Zone        TS        NR  Jacobian     Deriv CrossSect",
            )
        )
        for zone in group:
            state = states_by_zone[zone]
            lines.append(
                f"{state.zone} {state.counters.ts} {state.counters.nr} "
                f"{state.counters.jacobian} {state.counters.derivative} "
                f"{state.counters.cross_section}"
            )
        lines.extend(
            (
                "Timers Summary:",
                "        Total      1.000E-02",
            )
        )
    return "\n".join(lines) + "\n"


def _registered_case_executable(
    tmp_path: Path, case: RegressionCase, states: tuple[FinalState, ...]
) -> Path:
    """Return a fake XNet that writes fresh complete output in its work directory."""

    marker_text = "".join(
        f"{marker}\n" for marker in case.required_diagnostic_markers
    )
    diagnostic = base64.b64encode(
        (
            marker_text
            + _diagnostic_for_states(states, case.expected_diagnostic_groups)
        ).encode("utf-8")
    ).decode("ascii")
    body = (
        "from pathlib import Path\n"
        "import base64\n"
        f"Path('net_diag01').write_bytes(base64.b64decode({diagnostic!r}))\n"
        f"for name in {case.required_outputs[1:]!r}:\n"
        "    Path(name).write_bytes(b'fresh')"
    )
    return _make_executable(tmp_path, body)


def _registered_diagnostic_executable(
    tmp_path: Path, case: RegressionCase, diagnostic: str
) -> Path:
    encoded = base64.b64encode(diagnostic.encode("utf-8")).decode("ascii")
    body = (
        "from pathlib import Path\n"
        "import base64\n"
        f"Path('net_diag01').write_bytes(base64.b64decode({encoded!r}))\n"
        f"for name in {case.required_outputs[1:]!r}:\n"
        "    Path(name).write_bytes(b'fresh')"
    )
    return _make_executable(tmp_path, body)


def _run_registered_pytest(
    tmp_path: Path,
    case: RegressionCase,
    executable: Path,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            f"test/regression/test_regression.py::test_{case.name}",
            f"--xnet-executable={executable}",
            f"--basetemp={tmp_path / 'artifacts'}",
        ],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
        env=environment,
    )


def test_known_good_comparison_passes() -> None:
    states = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )
    compare_final_states(states, _matching_unit_reference())


def test_same_source_reference_reuses_policy_with_generated_values() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    policy = _matching_unit_reference()
    baseline = replace(state, temperature_gk=2.25, step=99)

    reference = reanchor_reference_to_states(
        policy,
        (baseline,),
        case_name="fake dense versus sparse",
    )

    assert reference.case_name == "fake dense versus sparse"
    assert reference.fields[1]["temperature_gk"] == replace(
        policy.fields[1]["temperature_gk"], value=2.25
    )
    assert reference.final_steps == {1: 99}
    assert reference.mass_fraction_tolerances == policy.mass_fraction_tolerances
    compare_final_states((baseline,), reference)
    with pytest.raises(ComparisonFailure, match="temperature_gk"):
        compare_final_states(
            (replace(baseline, temperature_gk=2.25001),), reference
        )


def test_same_source_reference_rejects_incomplete_composition() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    incomplete = dict(state.mass_fractions)
    incomplete.pop("zn60")

    with pytest.raises(SetupFailure, match="species for zone 1"):
        reanchor_reference_to_states(
            _matching_unit_reference(),
            (replace(state, mass_fractions=incomplete),),
            case_name="incomplete",
        )


def test_cross_solver_factor_scales_only_nonexact_difference_bounds() -> None:
    policy = _matching_unit_reference()
    fields = {1: dict(policy.fields[1])}
    fields[1]["density"] = Tolerance(4.0e6, 0.0, 0.0, exact=True)
    scaled = scale_comparison_tolerances(replace(policy, fields=fields), 2.0)

    assert scaled.fields[1]["temperature_gk"].atol == 2.0e-6
    assert scaled.fields[1]["temperature_gk"].rtol == 2.0e-6
    assert scaled.fields[1]["density"] == fields[1]["density"]
    assert scaled.mass_fraction_tolerances[1]["si28"].atol == 2.0e-8
    assert scaled.mass_fraction_sum_atols == policy.mass_fraction_sum_atols
    with pytest.raises(SetupFailure, match="finite and >= 1"):
        scale_comparison_tolerances(policy, 0.5)


def test_mass_fraction_expectation_comes_from_complete_reference() -> None:
    states = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )
    reference = _matching_unit_reference()
    changed_composition = {1: dict(reference.mass_fractions[1])}
    changed_composition[1]["si28"] = 0.16

    with pytest.raises(ComparisonFailure, match="si28 mass fraction"):
        compare_final_states(
            states, replace(reference, mass_fractions=changed_composition)
        )


def test_mass_fraction_tolerance_requires_composition_value(tmp_path: Path) -> None:
    reference_path = tmp_path / "inconsistent-reference.json"
    reference_path.write_text(
        json.dumps(
            {
                "case": "fake",
                "expected_zones": [1],
                "final_step": 42,
                "fields": {
                    name: {"value": value, "atol": 1e-6, "rtol": 1e-6}
                    for name, value in {
                        "target_time": 2.0,
                        "temperature_gk": 2.0,
                        "density": 4.0e6,
                        "electron_fraction": 0.5,
                    }.items()
                },
                "mass_fractions": {"c12": 1.0},
                "mass_fraction_selection": {"1": ["si28"]},
                "mass_fraction_tolerances": {
                    "si28": {"atol": 1e-8, "rtol": 1e-8}
                },
                "mass_fraction_sum_atol": 1e-8,
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(SetupFailure, match="no matching composition value: si28"):
        load_reference(reference_path)


def test_zone_specific_reference_values_are_expanded(tmp_path: Path) -> None:
    reference_path = tmp_path / "zone-specific-reference.json"
    reference_path.write_text(
        json.dumps(
            {
                "case": "fake",
                "expected_zones": [1, 2],
                "final_step": {"1": 42, "2": 43},
                "fields": {
                    name: {
                        "value": {"1": first, "2": second},
                        "atol": 1e-6,
                        "rtol": 1e-6,
                    }
                    for name, (first, second) in {
                        "target_time": (2.0, 1.0),
                        "temperature_gk": (2.0, 3.0),
                        "density": (4.0e6, 8.0e6),
                        "electron_fraction": (0.5, 0.5),
                    }.items()
                },
                "mass_fractions": {
                    "c12": {"1": 0.4, "2": 0.3},
                    "o16": {"1": 0.6, "2": 0.7},
                },
                "mass_fraction_selection": {
                    "1": ["c12"],
                    "2": ["c12"],
                },
                "mass_fraction_tolerances": {
                    "c12": {
                        "atol": {"1": 1e-8, "2": 2e-8},
                        "rtol": 1e-8,
                    }
                },
                "mass_fraction_sum_atol": {"1": 1e-8, "2": 2e-8},
            }
        ),
        encoding="utf-8",
    )

    reference = load_reference(reference_path)

    assert reference.final_steps == {1: 42, 2: 43}
    assert reference.fields[2]["temperature_gk"].value == 3.0
    assert reference.mass_fractions[2]["o16"] == 0.7
    assert reference.mass_fraction_tolerances[2]["c12"].atol == 2e-8
    assert reference.mass_fraction_sum_atols == {1: 1e-8, 2: 2e-8}


def test_zone_specific_reference_requires_every_zone(tmp_path: Path) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/tnsn_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["final_step"] = {"1": 2841}
    reference_path = tmp_path / "incomplete-zone-reference.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="must define zones"):
        load_reference(reference_path)


def test_comparison_schema_requires_explicit_nonpermissive_policies(
    tmp_path: Path,
) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/tnsn_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["fields"]["temperature_gk"].pop("exact")
    reference_path = tmp_path / "missing-exact-policy.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="exact"):
        load_reference(reference_path)


def test_registered_case_rejects_missing_comparison_schema_before_execution(
    tmp_path: Path,
) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document.pop("comparison_schema")
    reference_path = tmp_path / "missing-comparison-schema.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="must contain exactly"):
        run_and_compare(
            tmp_path / "not-run", replace(case, reference=reference_path),
            tmp_path / "work", timeout_seconds=1.0
        )


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda document: document.__setitem__(
                "comparison_schema", "xnet-comparison-v999"
            ),
            "unsupported comparison_schema",
        ),
        (
            lambda document: document["fields"]["temperature_gk"].update(
                {"exact": True, "atol": 1.0}
            ),
            "exact comparison with nonzero tolerance",
        ),
        (
            lambda document: document.pop("mass_fraction_printed_sum"),
            "mass_fraction_printed_sum is required",
        ),
        (
            lambda document: document["mass_fraction_printed_sum"].update(
                {"value": 0.5}
            ),
            "does not match the canonical composition sum",
        ),
        (
            lambda document: document["composition_norm_limits"].update(
                {"l1": -1.0}
            ),
            "negative tolerance",
        ),
        (
            lambda document: document["composition_norm_limits"].update(
                {"l2": 1.0}
            ),
            "one or both of l1 and linf",
        ),
    ],
)
def test_comparison_schema_rejects_malformed_policy(
    tmp_path: Path, mutate, message: str
) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/tnsn_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    mutate(document)
    reference_path = tmp_path / "malformed-comparison-policy.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match=message):
        load_reference(reference_path)


def test_comparison_norm_limit_round_trips_without_shrinking(tmp_path: Path) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/tnsn_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    reference_path = tmp_path / "round-trip-reference.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    reference = load_reference(reference_path)
    assert reference.composition_norm_limits is not None
    assert reference.composition_norm_limits[1].l1 == 2.0e-10
    assert reference.composition_norm_limits[1].linf == 2.0e-10


def test_composition_selection_rejects_missing_zone(tmp_path: Path) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/heat_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["mass_fraction_selection"].pop("6")
    reference_path = tmp_path / "missing-selection-zone.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="missing \\[6\\]"):
        load_reference(reference_path)


def test_composition_selection_rejects_unknown_zone(tmp_path: Path) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/heat_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["mass_fraction_selection"]["7"] = document[
        "mass_fraction_selection"
    ]["6"]
    reference_path = tmp_path / "unknown-selection-zone.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="unknown \\[7\\]"):
        load_reference(reference_path)


def test_composition_selection_rejects_duplicate_zone_key(tmp_path: Path) -> None:
    reference_path = tmp_path / "duplicate-selection-zone.json"
    reference_path.write_text(
        """{
  "expected_zones": [1],
  "mass_fraction_selection": {
    "1": ["c12"],
    "1": ["c12"]
  }
}
""",
        encoding="utf-8",
    )

    with pytest.raises(SetupFailure, match="duplicate JSON object key: 1"):
        load_reference(reference_path)


def test_composition_selection_rejects_duplicate_species(tmp_path: Path) -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["mass_fraction_selection"]["4"].append("o16")
    reference_path = tmp_path / "duplicate-selection-species.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="duplicate species: o16"):
        load_reference(reference_path)


def test_composition_selection_rejects_unknown_species(tmp_path: Path) -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["mass_fraction_selection"]["1"].append("xe999")
    reference_path = tmp_path / "unknown-selection-species.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="no matching composition value: xe999"):
        load_reference(reference_path)


def test_composition_selection_rejects_missing_policy_species(tmp_path: Path) -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["mass_fraction_selection"]["4"].remove("o16")
    reference_path = tmp_path / "missing-selection-species.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    reference = load_reference(reference_path)

    with pytest.raises(SetupFailure, match="per-zone composition policy for zone 4"):
        validate_reference_for_case(case, reference)


def test_composition_selection_requires_selected_species_tolerance(
    tmp_path: Path,
) -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["mass_fraction_tolerances"] = {
        "si28": document["mass_fraction_tolerances"]["all_selected"]
    }
    reference_path = tmp_path / "missing-selection-tolerance.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="no matching tolerance"):
        load_reference(reference_path)


@pytest.mark.parametrize("zone", [6.9, "6", True])
def test_reference_rejects_noninteger_zone_identifiers(
    tmp_path: Path, zone: object
) -> None:
    source = REPOSITORY_ROOT / "test/regression/cases/heat_alpha/reference/final_state.json"
    document = json.loads(source.read_text(encoding="utf-8"))
    document["expected_zones"][-1] = zone
    reference_path = tmp_path / "malformed-zone-reference.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="must be a list of integers"):
        load_reference(reference_path)


def test_value_within_tolerance_passes() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    policy = {1: dict(reference.fields[1])}
    policy[1]["temperature_gk"] = Tolerance(2.0, atol=0.1, rtol=0.0)
    compare_final_states(
        (replace(state, temperature_gk=2.05),), replace(reference, fields=policy)
    )


def test_value_outside_tolerance_fails() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    policy = {1: dict(reference.fields[1])}
    policy[1]["temperature_gk"] = Tolerance(2.0, atol=0.1, rtol=0.0)
    with pytest.raises(ComparisonFailure, match="temperature_gk"):
        compare_final_states(
            (replace(state, temperature_gk=2.2),), replace(reference, fields=policy)
        )


def test_exact_scalar_policy_rejects_any_difference() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    fields = {1: dict(reference.fields[1])}
    fields[1]["temperature_gk"] = Tolerance(2.0, 0.0, 0.0, exact=True)

    with pytest.raises(ComparisonFailure, match=r"allowed=0.000e\+00"):
        compare_final_states(
            (replace(state, temperature_gk=2.0 + 1.0e-12),),
            replace(reference, fields=fields),
        )


def test_relative_scalar_policy_passes_exactly_at_boundary() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    fields = {1: dict(reference.fields[1])}
    fields[1]["temperature_gk"] = Tolerance(2.0, 0.0, 0.25)

    compare_final_states(
        (replace(state, temperature_gk=2.5),), replace(reference, fields=fields)
    )


def test_combined_scalar_policy_honors_boundary_and_zero_reference() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    fields = {1: dict(reference.fields[1])}
    fields[1]["temperature_gk"] = Tolerance(2.0, 0.1, 0.2)
    bounded_reference = replace(reference, fields=fields)

    compare_final_states((replace(state, temperature_gk=2.5),), bounded_reference)
    with pytest.raises(ComparisonFailure, match="temperature_gk"):
        compare_final_states(
            (replace(state, temperature_gk=2.5 + 1.0e-12),),
            bounded_reference,
        )

    fields[1]["temperature_gk"] = Tolerance(0.0, 0.0, 0.2)
    with pytest.raises(ComparisonFailure, match="temperature_gk"):
        compare_final_states(
            (replace(state, temperature_gk=math.nextafter(0.0, math.inf)),),
            replace(reference, fields=fields),
        )


def test_nonexact_scalar_boundary_covers_only_binary_subtraction_roundoff() -> None:
    reference = load_reference(bdf_sn160_case(REPOSITORY_ROOT).reference)
    state = _state_from_reference(reference, 3)
    zone_reference = replace(reference, expected_zones=(3,))
    policy = reference.fields[3]["electron_fraction"]
    last_passing = float.fromhex("0x1.ff4369364f727p-2")
    first_failing = math.nextafter(last_passing, math.inf)

    # The decimal 1e-7 boundary is not exactly representable as a binary
    # subtraction.  Pin the adjacent binary floats so any larger decimal
    # guard, including 1e-13, moves this boundary and fails the test.
    assert last_passing.hex() == "0x1.ff4369364f727p-2"
    assert first_failing.hex() == "0x1.ff4369364f728p-2"
    assert abs(last_passing - policy.value) > policy.atol
    compare_final_states(
        (replace(state, electron_fraction=last_passing),), zone_reference
    )

    with pytest.raises(ComparisonFailure, match="electron_fraction"):
        compare_final_states(
            (replace(state, electron_fraction=first_failing),),
            zone_reference,
        )


def test_exact_selected_species_isolated_from_nonzero_normalization_limit() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = replace(
        _matching_unit_reference(),
        mass_fraction_tolerances={
            1: {"si28": ToleranceBounds(0.0, 0.0, exact=True)}
        },
        mass_fraction_sum_atols={1: 1.0},
    )
    changed = dict(state.mass_fractions)
    changed["si28"] = math.nextafter(changed["si28"], math.inf)

    with pytest.raises(ComparisonFailure, match="si28 mass fraction"):
        compare_final_states(
            (replace(state, mass_fractions=changed),), reference
        )


def test_zero_vector_limits_reject_nonzero_vector_with_normalization_slack() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = replace(
        _matching_unit_reference(),
        mass_fraction_tolerances={1: {}},
        composition_norm_limits={1: CompositionNormLimits(l1=0.0, linf=0.0)},
        mass_fraction_sum_atols={1: 1.0},
    )
    changed = dict(state.mass_fractions)
    amount = math.ulp(changed["c12"])
    changed["c12"] -= amount
    changed["ne20"] += amount

    with pytest.raises(ComparisonFailure) as failure:
        compare_final_states(
            (replace(state, mass_fractions=changed),), reference
        )
    assert "L1" in str(failure.value)
    assert "L-infinity" in str(failure.value)


def test_num_obs_001_old_bounds_fail_and_revised_composition_bounds_pass() -> None:
    torch_reference = load_reference(
        tnsn_torch47_case(REPOSITORY_ROOT).reference
    )
    torch_state = _state_from_reference(torch_reference, 1)
    torch_mass_fractions = dict(torch_state.mass_fractions)
    torch_mass_fractions["s31"] += 1.890e-9
    torch_mass_fractions["co55"] += 2.800e-9
    torch_mass_fractions["mn51"] -= 1.8275e-9
    torch_mass_fractions["p30"] += 0.2685e-9
    torch_observation = replace(
        torch_state, mass_fractions=torch_mass_fractions
    )
    torch_norms = calculate_composition_norms(
        (torch_observation,), torch_reference
    )[0]

    s31_difference = abs(
        torch_mass_fractions["s31"]
        - torch_reference.mass_fractions[1]["s31"]
    )
    co55_difference = abs(
        torch_mass_fractions["co55"]
        - torch_reference.mass_fractions[1]["co55"]
    )
    assert s31_difference > 1.0e-10
    assert co55_difference > 1.0e-10
    assert s31_difference <= torch_reference.mass_fraction_tolerances[1]["s31"].atol
    assert co55_difference <= torch_reference.mass_fraction_tolerances[1]["co55"].atol
    assert torch_norms.l1 == pytest.approx(6.786e-9, abs=1.0e-15)
    assert torch_norms.linf == pytest.approx(2.800e-9, abs=1.0e-15)
    assert torch_norms.l1 > 2.0e-10
    assert torch_norms.linf > 1.0e-10
    assert torch_norms.l1 <= torch_reference.composition_norm_limits[1].l1  # type: ignore[index,union-attr]
    assert torch_norms.linf <= torch_reference.composition_norm_limits[1].linf  # type: ignore[index,union-attr]
    torch_printed_difference = abs(
        math.fsum(torch_mass_fractions.values())
        - torch_reference.mass_fraction_printed_sum_tolerances[1].value  # type: ignore[index]
    )
    assert torch_printed_difference == pytest.approx(3.131e-9, abs=1.0e-15)
    assert 1.0e-10 < torch_printed_difference <= 4.7e-9
    compare_final_states((torch_observation,), torch_reference)
    old_torch_reference = replace(
        torch_reference,
        mass_fraction_tolerances={
            1: {
                species: replace(bounds, atol=1.0e-10)
                for species, bounds in torch_reference.mass_fraction_tolerances[1].items()
            }
        },
        composition_norm_limits={
            1: CompositionNormLimits(l1=2.0e-10, linf=1.0e-10)
        },
        mass_fraction_printed_sum_tolerances={
            1: replace(
                torch_reference.mass_fraction_printed_sum_tolerances[1],  # type: ignore[index]
                atol=1.0e-10,
            )
        },
    )
    with pytest.raises(ComparisonFailure) as old_torch_failure:
        compare_final_states((torch_observation,), old_torch_reference)
    old_torch_diagnostics = str(old_torch_failure.value)
    for diagnostic in (
        "s31 mass fraction",
        "co55 mass fraction",
        "L1",
        "L-infinity",
        "printed mass-fraction sum",
    ):
        assert diagnostic in old_torch_diagnostics

    bdf_reference = load_reference(bdf_sn160_case(REPOSITORY_ROOT).reference)
    old_bdf_reference = replace(
        bdf_reference,
        mass_fraction_printed_sum_tolerances={
            zone: replace(policy, atol=3.0e-7)
            for zone, policy in bdf_reference.mass_fraction_printed_sum_tolerances.items()  # type: ignore[union-attr]
        },
        mass_fraction_sum_atols={
            **bdf_reference.mass_fraction_sum_atols,
            3: 3.042424516643881e-5,
        },
    )
    for zone, actual_sum in (
        (1, 0.9999997928543),
        (3, 1.000030571079),
    ):
        state = _state_from_reference(bdf_reference, zone)
        mass_fractions = dict(state.mass_fractions)
        species = "fe56" if zone == 1 else "fe57"
        mass_fractions[species] += actual_sum - math.fsum(
            mass_fractions.values()
        )
        if zone == 3:
            state = replace(state, electron_fraction=0.49928059)
        observation = replace(state, mass_fractions=mass_fractions)
        printed_policy = bdf_reference.mass_fraction_printed_sum_tolerances[zone]  # type: ignore[index]
        assert abs(actual_sum - printed_policy.value) > 3.0e-7
        assert abs(actual_sum - printed_policy.value) <= printed_policy.atol
        if zone == 3:
            assert abs(actual_sum - 1.0) > 3.042424516643881e-5
            assert abs(actual_sum - 1.0) <= bdf_reference.mass_fraction_sum_atols[zone]
        compare_final_states(
            (observation,), replace(bdf_reference, expected_zones=(zone,))
        )
        with pytest.raises(ComparisonFailure) as old_bdf_failure:
            compare_final_states(
                (observation,),
                replace(old_bdf_reference, expected_zones=(zone,)),
            )
        old_bdf_diagnostics = str(old_bdf_failure.value)
        assert "printed mass-fraction sum" in old_bdf_diagnostics
        if zone == 3:
            assert "recomputed mass-fraction normalization" in old_bdf_diagnostics


def test_complete_vector_limits_are_independent_of_selected_species() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    unselected_change = dict(state.mass_fractions)
    unselected_change["ne20"] += 2.0e-5
    unselected_change["ar36"] -= 2.0e-5
    vector_reference = replace(
        reference,
        composition_norm_limits={1: CompositionNormLimits(l1=1.0e-5, linf=1.0e-5)},
    )

    with pytest.raises(ComparisonFailure, match="case fake zone 1 L1"):
        compare_final_states(
            (replace(state, mass_fractions=unselected_change),), vector_reference
        )

    selected_change = dict(state.mass_fractions)
    selected_change["si28"] += 2.0e-5
    selected_change["c12"] -= 2.0e-5
    permissive_vector_reference = replace(
        reference,
        composition_norm_limits={1: CompositionNormLimits(l1=1.0, linf=1.0)},
    )
    with pytest.raises(ComparisonFailure, match="si28 mass fraction"):
        compare_final_states(
            (replace(state, mass_fractions=selected_change),),
            permissive_vector_reference,
        )


def test_linf_limit_identifies_responsible_species_and_boundary() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    changed = dict(state.mass_fractions)
    changed["ne20"] += 0.125
    changed["c12"] -= 0.125
    limits = {1: CompositionNormLimits(linf=0.125)}
    compare_final_states(
        (replace(state, mass_fractions=changed),),
        replace(reference, composition_norm_limits=limits),
    )

    changed["ne20"] += 1.0e-6
    changed["c12"] -= 1.0e-6
    with pytest.raises(ComparisonFailure, match="species=c12"):
        compare_final_states(
            (replace(state, mass_fractions=changed),),
            replace(reference, composition_norm_limits=limits),
        )


def test_l1_and_linf_gates_fail_independently() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = replace(
        _matching_unit_reference(), mass_fraction_tolerances={1: {}}
    )

    broad = dict(state.mass_fractions)
    for donor, recipient in (("he4", "ne20"), ("c12", "mg24")):
        broad[donor] -= 0.02
        broad[recipient] += 0.02
    with pytest.raises(ComparisonFailure, match="case fake zone 1 L1"):
        compare_final_states(
            (replace(state, mass_fractions=broad),),
            replace(
                reference,
                composition_norm_limits={
                    1: CompositionNormLimits(l1=0.05, linf=0.03)
                },
            ),
        )

    localized = dict(state.mass_fractions)
    localized["he4"] -= 0.02
    localized["ne20"] += 0.02
    with pytest.raises(ComparisonFailure, match="case fake zone 1 L-infinity"):
        compare_final_states(
            (replace(state, mass_fractions=localized),),
            replace(
                reference,
                composition_norm_limits={
                    1: CompositionNormLimits(l1=0.05, linf=0.01)
                },
            ),
        )


def test_final_step_is_diagnostic_only() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]

    compare_final_states(
        (replace(state, step=142),), _matching_unit_reference()
    )


def test_solver_counters_are_diagnostic_only() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]

    compare_final_states(
        (replace(state, counters=SolverCounters(142, 143, 144, 145, 146)),),
        _matching_unit_reference(),
    )


def test_achieved_time_must_reach_target_time() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    with pytest.raises(ComparisonFailure, match="achieved_time"):
        compare_final_states((replace(state, time=1.9),), _matching_unit_reference())


def test_characterized_achieved_time_does_not_compound_target_tolerance() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    reference = _matching_unit_reference()
    fields = {1: dict(reference.fields[1])}
    fields[1]["target_time"] = Tolerance(2.0, atol=0.1, rtol=0.0)
    fields[1]["achieved_time"] = Tolerance(2.0, atol=0.1, rtol=0.0)

    with pytest.raises(ComparisonFailure, match="achieved_time"):
        compare_final_states(
            (replace(state, target_time=2.075, time=2.15),),
            replace(reference, fields=fields),
        )


def test_composition_norms_are_diagnostic_only() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]
    perturbed_mass_fractions = dict(state.mass_fractions)
    perturbed_mass_fractions["c12"] += 1e-4
    perturbed_mass_fractions["o16"] -= 1e-4
    perturbed = replace(state, mass_fractions=perturbed_mass_fractions)
    reference = _matching_unit_reference()

    diagnostics = calculate_composition_norms((perturbed,), reference)
    assert diagnostics == (
        CompositionNorms(
            zone=1,
            l1=pytest.approx(2e-4),
            l2=pytest.approx(2**0.5 * 1e-4),
            linf=pytest.approx(1e-4),
            linf_species="c12",
        ),
    )
    compare_final_states((perturbed,), reference)


def test_missing_final_output_is_a_parsing_failure() -> None:
    with pytest.raises(ParsingFailure, match="unexpected diagnostic structure"):
        parse_diagnostic(
            "Timers Summary:\n        Total 1.0E-02\n", (1,), ALPHA_SPECIES
        )


def test_malformed_output_is_a_parsing_failure() -> None:
    diagnostic = _fabricated_final_diagnostic().replace(
        "End     1    42", "End     1    malformed"
    )
    with pytest.raises(ParsingFailure, match="malformed End record"):
        parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)


def test_counter_zone_must_agree_with_end_zone() -> None:
    diagnostic = _fabricated_final_diagnostic().replace(
        "1        42        42", "2        42        42"
    )

    with pytest.raises(ParsingFailure, match="does not match its End record"):
        parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)


def test_grouped_parser_requires_declared_batch_structure() -> None:
    case = batch_alpha_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    states = _reference_states(reference)
    diagnostic = _diagnostic_for_states(states, case.expected_diagnostic_groups)
    assert tuple(
        state.zone
        for state in parse_diagnostic(
            diagnostic,
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    ) == case.expected_zones

    one_zone_groups = tuple((zone,) for zone in case.expected_zones)
    with pytest.raises(ParsingFailure, match="missing End record for zone 2"):
        parse_diagnostic(
            _diagnostic_for_states(states, one_zone_groups),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="expected End zone 5, found 6"):
        parse_diagnostic(
            _diagnostic_for_states(
                states, ((1, 2, 3, 4), (6, 5, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16))
            ),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="counter record does not match"):
        parse_diagnostic(
            diagnostic.replace(
                f"\n1 {states[0].counters.ts}",
                f"\n2 {states[0].counters.ts}",
                1,
            ),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="missing timer section for group 4"):
        parse_diagnostic(
            diagnostic.rsplit("Timers Summary:", 1)[0],
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="unexpected diagnostic structure before group 2"):
        parse_diagnostic(
            diagnostic.replace(
                "\nEnd 5 ",
                "\nCounters:  Zone        TS        NR  Jacobian     Deriv CrossSect\n"
                "99 1 1 1 1 1\nEnd 5 ",
                1,
            ),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="duplicate timer row Total in group 1"):
        parse_diagnostic(
            diagnostic.replace(
                "\nEnd 5 ",
                "\n        Total      1.000E-02\nEnd 5 ",
                1,
            ),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )
    with pytest.raises(ParsingFailure, match="unexpected diagnostic structure before group 2"):
        parse_diagnostic(
            diagnostic.replace(
                "\nEnd 5 ",
                "\nTimers Summary:\n        Total      1.000E-02\nEnd 5 ",
                1,
            ),
            case.expected_zones,
            case.expected_species,
            case.expected_diagnostic_groups,
        )


def test_grouped_parser_accepts_and_rejects_short_final_group() -> None:
    case = batch_alpha_case(REPOSITORY_ROOT)
    states = _reference_states(load_reference(case.reference))
    groups = ((1, 2, 3, 4, 5, 6), (7, 8, 9, 10, 11, 12), (13, 14, 15, 16))
    diagnostic = _diagnostic_for_states(states, groups)
    parse_diagnostic(diagnostic, case.expected_zones, case.expected_species, groups)

    mutations = (
        (_diagnostic_for_states(states, (*groups[:2], (13, 14, 15))), "missing End record for zone 16"),
        (_diagnostic_for_states(states, (*groups[:2], (13, 14, 15, 16, 1, 2))), "missing Counters record after group 3"),
        (diagnostic.replace("End 16 ", "End 17 "), "expected End zone 16, found 17"),
        ("\nTimers Summary:".join(diagnostic.rsplit("\nTimers Summary:", 1)[:-1]) + "\n16 1 1 1 2 2\nTimers Summary:" + diagnostic.rsplit("\nTimers Summary:", 1)[-1], "missing timer section"),
        (_diagnostic_for_states(states, (*groups[:2], (13, 14, 16, 15))), "expected End zone 15, found 16"),
        (diagnostic + "Timers Summary:\n        Total      1.000E-02\n", "unexpected extra diagnostic structure"),
    )
    for mutated, message in mutations:
        with pytest.raises(ParsingFailure, match=message):
            parse_diagnostic(mutated, case.expected_zones, case.expected_species, groups)


@pytest.mark.parametrize(
    "mutation",
    (
        "missing_group",
        "wrong_order",
        "wrong_membership",
        "cross_counter",
        "malformed_final",
        "singleton_groups",
    ),
)
def test_batch_alpha_registered_path_rejects_group_mutations(
    tmp_path: Path, mutation: str
) -> None:
    """Run structural failures through the registered external-process test."""

    case = batch_alpha_case(REPOSITORY_ROOT)
    states = _reference_states(load_reference(case.reference))
    groups = case.expected_diagnostic_groups
    if mutation == "missing_group":
        diagnostic = _diagnostic_for_states(states, groups[:3])
    elif mutation == "wrong_order":
        diagnostic = _diagnostic_for_states(
            states, (groups[0], (6, 5, 7, 8), groups[2], groups[3])
        )
    elif mutation == "wrong_membership":
        diagnostic = _diagnostic_for_states(
            states, (groups[0], (5, 6, 7, 9), groups[2], groups[3])
        )
    elif mutation == "cross_counter":
        diagnostic = _diagnostic_for_states(states, groups).replace(
            f"\n1 {states[0].counters.ts}",
            f"\n2 {states[0].counters.ts}",
            1,
        )
    elif mutation == "malformed_final":
        diagnostic = _diagnostic_for_states(states, groups).rsplit(
            "Timers Summary:", 1
        )[0]
    else:
        diagnostic = _diagnostic_for_states(
            states, tuple((zone,) for zone in case.expected_zones)
        )
    completed = _run_registered_pytest(
        tmp_path, case, _registered_diagnostic_executable(tmp_path, case, diagnostic)
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "parsing failure" in output


def test_parser_accepts_distinct_end_step_and_ts_counter() -> None:
    diagnostic = _fabricated_final_diagnostic().replace(
        "1        42        42", "1        43        42"
    )

    state = parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)[0]

    assert state.step == 42
    assert state.counters.ts == 43


def test_parser_preserves_source_labeled_solver_counters() -> None:
    state = parse_diagnostic(
        _fabricated_final_diagnostic(), (1,), ALPHA_SPECIES
    )[0]

    assert state.counters == SolverCounters(
        ts=42,
        nr=42,
        jacobian=42,
        derivative=43,
        cross_section=43,
    )


def test_parser_requires_exact_counter_heading() -> None:
    diagnostic = _fabricated_final_diagnostic().replace("CrossSect", "Rates")
    with pytest.raises(ParsingFailure, match="missing Counters record"):
        parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)


@pytest.mark.parametrize(
    "replacement",
    [
        "              1        42        42        42        43",
        "              1        42        42        42        43 malformed",
        "              1        42        42        42        -1        43",
    ],
)
def test_parser_rejects_malformed_counter_rows(replacement: str) -> None:
    diagnostic = _fabricated_final_diagnostic().replace(
        "              1        42        42        42        43        43",
        replacement,
    )

    with pytest.raises(ParsingFailure, match="malformed counter values"):
        parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)


@pytest.mark.parametrize("token", ["NaN", "+Inf", "-Infinity"])
def test_nonfinite_output_is_rejected(token: str) -> None:
    diagnostic = _fabricated_final_diagnostic().replace("4.0000000E+06", token)
    with pytest.raises(ParsingFailure, match="non-finite"):
        parse_diagnostic(diagnostic, (1,), ALPHA_SPECIES)


def test_timer_only_variation_is_excluded() -> None:
    first = parse_diagnostic(
        _fabricated_final_diagnostic(timer_total="1.000E-02"), (1,), ALPHA_SPECIES
    )
    second = parse_diagnostic(
        _fabricated_final_diagnostic(timer_total="9.999E+02"), (1,), ALPHA_SPECIES
    )
    assert first == second


def test_timer_exclusion_stops_before_unrecognized_content() -> None:
    normalized, count = strip_timer_sections(
        "Timers Summary:\n        Total 1.0E-02\nUnexpected diagnostic row\n"
    )
    assert count == 1
    assert normalized == "Unexpected diagnostic row"


def test_missing_executable_is_a_setup_failure(tmp_path: Path) -> None:
    with pytest.raises(SetupFailure, match="does not exist"):
        validate_executable(tmp_path / "missing-xnet")


def test_nonexecutable_file_is_a_setup_failure(tmp_path: Path) -> None:
    executable = tmp_path / "xnet"
    executable.write_text("not executable\n", encoding="utf-8")

    with pytest.raises(SetupFailure, match="not executable"):
        validate_executable(executable)


def test_nonzero_process_exit_is_an_execution_failure(tmp_path: Path) -> None:
    executable = _make_executable(
        tmp_path,
        "import sys\nprint('captured stdout')\n"
        "print('captured stderr', file=sys.stderr)\nraise SystemExit(7)",
    )
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    with pytest.raises(ExecutionFailure, match="nonzero status 7"):
        run_xnet(
            executable,
            _fake_case(tmp_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert (work_directory / "xnet.status.txt").read_text() == "return_code=7\n"
    assert (work_directory / "xnet.stdout.txt").read_text() == "captured stdout\n"
    assert (work_directory / "xnet.stderr.txt").read_text() == "captured stderr\n"


def test_signal_termination_is_an_execution_failure(tmp_path: Path) -> None:
    executable = _make_executable(
        tmp_path,
        "import os, signal\nos.kill(os.getpid(), signal.SIGTERM)",
    )
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    with pytest.raises(ExecutionFailure, match=f"signal {signal.SIGTERM}"):
        run_xnet(
            executable,
            _fake_case(tmp_path),
            work_directory,
            timeout_seconds=2.0,
        )


def test_timeout_is_an_execution_failure(tmp_path: Path) -> None:
    executable = _make_executable(tmp_path, "import time\ntime.sleep(2)")
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    with pytest.raises(ExecutionFailure, match="timed out"):
        run_xnet(
            executable,
            _fake_case(tmp_path),
            work_directory,
            timeout_seconds=0.05,
        )
    assert (work_directory / "xnet.status.txt").read_text() == "timeout\n"


def test_zero_exit_without_required_output_is_an_execution_failure(
    tmp_path: Path,
) -> None:
    executable = _make_executable(tmp_path, "raise SystemExit(0)")
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    with pytest.raises(ExecutionFailure, match="required fresh output"):
        run_xnet(
            executable,
            _fake_case(tmp_path),
            work_directory,
            timeout_seconds=2.0,
        )


@pytest.mark.parametrize(
    ("target", "empty"),
    [
        ("net_diag01", False),
        ("ev_fake_1", False),
        ("ts_fake_1", False),
        ("net_diag01", True),
        ("ev_fake_1", True),
        ("ts_fake_1", True),
    ],
)
def test_each_required_output_must_be_present_and_nonempty(
    tmp_path: Path, target: str, empty: bool
) -> None:
    required_outputs = ("net_diag01", "ev_fake_1", "ts_fake_1")
    writes = ["from pathlib import Path"]
    for filename in required_outputs:
        if filename == target and not empty:
            continue
        content = "b''" if filename == target else "b'fresh'"
        writes.append(f"Path({filename!r}).write_bytes({content})")
    executable = _make_executable(tmp_path, "\n".join(writes))
    work_directory = tmp_path / "work"
    work_directory.mkdir()

    with pytest.raises(ExecutionFailure, match=target):
        run_xnet(
            executable,
            _fake_case(tmp_path),
            work_directory,
            timeout_seconds=2.0,
        )


def test_nonempty_work_directory_is_a_setup_failure(tmp_path: Path) -> None:
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    (work_directory / "net_diag01").write_text("stale", encoding="utf-8")
    with pytest.raises(SetupFailure, match="refusing stale artifacts"):
        prepare_work_directory(tnsn_alpha_case(REPOSITORY_ROOT), work_directory)


def test_network_preprocessing_outputs_stay_in_work_directory(tmp_path: Path) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")
    local_network = work_directory / "Data_alpha"
    assert local_network.is_dir()
    assert not local_network.is_symlink()
    assert {path.name for path in local_network.iterdir()} == set(case.network_inputs)
    assert all((local_network / name).is_symlink() for name in case.network_inputs)

    generated = local_network / "nets3"
    generated.write_bytes(b"generated locally")
    assert generated.is_file()


def test_heat_alpha_stages_all_trajectories_and_expected_outputs(
    tmp_path: Path,
) -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")

    assert len(case.trajectories) == 6
    assert all(
        (work_directory / trajectory.name).is_symlink()
        for trajectory in case.trajectories
    )
    assert case.required_outputs == (
        "net_diag01",
        "ev_heat_alpha_1",
        "ts_heat_alpha_1",
        "ev_heat_alpha_2",
        "ts_heat_alpha_2",
        "ev_heat_alpha_3",
        "ts_heat_alpha_3",
        "ev_heat_alpha_4",
        "ts_heat_alpha_4",
        "ev_heat_alpha_5",
        "ts_heat_alpha_5",
        "ev_heat_alpha_6",
        "ts_heat_alpha_6",
    )


def test_torch47_definition_is_case_driven_and_stages_only_source_inputs(
    tmp_path: Path,
) -> None:
    case = tnsn_torch47_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")
    local_network = work_directory / "Data_torch47"

    sunet_species = tuple(
        line.strip().lower()
        for line in (case.network_data / "sunet")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    )

    assert len(case.expected_species) == 47
    assert case.expected_species == TORCH47_SPECIES == sunet_species
    assert len(set(case.expected_species)) == len(case.expected_species)
    assert {path.name for path in local_network.iterdir()} == set(
        case.network_inputs
    )
    assert all((local_network / name).is_symlink() for name in case.network_inputs)
    reference = load_reference(case.reference)
    assert comparison_species_for_zone(
        case.expected_species, reference.mass_fractions[1]
    ) == (
        "si28",
        "s31",
        "s32",
        "ar36",
        "ca40",
        "ti44",
        "cr48",
        "fe52",
        "co55",
        "ni56",
    )
    assert case.required_outputs == (
        "net_diag01",
        "ev_tnsn_torch47_1",
        "ts_tnsn_torch47_1",
    )


def test_sn160_definition_is_complete_and_stages_only_source_inputs(
    tmp_path: Path,
) -> None:
    case = heat_sn160_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")
    local_network = work_directory / "Data_SN160"
    sunet_species = tuple(
        line.strip().lower()
        for line in (case.network_data / "sunet")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    )

    assert len(case.expected_species) == 160
    assert case.expected_species == SN160_SPECIES == sunet_species
    assert len(set(case.expected_species)) == len(case.expected_species)
    assert {path.name for path in local_network.iterdir()} == set(
        case.network_inputs
    )
    assert all((local_network / name).is_symlink() for name in case.network_inputs)
    assert case.required_outputs == (
        "net_diag01",
        "ev_heat_sn160_1",
        "ts_heat_sn160_1",
        "ev_heat_sn160_2",
        "ts_heat_sn160_2",
        "ev_heat_sn160_3",
        "ts_heat_sn160_3",
        "ev_heat_sn160_4",
        "ts_heat_sn160_4",
        "ev_heat_sn160_5",
        "ts_heat_sn160_5",
        "ev_heat_sn160_6",
        "ts_heat_sn160_6",
    )


def test_bdf_sn160_definition_reuses_isolated_sn160_staging(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")
    local_network = work_directory / "Data_SN160"

    assert case.expected_species == SN160_SPECIES
    assert {path.name for path in local_network.iterdir()} == set(
        case.network_inputs
    )
    assert all((local_network / name).is_symlink() for name in case.network_inputs)
    assert case.required_outputs == (
        "net_diag01",
        "ev_bdf_sn160_1",
        "ts_bdf_sn160_1",
        "ev_bdf_sn160_2",
        "ts_bdf_sn160_2",
        "ev_bdf_sn160_3",
        "ts_bdf_sn160_3",
        "ev_bdf_sn160_4",
        "ts_bdf_sn160_4",
        "ev_bdf_sn160_5",
        "ts_bdf_sn160_5",
        "ev_bdf_sn160_6",
        "ts_bdf_sn160_6",
    )


def test_bdf_control_is_the_normalized_legacy_id_54_concatenation() -> None:
    settings = (REPOSITORY_ROOT / "test/test_settings_bdf").read_text(
        encoding="utf-8"
    )
    setup = (REPOSITORY_ROOT / "test/Test_Problems/setup_bdf_sn160").read_text(
        encoding="utf-8"
    )
    normalized = "\n".join(
        line.rstrip() for line in (settings + setup).splitlines()
    ) + "\n"
    normalized = normalized.replace("Test_Results/", "")
    normalized = normalized.replace("Test_Problems/", "")
    normalized = normalized.replace(
        "4         Blocking size for zone loop",
        "1         Blocking size for zone loop",
    )

    assert bdf_sn160_case(REPOSITORY_ROOT).control.read_text(
        encoding="utf-8"
    ) == normalized


def test_bdf_reference_records_end_steps_and_all_solver_counter_fields() -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)

    assert tuple(reference.final_steps) == case.expected_zones
    assert reference.solver_counters is not None
    assert tuple(reference.solver_counters) == case.expected_zones
    assert any(
        reference.solver_counters[zone].ts != reference.final_steps[zone]
        for zone in case.expected_zones
    )


def test_bdf_reference_rejects_incomplete_solver_counters(tmp_path: Path) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["solver_counters"].pop("cross_section")
    reference_path = tmp_path / "incomplete-bdf-counters.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="solver_counters must contain exactly"):
        load_reference(reference_path)


def test_bdf_run_rejects_missing_solver_counters_before_execution(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document.pop("solver_counters")
    reference_path = tmp_path / "missing-bdf-counters.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="solver_counters is required"):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


@pytest.mark.parametrize(
    "metadata_name",
    ["build", "legacy_provenance", "input_sha256", "python"],
)
def test_bdf_reference_rejects_missing_required_metadata(
    tmp_path: Path, metadata_name: str
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document.pop(metadata_name)
    reference_path = tmp_path / f"missing-{metadata_name}.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(SetupFailure, match="reference"):
        load_reference(reference_path)


def test_bdf_run_rejects_stale_input_hash_before_execution(tmp_path: Path) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    control_label = "test/regression/cases/bdf_sn160/control"
    document["input_sha256"][control_label] = "0" * 64
    reference_path = tmp_path / "stale-input-hash.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="reference input hash does not match"):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


def test_bdf_run_rejects_negative_reference_abundance_before_execution(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["mass_fractions"]["n"]["1"] = -1.0e-30
    reference_path = tmp_path / "negative-reference-abundance.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="mass_fractions contains negative"):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


def test_bdf_run_rejects_missing_reference_schema_before_execution(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document.pop("reference_schema")
    reference_path = tmp_path / "missing-reference-schema.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="reference schema does not match"):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


def test_bdf_run_rejects_mismatched_solver_provenance_before_execution(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["legacy_provenance"]["maintained_solver"] = "Bader-Deuflhard"
    reference_path = tmp_path / "wrong-solver-provenance.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="solver provenance does not match"):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


@pytest.mark.parametrize(
    ("metadata_path", "replacement"),
    [
        (
            ("baseline_status",),
            "characterization-only; independently validated scientific truth",
        ),
        (("generated_from_revision",), "0" * 40),
        (("generated_on",), "2026-08-04"),
        (("platform",), "contradictory platform"),
        (("compiler",), "contradictory compiler"),
        (("python",), "Python 0.0.0"),
        (("pytest",), "pytest 0.0.0"),
        (("build", "CMODE"), "DEBUG"),
        (
            ("legacy_provenance", "assembly"),
            ["test/test_settings_bdf", "test/Test_Problems/setup_heat_sn160"],
        ),
        (
            ("legacy_provenance", "normalized_changes"),
            ["contradictory normalization claim"],
        ),
    ],
)
def test_bdf_run_rejects_contradictory_characterization_metadata_before_execution(
    tmp_path: Path, metadata_path: tuple[str, ...], replacement: object
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    target = document
    for name in metadata_path[:-1]:
        target = target[name]
    target[metadata_path[-1]] = replacement
    reference_path = tmp_path / "contradictory-metadata.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")
    executable = _make_executable(tmp_path, "raise SystemExit(99)")
    work_directory = tmp_path / "work"

    with pytest.raises(
        SetupFailure,
        match="characterization metadata does not match the case definition",
    ):
        run_and_compare(
            executable,
            replace(case, reference=reference_path),
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


def test_case_rejects_sunet_species_mismatch(tmp_path: Path) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    network_data = tmp_path / "Data_alpha"
    network_data.mkdir()
    for filename in case.network_inputs:
        source = case.network_data / filename
        destination = network_data / filename
        if filename == "sunet":
            destination.write_text("p\nhe4\n", encoding="utf-8")
        else:
            destination.symlink_to(source.resolve())

    with pytest.raises(SetupFailure, match="network species input does not match"):
        prepare_work_directory(
            replace(case, network_data=network_data), tmp_path / "work"
        )


def test_case_rejects_duplicate_sunet_species(tmp_path: Path) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    network_data = tmp_path / "Data_alpha"
    network_data.mkdir()
    for filename in case.network_inputs:
        source = case.network_data / filename
        destination = network_data / filename
        if filename == "sunet":
            destination.write_text("he4\nhe4\n", encoding="utf-8")
        else:
            destination.symlink_to(source.resolve())

    with pytest.raises(SetupFailure, match="empty or contains duplicates"):
        prepare_work_directory(
            replace(case, network_data=network_data), tmp_path / "work"
        )


def test_sn160_parser_rejects_incomplete_and_duplicate_species() -> None:
    incomplete = _fabricated_diagnostic_with_species(SN160_SPECIES[:-1])
    with pytest.raises(ParsingFailure, match="incomplete abundance record"):
        parse_diagnostic(incomplete, (1,), SN160_SPECIES)

    duplicated_species = (SN160_SPECIES[0], *SN160_SPECIES[:-1])
    duplicated = _fabricated_diagnostic_with_species(duplicated_species)
    with pytest.raises(ParsingFailure, match="duplicate abundance"):
        parse_diagnostic(duplicated, (1,), SN160_SPECIES)


def test_torch47_parser_rejects_wrong_species_order() -> None:
    diagnostic = _fabricated_diagnostic_with_species(TORCH47_SPECIES)
    states = parse_diagnostic(diagnostic, (1,), TORCH47_SPECIES)
    assert tuple(states[0].mass_fractions) == TORCH47_SPECIES

    swapped_species = list(TORCH47_SPECIES)
    swapped_species[0], swapped_species[1] = swapped_species[1], swapped_species[0]
    swapped_diagnostic = _fabricated_diagnostic_with_species(tuple(swapped_species))
    with pytest.raises(ParsingFailure, match="unexpected species structure"):
        parse_diagnostic(swapped_diagnostic, (1,), TORCH47_SPECIES)


def test_torch47_reference_rejects_wrong_species_order(tmp_path: Path) -> None:
    case = tnsn_torch47_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    composition_items = list(document["mass_fractions"].items())
    composition_items[0], composition_items[1] = (
        composition_items[1],
        composition_items[0],
    )
    document["mass_fractions"] = dict(composition_items)
    reference_path = tmp_path / "reordered-reference.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    reference = load_reference(reference_path)
    with pytest.raises(SetupFailure, match="composition reference does not match"):
        validate_reference_for_case(case, reference)


def test_reference_rejects_wrong_case_association(tmp_path: Path) -> None:
    case = heat_sn160_case(REPOSITORY_ROOT)
    document = json.loads(case.reference.read_text(encoding="utf-8"))
    document["case"] = "bdf_sn160"
    reference_path = tmp_path / "wrong-case-reference.json"
    reference_path.write_text(json.dumps(document), encoding="utf-8")

    reference = load_reference(reference_path)
    with pytest.raises(SetupFailure, match="reference case does not match"):
        validate_reference_for_case(case, reference)


@pytest.mark.parametrize(
    ("case_factory", "wrong_reference_factory"),
    [
        (heat_sn160_case, bdf_sn160_case),
        (bdf_sn160_case, heat_sn160_case),
    ],
)
def test_sn160_references_cannot_be_swapped(
    tmp_path: Path, case_factory, wrong_reference_factory
) -> None:
    case = replace(
        case_factory(REPOSITORY_ROOT),
        reference=wrong_reference_factory(REPOSITORY_ROOT).reference,
    )
    executable = _make_executable(
        tmp_path,
        "from pathlib import Path\nPath('executed').write_text('unexpected')",
    )
    work_directory = tmp_path / "work"

    with pytest.raises(SetupFailure, match="reference case does not match"):
        run_and_compare(
            executable,
            case,
            work_directory,
            timeout_seconds=2.0,
        )
    assert not work_directory.exists()


def test_reference_rejects_duplicate_tolerance_declaration(tmp_path: Path) -> None:
    case = tnsn_torch47_case(REPOSITORY_ROOT)
    source = case.reference.read_text(encoding="utf-8")
    tolerance_start = source.index('"mass_fraction_tolerances": {')
    declaration_start = source.index('\n    "all_selected":', tolerance_start)
    duplicate = '\n    "all_selected": {"atol": 1.0, "rtol": 0.0, "exact": false},'
    reference_path = tmp_path / "duplicate-tolerance-reference.json"
    reference_path.write_text(
        source[:declaration_start] + duplicate + source[declaration_start:],
        encoding="utf-8",
    )

    with pytest.raises(SetupFailure, match="duplicate JSON object key: all_selected"):
        load_reference(reference_path)


def test_all_migrated_references_use_per_zone_comparison_species_policy() -> None:
    for case in (
        tnsn_alpha_case(REPOSITORY_ROOT),
        heat_alpha_case(REPOSITORY_ROOT),
        tnsn_torch47_case(REPOSITORY_ROOT),
        heat_sn160_case(REPOSITORY_ROOT),
        bdf_sn160_case(REPOSITORY_ROOT),
    ):
        reference = load_reference(case.reference)
        for zone in case.expected_zones:
            expected_selection = comparison_species_for_zone(
                case.expected_species, reference.mass_fractions[zone]
            )
            assert (
                tuple(reference.mass_fraction_tolerances[zone])
                == expected_selection
            )


def test_material_species_policy_does_not_require_silicon_products() -> None:
    cno_species = ("p", "he4", "c12", "n14", "o16", "ne20")
    mass_fractions = {
        "p": 1.0e-6,
        "he4": 5.0e-5,
        "c12": 0.2,
        "n14": 0.3,
        "o16": 0.4,
        "ne20": 0.099949,
    }

    assert comparison_species_for_zone(cno_species, mass_fractions) == (
        "c12",
        "n14",
        "o16",
        "ne20",
    )


def test_per_zone_material_threshold_is_inclusive_and_ordered() -> None:
    species = ("he4", "si28", "c12", "o16")

    assert comparison_species_for_zone(
        species,
        {"he4": 0.2, "si28": 0.0, "c12": 1.0e-4, "o16": 0.7999},
    ) == ("he4", "si28", "c12", "o16")


def test_comparison_anchors_are_retained_per_zone_when_present() -> None:
    assert comparison_species_for_zone(
        ("p", "si28"), {"p": 9.0e-5, "si28": 0.0}
    ) == ("si28",)


def test_heat_alpha_selected_zone_species_perturbation_fails() -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    state = _state_from_reference(reference, 4)
    mass_fractions = dict(state.mass_fractions)
    mass_fractions["o16"] += 1.0e-6
    mass_fractions["c12"] -= 1.0e-6

    zone_reference = replace(reference, expected_zones=(4,))
    with pytest.raises(ComparisonFailure, match="zone 4 o16 mass fraction"):
        compare_final_states(
            (replace(state, mass_fractions=mass_fractions),), zone_reference
        )


def test_complete_vector_limits_detect_unselected_zone_change() -> None:
    case = heat_alpha_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    state = _state_from_reference(reference, 3)
    mass_fractions = dict(state.mass_fractions)
    mass_fractions["o16"] += 1.0e-6
    mass_fractions["c12"] -= 1.0e-6
    perturbed = replace(state, mass_fractions=mass_fractions)

    assert "o16" not in reference.mass_fraction_tolerances[3]
    assert "o16" in reference.mass_fraction_tolerances[4]
    assert min(mass_fractions.values()) >= 0.0
    assert sum(mass_fractions.values()) == pytest.approx(
        sum(state.mass_fractions.values()), abs=1.0e-15
    )
    assert calculate_composition_norms((perturbed,), reference)[0].linf == pytest.approx(
        1.0e-6
    )
    with pytest.raises(ComparisonFailure, match="zone 3 L1"):
        compare_final_states((perturbed,), replace(reference, expected_zones=(3,)))


def test_diagnostic_only_species_still_receive_negative_fraction_check() -> None:
    reference = load_reference(heat_alpha_case(REPOSITORY_ROOT).reference)
    state = _state_from_reference(reference, 3)
    mass_fractions = dict(state.mass_fractions)
    transfer = mass_fractions["c12"] + 1.0e-6
    mass_fractions["c12"] -= transfer
    mass_fractions["o16"] += transfer

    zone_reference = replace(reference, expected_zones=(3,))
    with pytest.raises(ComparisonFailure, match="negative mass fractions: c12"):
        compare_final_states(
            (replace(state, mass_fractions=mass_fractions),), zone_reference
        )


def test_diagnostic_only_species_still_receive_composition_sum_check() -> None:
    reference = load_reference(heat_alpha_case(REPOSITORY_ROOT).reference)
    state = _state_from_reference(reference, 3)
    mass_fractions = dict(state.mass_fractions)
    mass_fractions["c12"] += 1.0e-6

    zone_reference = replace(reference, expected_zones=(3,))
    with pytest.raises(ComparisonFailure, match="mass-fraction sum"):
        compare_final_states(
            (replace(state, mass_fractions=mass_fractions),), zone_reference
        )


def test_printed_sum_gate_rejects_change_hidden_by_normalization_allowance() -> None:
    reference = load_reference(bdf_sn160_case(REPOSITORY_ROOT).reference)
    state = _state_from_reference(reference, 3)
    mass_fractions = dict(state.mass_fractions)
    for species in ("fe57", "ca42", "ti47"):
        mass_fractions[species] -= 1.5e-5

    norms = calculate_composition_norms(
        (replace(state, mass_fractions=mass_fractions),), reference
    )[0]
    assert norms.l1 <= reference.composition_norm_limits[3].l1  # type: ignore[index]
    assert norms.linf <= reference.composition_norm_limits[3].linf  # type: ignore[index]
    assert all(
        mass_fractions[species] == state.mass_fractions[species]
        for species in reference.mass_fraction_tolerances[3]
    )
    with pytest.raises(ComparisonFailure, match="printed mass-fraction sum"):
        compare_final_states(
            (replace(state, mass_fractions=mass_fractions),),
            replace(reference, expected_zones=(3,)),
        )


def test_torch47_network_specific_products_gate_comparison() -> None:
    case = tnsn_torch47_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    fields = reference.fields[1]
    mass_fractions = dict(reference.mass_fractions[1])
    mass_fractions["co55"] += 1.0e-4
    mass_fractions["s31"] -= 1.0e-4
    state = FinalState(
        zone=1,
        step=reference.final_steps[1],
        target_time=fields["target_time"].value,
        time=fields["target_time"].value,
        temperature_gk=fields["temperature_gk"].value,
        density=fields["density"].value,
        electron_fraction=fields["electron_fraction"].value,
        mass_fractions=mass_fractions,
        counters=SolverCounters(reference.final_steps[1], 0, 0, 0, 0),
    )

    with pytest.raises(ComparisonFailure, match="co55 mass fraction"):
        compare_final_states((state,), reference)


def test_duplicate_trajectory_basenames_are_a_setup_failure(tmp_path: Path) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    duplicated = replace(case, trajectories=(case.trajectories[0],) * 2)

    with pytest.raises(SetupFailure, match="trajectory definition"):
        prepare_work_directory(duplicated, tmp_path / "work")


def test_batch_alpha_stages_nested_prefix_inputs(tmp_path: Path) -> None:
    case = batch_alpha_case(REPOSITORY_ROOT)
    work_directory = prepare_work_directory(case, tmp_path / "work")
    for item in case.staged_inputs:
        staged = work_directory / item.destination
        assert staged.is_symlink()
        assert staged.resolve() == item.source.resolve()


def test_batch_alpha_control_is_the_normalized_legacy_id_61_concatenation() -> None:
    settings = (REPOSITORY_ROOT / "test/test_settings_batch").read_text(
        encoding="utf-8"
    )
    setup = (REPOSITORY_ROOT / "test/Test_Problems/setup_batch_alpha").read_text(
        encoding="utf-8"
    )
    normalized = "\n".join(line.rstrip() for line in (settings + setup).splitlines()) + "\n"
    assert batch_alpha_case(REPOSITORY_ROOT).control.read_text(encoding="utf-8") == normalized.replace(
        "Test_Results/", ""
    )


@pytest.mark.parametrize(
    ("index", "replacement", "message"),
    (
        (0, StagedInput(Path("missing-abundance"), Path("Data_alpha/ab_batch/ab_batch_01")), "declared staged input 1"),
        (16, StagedInput(Path("missing-trajectory"), Path("Test_Problems/th_batch/th_batch_01")), "declared staged input 17"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("../escape")), "invalid staged destination"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("/unsafe")), "invalid staged destination"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("Data_alpha/ab_batch/ab_batch_01")), "duplicate staged destination"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("net_diag01")), "staged destination is reserved"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("net_diag01/poison")), "staged destination is reserved"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("xnet.stdout.txt")), "staged destination is reserved"),
        (0, StagedInput(REPOSITORY_ROOT / "test/Data_alpha/ab_batch/ab_batch_01", Path("Data_alpha")), "overlapping staged destinations"),
    ),
)
def test_batch_alpha_staging_rejects_missing_or_unsafe_inputs(
    tmp_path: Path, index: int, replacement: StagedInput, message: str
) -> None:
    case = batch_alpha_case(REPOSITORY_ROOT)
    staged_inputs = list(case.staged_inputs)
    if message in ("duplicate staged destination", "staged destination is reserved"):
        staged_inputs.append(replacement)
    else:
        staged_inputs[index] = replacement
    with pytest.raises(SetupFailure, match=message):
        prepare_work_directory(
            replace(case, staged_inputs=tuple(staged_inputs)),
            tmp_path / "work",
        )


def test_missing_case_input_is_a_setup_failure(tmp_path: Path) -> None:
    case = replace(tnsn_alpha_case(REPOSITORY_ROOT), control=tmp_path / "missing-control")
    with pytest.raises(SetupFailure, match="complete control input"):
        prepare_work_directory(case, tmp_path / "work")


def test_missing_sn160_network_directory_is_a_setup_failure(tmp_path: Path) -> None:
    case = replace(
        heat_sn160_case(REPOSITORY_ROOT),
        network_data=tmp_path / "missing-Data_SN160",
    )
    with pytest.raises(SetupFailure, match="network data directory"):
        prepare_work_directory(case, tmp_path / "work")


def test_missing_sn160_abundance_is_a_setup_failure(tmp_path: Path) -> None:
    case = heat_sn160_case(REPOSITORY_ROOT)
    network_data = tmp_path / "Data_SN160"
    network_data.mkdir()
    for filename in case.network_inputs:
        if filename != "ab_co":
            (network_data / filename).symlink_to(
                (case.network_data / filename).resolve()
            )

    with pytest.raises(SetupFailure, match="required network source input"):
        prepare_work_directory(
            replace(case, network_data=network_data), tmp_path / "work"
        )


def test_missing_sn160_trajectory_is_a_setup_failure(tmp_path: Path) -> None:
    case = heat_sn160_case(REPOSITORY_ROOT)
    trajectories = (tmp_path / "missing-trajectory", *case.trajectories[1:])

    with pytest.raises(SetupFailure, match="thermodynamic trajectory 1"):
        prepare_work_directory(
            replace(case, trajectories=trajectories), tmp_path / "work"
        )


def test_missing_sn160_eos_table_is_a_setup_failure(tmp_path: Path) -> None:
    case = replace(
        heat_sn160_case(REPOSITORY_ROOT),
        helm_table=tmp_path / "missing-helm_table.dat",
    )
    with pytest.raises(SetupFailure, match="Helmholtz EOS table"):
        prepare_work_directory(case, tmp_path / "work")


def test_failing_regression_makes_pytest_return_nonzero(tmp_path: Path) -> None:
    executable = _make_executable(tmp_path, "raise SystemExit(9)")
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "test/regression/test_regression.py",
            f"--xnet-executable={executable}",
        ],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
        env=environment,
    )
    assert completed.returncode != 0
    assert "execution failure" in completed.stdout + completed.stderr


CASE_FACTORIES = (
    tnsn_alpha_case,
    heat_alpha_case,
    tnsn_torch47_case,
    heat_sn160_case,
    bdf_sn160_case,
    nse_sn160_case,
    batch_alpha_case,
)


def _reference_states(reference: CharacterizationReference) -> tuple[FinalState, ...]:
    return tuple(_state_from_reference(reference, zone) for zone in reference.expected_zones)


def _move_mass_fraction(
    state: FinalState, donor: str, recipient: str, amount: float
) -> FinalState:
    mass_fractions = dict(state.mass_fractions)
    assert donor != recipient
    assert mass_fractions[donor] >= amount
    mass_fractions[donor] -= amount
    mass_fractions[recipient] += amount
    return replace(state, mass_fractions=mass_fractions)


def test_registered_num_obs_001_composition_observations_pass(
    tmp_path: Path,
) -> None:
    torch_case = tnsn_torch47_case(REPOSITORY_ROOT)
    torch_reference = load_reference(torch_case.reference)
    torch_states = list(_reference_states(torch_reference))
    torch_mass_fractions = dict(torch_states[0].mass_fractions)
    torch_mass_fractions["s31"] += 1.890e-9
    torch_mass_fractions["co55"] += 2.800e-9
    torch_mass_fractions["mn51"] -= 1.8275e-9
    torch_mass_fractions["p30"] += 0.2685e-9
    torch_states[0] = replace(
        torch_states[0], mass_fractions=torch_mass_fractions
    )
    torch_executable = _registered_case_executable(
        tmp_path, torch_case, tuple(torch_states)
    )
    torch_completed = _run_registered_pytest(
        tmp_path, torch_case, torch_executable
    )
    assert torch_completed.returncode == 0, (
        torch_completed.stdout + torch_completed.stderr
    )
    torch_artifact = next(
        (tmp_path / "artifacts").glob(
            "test_tnsn_torch47*/tnsn_torch47/composition_error_norms.json"
        )
    )
    torch_diagnostics = json.loads(torch_artifact.read_text(encoding="utf-8"))
    assert torch_diagnostics["zones"] == [
        {
            "zone": 1,
            "l1": pytest.approx(6.786e-9, abs=1.0e-15),
            "l2": pytest.approx(
                math.sqrt(
                    (1.890e-9) ** 2
                    + (2.800e-9) ** 2
                    + (1.8275e-9) ** 2
                    + (0.2685e-9) ** 2
                ),
                abs=1.0e-15,
            ),
            "linf": pytest.approx(2.800e-9, abs=1.0e-15),
            "linf_species": "co55",
            "l1_limit": 1.1e-8,
            "linf_limit": 4.3e-9,
        }
    ]

    bdf_case = bdf_sn160_case(REPOSITORY_ROOT)
    bdf_reference = load_reference(bdf_case.reference)
    bdf_states = list(_reference_states(bdf_reference))
    for state_index, actual_sum, species in (
        (0, 0.9999997928543, "fe56"),
        (2, 1.000030571079, "fe57"),
    ):
        state = bdf_states[state_index]
        mass_fractions = dict(state.mass_fractions)
        mass_fractions[species] += actual_sum - math.fsum(
            mass_fractions.values()
        )
        bdf_states[state_index] = replace(
            state,
            electron_fraction=(
                0.49928059
                if state.zone == 3
                else state.electron_fraction
            ),
            mass_fractions=mass_fractions,
        )
    bdf_executable = _registered_case_executable(
        tmp_path, bdf_case, tuple(bdf_states)
    )
    bdf_completed = _run_registered_pytest(
        tmp_path, bdf_case, bdf_executable
    )
    assert bdf_completed.returncode == 0, (
        bdf_completed.stdout + bdf_completed.stderr
    )


def test_registered_tnsn_alpha_sci_001_counterexample_fails_selected_and_vector_gates(
    tmp_path: Path,
) -> None:
    case = tnsn_alpha_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    states = []
    for state in _reference_states(reference):
        mass_fractions = dict(state.mass_fractions)
        mass_fractions["si28"] -= 1.0e-8
        mass_fractions["s32"] += 1.0e-8
        states.append(replace(state, mass_fractions=mass_fractions))

    executable = _registered_case_executable(tmp_path, case, tuple(states))
    completed = _run_registered_pytest(tmp_path, case, executable)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "si28 mass fraction" in output
    assert "s32 mass fraction" in output
    assert "L1" in output
    assert "L-infinity" in output


def test_registered_printed_sum_mutation_fails_while_other_gates_pass(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    states = list(_reference_states(reference))
    state = states[2]
    mass_fractions = dict(state.mass_fractions)
    for species in ("fe57", "ca42", "ti47"):
        mass_fractions[species] -= 1.5e-5
    states[2] = replace(state, mass_fractions=mass_fractions)

    executable = _registered_case_executable(tmp_path, case, tuple(states))
    completed = _run_registered_pytest(tmp_path, case, executable)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "printed mass-fraction sum" in output
    assert "recomputed mass-fraction normalization" not in output
    assert "L1" not in output
    assert "L-infinity" not in output


def test_registered_normalization_mutation_fails_structural_gate(
    tmp_path: Path,
) -> None:
    case = bdf_sn160_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    states = list(_reference_states(reference))
    state = states[2]
    mass_fractions = dict(state.mass_fractions)
    for species in ("fe57", "ca42", "ti47"):
        assert species not in reference.mass_fraction_tolerances[3]
        mass_fractions[species] += 4.0e-7
    states[2] = replace(state, mass_fractions=mass_fractions)

    executable = _registered_case_executable(tmp_path, case, tuple(states))
    completed = _run_registered_pytest(tmp_path, case, executable)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "printed mass-fraction sum" in output
    assert "recomputed mass-fraction normalization" in output
    assert "L1" not in output
    assert "L-infinity" not in output


def test_nse_case_rejects_missing_initialization_marker(tmp_path: Path) -> None:
    case = nse_sn160_case(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    executable = _registered_diagnostic_executable(
        tmp_path,
        case,
        _diagnostic_for_states(
            _reference_states(reference), case.expected_diagnostic_groups
        ),
    )

    completed = _run_registered_pytest(tmp_path, case, executable)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "parsing failure" in output
    assert "Initial abundances from NSE" in output


@pytest.mark.parametrize("case_factory", CASE_FACTORIES, ids=lambda item: item.__name__)
@pytest.mark.parametrize(
    ("mutation", "expected_diagnostic"),
    (
        ("temperature", "temperature_gk"),
        ("selected_species", "mass fraction"),
        ("localized_vector", "L-infinity"),
        ("broad_vector", "L1"),
    ),
)
def test_registered_cases_reject_controlled_policy_violations(
    tmp_path: Path,
    case_factory,
    mutation: str,
    expected_diagnostic: str,
) -> None:
    """Exercise each policy gate through the normal registered pytest entry point."""

    case = case_factory(REPOSITORY_ROOT)
    reference = load_reference(case.reference)
    states = list(_reference_states(reference))
    zone = reference.expected_zones[0]
    state_index = 0
    state = states[state_index]

    if mutation == "temperature":
        policy = reference.fields[zone]["temperature_gk"]
        amount = 1.0e-6 if policy.exact else 3.0 * policy.atol
        states[state_index] = replace(
            state, temperature_gk=state.temperature_gk + amount
        )
    elif mutation == "selected_species":
        selected = reference.mass_fraction_tolerances[zone]
        recipient = max(selected, key=state.mass_fractions.__getitem__)
        donor = max(
            (species for species in state.mass_fractions if species != recipient),
            key=state.mass_fractions.__getitem__,
        )
        bounds = selected[recipient]
        amount = (
            1.0e-6
            if bounds.exact
            else 3.0
            * (
                bounds.atol
                + bounds.rtol * abs(state.mass_fractions[recipient])
            )
        )
        states[state_index] = _move_mass_fraction(state, donor, recipient, amount)
    elif mutation == "localized_vector":
        assert reference.composition_norm_limits is not None
        limits = reference.composition_norm_limits[zone]
        assert limits.linf is not None
        donor = max(state.mass_fractions, key=state.mass_fractions.__getitem__)
        recipient = next(species for species in state.mass_fractions if species != donor)
        states[state_index] = _move_mass_fraction(
            state, donor, recipient, max(1.0e-6, 3.0 * limits.linf)
        )
    else:
        assert reference.composition_norm_limits is not None
        limits = reference.composition_norm_limits[zone]
        assert limits.l1 is not None
        ordered = sorted(
            state.mass_fractions, key=state.mass_fractions.__getitem__
        )
        donors = ordered[-4:]
        recipients = ordered[:4]
        amount = max(1.0e-6, limits.l1 / 2.0)
        changed = state
        for donor, recipient in zip(donors, recipients, strict=True):
            changed = _move_mass_fraction(changed, donor, recipient, amount)
        states[state_index] = changed

    executable = _registered_case_executable(tmp_path, case, tuple(states))
    completed = _run_registered_pytest(tmp_path, case, executable)
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0, output
    assert "comparison failure" in output
    assert expected_diagnostic in output


@pytest.mark.parametrize("case_factory", CASE_FACTORIES, ids=lambda item: item.__name__)
def test_registered_execution_cannot_modify_reference(
    tmp_path: Path, case_factory
) -> None:
    case = case_factory(REPOSITORY_ROOT)
    before = case.reference.read_bytes()
    reference = load_reference(case.reference)
    executable = _registered_case_executable(
        tmp_path, case, _reference_states(reference)
    )

    completed = _run_registered_pytest(tmp_path, case, executable)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert case.reference.read_bytes() == before
