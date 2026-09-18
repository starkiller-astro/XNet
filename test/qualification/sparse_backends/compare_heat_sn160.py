#!/usr/bin/env python3
"""Compare same-source dense and qualified sparse heat_sn160 endpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT / "test" / "regression"))

from xnet_regression import (  # noqa: E402
    calculate_composition_norms,
    compare_final_states,
    heat_sn160_case,
    load_reference,
    parse_diagnostic,
    prepare_work_directory,
    reanchor_reference_to_states,
    RegressionFailure,
    run_xnet,
    scale_comparison_tolerances,
    validate_reference_for_case,
)

CROSS_SOLVER_TOLERANCE_FACTOR = 2.0


def _run_endpoint(executable: Path, work_directory: Path, timeout: float):
    case = heat_sn160_case(REPOSITORY_ROOT)
    prepared = prepare_work_directory(case, work_directory)
    result = run_xnet(
        executable,
        case,
        prepared,
        timeout_seconds=timeout,
    )
    diagnostic = (prepared / "net_diag01").read_text(encoding="utf-8")
    states = parse_diagnostic(
        diagnostic,
        case.expected_zones,
        case.expected_species,
        case.expected_diagnostic_groups,
    )
    return result, states


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _executable_hashes(dense_executable: Path, sparse_executable: Path):
    dense_hash = _sha256(dense_executable)
    sparse_hash = _sha256(sparse_executable)
    if dense_hash == sparse_hash:
        raise RegressionFailure(
            "dense and sparse qualification executables have identical SHA-256 "
            f"hashes ({dense_hash}); distinct backend builds are required"
        )
    return dense_hash, sparse_hash


def _field_value(state, name: str) -> float:
    return state.time if name == "achieved_time" else getattr(state, name)


def _comparison_report(provider: str, dense_states, sparse_states, reference):
    dense_by_zone = {state.zone: state for state in dense_states}
    sparse_by_zone = {state.zone: state for state in sparse_states}
    norms = calculate_composition_norms(sparse_states, reference)
    zones = []
    for norm in norms:
        zone = norm.zone
        dense = dense_by_zone[zone]
        sparse = sparse_by_zone[zone]
        fields = {}
        for name, policy in reference.fields[zone].items():
            absolute_difference = abs(_field_value(sparse, name) - policy.value)
            allowed = 0.0 if policy.exact else policy.atol + policy.rtol * abs(policy.value)
            fields[name] = {
                "dense": policy.value,
                "sparse": _field_value(sparse, name),
                "absolute_difference": absolute_difference,
                "allowed": allowed,
            }

        selected = reference.mass_fraction_tolerances[zone]
        selected_differences = {
            species: abs(
                sparse.mass_fractions[species] - dense.mass_fractions[species]
            )
            for species in selected
        }
        maximum_species = max(selected_differences, key=selected_differences.__getitem__)
        zones.append(
            {
                "zone": zone,
                "dense_step": dense.step,
                "sparse_step": sparse.step,
                "fields": fields,
                "composition": {
                    "species_count": len(sparse.mass_fractions),
                    "selected_species_count": len(selected),
                    "l1": norm.l1,
                    "l2": norm.l2,
                    "linf": norm.linf,
                    "linf_species": norm.linf_species,
                    "maximum_selected_absolute_difference": selected_differences[
                        maximum_species
                    ],
                    "maximum_selected_species": maximum_species,
                    "dense_printed_sum": math.fsum(dense.mass_fractions.values()),
                    "sparse_printed_sum": math.fsum(sparse.mass_fractions.values()),
                },
            }
        )
    return {
        "schema": "xnet-sparse-qualification-v1",
        "provider": provider,
        "case": "heat_sn160",
        "comparison_policy": (
            "xnet-comparison-v1 tolerances reanchored to the same-source dense run; "
            f"non-exact difference bounds multiplied by {CROSS_SOLVER_TOLERANCE_FACTOR:g}"
        ),
        "zones": zones,
    }


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provider", required=True, choices=("ma48", "pardiso-mkl"))
    parser.add_argument("--dense-executable", required=True, type=Path)
    parser.add_argument("--sparse-executable", required=True, type=Path)
    parser.add_argument("--work-directory", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=180.0)
    return parser.parse_args()


def main() -> int:
    arguments = _arguments()
    dense_hash, sparse_hash = _executable_hashes(
        arguments.dense_executable,
        arguments.sparse_executable,
    )
    case = heat_sn160_case(REPOSITORY_ROOT)
    policy = load_reference(case.reference)
    validate_reference_for_case(case, policy)

    dense_result, dense_states = _run_endpoint(
        arguments.dense_executable,
        arguments.work_directory / "dense",
        arguments.timeout,
    )
    sparse_result, sparse_states = _run_endpoint(
        arguments.sparse_executable,
        arguments.work_directory / arguments.provider,
        arguments.timeout,
    )
    comparison_reference = reanchor_reference_to_states(
        policy,
        dense_states,
        case_name=f"heat_sn160 dense versus {arguments.provider}",
    )
    comparison_reference = scale_comparison_tolerances(
        comparison_reference,
        CROSS_SOLVER_TOLERANCE_FACTOR,
    )
    compare_final_states(sparse_states, comparison_reference)

    report = _comparison_report(
        arguments.provider,
        dense_states,
        sparse_states,
        comparison_reference,
    )
    report["executables"] = {
        "dense": {
            "path": str(dense_result.executable),
            "sha256": dense_hash,
            "return_code": dense_result.return_code,
        },
        "sparse": {
            "path": str(sparse_result.executable),
            "sha256": sparse_hash,
            "return_code": sparse_result.return_code,
        },
    }
    report_path = arguments.work_directory / "qualification-report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    maximum_linf = max(zone["composition"]["linf"] for zone in report["zones"])
    print(
        f"{arguments.provider} heat_sn160 agrees with the same-source dense run; "
        f"maximum full-composition L-infinity difference={maximum_linf:.3e}"
    )
    print(f"qualification report: {report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, UnicodeError, RegressionFailure) as error:
        print(f"sparse qualification failed: {error}", file=sys.stderr)
        raise SystemExit(1) from None
