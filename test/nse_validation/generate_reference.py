#!/usr/bin/env python3
"""Generate the frozen independent NSE reference outside ordinary tests."""

from __future__ import annotations

import argparse
import decimal
import hashlib
import itertools
import json
import platform
from decimal import Decimal, ROUND_CEILING, localcontext
from pathlib import Path
from typing import Any

from extract_inputs import REPOSITORY_ROOT, canonical_bytes, extract, sha256
from reference_solver import (
    State,
    build_species,
    composition_norms,
    evaluate,
    solve,
)


GENERATOR_VERSION = "xnet-independent-nse-reference-v3"
REFERENCE_PRECISION = 50
PRECISION_CHECKS = (35, 65)
REFERENCE_RESIDUAL_LIMIT = Decimal("1e-25")
ROUTE_DIFFERENCE_LIMIT = Decimal("1e-24")
XNET_SOLVER_RESIDUAL = Decimal("1e-8")
START_OFFSETS = (("0", "0"), ("0.1", "-0.1"), ("2", "-2"), ("-2", "2"))
STATES = (
    ("symmetric", "1e9", "7", "0.50"),
    ("neutron_rich", "1e9", "7", "0.45"),
    ("proton_rich_low_density", "1e7", "6.5", "0.55"),
)


def binary64_string(value: Decimal) -> str:
    return format(float(value), ".17e")


def round_up_significant(value: Decimal, digits: int) -> Decimal:
    if value <= 0:
        return Decimal(0)
    quantum = Decimal(1).scaleb(value.adjusted() - digits + 1)
    return value.quantize(quantum, rounding=ROUND_CEILING)


def floating_point_budget(species_count: int) -> dict[str, Any]:
    # The inspected XNet partition-factor and abundance expressions require
    # fewer than 32 scalar arithmetic/transcendental evaluations per species.
    # Treat each as one correctly-rounded binary64 operation, then include the
    # complete serial accumulation.  This is an explicit portability
    # assumption, not an observed XNet discrepancy.
    operation_count = 32 * species_count + species_count - 1
    unit_roundoff = Decimal(2) ** -53
    gamma = Decimal(operation_count) * unit_roundoff / (
        Decimal(1) - Decimal(operation_count) * unit_roundoff
    )
    accepted = round_up_significant(gamma, 1)
    return {
        "binary64_unit_roundoff": str(unit_roundoff),
        "operation_equivalent_count": operation_count,
        "gamma_n": str(gamma),
        "accepted_composition_budget": str(accepted),
        "derivation": (
            "gamma_n = n*u/(1-n*u), with an inspected upper bound of 32 "
            "correctly-rounded operation-equivalents per species for the "
            "partition-factor and abundance expressions, plus a serial "
            "489-term accumulation; the result is rounded upward to one "
            "significant digit."
        ),
    }


def jacobian_condition(
    manifest: dict[str, Any], state: State, eta_n: Decimal, eta_p: Decimal
) -> Decimal:
    species = build_species(manifest, state, REFERENCE_PRECISION)
    jacobian = evaluate(species, eta_n, eta_p, state.ye).jacobian
    a, b = jacobian[0]
    c, d = jacobian[1]
    determinant = a * d - b * c
    inverse = ((d / determinant, -b / determinant), (-c / determinant, a / determinant))
    norm = max(abs(a) + abs(b), abs(c) + abs(d))
    inverse_norm = max(
        abs(inverse[0][0]) + abs(inverse[0][1]),
        abs(inverse[1][0]) + abs(inverse[1][1]),
    )
    return norm * inverse_norm


def solve_from_starts(
    manifest: dict[str, Any], state: State, route: str
) -> tuple[Any, list[dict[str, Any]]]:
    solutions = []
    diagnostics = []
    for offset in START_OFFSETS:
        solution = solve(
            manifest,
            state,
            REFERENCE_PRECISION,
            route=route,
            offset=offset,
        )
        solutions.append(solution)
        diagnostics.append(
            {
                "offset": list(offset),
                "iterations": solution.iterations,
                "mass_residual": str(solution.mass_residual),
                "charge_residual": str(solution.charge_residual),
                "xnet_basis_charge_residual": str(
                    solution.xnet_charge_residual
                ),
            }
        )
    primary = solutions[0]
    differences = [
        composition_norms(primary.composition, item.composition)
        for item in solutions[1:]
    ]
    if any(max(pair) > ROUTE_DIFFERENCE_LIMIT for pair in differences):
        raise RuntimeError(f"{state.state_id} {route} starting points disagree")
    return primary, diagnostics


def constraint_boundary(
    manifest: dict[str, Any], state: State, nominal: Any, residual_gate: Decimal
) -> dict[str, Any]:
    results = []
    for mass_sign, charge_sign in itertools.product((-1, 0, 1), repeat=2):
        if mass_sign == 0 and charge_sign == 0:
            continue
        mass_target = Decimal(mass_sign) * residual_gate
        charge_target = Decimal(charge_sign) * residual_gate
        shifted = solve(
            manifest,
            state,
            REFERENCE_PRECISION,
            constraint_target=(str(mass_target), str(charge_target)),
        )
        l1, linf = composition_norms(nominal.composition, shifted.composition)
        results.append(
            {
                "mass_constraint_target": str(mass_target),
                "charge_constraint_target": str(charge_target),
                "l1": str(l1),
                "linf": str(linf),
            }
        )
    return {
        "evaluated_targets": results,
        "maximum_l1": str(max(Decimal(item["l1"]) for item in results)),
        "maximum_linf": str(max(Decimal(item["linf"]) for item in results)),
    }


def make_state_record(
    manifest: dict[str, Any], state: State, fp_budget: Decimal
) -> dict[str, Any]:
    analytic, analytic_starts = solve_from_starts(
        manifest, state, "analytic-newton"
    )
    numeric, numeric_starts = solve_from_starts(manifest, state, "numeric-newton")
    route_l1, route_linf = composition_norms(
        analytic.composition, numeric.composition
    )
    if max(route_l1, route_linf) > ROUTE_DIFFERENCE_LIMIT:
        raise RuntimeError(
            f"analytic/numerical Jacobian variants disagree for {state.state_id}"
        )
    if max(
        abs(analytic.mass_residual),
        abs(analytic.charge_residual),
        abs(analytic.xnet_charge_residual),
    ) > REFERENCE_RESIDUAL_LIMIT:
        raise RuntimeError(f"reference residual is too large for {state.state_id}")

    precision_checks = []
    for precision in PRECISION_CHECKS:
        checked = solve(manifest, state, precision)
        l1, linf = composition_norms(analytic.composition, checked.composition)
        binary64_identical = all(
            float(left).hex() == float(right).hex()
            for left, right in zip(
                analytic.composition, checked.composition, strict=True
            )
        )
        if not binary64_identical:
            raise RuntimeError(
                f"{precision}-digit result changes stored values for {state.state_id}"
            )
        precision_checks.append(
            {
                "precision_decimal_digits": precision,
                "l1": str(l1),
                "linf": str(linf),
                "stored_binary64_identical": binary64_identical,
            }
        )

    stored = [binary64_string(value) for value in analytic.composition]
    stored_decimals = [Decimal.from_float(float(value)) for value in stored]
    serialization_l1, serialization_linf = composition_norms(
        analytic.composition, stored_decimals
    )
    residual_gate = XNET_SOLVER_RESIDUAL + fp_budget
    boundary = constraint_boundary(manifest, state, analytic, residual_gate)
    l1_budget = (
        Decimal(boundary["maximum_l1"])
        + route_l1
        + serialization_l1
        + fp_budget
    )
    linf_budget = (
        Decimal(boundary["maximum_linf"])
        + route_linf
        + serialization_linf
        + fp_budget
    )
    l1_tolerance = round_up_significant(l1_budget, 4)
    linf_tolerance = round_up_significant(linf_budget, 4)
    reconstructed_ye_budget = (
        residual_gate + abs(state.ye) * residual_gate + fp_budget
    )
    reconstructed_ye_tolerance = round_up_significant(
        reconstructed_ye_budget, 5
    )

    dominant = sorted(
        (
            (value, item["name"])
            for item, value in zip(
                manifest["species"], analytic.composition, strict=True
            )
        ),
        reverse=True,
    )[:5]
    return {
        "id": state.state_id,
        "inputs": {
            "rho_g_cm3": binary64_string(state.rho),
            "rho_binary64_hex": float(state.rho).hex(),
            "t9_gk": binary64_string(state.t9),
            "t9_binary64_hex": float(state.t9).hex(),
            "ye": binary64_string(state.ye),
            "ye_binary64_hex": float(state.ye).hex(),
        },
        "reference": {
            "eta_n": str(analytic.eta_n),
            "eta_p": str(analytic.eta_p),
            "mass_residual": str(analytic.mass_residual),
            "charge_residual": str(analytic.charge_residual),
            "xnet_basis_charge_residual": str(
                analytic.xnet_charge_residual
            ),
            "jacobian_condition_inf": str(
                jacobian_condition(
                    manifest, state, analytic.eta_n, analytic.eta_p
                )
            ),
            "analytic_starts": analytic_starts,
            "numeric_starts": numeric_starts,
            "route_difference_l1": str(route_l1),
            "route_difference_linf": str(route_linf),
            "precision_checks": precision_checks,
            "serialization_l1": str(serialization_l1),
            "serialization_linf": str(serialization_linf),
            "constraint_boundary": boundary,
        },
        "tolerances": {
            "mass_normalization_absolute": str(residual_gate),
            "xnet_charge_residual_absolute": str(residual_gate),
            "reconstructed_ye_absolute": str(reconstructed_ye_tolerance),
            "composition_l1_absolute": str(l1_tolerance),
            "composition_linf_absolute": str(linf_tolerance),
            "finite": "required exactly",
            "nonnegative": "required exactly",
            "species_identity_and_order": "required exactly",
            "derivation": {
                "composition_l1_unrounded_budget": str(l1_budget),
                "composition_linf_unrounded_budget": str(linf_budget),
                "terms": (
                    "maximum composition displacement at the four corners and four "
                    "edge midpoints sampled on the accepted mass/XNet-charge residual "
                    "box + analytic/numerical-Jacobian difference + binary64 "
                    "serialization difference + binary64 operation budget; rounded "
                    "upward to four significant digits"
                ),
                "reconstructed_ye": (
                    "|sum(qX)-Ye| <= |sum((q-Ye)X)| + "
                    "|Ye|*|sum(X)-1| plus the binary64 operation budget"
                ),
            },
        },
        "dominant_species": [
            {"name": name, "mass_fraction": binary64_string(value)}
            for value, name in dominant
        ],
        "composition": [
            {
                "name": item["name"],
                "a": item["a"],
                "z": item["z"],
                "n": item["n"],
                "mass_fraction": value,
            }
            for item, value in zip(manifest["species"], stored, strict=True)
        ],
    }


def scientific_dataset_hash(payload: dict[str, Any]) -> str:
    scientific = {
        key: value
        for key, value in payload.items()
        if key not in ("scientific_dataset_sha256", "reference_data_sha256")
    }
    # This independently checked identity is provenance for the unchanged
    # numerical dataset, so adding it must not relabel the accepted results.
    scientific["network"] = {
        key: value
        for key, value in scientific["network"].items()
        if key != "scientific_input_sha256"
    }
    return hashlib.sha256(canonical_bytes(scientific)).hexdigest()


def reference_data_text(payload: dict[str, Any]) -> str:
    lines = [
        "XNET_NSE_REFERENCE_V2",
        f"{payload['network']['species_count']} {len(payload['states'])}",
        payload["network"]["order_sha256"],
    ]
    for state in payload["states"]:
        inputs = state["inputs"]
        tolerances = state["tolerances"]
        lines.extend(
            (
                f"STATE {state['id']}",
                f"{inputs['rho_g_cm3']} {inputs['t9_gk']} {inputs['ye']}",
                " ".join(
                    tolerances[key]
                    for key in (
                        "mass_normalization_absolute",
                        "xnet_charge_residual_absolute",
                        "reconstructed_ye_absolute",
                        "composition_l1_absolute",
                        "composition_linf_absolute",
                    )
                ),
            )
        )
        lines.extend(
            f"{item['name']:>5} {item['a']} {item['z']} {item['n']} "
            f"{item['mass_fraction']}"
            for item in state["composition"]
        )
    return "\n".join(lines) + "\n"


def generate(network_directory: Path) -> tuple[dict[str, Any], str]:
    manifest = extract(network_directory, REPOSITORY_ROOT)
    if manifest["network"]["species_count"] != 489:
        raise RuntimeError("the frozen NSE network must contain 489 species")
    fp = floating_point_budget(manifest["network"]["species_count"])
    with localcontext() as context:
        context.prec = REFERENCE_PRECISION + 30
        state_records = [
            make_state_record(
                manifest,
                State.from_strings(state_id, rho, t9, ye),
                Decimal(fp["accepted_composition_budget"]),
            )
            for state_id, rho, t9, ye in STATES
        ]
    directory = Path(__file__).resolve().parent
    generator_files = (
        directory / "extract_inputs.py",
        directory / "reference_solver.py",
        directory / "generate_reference.py",
    )
    payload: dict[str, Any] = {
        "schema": "xnet-independent-nse-reference-v2",
        "issue": "https://github.com/jaharris87/XNet/issues/41",
        "xnet_base_commit": manifest["xnet_base_commit"],
        "scientific_authority": {
            "primary": (
                "Seitenzahl et al. 2009, DOI 10.1016/j.adt.2008.08.001, "
                "equations (2)-(10)"
            ),
            "corroborating_formulations": [
                "Hix and Meyer 2006, DOI 10.1016/j.nuclphysa.2004.10.009, equations (21)-(22)",
                "Lippuner and Roberts 2017, DOI 10.3847/1538-4365/aa94cb, Appendix B",
            ],
        },
        "generator": {
            "version": GENERATOR_VERSION,
            "language": "Python standard library only",
            "arithmetic": f"Decimal, {REFERENCE_PRECISION} requested decimal digits",
            "runtime": {
                "python_implementation": platform.python_implementation(),
                "python_version": platform.python_version(),
                "decimal_module_version": decimal.__version__,
                "libmpdec_version": decimal.__libmpdec_version__,
            },
            "residual_limit": str(REFERENCE_RESIDUAL_LIMIT),
            "route_difference_limit": str(ROUTE_DIFFERENCE_LIMIT),
            "files": {
                path.name: sha256(path)
                for path in generator_files
            },
        },
        "network": {
            **manifest["network"],
            "source_selection": "test/build_net/sunet.torch489",
            "source_selection_sha256": sha256(
                REPOSITORY_ROOT / "test/build_net/sunet.torch489"
            ),
            "build_input": "test/nse_validation/network/build_input.namelist",
            "build_input_sha256": sha256(
                directory / "network" / "build_input.namelist"
            ),
            "generated_with": (
                "test/build_net using winvne_JINAv22 and the mass sources "
                "recorded below; weak and neutrino-rate output disabled"
            ),
            "network_builder_source_hashes": {
                path: sha256(REPOSITORY_ROOT / path)
                for path in (
                    "test/build_net/net_module.f90",
                    "test/build_net/partf_module.f90",
                    "test/build_net/ffn_module.f90",
                    "test/build_net/nnu_module.f90",
                    "test/build_net/file_module.f90",
                    "test/build_net/reaclib_reader.f90",
                )
            },
            "canonical_input_sha256": manifest["canonical_input_sha256"],
            "scientific_input_sha256": manifest["scientific_input_sha256"],
        },
        "conventions": manifest["conventions"],
        "constants": manifest["constants"],
        "raw_data_provenance": manifest["raw_data_provenance"],
        "source_hashes": manifest["source_hashes"],
        "floating_point_budget": fp,
        "states": state_records,
        "limitations": [
            "The result is equilibrium over exactly this finite 489-species set.",
            "Screening and every other nonideal correction are disabled.",
            "The test validates a static NSE composition, not the timescale for reaching NSE.",
            "It does not validate reaction rates, weak evolution, or screened NSE.",
        ],
    }
    payload["scientific_dataset_sha256"] = scientific_dataset_hash(payload)
    data_text = reference_data_text(payload)
    payload["reference_data_sha256"] = hashlib.sha256(data_text.encode("ascii")).hexdigest()
    return payload, data_text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("network_directory", type=Path)
    parser.add_argument("json_output", type=Path)
    parser.add_argument("fortran_data_output", type=Path)
    args = parser.parse_args()
    payload, data_text = generate(args.network_directory.resolve())
    args.json_output.write_bytes(canonical_bytes(payload))
    args.fortran_data_output.write_text(data_text, encoding="ascii")
    print(payload["scientific_dataset_sha256"])
    print(payload["reference_data_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
