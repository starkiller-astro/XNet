#!/usr/bin/env python3
"""Focused tests for the independent reference and its retained inputs."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import tempfile
import unittest
from unittest.mock import patch
from decimal import Decimal, localcontext
from pathlib import Path

from extract_inputs import REPOSITORY_ROOT, extract, scientific_input_sha256, sha256
from generate_reference import scientific_dataset_hash
from preflight_states import analyze, compare_pair
from reference_solver import (
    Species,
    State,
    build_species,
    composition_norms,
    evaluate,
    solve,
)


DIRECTORY = Path(__file__).resolve().parent
NETWORK_DIRECTORY = DIRECTORY / "network"
REFERENCE_JSON = DIRECTORY / "reference.json"
REFERENCE_DATA = DIRECTORY / "reference.dat"


class ReferenceToolingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads(REFERENCE_JSON.read_text(encoding="utf-8"))
        cls.manifest = extract(NETWORK_DIRECTORY, REPOSITORY_ROOT)

    def test_retained_network_and_provenance(self) -> None:
        selected = REPOSITORY_ROOT / "test/build_net/sunet.torch489"
        self.assertEqual(
            (NETWORK_DIRECTORY / "sunet").read_bytes(), selected.read_bytes()
        )
        self.assertEqual(
            self.payload["network"]["scientific_input_sha256"],
            self.manifest["scientific_input_sha256"],
        )
        self.assertEqual(
            self.payload["network"]["sunet_sha256"],
            sha256(NETWORK_DIRECTORY / "sunet"),
        )
        self.assertEqual(
            self.payload["network"]["netwinv_sha256"],
            sha256(NETWORK_DIRECTORY / "netwinv"),
        )
        self.assertEqual(
            self.payload["network"]["build_input_sha256"],
            sha256(NETWORK_DIRECTORY / "build_input.namelist"),
        )
        historical_hashes = (
            self.payload["generator"]["files"],
            self.payload["network"]["network_builder_source_hashes"],
            self.payload["source_hashes"],
        )
        for hashes in historical_hashes:
            for value in hashes.values():
                self.assertRegex(value, r"^[0-9a-f]{64}$")
        self.assertEqual(
            self.payload["generator"]["runtime"],
            {
                "python_implementation": "CPython",
                "python_version": "3.13.0",
                "decimal_module_version": "1.70",
                "libmpdec_version": "4.0.1",
            },
        )
        archive = self.payload["raw_data_provenance"]["archive"]
        self.assertEqual(
            archive["initial_database_commit"],
            "77141ca2a3dfc9fa9fd52ef0fcf39a49d74c08e1",
        )
        self.assertEqual(
            archive["mass_reac1_git_blob_sha1"],
            "acd42b416e52edea020990d68631fb2b8063275d",
        )
        self.assertEqual(
            archive["winvne_JINAv22_git_blob_sha1"],
            "5bac4be2a2095bfd6a6fe3b1282bea4bcd6040ea",
        )

    def test_complete_dataset_identity_and_hashes(self) -> None:
        expected_names = [item["name"] for item in self.manifest["species"]]
        self.assertEqual(len(expected_names), 489)
        self.assertEqual(len(set(expected_names)), len(expected_names))
        for state in self.payload["states"]:
            names = [item["name"] for item in state["composition"]]
            self.assertEqual(names, expected_names)
            self.assertEqual(len(set(names)), len(names))
            self.assertEqual(
                [
                    (item["a"], item["z"], item["n"])
                    for item in state["composition"]
                ],
                [
                    (item["a"], item["z"], item["n"])
                    for item in self.manifest["species"]
                ],
            )
            numeric_tolerances = [
                float(state["tolerances"][name])
                for name in (
                    "mass_normalization_absolute",
                    "xnet_charge_residual_absolute",
                    "reconstructed_ye_absolute",
                    "composition_l1_absolute",
                    "composition_linf_absolute",
                )
            ]
            self.assertTrue(
                all(math.isfinite(value) and value > 0.0 for value in numeric_tolerances)
            )
            values = [float(item["mass_fraction"]) for item in state["composition"]]
            self.assertTrue(all(math.isfinite(value) for value in values))
            self.assertTrue(all(value >= 0.0 for value in values))
            self.assertLess(abs(math.fsum(values) - 1.0), 5.0e-15)
            reconstructed_ye = math.fsum(
                value * item["z"] / item["a"]
                for value, item in zip(
                    values, self.manifest["species"], strict=True
                )
            )
            self.assertLess(
                abs(reconstructed_ye - float(state["inputs"]["ye"])), 5.0e-15
            )
            reference = state["reference"]
            self.assertLess(abs(Decimal(reference["mass_residual"])), Decimal("1e-25"))
            self.assertLess(abs(Decimal(reference["charge_residual"])), Decimal("1e-25"))
            self.assertLess(
                abs(Decimal(reference["xnet_basis_charge_residual"])),
                Decimal("1e-25"),
            )
            self.assertGreater(Decimal(reference["jacobian_condition_inf"]), 0)
            self.assertLess(Decimal(reference["route_difference_l1"]), Decimal("1e-24"))
            self.assertLess(Decimal(reference["route_difference_linf"]), Decimal("1e-24"))
        self.assertEqual(
            self.payload["scientific_dataset_sha256"],
            scientific_dataset_hash(self.payload),
        )
        self.assertEqual(
            self.payload["reference_data_sha256"],
            hashlib.sha256(REFERENCE_DATA.read_bytes()).hexdigest(),
        )

    def test_source_byte_change_does_not_change_scientific_input_identity(self) -> None:
        util_path = REPOSITORY_ROOT / "source/xnet_util.F90"

        def changed_source_hash(path: Path) -> str:
            if path.resolve() == util_path.resolve():
                return "0" * 64
            return sha256(path)

        with patch("extract_inputs.sha256", side_effect=changed_source_hash):
            changed = extract(NETWORK_DIRECTORY, REPOSITORY_ROOT)
        self.assertNotEqual(
            self.manifest["source_hashes"]["source/xnet_util.F90"],
            changed["source_hashes"]["source/xnet_util.F90"],
        )
        self.assertNotEqual(
            self.manifest["canonical_input_sha256"],
            changed["canonical_input_sha256"],
        )
        self.assertEqual(
            self.manifest["scientific_input_sha256"],
            changed["scientific_input_sha256"],
        )
        equivalent_literal = copy.deepcopy(self.manifest)
        equivalent_literal["constants"]["bok"]["source_lexeme"] = "0.086173303"
        self.assertEqual(
            self.manifest["scientific_input_sha256"],
            scientific_input_sha256(equivalent_literal),
        )

    def test_scientific_identity_ignores_reconciliation_provenance(self) -> None:
        unreconciled = extract(
            NETWORK_DIRECTORY, REPOSITORY_ROOT, reconcile_provenance=False
        )
        self.assertNotEqual(
            self.manifest["raw_data_provenance"],
            unreconciled["raw_data_provenance"],
        )
        self.assertEqual(
            self.manifest["scientific_input_sha256"],
            unreconciled["scientific_input_sha256"],
        )

    def test_scientific_identity_rejects_extracted_constant_change(self) -> None:
        mutated = copy.deepcopy(self.manifest)
        mutated["constants"]["bok"]["hex"] = "0x1.0p-3"
        self.assertNotEqual(
            self.payload["network"]["scientific_input_sha256"],
            scientific_input_sha256(mutated),
        )

    def test_reference_integrity_rejects_data_and_result_changes(self) -> None:
        mutated = copy.deepcopy(self.payload)
        mutated["states"][0]["reference"]["mass_residual"] = "9.99"
        self.assertNotEqual(
            self.payload["scientific_dataset_sha256"],
            scientific_dataset_hash(mutated),
        )
        mutated_data = REFERENCE_DATA.read_bytes() + b"\n"
        self.assertNotEqual(
            self.payload["reference_data_sha256"],
            hashlib.sha256(mutated_data).hexdigest(),
        )

    def test_fortran_data_has_complete_stable_order(self) -> None:
        lines = REFERENCE_DATA.read_text(encoding="ascii").splitlines()
        self.assertEqual(lines[0], "XNET_NSE_REFERENCE_V2")
        species_count, state_count = (int(value) for value in lines[1].split())
        self.assertEqual((species_count, state_count), (489, 3))
        self.assertEqual(lines[2], self.manifest["network"]["order_sha256"])
        cursor = 3
        expected_names = [item["name"] for item in self.manifest["species"]]
        for expected_state in self.payload["states"]:
            self.assertEqual(lines[cursor], f"STATE {expected_state['id']}")
            self.assertEqual(
                lines[cursor + 1].split(),
                [
                    expected_state["inputs"]["rho_g_cm3"],
                    expected_state["inputs"]["t9_gk"],
                    expected_state["inputs"]["ye"],
                ],
            )
            self.assertEqual(
                lines[cursor + 2].split(),
                [
                    expected_state["tolerances"][name]
                    for name in (
                        "mass_normalization_absolute",
                        "xnet_charge_residual_absolute",
                        "reconstructed_ye_absolute",
                        "composition_l1_absolute",
                        "composition_linf_absolute",
                    )
                ],
            )
            cursor += 3
            names = []
            identities = []
            values = []
            for _ in range(species_count):
                name, aa, zz, nn, value = lines[cursor].split()
                names.append(name)
                identities.append((int(aa), int(zz), int(nn)))
                values.append(value)
                cursor += 1
            self.assertEqual(names, expected_names)
            self.assertEqual(
                identities,
                [
                    (item["a"], item["z"], item["n"])
                    for item in self.manifest["species"]
                ],
            )
            self.assertEqual(
                values,
                [item["mass_fraction"] for item in expected_state["composition"]],
            )
        self.assertEqual(cursor, len(lines))

    def test_two_species_constraints_have_closed_form_solution(self) -> None:
        with localcontext() as context:
            context.prec = 80
            ye = Decimal("0.4")
            species = (
                Species("n", 1, 0, 1, Decimal(0)),
                Species("p", 1, 1, 0, Decimal(0)),
            )
            result = evaluate(species, (Decimal(1) - ye).ln(), ye.ln(), ye)
        self.assertLess(max(abs(value) for value in result.residual), Decimal("1e-70"))
        self.assertLess(
            abs(result.normalized_composition[0] - (Decimal(1) - ye)),
            Decimal("1e-70"),
        )
        self.assertLess(
            abs(result.normalized_composition[1] - ye), Decimal("1e-70")
        )

    def test_analytic_jacobian_matches_independent_finite_difference(self) -> None:
        record = self.payload["states"][0]
        state = State.from_strings(
            record["id"],
            record["inputs"]["rho_g_cm3"],
            record["inputs"]["t9_gk"],
            record["inputs"]["ye"],
        )
        with localcontext() as context:
            context.prec = 80
            species = build_species(self.manifest, state, 60)
            eta_n = Decimal(record["reference"]["eta_n"])
            eta_p = Decimal(record["reference"]["eta_p"])
            analytic = evaluate(species, eta_n, eta_p, state.ye).jacobian
            step = Decimal("1e-20")
            n_plus = evaluate(species, eta_n + step, eta_p, state.ye).residual
            n_minus = evaluate(species, eta_n - step, eta_p, state.ye).residual
            p_plus = evaluate(species, eta_n, eta_p + step, state.ye).residual
            p_minus = evaluate(species, eta_n, eta_p - step, state.ye).residual
            numeric = (
                (
                    (n_plus[0] - n_minus[0]) / (2 * step),
                    (p_plus[0] - p_minus[0]) / (2 * step),
                ),
                (
                    (n_plus[1] - n_minus[1]) / (2 * step),
                    (p_plus[1] - p_minus[1]) / (2 * step),
                ),
            )
        for analytic_row, numeric_row in zip(analytic, numeric, strict=True):
            for analytic_value, numeric_value in zip(
                analytic_row, numeric_row, strict=True
            ):
                self.assertLess(
                    abs(analytic_value - numeric_value), Decimal("1e-35")
                )

    def test_mass_input_perturbation_changes_independent_result(self) -> None:
        record = self.payload["states"][0]
        state = State.from_strings(
            record["id"],
            record["inputs"]["rho_g_cm3"],
            record["inputs"]["t9_gk"],
            record["inputs"]["ye"],
        )
        nominal = solve(self.manifest, state, 35)
        perturbed_manifest = copy.deepcopy(self.manifest)
        target = next(
            item for item in perturbed_manifest["species"] if item["name"] == "co55"
        )
        changed = float.fromhex(target["binding_energy_mev"]["hex"]) + 0.010
        target["binding_energy_mev"] = {
            "decimal": format(changed, ".17g"),
            "hex": changed.hex(),
        }
        perturbed = solve(perturbed_manifest, state, 35)
        l1, linf = composition_norms(nominal.composition, perturbed.composition)
        self.assertGreater(
            l1, Decimal(record["tolerances"]["composition_l1_absolute"])
        )
        self.assertGreater(
            linf, Decimal(record["tolerances"]["composition_linf_absolute"])
        )

    def test_precision_checks_are_binary64_stable(self) -> None:
        for state in self.payload["states"]:
            checks = state["reference"]["precision_checks"]
            self.assertEqual(
                [item["precision_decimal_digits"] for item in checks], [35, 65]
            )
            self.assertTrue(
                all(item["stored_binary64_identical"] for item in checks)
            )

    def test_preflight_rejects_different_shared_scientific_inputs(self) -> None:
        left = copy.deepcopy(self.manifest)
        right = copy.deepcopy(self.manifest)
        right["species"][0]["mass_excess_mev"]["hex"] = "0x1.0p+99"
        with tempfile.TemporaryDirectory() as temporary_directory:
            left_path = Path(temporary_directory) / "left.json"
            right_path = Path(temporary_directory) / "right.json"
            left_path.write_text(json.dumps(left), encoding="utf-8")
            right_path.write_text(json.dumps(right), encoding="utf-8")
            result = {
                "network": left["network"],
                "manifest": str(left_path),
                "manifest_sha256": hashlib.sha256(left_path.read_bytes()).hexdigest(),
                "_manifest_snapshot": left,
                "states": [],
            }
            changed = {
                "network": right["network"],
                "manifest": str(right_path),
                "manifest_sha256": hashlib.sha256(right_path.read_bytes()).hexdigest(),
                "_manifest_snapshot": right,
                "states": [],
            }
            with self.assertRaisesRegex(
                RuntimeError, "scientific inputs differ for shared species n"
            ):
                compare_pair(result, changed)

    def test_preflight_rejects_manifest_rewritten_after_analysis(self) -> None:
        manifest = copy.deepcopy(self.manifest)
        with tempfile.TemporaryDirectory() as temporary_directory:
            manifest_path = Path(temporary_directory) / "manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with patch("preflight_states.CANDIDATE_STATES", ()):
                result = analyze(manifest_path, 35)
            manifest["species"][0]["mass_excess_mev"]["hex"] = "0x1.0p+99"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(
                RuntimeError, "preflight manifest changed after analysis"
            ):
                compare_pair(result, result)


if __name__ == "__main__":
    unittest.main()
