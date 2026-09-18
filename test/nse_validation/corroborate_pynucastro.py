#!/usr/bin/env python3
"""Optional secondary check with pynucastro's independently maintained solver.

This does not generate retained reference data and is never run by the ordinary
test target.  It configures pynucastro with the exact finite XNet inputs so that
the comparison is not obscured by its different default masses, translational
masses, constants, spins, or partition functions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

from extract_inputs import extract


PINNED_COMMIT = "a7268f86f42c556172578ca53293cf00b6c539ab"
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


class XNetLogLinearPartition:
    """Expose XNet's normalized, log-linear partition interpolation."""

    def __init__(self, t9: list[float], factors: list[float]):
        self.t9 = t9
        self.logs = [math.log(value) for value in factors]

    def eval(self, temperature_k: float) -> float:
        value = temperature_k * 1.0e-9
        if value <= self.t9[0]:
            return self.logs[0]
        if value >= self.t9[-1]:
            return self.logs[-1]
        upper = next(index for index, point in enumerate(self.t9) if point >= value)
        lower = upper - 1
        weight = (value - self.t9[lower]) / (self.t9[upper] - self.t9[lower])
        return (1.0 - weight) * self.logs[lower] + weight * self.logs[upper]


def verify_checkout(path: Path) -> None:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    actual = result.stdout.strip()
    if actual != PINNED_COMMIT:
        raise RuntimeError(f"pynucastro source is {actual}, expected {PINNED_COMMIT}")


def verify_imported_sources(checkout: Path, package_file: str) -> None:
    installed_root = Path(package_file).resolve().parent.parent
    for relative in (
        "pynucastro/constants/constants.py",
        "pynucastro/networks/nse_network.py",
        "pynucastro/networks/rate_collection.py",
        "pynucastro/nucdata/nucleus.py",
    ):
        expected_source = subprocess.run(
            ["git", "-C", str(checkout), "show", f"{PINNED_COMMIT}:{relative}"],
            check=True,
            capture_output=True,
        ).stdout
        expected = hashlib.sha256(expected_source).digest()
        actual = hashlib.sha256((installed_root / relative).read_bytes()).digest()
        if actual != expected:
            raise RuntimeError(
                f"imported pynucastro source differs from pinned commit for {relative}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pynucastro_source",
        type=Path,
        help="checkout used to install pynucastro; its exact commit is verified",
    )
    args = parser.parse_args()
    verify_checkout(args.pynucastro_source.resolve())

    import pynucastro as pyna
    from pynucastro.constants import constants as pyna_constants

    verify_imported_sources(args.pynucastro_source.resolve(), pyna.__file__)

    network_directory = Path(__file__).with_name("network")
    reference = json.loads(Path(__file__).with_name("reference.json").read_text())
    scientific_inputs = extract(network_directory, REPOSITORY_ROOT)
    constants = {
        name: float(item["decimal"])
        for name, item in reference["constants"].items()
    }

    # pynucastro otherwise uses actual nuclear translational masses and current
    # scipy constants.  Reconcile those choices to the finite XNet problem.
    pyna_constants.m_u_C18 = 1.0 / constants["avn"]
    pyna_constants.k = constants["bok"] * constants["epmev"] / 1.0e9
    pyna_constants.k_MeV = constants["bok"] / 1.0e9
    pyna_constants.hbar = constants["hbar"] * constants["epmev"]
    pyna_constants.erg2MeV = 1.0 / constants["epmev"]

    names = [item["name"] for item in scientific_inputs["species"]]
    nuclei = [pyna.Nucleus.from_cache(name) for name in names]
    t9_grid = [
        float(point["decimal"])
        for point in scientific_inputs["temperature_grid_gk"]
    ]
    for nucleus, data in zip(nuclei, scientific_inputs["species"], strict=True):
        if nucleus.raw != data["name"]:
            raise RuntimeError(f"species mismatch for {data['name']}")
        nucleus.A_nuc = float(data["a"])
        nucleus.nucbind = float(data["binding_energy_mev"]["decimal"]) / data["a"]
        nucleus.spin_states = float(data["ground_state_degeneracy"]["decimal"])
        factors = [float(item["decimal"]) for item in data["partition_factors"]]
        nucleus.partition_function = XNetLogLinearPartition(t9_grid, factors)

    network = pyna.NSENetwork(inert_nuclei=nuclei, use_unreliable_spins=True)
    print(f"pynucastro commit={PINNED_COMMIT} species={len(nuclei)} inputs=reconciled")
    for state in reference["states"]:
        inputs = state["inputs"]
        composition = network.get_comp_nse(
            float(inputs["rho_g_cm3"]),
            float(inputs["t9_gk"]) * 1.0e9,
            float(inputs["ye"]),
            tol=1.0e-11,
            use_coulomb_corr=False,
        )
        expected = {
            entry["name"]: float(entry["mass_fraction"])
            for entry in state["composition"]
        }
        differences = [
            abs(composition.X[nucleus] - expected[name])
            for name, nucleus in zip(names, nuclei, strict=True)
        ]
        mass = sum(composition.X.values()) - 1.0
        charge = sum(
            nucleus.Z / nucleus.A * composition.X[nucleus] for nucleus in nuclei
        ) - float(inputs["ye"])
        if not all(math.isfinite(value) for value in (*differences, mass, charge)):
            raise RuntimeError(f"nonfinite pynucastro result for {state['id']}")
        print(
            f"{state['id']} L1={sum(differences):.9e} "
            f"Linf={max(differences):.9e} mass={mass:.3e} charge={charge:.3e}"
        )


if __name__ == "__main__":
    main()
