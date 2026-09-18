#!/usr/bin/env python3
"""Extract XNet nuclear inputs without importing or executing XNet NSE code.

This program is data-compatibility tooling, not an NSE solver.  It records the
exact finite species set, source decimal values, binary64 values consumed by
the tracked 64-bit-real build configurations, and the source-file hashes used
to define one independent equilibrium problem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CONSTANT_NAMES = ("pi", "hbar", "avn", "bok", "epmev")
MASS_SOURCE_FILES = {
    "ame03": "mass_ame03.dat",
    "ame03extrap": "mass_ame03extrap.dat",
    "ame11": "mass_ame11.dat",
    "ame11extrap": "mass_ame11extrap.dat",
    "reac1": "mass_reac1.dat",
    "frdm": "mass_frdm.dat",
}
BUILD_NET_ARCHIVE = {
    "repository": "https://github.com/jaharris87/build_net",
    "initial_database_commit": "77141ca2a3dfc9fa9fd52ef0fcf39a49d74c08e1",
    "initial_database_commit_date": "2017-01-24",
    "xnet_subtree_import_commit": "90e9363d5f9443a8ad2d5e986232c1f60bb5b96a",
    "mass_reac1_git_blob_sha1": "acd42b416e52edea020990d68631fb2b8063275d",
    "winvne_JINAv22_git_blob_sha1": "5bac4be2a2095bfd6a6fe3b1282bea4bcd6040ea",
}


class ExtractionError(RuntimeError):
    """Raised when a source input cannot be reconciled unambiguously."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binary64(value: str | Decimal | float) -> dict[str, str]:
    parsed = float(value)
    return {"decimal": format(parsed, ".17g"), "hex": parsed.hex()}


def canonical_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def scientific_input_sha256(payload: dict[str, Any]) -> str:
    """Identify the retained calculation inputs, independent of source bytes."""
    constants = {
        name: {"hex": record["hex"], "units": record["units"]}
        for name, record in payload["constants"].items()
    }
    species = [
        {
            "index": record["index"],
            "name": record["name"],
            "a": record["a"],
            "z": record["z"],
            "n": record["n"],
            "spin_hex": record["spin"]["hex"],
            "ground_state_degeneracy_hex": record[
                "ground_state_degeneracy"
            ]["hex"],
            "mass_excess_mev_hex": record["mass_excess_mev"]["hex"],
            "partition_factor_hexes": [
                value["hex"] for value in record["partition_factors"]
            ],
            "binding_energy_mev_hex": record["binding_energy_mev"]["hex"],
            "translational_mass_g_hex": record["translational_mass_g"]["hex"],
        }
        for record in payload["species"]
    ]
    scientific_inputs = {
        "schema": payload["schema"],
        "network": {
            key: payload["network"][key]
            for key in (
                "species_count",
                "sunet_sha256",
                "netwinv_sha256",
                "order_sha256",
            )
        },
        "conventions": payload["conventions"],
        "constants": constants,
        "temperature_grid_gk_hexes": [
            value["hex"] for value in payload["temperature_grid_gk"]
        ],
        "species": species,
    }
    return hashlib.sha256(canonical_bytes(scientific_inputs)).hexdigest()


def parse_constants(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8")
    constants: dict[str, dict[str, str]] = {}
    for name in CONSTANT_NAMES:
        match = re.search(
            rf"Real\(dp\),\s*Parameter\s*::\s*{name}\s*=\s*([^\s!]+)", text
        )
        if match is None:
            raise ExtractionError(f"missing XNet constant {name} in {path}")
        lexeme = match.group(1)
        constants[name] = {
            "source_lexeme": lexeme,
            "units": {
                "pi": "dimensionless",
                "hbar": "MeV s",
                "avn": "mol^-1",
                "bok": "MeV GK^-1",
                "epmev": "erg MeV^-1",
            }[name],
            **binary64(lexeme),
        }
    return constants


def parse_temperature_grid(line: str) -> list[dict[str, str | int]]:
    if len(line) < 72:
        raise ExtractionError("netwinv partition-temperature grid is truncated")
    integers = [int(line[index : index + 3]) for index in range(0, 72, 3)]
    if len(integers) != 24:
        raise ExtractionError("netwinv partition-temperature grid must have 24 nodes")
    result = []
    for index, value in enumerate(integers):
        scale = Decimal("0.1") if index == 23 else Decimal("0.01")
        decimal_value = Decimal(value) * scale
        result.append(
            {
                "encoded_integer": value,
                "source_decimal": str(decimal_value),
                **binary64(decimal_value),
            }
        )
    return result


def parse_netwinv(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ExtractionError(f"{path} is truncated")
    count = int(lines[0])
    grid = parse_temperature_grid(lines[1])
    original_names = lines[2 : 2 + count]
    if len(original_names) != count or any(len(name) != 5 for name in original_names):
        raise ExtractionError(f"{path} has an invalid five-character species list")
    names = [name.strip().lower() for name in original_names]
    if len(set(names)) != count:
        raise ExtractionError(f"{path} contains duplicate normalized species names")

    cursor = 2 + count
    records: list[dict[str, Any]] = []
    for index, (original_name, name) in enumerate(zip(original_names, names, strict=True)):
        if cursor + 3 >= len(lines):
            raise ExtractionError(f"{path} is truncated at species {name}")
        fields = lines[cursor].split()
        partition_lexemes = " ".join(lines[cursor + 1 : cursor + 4]).split()
        cursor += 4
        if len(fields) != 6 or len(partition_lexemes) != 24:
            raise ExtractionError(f"malformed netwinv record for {name}")
        if fields[0].lower() != name:
            raise ExtractionError(f"netwinv header/record mismatch for {name}")
        a_lexeme, z_lexeme, n_lexeme, spin_lexeme, mex_lexeme = fields[1:]
        a = int(Decimal(a_lexeme))
        z = int(z_lexeme)
        n = int(n_lexeme)
        if a != z + n:
            raise ExtractionError(f"A != Z + N for {name}")
        records.append(
            {
                "index": index + 1,
                "name": name,
                "xnet_name": original_name,
                "a": a,
                "z": z,
                "n": n,
                "spin": {"source_lexeme": spin_lexeme, **binary64(spin_lexeme)},
                "ground_state_degeneracy": binary64(
                    2.0 * float(spin_lexeme) + 1.0
                ),
                "mass_excess_mev": {
                    "source_lexeme": mex_lexeme,
                    **binary64(mex_lexeme),
                },
                "partition_factors": [
                    {"source_lexeme": value, **binary64(value)}
                    for value in partition_lexemes
                ],
            }
        )
    if any(line.strip() for line in lines[cursor:]):
        raise ExtractionError(f"{path} has unexpected trailing records")
    return {"count": count, "temperature_grid_gk": grid, "species": records}


def parse_master_partition(path: Path) -> dict[str, dict[str, Any]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    maximum_count = int(lines[0])
    names: list[str] = []
    cursor = 2
    while cursor < len(lines):
        name = lines[cursor].strip().lower()
        cursor += 1
        if names and name == names[-1]:
            break
        names.append(name)
        if len(names) > maximum_count:
            raise ExtractionError("master partition species list exceeds its header bound")
    else:
        raise ExtractionError("master partition species-list sentinel is missing")
    records: dict[str, dict[str, Any]] = {}
    for name in names:
        fields = lines[cursor].split()
        partition = " ".join(lines[cursor + 1 : cursor + 4]).split()
        cursor += 4
        if len(fields) != 7 or len(partition) != 24 or fields[0].lower() != name:
            raise ExtractionError(f"malformed master partition record for {name}")
        records[name] = {
            "a": int(Decimal(fields[1])),
            "z": int(fields[2]),
            "n": int(fields[3]),
            "spin": fields[4],
            "original_mass_excess_mev": fields[5],
            "mass_source": fields[6].lower(),
            "partition": partition,
        }
    return records


def parse_mass_source(path: Path) -> dict[tuple[int, int], Decimal]:
    values: dict[tuple[int, int], Decimal] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(" "):
            continue
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            z = int(fields[0])
            a = int(fields[1])
            mass_excess_kev = Decimal(fields[3])
        except (ValueError, ArithmeticError):
            continue
        values[(z, a)] = mass_excess_kev * Decimal("0.001")
    return values


def reconcile_raw_provenance(
    parsed: dict[str, Any], repository_root: Path
) -> dict[str, Any]:
    build_data = repository_root / "test" / "build_net"
    partition_path = build_data / "partf_data" / "winvne_JINAv22"
    master = parse_master_partition(partition_path)
    mass_tables = {
        label: parse_mass_source(build_data / "mass_data" / filename)
        for label, filename in MASS_SOURCE_FILES.items()
    }
    source_counts: dict[str, int] = {}
    for species in parsed["species"]:
        name = species["name"]
        if name not in master:
            raise ExtractionError(f"{name} is absent from winvne_JINAv22")
        source = master[name]
        if (source["a"], source["z"], source["n"]) != (
            species["a"],
            species["z"],
            species["n"],
        ):
            raise ExtractionError(f"raw partition identity differs for {name}")
        if float(source["spin"]) != float(species["spin"]["source_lexeme"]):
            raise ExtractionError(f"raw partition spin differs for {name}")
        if [float(value) for value in source["partition"]] != [
            float(value["source_lexeme"]) for value in species["partition_factors"]
        ]:
            raise ExtractionError(f"raw partition factors differ for {name}")
        label = source["mass_source"]
        if label in mass_tables:
            try:
                source_mass = mass_tables[label][(species["z"], species["a"])]
            except KeyError as error:
                raise ExtractionError(
                    f"{label} has no mass entry for {name}"
                ) from error
        else:
            source_mass = Decimal(source["original_mass_excess_mev"])
        stored_mass = Decimal(species["mass_excess_mev"]["source_lexeme"])
        if abs(source_mass - stored_mass) > Decimal("5.0001e-9"):
            raise ExtractionError(
                f"stored mass for {name} is not the source value rounded to eight decimals"
            )
        species["mass_source"] = label
        species["mass_source_value_mev"] = str(source_mass)
        source_counts[label] = source_counts.get(label, 0) + 1
    return {
        "partition_source": "test/build_net/partf_data/winvne_JINAv22",
        "partition_source_sha256": sha256(partition_path),
        "archive": {
            **BUILD_NET_ARCHIVE,
            "partition_snapshot": (
                "JINA REACLIB V2.2 snapshot dated 2016-11-14; the exact retained "
                "bytes are fixed by the commit, Git blob, and SHA-256 recorded here"
            ),
            "partition_archive_url": (
                "https://reaclib.jinaweb.org/library.php?action=viewsnapshots"
            ),
            "mass_reac1_origin": (
                "The retained file header identifies the JINA Nuclide Database "
                "evaluation label reac1. The public build_net archive fixes its "
                "exact bytes, but the original JINA per-record publication or "
                "snapshot identifier is not recoverable from retained metadata."
            ),
        },
        "mass_source_counts": dict(sorted(source_counts.items())),
        "mass_source_files": {
            label: {
                "path": f"test/build_net/mass_data/{filename}",
                "sha256": sha256(build_data / "mass_data" / filename),
            }
            for label, filename in MASS_SOURCE_FILES.items()
        },
        "reconciliation": (
            "Every selected spin and partition factor matches winvne_JINAv22; "
            "every stored mass excess matches the source selected by that file "
            "after the build_net f15.8 MeV rounding."
        ),
    }


def extract(
    network_directory: Path,
    repository_root: Path,
    *,
    reconcile_provenance: bool = True,
) -> dict[str, Any]:
    network_directory = network_directory.resolve()
    sunet_path = network_directory / "sunet"
    netwinv_path = network_directory / "netwinv"
    parsed = parse_netwinv(netwinv_path)
    sunet_names = [
        line.strip().lower()
        for line in sunet_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    parsed_names = [species["name"] for species in parsed["species"]]
    if sunet_names != parsed_names:
        raise ExtractionError("sunet and netwinv species/order differ")

    constants_path = repository_root / "source" / "xnet_constants.F90"
    constants = parse_constants(constants_path)
    avn = float.fromhex(constants["avn"]["hex"])
    neutron = next(item for item in parsed["species"] if item["z"] == 0 and item["n"] == 1)
    proton = next(item for item in parsed["species"] if item["z"] == 1 and item["n"] == 0)
    mex_n = float.fromhex(neutron["mass_excess_mev"]["hex"])
    mex_p = float.fromhex(proton["mass_excess_mev"]["hex"])
    for species in parsed["species"]:
        mex = float.fromhex(species["mass_excess_mev"]["hex"])
        binding = species["n"] * mex_n + species["z"] * mex_p - mex
        species["binding_energy_mev"] = binary64(binding)
        species["translational_mass_g"] = binary64(species["a"] / avn)

    try:
        relative_network = network_directory.relative_to(repository_root)
    except ValueError:
        relative_network = Path("external-preflight") / network_directory.name
    payload: dict[str, Any] = {
        "schema": "xnet-independent-nse-input-v1",
        "xnet_base_commit": "65271bcbeea430534c1adc92748ef13bea10c228",
        "network": {
            "directory": relative_network.as_posix(),
            "species_count": parsed["count"],
            "sunet_sha256": sha256(sunet_path),
            "netwinv_sha256": sha256(netwinv_path),
            "order_sha256": hashlib.sha256(
                ("\n".join(parsed_names) + "\n").encode("ascii")
            ).hexdigest(),
        },
        "conventions": {
            "abundance": "X_i is mass fraction; Y_i = X_i / A_i",
            "density": "g cm^-3 with translational mass m_i = A_i / N_A",
            "electron_fraction": "Ye = sum_i (Z_i/A_i) X_i",
            "free_neutrons_and_protons": "ordinary finite-network species",
            "partition_function": (
                "(2J_i+1) times the JINA REACLIB normalized partition factor; "
                "linear interpolation in ln(G) versus T9 with endpoint clamping"
            ),
            "screening": "disabled exactly; Coulomb exponent h_i = 0",
            "other_nonideal_terms": "none in the validated xnet_nse abundance equation",
            "real_kind": (
                "binary64; tracked compiler configurations promote default real to 64 bits"
            ),
        },
        "constants": constants,
        "temperature_grid_gk": parsed["temperature_grid_gk"],
        "species": parsed["species"],
        "raw_data_provenance": (
            reconcile_raw_provenance(parsed, repository_root)
            if reconcile_provenance
            else {
                "reconciliation": (
                    "Skipped for exploratory network-boundary preflight; direct "
                    "netwinv values remain recorded exactly."
                )
            }
        ),
        "source_hashes": {
            path: sha256(repository_root / path)
            for path in (
                "source/xnet_constants.F90",
                "source/xnet_data.F90",
                "source/xnet_nse.F90",
                "source/xnet_util.F90",
            )
        },
    }
    payload["canonical_input_sha256"] = hashlib.sha256(canonical_bytes(payload)).hexdigest()
    payload["scientific_input_sha256"] = scientific_input_sha256(payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "network_directory",
        type=Path,
        help="tracked XNet Data_* directory containing sunet and netwinv",
    )
    parser.add_argument("output", type=Path, help="canonical JSON manifest to create")
    parser.add_argument(
        "--skip-raw-reconciliation",
        action="store_true",
        help="allow exploratory preflight when historical raw inputs are incomplete",
    )
    args = parser.parse_args()
    payload = extract(
        args.network_directory,
        REPOSITORY_ROOT,
        reconcile_provenance=not args.skip_raw_reconciliation,
    )
    args.output.write_bytes(canonical_bytes(payload))
    print(payload["canonical_input_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
