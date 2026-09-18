#!/usr/bin/env python3

"""Exercise build_net with tiny sources and an isolated production smoke."""

from __future__ import annotations

import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


class ContractError(RuntimeError):
    pass


MASTER_SPECIES = (
    ("n", 1.0, 0, 1, 0.5, 8.07131710),
    ("p", 1.0, 1, 0, 0.5, 7.28897060),
    ("he4", 4.0, 2, 2, 0.0, 2.42491560),
    ("c12", 12.0, 6, 6, 0.0, 0.0),
    ("o16", 16.0, 8, 8, 0.0, -4.73700140),
    ("ne20", 20.0, 10, 10, 0.0, -7.04193060),
)
EXPECTED_SPECIES = ("n", "p", "he4", "c12", "o16")
REQUESTED_NAMES = ("N01", "H1", "HE4", "C12", "O16")
EXPECTED_MASS = {row[0]: row[5] for row in MASTER_SPECIES}
FORWARD_Q = EXPECTED_MASS["he4"] + EXPECTED_MASS["c12"] - EXPECTED_MASS["o16"]
WEAK_Q = EXPECTED_MASS["p"] - EXPECTED_MASS["n"]
GENERATED_ARTIFACTS = (
    "nuc_data",
    "nets3",
    "nets4",
    "ab_blank",
    "match_data",
    "match_read",
    "sparse_ind",
    "matr_shape",
    "net_desc",
    "net_diag",
)


def fail(message: str) -> None:
    raise ContractError(message)


def fixed_rate_header(
    chapter: int,
    names: tuple[str, ...] = (),
    descriptor: str = "",
    resonance: str = "",
    reverse: str = "",
    q_value: float = 0.0,
) -> str:
    prefix = f"{chapter:1d}" + " " * 4 if chapter < 10 else f"{chapter:2d}" + " " * 3
    participants = tuple(names) + ("",) * (6 - len(names))
    return (
        prefix
        + "".join(f"{name:>5}" for name in participants)
        + " " * 8
        + f"{descriptor:<4.4}{resonance:<1.1}{reverse:<1.1}"
        + " " * 3
        + f"{q_value:12.5E}"
    )


def fixed_coefficients(values: tuple[float, ...]) -> list[str]:
    return [
        "".join(f"{value:13.6E}" for value in values[:4]),
        "".join(f"{value:13.6E}" for value in values[4:]),
    ]


def write_reaclib(path: Path) -> None:
    zero = (0.0,) * 7
    slow = (-100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    lines: list[str] = []
    for chapter in range(1, 12):
        lines.append(fixed_rate_header(chapter))
        lines.extend(fixed_coefficients(zero))
        if chapter == 2:
            lines.append(
                fixed_rate_header(0, ("o16", "he4", "c12"), "syn1", "", "v", -FORWARD_Q)
            )
            lines.extend(fixed_coefficients(slow))
        if chapter == 4:
            lines.append(fixed_rate_header(0, ("he4", "c12", "o16"), "syn1", "", "", FORWARD_Q))
            lines.extend(fixed_coefficients(slow))
            lines.append(fixed_rate_header(0, ("he4", "o16", "ne20"), "drop", "", "", 9.0))
            lines.extend(fixed_coefficients(slow))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_partition_source(path: Path) -> None:
    lines = [
        f"{len(MASTER_SPECIES):5d}",
        "010015020030040050060070080090100150200250300350400450500600700800900100",
    ]
    lines.extend(f"{name:>5}" for name, *_ in MASTER_SPECIES)
    lines.append(f"{MASTER_SPECIES[-1][0]:>5}")
    for index, (name, aa, zz, nn, spin, mass) in enumerate(MASTER_SPECIES, start=1):
        raw_mass = mass + 0.125 * index
        lines.append(
            f"{name:>5}{aa:12.3f}{zz:4d}{nn:4d}{spin:6.1f}{raw_mass:10.3f} {'ame11':<11}"
        )
        values = (1.0,) * 24
        for offset in range(0, 24, 8):
            lines.append("".join(f"{value:12.5E}" for value in values[offset : offset + 8]))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_mass_source(path: Path) -> None:
    lines = [f"synthetic mass header {index:02d}" for index in range(1, 16)]
    first_data_index = len(lines)
    lines.extend(
        f"{zz:d} {int(aa):d} synthetic {mass * 1.0e3:.8f}"
        for _, aa, zz, _, _, mass in MASTER_SPECIES
    )
    # A missing uncertainty after a valid mass must not discard that mass.
    lines[first_data_index] += " #"
    # The retained JINA-derived tables use '#' for unavailable evaluations.
    # An unselected unavailable record must not make network construction fail.
    lines.append("75 203 synthetic #")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def weak_values_line(prefix: str) -> str:
    values = (-10.0, -10.0, -2.0, -10.0, -10.0, -2.0)
    return (
        prefix
        + f"{values[0]:8.3f} {values[1]:8.3f}  {values[2]:8.3f}"
        + f"  {values[3]:8.3f}  {values[4]:8.3f}  {values[5]:8.3f}"
    )


def write_weak_sources(directory: Path) -> None:
    directory.mkdir()
    element_lines = [f"{index:3d} {'xx':>2}" for index in range(119)]
    (directory / "element_list.txt").write_text("\n".join(element_lines) + "\n", encoding="ascii")
    lines = ["synthetic public weak source"]
    lines.append(weak_values_line(f" {'1':<3} {1:3d}" + " " * 33))
    lines.extend(weak_values_line(" " * 41) for _ in range(142))
    (directory / "updated_rate_table.txt").write_text("\n".join(lines) + "\n", encoding="ascii")


def namelist_text(weak_enabled: bool) -> str:
    weak_flag = ".true." if weak_enabled else ".false."
    return f"""&file_input
  new_data_dir = './out',
/
&net_input
  sunet_fname = 'sunet.requested',
  netsu_data_dir = './reaclib_data',
  netsu_in_fname = 'synthetic_reaclib',
  reaclib_ver = 1,
  no910 = .false.,
  netsu_out_fname = 'netsu',
/
&ffn_input
  netweak_flag = {weak_flag},
  netweak_data_dir = './weak_data',
  netweak_in_fname = 'updated_rate_table.txt',
  netweak_out_fname = 'netweak',
  element_list_fname = 'element_list.txt',
/
&partf_input
  netwinv_data_dir = './partf_data',
  netwinv_in_fname = 'synthetic_netwinv',
  mass_data_dir = './mass_data',
  ame03_fname = 'mass_ame03.dat',
  ame03extrap_fname = 'mass_ame03extrap.dat',
  ame11_fname = 'mass_ame11.dat',
  ame11extrap_fname = 'mass_ame11extrap.dat',
  reac1_fname = 'mass_reac1.dat',
  frdm_fname = 'mass_frdm.dat',
  netwinv_out_fname = 'netwinv',
/
&nnu_input
  netneutr_flag = .false.,
  netneutr_data_dir = './private-neutrino-data-not-present',
  netneutr_in_fname = 'neutrino.data',
  netneutr_out_fname = 'netneutr',
/
"""


def prepare_case(
    directory: Path,
    *,
    requested: tuple[str, ...] = REQUESTED_NAMES,
    weak_enabled: bool = True,
    malformed_namelist: bool = False,
    malformed_reaclib: bool = False,
    malformed_reaclib_late: bool = False,
    malformed_mass_hash_record: bool = False,
    missing_mass: bool = False,
    selected_unavailable_mass: bool = False,
    missing_weak: bool = False,
) -> None:
    directory.mkdir(parents=True)
    (directory / "out").mkdir()
    (directory / "reaclib_data").mkdir()
    (directory / "partf_data").mkdir()
    (directory / "mass_data").mkdir()
    (directory / "sunet.requested").write_text(
        "\n".join(requested) + "\n", encoding="ascii"
    )
    write_reaclib(directory / "reaclib_data" / "synthetic_reaclib")
    if malformed_reaclib:
        (directory / "reaclib_data" / "synthetic_reaclib").write_text(
            "not a REACLIB record\n", encoding="ascii"
        )
    if malformed_reaclib_late:
        reaclib_path = directory / "reaclib_data" / "synthetic_reaclib"
        lines = reaclib_path.read_text(encoding="ascii").splitlines()
        lines[3] = "not a later REACLIB header"
        reaclib_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    write_partition_source(directory / "partf_data" / "synthetic_netwinv")
    mass_names = (
        "mass_ame03.dat",
        "mass_ame03extrap.dat",
        "mass_ame11.dat",
        "mass_ame11extrap.dat",
        "mass_reac1.dat",
        "mass_frdm.dat",
    )
    for name in mass_names:
        write_mass_source(directory / "mass_data" / name)
    ame11_path = directory / "mass_data" / "mass_ame11.dat"
    if malformed_mass_hash_record:
        lines = ame11_path.read_text(encoding="ascii").splitlines()
        lines.append("this is not a mass record #")
        ame11_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    if selected_unavailable_mass:
        lines = ame11_path.read_text(encoding="ascii").splitlines()
        for index, line in enumerate(lines):
            if line.split()[:3] == ["6", "12", "synthetic"]:
                lines[index] = "6 12 synthetic #"
                break
        else:
            fail("could not locate selected c12 mass fixture")
        ame11_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    if missing_mass:
        ame11_path.unlink()
    if weak_enabled:
        write_weak_sources(directory / "weak_data")
        if missing_weak:
            (directory / "weak_data" / "updated_rate_table.txt").unlink()
    content = namelist_text(weak_enabled)
    if malformed_namelist:
        content = content.replace("  new_data_dir = './out'", "  unknown_file_control = 1")
    (directory / "input.namelist").write_text(content, encoding="ascii")


def run_process(
    command: list[str],
    directory: Path,
    *,
    input_text: str | None = None,
    timeout: float = 30.0,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        cwd=directory,
        input=input_text,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    stem = Path(command[0]).name
    (directory / f"{stem}.stdout.txt").write_text(result.stdout, encoding="utf-8")
    (directory / f"{stem}.stderr.txt").write_text(result.stderr, encoding="utf-8")
    return result


def require_success(label: str, result: subprocess.CompletedProcess[str]) -> None:
    if result.returncode != 0:
        fail(f"{label} returned {result.returncode}: {result.stderr or result.stdout}")


def require_failure(label: str, result: subprocess.CompletedProcess[str]) -> None:
    if result.returncode == 0:
        fail(f"{label} unexpectedly returned zero")
    if "error" not in (result.stdout + result.stderr).lower():
        fail(f"{label} failed without a clear error diagnostic")


def require_failure_diagnostic(
    label: str, result: subprocess.CompletedProcess[str], diagnostic: str
) -> None:
    require_failure(label, result)
    if diagnostic.lower() not in (result.stdout + result.stderr).lower():
        fail(f"{label} failed without expected diagnostic: {diagnostic}")


def require_reader_rejection(
    label: str, result: subprocess.CompletedProcess[str], diagnostic: str
) -> None:
    output = result.stdout + result.stderr
    if result.returncode == 0:
        fail(f"{label} unexpectedly reached dependent reader use")
    if diagnostic.lower() not in output.lower():
        fail(f"{label} failed without expected early diagnostic: {diagnostic}")
    termination = f"parallel_abort(): {diagnostic}"
    if termination.lower() not in output.lower():
        fail(f"{label} did not use the expected reader termination path")
    if "build_net production reader semantics passed" in output:
        fail(f"{label} reported rejection after dependent reader use")


def floats_from_fixed(line: str, width: int) -> tuple[float, ...]:
    values = []
    for offset in range(0, len(line), width):
        field = line[offset : offset + width].strip()
        if field:
            values.append(float(field.replace("D", "E")))
    return tuple(values)


def parse_netsu(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines:
        fail("netsu is empty")
    header = tuple(int(value) for value in lines[0].split())
    if len(header) != 2:
        fail("netsu count header is malformed")
    records = []
    chapter = 0
    index = 1
    while index < len(lines):
        if index + 2 >= len(lines):
            fail("netsu record is truncated")
        line = lines[index]
        try:
            marker = int(line[:2])
            q_value = float(line[52:64].replace("D", "E"))
        except ValueError as error:
            raise ContractError(f"netsu header is malformed at line {index + 1}") from error
        names = tuple(line[5 + 5 * item : 10 + 5 * item].strip() for item in range(6))
        descriptor = line[43:47].strip()
        resonance = line[47:48]
        reverse = line[48:49]
        coefficients = floats_from_fixed(lines[index + 1], 13) + floats_from_fixed(lines[index + 2], 13)
        if len(coefficients) != 7:
            fail(f"netsu coefficient count is wrong at line {index + 2}")
        if marker > 0:
            chapter = marker
        records.append((chapter, marker, names, descriptor, resonance, reverse, q_value, coefficients))
        index += 3
    return {"counts": header, "records": tuple(records)}


def parse_netwinv(path: Path) -> dict[str, object]:
    lines = path.read_text(encoding="ascii").splitlines()
    try:
        count = int(lines[0])
    except (IndexError, ValueError) as error:
        raise ContractError("netwinv header is malformed") from error
    species = tuple(line.strip() for line in lines[2 : 2 + count])
    records = []
    index = 2 + count
    for expected_name in species:
        if index + 3 >= len(lines):
            fail("netwinv nuclear record is truncated")
        fields = lines[index].split()
        if len(fields) != 6:
            fail(f"netwinv nuclear record is malformed for {expected_name}")
        name, aa, zz, nn, spin, mass = fields
        values = tuple(
            float(value.replace("D", "E"))
            for line in lines[index + 1 : index + 4]
            for value in line.split()
        )
        if len(values) != 24:
            fail(f"netwinv partition count is wrong for {expected_name}")
        records.append((name, float(aa), int(zz), int(nn), float(spin), float(mass), values))
        index += 4
    if index != len(lines):
        fail("netwinv has unexpected trailing records")
    return {"count": count, "grid": lines[1], "species": species, "records": tuple(records)}


def parse_netweak(path: Path) -> tuple[tuple[object, ...], ...]:
    lines = path.read_text(encoding="ascii").splitlines()
    records = []
    index = 0
    while index < len(lines):
        header = lines[index]
        if len(header) < 64:
            fail("netweak header is malformed")
        names = (header[5:10].strip(), header[10:15].strip())
        descriptor = header[43:46].strip()
        try:
            q_value = float(header[52:64].replace("D", "E"))
        except ValueError as error:
            raise ContractError("netweak Q value is malformed") from error
        index += 1
        values = []
        while index < len(lines) and len(values) < 286:
            try:
                values.extend(float(value) for value in lines[index].split())
            except ValueError as error:
                raise ContractError("netweak table is malformed") from error
            index += 1
        if len(values) != 286:
            fail("netweak table is truncated")
        records.append((names, descriptor, q_value, tuple(values)))
    return tuple(records)


def close(actual: float, expected: float, tolerance: float = 5.0e-6) -> bool:
    return math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance)


def verify_output(directory: Path, *, weak_enabled: bool) -> tuple[object, ...]:
    sunet = tuple(line.strip() for line in (directory / "sunet").read_text(encoding="ascii").splitlines())
    if sunet != EXPECTED_SPECIES:
        fail(f"sunet species/order {sunet} != {EXPECTED_SPECIES}")
    if (directory / "netneutr").exists():
        fail("netneutr was created while private neutrino data was disabled")

    partition = parse_netwinv(directory / "netwinv")
    if partition["species"] != EXPECTED_SPECIES:
        fail("netwinv species/order is inconsistent with sunet")
    if partition["count"] != len(EXPECTED_SPECIES):
        fail("netwinv count is inconsistent with sunet")
    for record, expected in zip(partition["records"], MASTER_SPECIES):
        name, aa, zz, nn, spin, mass, values = record
        expected_name, expected_aa, expected_zz, expected_nn, expected_spin, expected_mass = expected
        if name != expected_name or name not in EXPECTED_SPECIES:
            fail(f"unexpected netwinv species record {name}")
        if (aa, zz, nn, spin) != (expected_aa, expected_zz, expected_nn, expected_spin):
            fail(f"wrong nuclear identity for {name}")
        if not close(mass, expected_mass, 5.0e-8):
            fail(f"wrong copied mass for {name}: {mass}")
        if values != (1.0,) * 24:
            fail(f"wrong partition values for {name}")

    rates = parse_netsu(directory / "netsu")
    weak_records = parse_netweak(directory / "netweak")
    expected_weak_count = 2 if weak_enabled else 0
    if rates["counts"] != (expected_weak_count, 0):
        fail(f"wrong netsu weak/neutrino counts: {rates['counts']}")
    if len(weak_records) != expected_weak_count:
        fail("netweak record count is inconsistent with netsu")

    selected_rates = []
    for chapter, marker, names, descriptor, resonance, reverse, q_value, coefficients in rates["records"]:
        participants = tuple(name for name in names if name)
        if marker == 0 and any(name not in EXPECTED_SPECIES for name in participants):
            fail(f"retained reaction has an absent participant: {participants}")
        if marker == 0:
            selected_rates.append((chapter, participants, descriptor, resonance, reverse, q_value, coefficients))

    strong = {(record[0], record[1], record[2]): record for record in selected_rates if record[2] == "syn1"}
    reverse_key = (2, ("o16", "he4", "c12"), "syn1")
    forward_key = (4, ("he4", "c12", "o16"), "syn1")
    if set(strong) != {reverse_key, forward_key}:
        fail(f"wrong retained forward/reverse reactions: {tuple(strong)}")
    if (
        strong[reverse_key][3].strip()
        or strong[forward_key][3].strip()
        or strong[reverse_key][4] != "v"
        or strong[forward_key][4].strip()
    ):
        fail("forward/reverse flags are wrong")
    if not close(strong[reverse_key][5], -FORWARD_Q) or not close(strong[forward_key][5], FORWARD_Q):
        fail("strong reaction Q values are inconsistent with copied masses")
    if any(record[2] == "drop" or "ne20" in record[1] for record in selected_rates):
        fail("reaction with an unrequested participant was retained")

    weak_rates = [record for record in selected_rates if record[2] == "ffn"]
    if len(weak_rates) != expected_weak_count:
        fail("netsu weak record count is inconsistent with netweak")
    if weak_enabled:
        expected_pairs = {("p", "n"): WEAK_Q, ("n", "p"): -WEAK_Q}
        for names, descriptor, q_value, values in weak_records:
            if descriptor != "ecr" or names not in expected_pairs:
                fail(f"unexpected netweak record {names}")
            if not close(q_value, expected_pairs[names]):
                fail(f"weak Q value is inconsistent with copied masses for {names}")
            if len(values) != 286 or not all(math.isfinite(value) for value in values):
                fail(f"invalid netweak table for {names}")

    return (
        sunet,
        partition,
        rates,
        weak_records,
    )


def expect_verifier_failure(label: str, directory: Path, *, weak_enabled: bool) -> None:
    try:
        verify_output(directory, weak_enabled=weak_enabled)
    except ContractError:
        return
    fail(f"semantic verifier accepted controlled corruption: {label}")


def copy_output(source: Path, destination: Path) -> None:
    shutil.copytree(source, destination)


def mutate_reaction_participant(directory: Path) -> None:
    path = directory / "netsu"
    lines = path.read_text(encoding="ascii").splitlines()
    for index, line in enumerate(lines):
        if line[43:47].strip() == "syn1" and line[5:10].strip() == "he4":
            lines[index] = line[:15] + f"{'ne20':>5}" + line[20:]
            path.write_text("\n".join(lines) + "\n", encoding="ascii")
            return
    fail("could not locate forward reaction for mutation")


def mutate_species_order(directory: Path) -> None:
    path = directory / "netwinv"
    lines = path.read_text(encoding="ascii").splitlines()
    lines[2], lines[3] = lines[3], lines[2]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def mutate_species_tail_order(directory: Path) -> None:
    path = directory / "netwinv"
    lines = path.read_text(encoding="ascii").splitlines()
    lines[int(lines[0]) + 1] = f"{'ne20':>5}"
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def truncate_netwinv(directory: Path, *, after_names: bool) -> None:
    path = directory / "netwinv"
    lines = path.read_text(encoding="ascii").splitlines()
    retained = 2 + int(lines[0]) if after_names else 1
    path.write_text("\n".join(lines[:retained]) + "\n", encoding="ascii")


def mutate_sunet_count(directory: Path) -> None:
    path = directory / "sunet"
    lines = path.read_text(encoding="ascii").splitlines()
    lines.append("ne20")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    truncate_netwinv(directory, after_names=False)


def mutate_sunet_order(directory: Path) -> None:
    path = directory / "sunet"
    lines = path.read_text(encoding="ascii").splitlines()
    lines[0], lines[1] = lines[1], lines[0]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    truncate_netwinv(directory, after_names=True)


def mutate_sunet_tail_order(directory: Path) -> None:
    path = directory / "sunet"
    lines = path.read_text(encoding="ascii").splitlines()
    lines[-1] = "ne20"
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    truncate_netwinv(directory, after_names=True)


def mutate_sunet_padding(directory: Path) -> None:
    path = directory / "sunet"
    lines = path.read_text(encoding="ascii").splitlines()
    path.write_text("\n".join(line.strip().ljust(5) for line in lines) + "\n", encoding="ascii")


def mutate_netwinv_count(directory: Path) -> None:
    path = directory / "netwinv"
    lines = path.read_text(encoding="ascii").splitlines()
    lines[0] = f"{int(lines[0]) + 1:5d}"
    path.write_text(lines[0] + "\n", encoding="ascii")


def mutate_netwinv_order(directory: Path) -> None:
    mutate_species_order(directory)
    truncate_netwinv(directory, after_names=True)


def mutate_netwinv_tail_order(directory: Path) -> None:
    mutate_species_tail_order(directory)
    truncate_netwinv(directory, after_names=True)


def mutate_mass(directory: Path) -> None:
    path = directory / "netwinv"
    lines = path.read_text(encoding="ascii").splitlines()
    count = int(lines[0])
    index = 2 + count
    for _ in range(count):
        if lines[index][:5].strip() == "c12":
            lines[index] = lines[index][:31] + f"{99.0:15.8f}"
            path.write_text("\n".join(lines) + "\n", encoding="ascii")
            return
        index += 4
    fail("could not locate c12 mass for mutation")


def mutate_q_value(directory: Path) -> None:
    path = directory / "netsu"
    lines = path.read_text(encoding="ascii").splitlines()
    for index, line in enumerate(lines):
        if line[43:47].strip() == "syn1" and line[5:10].strip() == "he4":
            lines[index] = line[:52] + f"{FORWARD_Q + 1.0:12.5E}" + line[64:]
            path.write_text("\n".join(lines) + "\n", encoding="ascii")
            return
    fail("could not locate forward Q value for mutation")


def write_smoke_inputs(directory: Path, helm_table: Path) -> None:
    (directory / "initial_abundances").write_text(
        "synthetic he4 state\n"
        "n 0.0 p 0.0 he4 2.5e-1 c12 0.0\n"
        "o16 0.0\n",
        encoding="ascii",
    )
    (directory / "thermo").write_text(
        "two-point constant one-zone history\n"
        "0.0\n"
        "1.0e-10\n"
        "-1.0e-10\n"
        "0.0 1.0e-2 1.0e4 0.5\n"
        "1.0e-10 1.0e-2 1.0e4 0.5\n",
        encoding="ascii",
    )
    species_line = "".join(f"{name:>5}" for name in EXPECTED_SPECIES)
    count_line = "# Species to output in ASCII output (format 14a5):" + f"{len(EXPECTED_SPECIES):4d}"
    control = f"""## Problem Description
build_net interoperability smoke
one zone and one short interval
structural completion only
## Job Controls
1
1
1
0
0
## Neutrinos
0
## NSE Initial Conditions
11.0
## Integration Controls
1
20
5
4
0
1.0e-1
1.0e-7
1.0e-6
1.0e-4
1.0e-30
2.0
## Self-heating Controls
0
1.0e-2
1.0e-4
## Zone Batching Controls
1
## Output Controls
0
0
# ASCII output filename root, network will append zone number
ev_build_net_
# Binary output filename root, network will append zone number
ts_build_net_
{count_line}
{species_line}
## Input Controls
# Nuclear Data Directory
.
# Initial Abundance and Thermodynamic Trajectory Files
initial_abundances
thermo
"""
    (directory / "control").write_text(control, encoding="ascii")
    (directory / "helm_table.dat").symlink_to(helm_table)


def verify_smoke(directory: Path, result: subprocess.CompletedProcess[str]) -> None:
    require_success("one-zone XNet smoke", result)
    diagnostics = sorted(directory.glob("net_diag[0-9]*"))
    if len(diagnostics) != 1 or diagnostics[0].stat().st_size == 0:
        fail("one-zone XNet smoke did not produce one nonempty diagnostic")
    text = diagnostics[0].read_text(encoding="ascii")
    match = re.search(
        r"^End\s+(\d+)\s+(\d+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)",
        text,
        flags=re.MULTILINE,
    )
    if match is None:
        fail("one-zone XNet smoke has no final End record")
    zone, steps = int(match.group(1)), int(match.group(2))
    target, achieved = float(match.group(3)), float(match.group(4))
    if zone != 1 or steps < 1:
        fail(f"one-zone XNet smoke did not evolve: zone={zone}, steps={steps}")
    if not close(target, 1.0e-10, 1.0e-18) or not close(achieved, target, 1.0e-18):
        fail(f"one-zone XNet smoke did not reach its target: {achieved} of {target}")
    if "Counters:" not in text:
        fail("one-zone XNet smoke has no completion counters")


def main(argv: list[str]) -> int:
    if len(argv) != 8:
        print(
            "usage: run_build_net_contracts.py WORK_DIR BUILD_NET NET_SETUP "
            "READER_CHECK INPUT_MUTATOR XNET HELM_TABLE",
            file=sys.stderr,
        )
        return 2
    work_dir = Path(argv[1]).resolve()
    build_net = Path(argv[2]).resolve()
    net_setup = Path(argv[3]).resolve()
    reader_check = Path(argv[4]).resolve()
    input_mutator = Path(argv[5]).resolve()
    xnet = Path(argv[6]).resolve()
    helm_table = Path(argv[7]).resolve()
    if "/build/" not in str(work_dir) or work_dir.name != "build-net-work":
        fail(f"refusing unexpected build_net work directory: {work_dir}")
    for required in (build_net, net_setup, reader_check, input_mutator, xnet, helm_table):
        if not required.is_file():
            fail(f"required executable/input does not exist: {required}")
    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True)

    positive_semantics = []
    for name in ("positive-a", "positive-b"):
        case = work_dir / name
        prepare_case(case)
        result = run_process([str(build_net)], case)
        require_success(name, result)
        positive_semantics.append(verify_output(case / "out", weak_enabled=True))
    if positive_semantics[0] != positive_semantics[1]:
        fail("two identical build_net inputs produced different semantic output")

    weak_off = work_dir / "weak-off"
    prepare_case(weak_off, weak_enabled=False)
    require_success("weak-data-off build", run_process([str(build_net)], weak_off))
    verify_output(weak_off / "out", weak_enabled=False)

    negative_cases = (
        ("duplicate-request", {"requested": ("N01", "H1", "p", "HE4", "C12", "O16")}),
        ("blank-request", {"requested": ("N01", "H1", "", "C12", "O16")}),
        ("unavailable-request", {"requested": ("N01", "H1", "HE4", "C12", "FE56")}),
        ("unknown-request", {"requested": ("N01", "H1", "HE4", "C12", "XY99")}),
        ("malformed-namelist", {"malformed_namelist": True}),
        ("malformed-reaclib-initial", {"malformed_reaclib": True}),
        ("malformed-reaclib-late", {"malformed_reaclib_late": True}),
        ("malformed-mass-hash-record", {"malformed_mass_hash_record": True}),
        ("missing-required-mass", {"missing_mass": True}),
        ("missing-enabled-weak", {"missing_weak": True}),
    )
    for name, options in negative_cases:
        case = work_dir / name
        prepare_case(case, **options)
        require_failure(name, run_process([str(build_net)], case))

    selected_missing = work_dir / "selected-unavailable-mass"
    prepare_case(selected_missing, selected_unavailable_mass=True)
    require_failure_diagnostic(
        "selected-unavailable-mass",
        run_process([str(build_net)], selected_missing),
        "Mass source has no value for c12",
    )

    positive_output = work_dir / "positive-a" / "out"
    mutations = (
        ("absent reaction participant", mutate_reaction_participant),
        ("inconsistent species order", mutate_species_order),
        ("corrupt copied mass", mutate_mass),
        ("corrupt reaction Q value", mutate_q_value),
    )
    for index, (label, mutation) in enumerate(mutations, start=1):
        destination = work_dir / f"mutation-{index}"
        copy_output(positive_output, destination)
        mutation(destination)
        expect_verifier_failure(label, destination, weak_enabled=True)

    downstream = work_dir / "downstream"
    copy_output(positive_output, downstream)
    setup_result = run_process([str(net_setup)], downstream, input_text="build_net synthetic fixture\n")
    require_success("build_net to net_setup", setup_result)
    missing = [name for name in GENERATED_ARTIFACTS if not (downstream / name).is_file() or (downstream / name).stat().st_size == 0]
    if missing:
        fail(f"net_setup did not create required artifacts: {', '.join(missing)}")
    reader_result = run_process([str(reader_check), str(downstream)], downstream)
    require_success("build_net production readers", reader_result)

    padding_case = work_dir / "reader-sunet-padding"
    copy_output(downstream, padding_case)
    mutate_sunet_padding(padding_case)
    padding_result = run_process([str(reader_check), str(padding_case)], padding_case)
    require_success("sunet padding preservation", padding_result)

    ascii_mutations = (
        (
            "sunet-count",
            mutate_sunet_count,
            "sunet and netwinv nuclei counts do not match",
        ),
        (
            "sunet-order",
            mutate_sunet_order,
            "sunet and netwinv nuclei order does not match for inuc=1",
        ),
        (
            "sunet-order-tail",
            mutate_sunet_tail_order,
            "sunet and netwinv nuclei order does not match for inuc=5",
        ),
        (
            "netwinv-count",
            mutate_netwinv_count,
            "sunet and netwinv nuclei counts do not match",
        ),
        (
            "netwinv-order",
            mutate_netwinv_order,
            "sunet and netwinv nuclei order does not match for inuc=1",
        ),
        (
            "netwinv-order-tail",
            mutate_netwinv_tail_order,
            "sunet and netwinv nuclei order does not match for inuc=5",
        ),
    )
    for name, mutation, diagnostic in ascii_mutations:
        case = work_dir / f"reader-{name}"
        copy_output(downstream, case)
        mutation(case)
        result = run_process([str(reader_check), str(case)], case)
        require_reader_rejection(name, result, diagnostic)

    binary_mutations = (
        ("nets4-count", "nets4 nuclei count does not match nuclear data"),
        ("nets4-order", "nets4 nuclei order does not match nuclear data for inuc=1"),
        (
            "nets4-order-tail",
            "nets4 nuclei order does not match nuclear data for inuc=5",
        ),
        ("match-header", "Error reading match_data file"),
        (
            "match-count-1",
            "match_data reaction count does not match nets4 for group=1",
        ),
        (
            "match-count-2",
            "match_data reaction count does not match nets4 for group=2",
        ),
        (
            "match-count-3",
            "match_data reaction count does not match nets4 for group=3",
        ),
        (
            "match-count-4",
            "match_data reaction count does not match nets4 for group=4",
        ),
    )
    for mode, diagnostic in binary_mutations:
        case = work_dir / f"reader-{mode}"
        copy_output(downstream, case)
        mutation_result = run_process([str(input_mutator), mode, str(case)], case)
        require_success(f"{mode} mutation", mutation_result)
        result = run_process([str(reader_check), str(case)], case)
        require_reader_rejection(mode, result, diagnostic)

    write_smoke_inputs(downstream, helm_table)
    smoke_result = run_process([str(xnet)], downstream, timeout=60.0)
    verify_smoke(downstream, smoke_result)

    print("build_net construction, preprocessing, and one-zone smoke contracts passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv))
    except (ContractError, OSError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
