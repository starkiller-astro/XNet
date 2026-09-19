"""Focused effectiveness tests for the manual Frontier qualification package."""

from dataclasses import replace
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


FRONTIER_DIRECTORY = Path(__file__).with_name("frontier")
sys.path.insert(0, str(FRONTIER_DIRECTORY))

import frontier_qualification as frontier_module  # noqa: E402
from frontier_qualification import (  # noqa: E402
    CPU_BUILD_VARIABLES,
    GPU_BUILD_VARIABLES,
    MANIFEST_SCHEMA,
    FrontierFailure,
    compare_ascii_endpoints,
    compare_endpoint_states,
    inventory_regular_files,
    load_policy,
    parse_linalg_probe,
    verify_source_binding,
    validate_manifest,
)
from submit_frontier import (  # noqa: E402
    _stage_source,
    classify_submission_failure,
    collect_submission_diagnostics,
    failed_submission_manifest,
    finalize_submission_manifest,
)
from parallel_zones import AsciiEndpoint  # noqa: E402
from xnet_regression import ALPHA_SPECIES, FinalState, SolverCounters  # noqa: E402


POLICY = FRONTIER_DIRECTORY / "comparison_policy.json"
HASH = "a" * 64
SOURCE_SHA = "b" * 40
TEST_POLICY = load_policy(POLICY)


def _state(zone: int) -> FinalState:
    fractions = {species: 0.0 for species in ALPHA_SPECIES}
    fractions.update({"si28": 0.4, "s32": 0.3, "ni56": 0.3})
    return FinalState(
        zone=zone,
        step=10,
        target_time=1.0,
        time=1.0,
        temperature_gk=5.0,
        density=1.0e8,
        electron_fraction=0.5,
        mass_fractions=fractions,
        counters=SolverCounters(10, 20, 20, 21, 21),
    )


def _artifact(path: str) -> dict[str, object]:
    return {"path": path, "size": 1, "sha256": HASH}


def _bounded(value: float, bounds: object) -> dict[str, float]:
    allowed = bounds.atol + bounds.rtol * abs(value)
    return {
        "observed": value,
        "expected": value,
        "difference": 0.0,
        "allowed": allowed,
    }


def _endpoint_evidence(zones: range) -> dict[str, object]:
    return {
        "status": "passed",
        "policy_status": "accepted after review of Frontier job 5230609",
        "maximum_selected_fraction_of_allowed": 0.0,
        "zones": [
            {
                "zone": zone,
                "scalar_differences": {
                    "achieved_time": _bounded(
                        1.0, TEST_POLICY.scalar_fields["achieved_time"]
                    ),
                    "temperature_gk": _bounded(
                        5.0, TEST_POLICY.scalar_fields["temperature_gk"]
                    ),
                    "density": _bounded(
                        1.0e8, TEST_POLICY.scalar_fields["density"]
                    ),
                    "electron_fraction": _bounded(
                        0.5, TEST_POLICY.scalar_fields["electron_fraction"]
                    ),
                },
                "selected_species": [
                    name
                    for name, value in _state(zone).mass_fractions.items()
                    if name in TEST_POLICY.anchors
                    or value >= TEST_POLICY.material_threshold
                ],
                "selected_differences": {
                    name: _bounded(value, TEST_POLICY.selected)
                    for name, value in _state(zone).mass_fractions.items()
                    if name in TEST_POLICY.anchors
                    or value >= TEST_POLICY.material_threshold
                },
                "observed_mass_fractions": dict(_state(zone).mass_fractions),
                "reference_mass_fractions": dict(_state(zone).mass_fractions),
                "composition_l1": 0.0,
                "composition_linf": 0.0,
                "composition_linf_species": next(
                    iter(_state(zone).mass_fractions)
                ),
            }
            for zone in zones
        ],
    }


def _ascii_evidence(zones: range, *, nonzero_neutrino: bool = False) -> dict[str, object]:
    return {
        "status": "passed",
        "maximum_fraction_of_allowed": 0.0,
        "reported_only_fields": ["energy_generation_rate"],
        "zones": [
            {
                "zone": zone,
                "field_differences": {
                    "neutrino_loss_rate": {
                        "comparison": "bounded",
                        **_bounded(
                            1.0 if nonzero_neutrino else 0.0,
                            TEST_POLICY.ascii_fields["neutrino_loss_rate"],
                        ),
                    },
                    "timestep": {
                        "comparison": "bounded",
                        **_bounded(
                            1.0, TEST_POLICY.ascii_fields["timestep"]
                        ),
                    },
                    "energy_generation_rate": {
                        "comparison": "reported_only",
                        "observed": 1.0,
                        "expected": 1.0,
                        "difference": 0.0,
                    },
                },
            }
            for zone in zones
        ],
    }


def _executable(artifact: str, link_evidence: str) -> dict[str, object]:
    return {
        "artifact": artifact,
        "size": 1,
        "sha256": HASH,
        "link_evidence": link_evidence,
        "linked_libraries": ["libc.so.6"],
    }


def _build(label: str, variables: dict[str, str], targets: tuple[str, ...]) -> dict[str, object]:
    resolved = {**variables, "FC": "ftn", "LDR": "ftn", "FFLAGS": "-O2", "LDFLAGS": ""}
    resolved["CRAY_OMP_PREPROCESS"] = "yes" if label == "gpu" else "no"
    resolved["XNET_CPP"] = "cpp"
    return {
        "status": "passed",
        "variables": variables,
        "resolved_variables": f"build/{label}/resolved-variables.stdout.txt",
        "resolved_variable_values": resolved,
        "timeout_seconds": 1800.0,
        "runtime_seconds": 1.0,
        "executables": {
            target: _executable(
                f"bin/{target}-{label}", f"build/{label}/{target}-link.stdout.txt"
            )
            for target in targets
        },
    }


def _manifest() -> dict[str, object]:
    required_artifacts = [
        "source-verification.json",
        "environment/compiler.stdout.txt",
        "environment/preprocessor.stdout.txt",
        "environment/rocm.stdout.txt",
        "environment/gpu.stdout.txt",
        "build/cpu/resolved-variables.stdout.txt",
        "build/gpu/resolved-variables.stdout.txt",
        "bin/xnet-cpu",
        "build/cpu/xnet-link.stdout.txt",
        "bin/xnet-gpu",
        "build/gpu/xnet-link.stdout.txt",
        "bin/frontier_gpu_linalg_probe-gpu",
        "build/gpu/frontier_gpu_linalg_probe-link.stdout.txt",
        "runs/gpu-linalg/probe.stdout.txt",
        "runs/partial-batch/cpu/cpu/output",
        "runs/partial-batch/gpu/gpu/output",
        "runs/heat-sn160/cpu/cpu/heat-output",
        "runs/heat-sn160/gpu/gpu/heat-output",
    ]
    command_evidence = {
        "artifact": "",
        "sha256": HASH,
        "runtime_seconds": 1.0,
        "summary": "version 1",
    }
    return {
        "schema": MANIFEST_SCHEMA,
        "status": "passed",
        "failure": None,
        "source": {
            "sha": SOURCE_SHA,
            "worktree_clean": True,
            "archive_sha256": HASH,
            "archive_commit_sha": SOURCE_SHA,
            "tree_sha256": HASH,
            "verified_before_build": True,
        },
        "environment": {
            "modules": [
                "PrgEnv-cray/1",
                "rocm/1",
                "craype-accel-amd-gfx90a",
                "hipfort/1",
                "cray-python/3.12.12",
            ],
            "compiler": {**command_evidence, "artifact": "environment/compiler.stdout.txt"},
            "preprocessor": {**command_evidence, "artifact": "environment/preprocessor.stdout.txt"},
            "rocm": {**command_evidence, "artifact": "environment/rocm.stdout.txt"},
            "gpu": {**command_evidence, "artifact": "environment/gpu.stdout.txt"},
            "gpu_model": "AMD Instinct MI250X",
            "rocm_version": "rocm-1",
            "hipfort_module": "hipfort/1",
        },
        "slurm": {
            "job_id": "123",
            "account_supplied": True,
            "partition": "batch",
            "qos_supplied": False,
            "reservation_supplied": False,
            "nodes": 1,
            "tasks": 1,
            "cpus_per_task": 7,
            "gpus_per_task": 1,
            "time_limit": "00:20:00",
        },
        "builds": {
            "cpu": _build("cpu", CPU_BUILD_VARIABLES, ("xnet",)),
            "gpu": _build(
                "gpu", GPU_BUILD_VARIABLES, ("xnet", "frontier_gpu_linalg_probe")
            ),
        },
        "inputs": [_artifact(f"inputs/input-{index:02d}") for index in range(38)],
        "checks": {
            "gpu_linalg": {
                "status": "passed",
                "device_count": 1,
                "offloaded": True,
                "data_present": True,
                "residual_limit": 1.0e-12,
                "batches": [
                    {"batch": batch, "info": 0, "relative_residual": 0.0}
                    for batch in (1, 2)
                ],
                "timeout_seconds": 60.0,
                "runtime_seconds": 1.0,
                "output_inventory": [_artifact("probe.stdout.txt")],
            },
            "partial_batch": {
                "status": "passed",
                "fixture": "ten distinguishable zones, nzbatchmx=4",
                "zones": list(range(1, 11)),
                "inactive_final_batch_lanes": 2,
                "timeout_seconds": 180.0,
                "cpu_runtime_seconds": 1.0,
                "gpu_runtime_seconds": 1.0,
                "endpoint_comparison": _endpoint_evidence(range(1, 11)),
                "ascii_comparison": _ascii_evidence(range(1, 11)),
                "cpu_output_inventory": [_artifact("cpu/output")],
                "gpu_output_inventory": [_artifact("gpu/output")],
            },
            "heat_sn160": {
                "status": "passed",
                "case": "heat_sn160",
                "zones": list(range(1, 7)),
                "timeout_seconds": 600.0,
                "cpu_runtime_seconds": 1.0,
                "gpu_runtime_seconds": 1.0,
                "endpoint_comparison": _endpoint_evidence(range(1, 7)),
                "ascii_comparison": _ascii_evidence(
                    range(1, 7), nonzero_neutrino=True
                ),
                "nonzero_neutrino_loss_zones": list(range(1, 7)),
                "cpu_output_inventory": [_artifact("cpu/heat-output")],
                "gpu_output_inventory": [_artifact("gpu/heat-output")],
            },
        },
        "artifact_inventory": [_artifact(path) for path in required_artifacts],
        "started_at_utc": "2026-08-10T00:00:00+00:00",
        "finished_at_utc": "2026-08-10T00:01:00+00:00",
        "runtime_seconds": 60.0,
    }


def test_policy_accepts_identity_and_rejects_material_endpoint_perturbation() -> None:
    policy = load_policy(POLICY)
    reference = (_state(1),)
    result = compare_endpoint_states(reference, reference, policy, "identity")
    assert result["status"] == "passed"

    fractions = dict(reference[0].mass_fractions)
    fractions["si28"] += 1.0e-3
    fractions["s32"] -= 1.0e-3
    perturbed = (replace(reference[0], mass_fractions=fractions),)
    with pytest.raises(FrontierFailure, match="endpoint comparison failed"):
        compare_endpoint_states(perturbed, reference, policy, "perturbed")


def test_ascii_policy_reports_energy_and_rejects_bounded_field_difference() -> None:
    policy = load_policy(POLICY)
    reference = (AsciiEndpoint(1, 1.0, 2.0, 3.0),)
    energy_difference = (replace(reference[0], energy_generation_rate=2.0),)
    result = compare_ascii_endpoints(energy_difference, reference, policy, "reported")
    energy = result["zones"][0]["field_differences"]["energy_generation_rate"]
    assert energy["comparison"] == "reported_only"
    assert energy["observed"] == 2.0
    assert energy["expected"] == 1.0
    assert energy["difference"] == 1.0
    assert result["reported_only_fields"] == ["energy_generation_rate"]

    perturbed = (replace(reference[0], timestep=4.0),)
    with pytest.raises(
        FrontierFailure,
        match=r"timestep: observed 4\.00000000e\+00, expected 3\.00000000e\+00",
    ):
        compare_ascii_endpoints(perturbed, reference, policy, "perturbed")


def test_linalg_report_requires_device_offload_success_and_small_residuals() -> None:
    report = "\n".join(
        (
            "XNET_GPU_LINALG device_count    1",
            "XNET_GPU_LINALG device    0",
            "XNET_GPU_LINALG offloaded T",
            "XNET_GPU_LINALG data_present T",
            "XNET_GPU_LINALG batch    1 info    0 relative_residual   1.0E-16",
            "XNET_GPU_LINALG batch    2 info    0 relative_residual   2.0E-16",
            "XNET_GPU_LINALG status passed",
        )
    )
    assert parse_linalg_probe(report)["status"] == "passed"
    with pytest.raises(FrontierFailure, match="stayed on the host"):
        parse_linalg_probe(report.replace("offloaded T", "offloaded F"))
    with pytest.raises(FrontierFailure, match="solve data is not present"):
        parse_linalg_probe(report.replace("data_present T", "data_present F"))
    with pytest.raises(FrontierFailure, match="batch 2"):
        parse_linalg_probe(report.replace("2.0E-16", "2.0E-2"))
    with pytest.raises(FrontierFailure, match="batch 1"):
        parse_linalg_probe(report.replace("batch    1 info    0", "batch    1 info    3"))


def test_inventory_hashes_all_regular_outputs_but_not_staged_links(tmp_path: Path) -> None:
    (tmp_path / "nested").mkdir()
    (tmp_path / "nested" / "output").write_text("result\n", encoding="utf-8")
    (tmp_path / "link").symlink_to(tmp_path / "nested" / "output")
    inventory = inventory_regular_files(tmp_path)
    assert [item["path"] for item in inventory] == ["nested/output"]
    assert inventory[0]["size"] == 7


def test_manifest_validator_requires_complete_success_evidence() -> None:
    manifest = _manifest()
    validate_manifest(manifest)
    manifest["checks"]["partial_batch"]["zones"] = list(range(1, 10))
    with pytest.raises(FrontierFailure, match="partial-batch evidence"):
        validate_manifest(manifest)

    manifest = _manifest()
    del manifest["checks"]["gpu_linalg"]["data_present"]
    with pytest.raises(FrontierFailure, match="linear-algebra evidence"):
        validate_manifest(manifest)

    manifest = _manifest()
    manifest["slurm"]["time_limit"] = "unknown"
    with pytest.raises(FrontierFailure, match="Slurm evidence"):
        validate_manifest(manifest)


def test_manifest_validator_rejects_controlled_false_success_mutants() -> None:
    mutations = (
        lambda manifest: manifest["builds"]["cpu"].update({"executables": {}}),
        lambda manifest: manifest["checks"]["gpu_linalg"].update({"device_count": 0}),
        lambda manifest: manifest["checks"]["gpu_linalg"].update(
            {"residual_limit": 1.0}
        ),
        lambda manifest: manifest["checks"]["gpu_linalg"]["batches"][0].update(
            {"info": 2}
        ),
        lambda manifest: manifest["checks"]["gpu_linalg"]["batches"][1].update(
            {"relative_residual": 1.0}
        ),
        lambda manifest: manifest["checks"]["partial_batch"].update(
            {"status": "failed"}
        ),
        lambda manifest: manifest["checks"]["partial_batch"].update(
            {"inactive_final_batch_lanes": 1}
        ),
        lambda manifest: manifest["checks"]["partial_batch"][
            "endpoint_comparison"
        ].update({"zones": []}),
        lambda manifest: manifest["checks"]["heat_sn160"].update(
            {"nonzero_neutrino_loss_zones": []}
        ),
        lambda manifest: manifest["checks"]["heat_sn160"].update(
            {"status": "failed"}
        ),
        lambda manifest: manifest["checks"]["heat_sn160"]["ascii_comparison"][
            "zones"
        ][0]["field_differences"]["neutrino_loss_rate"].update(
            {"observed": 100.0, "expected": 1.0, "difference": 99.0, "allowed": 99.0}
        ),
        lambda manifest: manifest["checks"]["heat_sn160"]["ascii_comparison"].update(
            {"maximum_fraction_of_allowed": 0.5}
        ),
        lambda manifest: manifest["checks"]["partial_batch"]["endpoint_comparison"][
            "zones"
        ][0]["scalar_differences"]["density"].update(
            {"observed": 1.0e99, "expected": 1.0e8, "difference": 1.0e99, "allowed": 1.0e99}
        ),
        lambda manifest: manifest["checks"]["partial_batch"]["endpoint_comparison"].update(
            {"maximum_selected_fraction_of_allowed": 0.5}
        ),
        lambda manifest: manifest["checks"]["partial_batch"]["endpoint_comparison"][
            "zones"
        ][0].update({"selected_species": ["bogus"]}),
        lambda manifest: manifest["environment"]["modules"].pop(),
        lambda manifest: manifest["environment"]["compiler"].update(
            {"sha256": "c" * 64}
        ),
        lambda manifest: manifest["builds"]["cpu"]["executables"]["xnet"].update(
            {"size": 2}
        ),
        lambda manifest: manifest["artifact_inventory"][0].update(
            {"path": "../outside"}
        ),
        lambda manifest: manifest["artifact_inventory"][0].update({"size": -1}),
        lambda manifest: manifest["artifact_inventory"].append(
            dict(manifest["artifact_inventory"][0])
        ),
        lambda manifest: manifest["artifact_inventory"].pop(0),
    )
    for mutate in mutations:
        manifest = _manifest()
        mutate(manifest)
        with pytest.raises(FrontierFailure):
            validate_manifest(manifest)


def test_staged_archive_is_verified_against_commit_and_extracted_tree(
    tmp_path: Path,
) -> None:
    repository = tmp_path / "repository"
    repository.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
    (repository / "source.txt").write_text("candidate\n", encoding="utf-8")
    (repository / "nested").mkdir()
    (repository / "nested" / "input.txt").write_text("input\n", encoding="utf-8")
    (repository / "input-link").symlink_to("nested/input.txt")
    subprocess.run(["git", "add", "."], cwd=repository, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=Qualification Test",
            "-c",
            "user.email=qualification@example.invalid",
            "commit",
            "-q",
            "-m",
            "candidate",
        ],
        cwd=repository,
        check=True,
    )
    source_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    artifacts = tmp_path / "artifacts"
    artifacts.mkdir()
    archive_sha256 = _stage_source(repository, artifacts)
    assert not (artifacts / "source").exists()
    extracted = artifacts / "source"
    extracted.mkdir()
    subprocess.run(
        ["tar", "-xf", str(artifacts / "source.tar"), "-C", str(extracted)],
        check=True,
    )

    evidence = verify_source_binding(
        extracted, artifacts, source_sha, archive_sha256
    )
    assert evidence["verified_before_build"] is True
    assert evidence["archive_commit_sha"] == source_sha

    with pytest.raises(FrontierFailure, match="archive SHA-256 differs"):
        verify_source_binding(extracted, artifacts, source_sha, "c" * 64)
    with pytest.raises(FrontierFailure, match="archive commit differs"):
        verify_source_binding(extracted, artifacts, "c" * 40, archive_sha256)

    (extracted / "source.txt").write_text("mutated while queued\n", encoding="utf-8")
    with pytest.raises(FrontierFailure, match="extracted source tree differs"):
        verify_source_binding(extracted, artifacts, source_sha, archive_sha256)


def test_submission_finalizes_redacted_resources_and_complete_inventory(
    tmp_path: Path,
) -> None:
    (tmp_path / "qualification_manifest.json").write_text("intermediate\n")
    (tmp_path / "submission.json").write_text("submitted\n")
    (tmp_path / "slurm.stdout.txt").write_text("complete\n")
    source = tmp_path / "source"
    source.mkdir()
    (source / "ignored.o").write_text("object\n")

    finalized = finalize_submission_manifest(
        _manifest(),
        tmp_path,
        job_id="456",
        partition="batch",
        qos_supplied=False,
        reservation_supplied=True,
        cpus_per_task=7,
        time_limit="00:20:00",
    )

    assert finalized["slurm"] == {
        "job_id": "456",
        "account_supplied": True,
        "partition": "batch",
        "qos_supplied": False,
        "reservation_supplied": True,
        "nodes": 1,
        "tasks": 1,
        "cpus_per_task": 7,
        "gpus_per_task": 1,
        "time_limit": "00:20:00",
    }
    assert [item["path"] for item in finalized["artifact_inventory"]] == [
        "slurm.stdout.txt",
        "submission.json",
    ]


def test_failed_manifest_retains_classification_without_false_success() -> None:
    manifest = _manifest()
    manifest.update(
        {
            "status": "failed",
            "failure": {
                "category": "allocation",
                "phase": "slurm-step",
                "message": "GPU unavailable",
            },
            "environment": {},
            "builds": {},
            "inputs": [],
            "checks": {},
        }
    )
    validate_manifest(manifest, require_pass=False)
    with pytest.raises(FrontierFailure, match="status is not passed"):
        validate_manifest(manifest)

    manifest["failure"]["category"] = "invented"
    with pytest.raises(FrontierFailure, match="failure classification"):
        validate_manifest(manifest, require_pass=False)


@pytest.mark.parametrize(
    ("message", "category"),
    (
        ("Invalid account or account/partition combination", "allocation"),
        ("source archive SHA-256 differs before extraction", "source"),
        ("Unable to contact slurm controller", "facility"),
        ("JOB CANCELLED AT DEADLINE", "queue"),
        ("unrecognized option", "submission"),
    ),
)
def test_submission_failures_are_classified(message: str, category: str) -> None:
    assert classify_submission_failure(message) == category


def test_pre_run_failure_uses_slurm_stream_and_writes_structured_manifest(
    tmp_path: Path,
) -> None:
    (tmp_path / "slurm.stderr.txt").write_text(
        "source archive SHA-256 differs before extraction\n", encoding="utf-8"
    )
    completed = subprocess.CompletedProcess(
        ["sbatch"], returncode=2, stdout="123\n", stderr=""
    )
    diagnostics = collect_submission_diagnostics(completed, tmp_path)
    category = classify_submission_failure(diagnostics)
    assert category == "source"

    manifest = failed_submission_manifest(
        tmp_path,
        source_sha=SOURCE_SHA,
        archive_sha256=HASH,
        category=category,
        message="runner did not start",
        job_id="123",
        partition="batch",
        qos_supplied=False,
        reservation_supplied=False,
        cpus_per_task=7,
        time_limit="00:20:00",
    )
    validate_manifest(manifest, require_pass=False)
    assert manifest["failure"]["category"] == "source"


def test_manifest_schema_file_is_versioned_and_matches_runner() -> None:
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads(
        (FRONTIER_DIRECTORY / "manifest.schema.json").read_text(encoding="utf-8")
    )
    assert schema["properties"]["schema"]["const"] == _manifest()["schema"]
    jsonschema.Draft202012Validator(schema).validate(_manifest())

    failed = _manifest()
    failed.update(
        {
            "status": "failed",
            "failure": {
                "category": "source",
                "phase": "pre-run",
                "message": "archive mismatch",
            },
            "environment": {},
            "builds": {},
            "inputs": [],
            "checks": {},
        }
    )
    jsonschema.Draft202012Validator(schema).validate(failed)

    contradictory_failed = dict(failed)
    contradictory_failed["failure"] = None
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(contradictory_failed)

    contradictory_passed = _manifest()
    contradictory_passed["failure"] = {
        "category": "test",
        "phase": "manifest",
        "message": "contradiction",
    }
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(contradictory_passed)


def test_cray_wrapper_preserves_fortran_and_expands_variadic_macros(
    tmp_path: Path,
) -> None:
    if shutil.which("cpp") is None:
        pytest.skip("system C preprocessor is unavailable")
    source = tmp_path / "probe.F90"
    source.write_text(
        '#include "xnet_macros.fh"\n'
        "Integer Function probe()\n"
        "Implicit None\n"
        "!XDIR declare target\n"
        'character(len=*), parameter :: joined = "a" // "b"\n'
        "integer :: first, second\n"
        "!XDIR XENTER_DATA XASYNC(1) &\n"
        "!XDIR XCOPYIN(probe)\n"
        "!XDIR XUPDATE XWAIT(1) &\n"
        "!XDIR XHOST(probe)\n"
        "!XDIR XWAIT(1)\n"
        "!XDIR parallel &\n"
        "!XDIR XPRESENT(probe) &\n"
        "!XDIR XPRIVATE(first,second)\n"
        "!XDIR XLOOP(1) XASYNC(1) &\n"
        "!XDIR XPRESENT(probe)\n"
        "Do first = 1, 1\n"
        "EndDo\n"
        "!XDIR XLOOP_SERIAL(1)\n"
        "Do second = 1, 1\n"
        "EndDo\n"
        "probe = 0\n"
        "End Function probe\n",
        encoding="utf-8",
    )
    capture = tmp_path / "preprocessed.f90"
    compiler = tmp_path / "capture-compiler"
    compiler.write_text(
        "#!/bin/bash\n"
        "for argument in \"$@\"; do source_file=$argument; done\n"
        "cp \"${source_file}\" \"${XNET_CAPTURE}\"\n"
        "printf '%s\\n' \"${source_file##*/}\" > \"${XNET_CAPTURE_NAME}\"\n",
        encoding="utf-8",
    )
    compiler.chmod(0o755)
    environment = os.environ.copy()
    environment.update(
        {
            "XNET_CRAY_FTN": str(compiler),
            "XNET_CAPTURE": str(capture),
            "XNET_CAPTURE_NAME": str(tmp_path / "source-name.txt"),
        }
    )
    wrapper = FRONTIER_DIRECTORY.parents[2] / "make" / "crayftn_cpp.sh"
    completed = subprocess.run(
        [
            str(wrapper),
            "-DXNET_OMP_OL",
            f"-I{FRONTIER_DIRECTORY.parents[2] / 'source'}",
            "-eZ",
            "-c",
            str(source),
            "-o",
            str(tmp_path / "probe.o"),
        ],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert (tmp_path / "source-name.txt").read_text(encoding="utf-8") == "probe.f90\n"
    preprocessed = capture.read_text(encoding="utf-8")
    normalized = "\n".join(" ".join(line.split()) for line in preprocessed.splitlines())
    assert 'joined = "a" // "b"' in preprocessed
    assert preprocessed.index("Implicit None") < preprocessed.index("!$omp declare target")
    assert "!$omp target enter data &" in normalized
    assert "!$omp target update &" in normalized
    assert "!$omp from(probe)" in preprocessed
    assert "!$omp parallel &" in normalized
    assert "!$omp private(first,second)" in preprocessed
    assert "!$omp target teams distribute parallel do simd collapse(1)" in normalized
    assert "!present" not in preprocessed
    assert "nowait" not in preprocessed
    assert "barrier" not in preprocessed
    assert "\n!$omp\n" not in preprocessed
    assert "!$omp nothing" not in preprocessed


def test_cray_gpu_preprocessing_remains_selected_after_mpi_compiler_override(
    tmp_path: Path,
) -> None:
    repository = FRONTIER_DIRECTORY.parents[2]
    hipfort = tmp_path / "hipfort"
    (hipfort / "include" / "hipfort" / "amdgcn").mkdir(parents=True)
    rocm = tmp_path / "rocm"
    rocm.mkdir()
    completed = subprocess.run(
        [
            "make",
            "-C",
            str(repository),
            "--no-print-directory",
            "MACHINE=frontier",
            "PE_ENV=CRAY",
            f"BUILD_DIR={tmp_path / 'frontier-wrapper'}",
            "GPU_MODE=ON",
            "GPU_BACKEND=HIP",
            "GPU_LAPACK_VER=ROCM",
            f"HIPFORT_DIR={hipfort}",
            f"ROCM_DIR={rocm}",
            "OPENMP_OL_MODE=ON",
            "MPI_MODE=ON",
            "print-FC",
            "print-CRAY_OMP_PREPROCESS",
            "print-XNET_CPP",
            "print-LDR",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "FC = ftn" in completed.stdout
    assert "CRAY_OMP_PREPROCESS = yes" in completed.stdout
    assert "XNET_CPP = cpp" in completed.stdout
    assert "LDR = ftn" in completed.stdout


@pytest.mark.parametrize(
    ("label", "variables", "targets"),
    (
        ("cpu", CPU_BUILD_VARIABLES, ("xnet",)),
        ("gpu", GPU_BUILD_VARIABLES, ("xnet", "frontier_gpu_linalg_probe")),
    ),
)
def test_frontier_build_uses_selected_build_name_and_separate_evidence_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    label: str,
    variables: dict[str, str],
    targets: tuple[str, ...],
) -> None:
    source_root = tmp_path / "source-tree"
    artifact_root = tmp_path / "qualification-output"
    (source_root / "source").mkdir(parents=True)
    artifact_root.mkdir()
    for target in targets:
        executable = source_root / "build" / variables["BUILD_NAME"] / "bin" / target
        executable.parent.mkdir(parents=True, exist_ok=True)
        executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        executable.chmod(0o755)

    calls: list[tuple[list[str], Path, Path, str]] = []

    def record_command(
        command: list[str],
        working_directory: Path,
        evidence_directory: Path,
        command_label: str,
        timeout_seconds: float,
    ) -> tuple[subprocess.CompletedProcess[str], float]:
        del timeout_seconds
        calls.append((command, working_directory, evidence_directory, command_label))
        if command_label == "resolved-variables":
            names = [
                argument.removeprefix("print-")
                for argument in command
                if argument.startswith("print-")
            ]
            stdout = "".join(f"{name} = checked-{name}\n" for name in names)
        elif command_label.endswith("-link"):
            stdout = "/lib/libc.so\n"
        else:
            stdout = ""
        return subprocess.CompletedProcess(command, 0, stdout, ""), 0.25

    monkeypatch.setattr(frontier_module, "_record_command", record_command)
    result = frontier_module._build_configuration(
        source_root,
        artifact_root,
        label,
        variables,
        jobs=8,
        targets=targets,
    )

    expected_evidence = artifact_root / "build" / label
    make_calls = [call for call in calls if call[0][0] == "make"]
    assert [call[3] for call in make_calls] == ["clean", "build", "resolved-variables"]
    assert all(f"BUILD_NAME={variables['BUILD_NAME']}" in call[0] for call in make_calls)
    assert all(call[1] == source_root for call in calls)
    assert all(call[2] == expected_evidence for call in calls)
    for target in targets:
        assert (artifact_root / "bin" / f"{target}-{label}").is_file()
    assert result["status"] == "passed"


def test_accelerator_routine_directives_follow_ordered_specification_statements() -> None:
    repository = FRONTIER_DIRECTORY.parents[2]
    source_files = list((repository / "source").glob("*.F90"))
    source_files.extend((repository / "tools" / "starkiller-helmholtz").glob("*.F90"))

    for source_file in source_files:
        lines = source_file.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            if "XROUTINE_SEQ" not in line and "XROUTINE_VECTOR" not in line:
                continue
            for following in lines[index + 1 :]:
                statement = following.strip()
                if not statement or statement.startswith("!") or statement.startswith("#"):
                    continue
                assert not statement.lower().startswith(("use ", "implicit none")), (
                    f"{source_file}:{index + 1}: accelerator routine directive "
                    f"precedes {statement!r}"
                )
                break


def test_accelerator_clauses_use_backend_macros() -> None:
    repository = FRONTIER_DIRECTORY.parents[2]
    source_files = list((repository / "source").glob("*.F90"))
    source_files.extend((repository / "tools" / "starkiller-helmholtz").glob("*.F90"))

    for source_file in source_files:
        for line_number, line in enumerate(
            source_file.read_text(encoding="utf-8").splitlines(), start=1
        ):
            directive = line.strip()
            if not directive.startswith("!XDIR"):
                continue
            clause = directive.removeprefix("!XDIR").lstrip()
            assert not clause.startswith(("ASYNC(", "HOST(", "PRIVATE(")), (
                f"{source_file}:{line_number}: accelerator clause bypasses its "
                "backend macro"
            )


def test_openmp_device_pointer_helpers_query_mapped_addresses() -> None:
    if shutil.which("cpp") is None:
        pytest.skip("system C preprocessor is unavailable")
    repository = FRONTIER_DIRECTORY.parents[2]
    completed = subprocess.run(
        [
            "cpp",
            "-P",
            "-C",
            "-nostdinc",
            "-DXNET_OMP_OL",
            f"-I{repository / 'source'}",
            str(repository / "source" / "xnet_gpu.F90"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.count("omp_get_mapped_ptr( C_LOC( a )") == 3
    assert "use_device_ptr" not in completed.stdout


def test_supported_gpu_batched_factor_and_solve_paths() -> None:
    if shutil.which("cpp") is None:
        pytest.skip("system C preprocessor is unavailable")
    repository = FRONTIER_DIRECTORY.parents[2]
    configurations = (
        (
            ("XNET_CUDA", "XNET_OACC", "XNET_LA_CUBLAS"),
            "cublasDgetrfBatched",
            "cublasDgetrsBatched",
        ),
        (
            ("XNET_HIP", "XNET_OMP_OL", "XNET_LA_ROCM"),
            "hipblasDgetrfStridedBatched",
            "hipblasDgetrsStridedBatched",
        ),
        (
            ("XNET_SYCL", "XNET_OMP_OL", "XNET_LA_ONEMKL"),
            "DGETRF_BATCH_STRIDED",
            "DGETRS_BATCH_STRIDED",
        ),
        (
            ("XNET_CUDA", "XNET_OACC", "XNET_LA_MAGMA"),
            "magma_dgetrf_batched",
            "magma_dgetrs_batched",
        ),
    )

    for macros, factor_call, solve_call in configurations:
        command = ["cpp", "-P", "-C", "-nostdinc", "-DXNET_GPU"]
        command.extend(f"-D{macro}" for macro in macros)
        command.extend(
            [
                f"-I{repository / 'source'}",
                str(repository / "source" / "xnet_linalg.F90"),
            ]
        )
        completed = subprocess.run(
            command, capture_output=True, text=True, check=False
        )
        assert completed.returncode == 0, completed.stderr

        batched_solve = completed.stdout.split(
            "Subroutine LinearSolveBatched(", 1
        )[1]
        batched_solve = batched_solve.split(
            "End Subroutine LinearSolveBatched", 1
        )[0]
        assert "Call LinearSolveBatched_GPU" in batched_solve

        gpu_solve = completed.stdout.split("Subroutine LinearSolveBatched_GPU", 1)[1]
        gpu_solve = gpu_solve.split("End Subroutine LinearSolveBatched_GPU", 1)[0]
        assert "Call LUDecompBatched_GPU" in gpu_solve
        assert "Call LUBksubBatched_GPU" in gpu_solve

        factor = completed.stdout.split("Subroutine LUDecompBatched_GPU", 1)[1]
        factor = factor.split("End Subroutine LUDecompBatched_GPU", 1)[0]
        assert factor_call in factor

        bksub = completed.stdout.split("Subroutine LUBksubBatched_GPU", 1)[1]
        bksub = bksub.split("End Subroutine LUBksubBatched_GPU", 1)[0]
        assert solve_call in bksub

        if "XNET_LA_CUBLAS" in macros:
            assert "C_NULL_PTR" in factor
            assert "C_NULL_PTR" in bksub
        else:
            assert "If ( .not. lpiv ) Call xnet_terminate" in factor
            assert "If ( .not. lpiv ) Call xnet_terminate" in bksub

        if "XNET_LA_ROCM" in macros:
            assert "da_base = dev_ptr( a(1,1) )" in factor
            assert "dipiv_base = dev_ptr( ipiv(1) )" in factor
            assert "stridea = Int( lda * n, C_INT64_T )" in factor
            assert "strideipiv = Int( n, C_INT64_T )" in factor
            assert "db_base = dev_ptr( b(1,1) )" in bksub
            assert "strideb = Int( ldb * nrhs, C_INT64_T )" in bksub
            assert "Call stream_sync( stream )" in bksub

        jacobian_command = command[:-1] + [
            str(repository / "source" / "xnet_jacobian_dense.F90")
        ]
        jacobian = subprocess.run(
            jacobian_command, capture_output=True, text=True, check=False
        )
        assert jacobian.returncode == 0, jacobian.stderr
        assert "LinearSolveBatched_GPU" in jacobian.stdout
        assert "LUDecompBatched_GPU" in jacobian.stdout
        assert "LUBksubBatched_GPU" in jacobian.stdout
        assert "StridedBatched_GPU" not in jacobian.stdout
        assert "ROCM_Strided" not in jacobian.stdout

    source_text = (repository / "source" / "xnet_linalg.F90").read_text(
        encoding="utf-8"
    )
    assert "StridedBatched_GPU" not in source_text
    assert "ROCM_Strided" not in source_text
    jacobian_text = (
        repository / "source" / "xnet_jacobian_dense.F90"
    ).read_text(encoding="utf-8")
    assert "XNET_LA_ROCM" not in jacobian_text
    assert "XNET_LA_ONEMKL" not in jacobian_text


def test_timestep_output_updates_all_device_computed_fields() -> None:
    if shutil.which("cpp") is None:
        pytest.skip("system C preprocessor is unavailable")
    repository = FRONTIER_DIRECTORY.parents[2]
    for backend, expected_update in (
        ("XNET_OMP_OL", "!$omp from(t,t9,rho,tdel,edot,sqnu,y,kmon)"),
        ("XNET_OACC", "!$acc host(t,t9,rho,tdel,edot,sqnu,y,kmon)"),
    ):
        completed = subprocess.run(
            [
                "cpp",
                "-P",
                "-C",
                "-nostdinc",
                "-DXNET_GPU",
                f"-D{backend}",
                f"-I{repository / 'source'}",
                str(repository / "source" / "xnet_output.F90"),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
        timestep_output = completed.stdout.split("Subroutine ts_output", 1)[1]
        timestep_output = timestep_output.split("End Subroutine ts_output", 1)[0]
        assert expected_update in timestep_output


def test_helmholtz_allocatable_lifetime_matches_accelerator_model() -> None:
    if shutil.which("cpp") is None:
        pytest.skip("system C preprocessor is unavailable")
    repository = FRONTIER_DIRECTORY.parents[2]
    completed = subprocess.run(
        [
            "cpp",
            "-P",
            "-C",
            "-nostdinc",
            "-DXNET_GPU",
            "-DXNET_OMP_OL",
            f"-I{repository / 'source'}",
            str(repository / "tools" / "starkiller-helmholtz" / "actual_eos.F90"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    declarations = completed.stdout.split("contains", 1)[0]
    assert "!$omp declare target link(itmax, jtmax, d, t)" in declarations
    assert "!$omp declare target link(ttol, dtol)" in declarations

    omp_eos_type = subprocess.run(
        [
            "cpp",
            "-P",
            "-C",
            "-nostdinc",
            "-DXNET_GPU",
            "-DXNET_OMP_OL",
            f"-I{repository / 'source'}",
            str(repository / "tools" / "starkiller-helmholtz" / "eos_type.F90"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert omp_eos_type.returncode == 0, omp_eos_type.stderr
    declarations = omp_eos_type.stdout.split("contains", 1)[0]
    assert (
        "!$omp declare target link(mintemp, maxtemp, mindens, maxdens)"
        in declarations
    )

    initialization = completed.stdout.split("subroutine actual_eos_init", 1)[1]
    initialization = initialization.split("end subroutine actual_eos_init", 1)[0]
    assert "!$omp target enter data" in initialization
    assert "map(to:itmax, jtmax, d, t)" in initialization
    assert "map(to:ttol, dtol)" in initialization
    assert "!$omp target update" not in initialization
    assert "map(alloc:" not in initialization
    assert "always" not in initialization

    finalization = completed.stdout.split("subroutine actual_eos_finalize", 1)[1]
    finalization = finalization.split("end subroutine actual_eos_finalize", 1)[0]
    assert "!$omp target exit data" in finalization
    assert "map(release:itmax, jtmax, d, t)" in finalization
    assert "map(release:ttol, dtol)" in finalization

    openacc = subprocess.run(
        [
            "cpp",
            "-P",
            "-C",
            "-nostdinc",
            "-DXNET_GPU",
            "-DXNET_OACC",
            f"-I{repository / 'source'}",
            str(repository / "tools" / "starkiller-helmholtz" / "actual_eos.F90"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert openacc.returncode == 0, openacc.stderr
    declarations = openacc.stdout.split("contains", 1)[0]
    assert "!$acc declare create(itmax, jtmax, d, t)" in declarations
    assert "!$acc declare create(ttol, dtol)" in declarations

    openacc_eos_type = subprocess.run(
        [
            "cpp",
            "-P",
            "-C",
            "-nostdinc",
            "-DXNET_GPU",
            "-DXNET_OACC",
            f"-I{repository / 'source'}",
            str(repository / "tools" / "starkiller-helmholtz" / "eos_type.F90"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert openacc_eos_type.returncode == 0, openacc_eos_type.stderr
    declarations = openacc_eos_type.stdout.split("contains", 1)[0]
    assert (
        "!$acc declare create(mintemp, maxtemp, mindens, maxdens)" in declarations
    )

    initialization = openacc.stdout.split("subroutine actual_eos_init", 1)[1]
    initialization = initialization.split("end subroutine actual_eos_init", 1)[0]
    assert "!$acc update" in initialization
    assert "!$acc device(itmax, jtmax, d, t)" in initialization
    assert "!$acc device(ttol, dtol)" in initialization
    assert "!$acc enter data" not in initialization

    finalization = openacc.stdout.split("subroutine actual_eos_finalize", 1)[1]
    finalization = finalization.split("end subroutine actual_eos_finalize", 1)[0]
    assert "!$acc exit data" not in finalization
