#!/usr/bin/env python3
"""Stage an exact clean XNet commit and submit the Frontier qualification."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Sequence

from frontier_qualification import (
    FrontierFailure,
    MANIFEST_SCHEMA,
    inventory_regular_files,
    validate_manifest,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
TIME_PATTERN = re.compile(r"^(?:\d+-)?\d{1,2}:\d{2}:\d{2}$")
REQUIRED_MODULE_MARKERS = (
    "PrgEnv-cray",
    "rocm",
    "craype-accel-amd-gfx90a",
    "hipfort",
    "cray-python",
)


class SubmissionFailure(RuntimeError):
    def __init__(self, category: str, message: str):
        super().__init__(message)
        self.category = category


def _run(command: Sequence[str], directory: Path) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            cwd=directory,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise SubmissionFailure("environment", f"could not run {command[0]}: {error}") from error


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _preflight(source_root: Path) -> str:
    status = _run(
        ["git", "status", "--porcelain", "--untracked-files=all"], source_root
    )
    if status.returncode != 0:
        raise SubmissionFailure(
            "source", f"could not inspect worktree: {status.stderr.strip()}"
        )
    if status.stdout:
        raise SubmissionFailure(
            "source", "worktree is not clean; commit the exact review candidate before submission"
        )
    revision = _run(["git", "rev-parse", "HEAD"], source_root)
    if revision.returncode != 0:
        raise SubmissionFailure(
            "source", f"could not resolve source SHA: {revision.stderr.strip()}"
        )
    modules = tuple(filter(None, os.environ.get("LOADEDMODULES", "").split(":")))
    missing = [
        marker
        for marker in REQUIRED_MODULE_MARKERS
        if not any(marker in item for item in modules)
    ]
    if missing:
        raise SubmissionFailure(
            "environment", "load required Frontier modules: " + ", ".join(missing)
        )
    for variable in ("ROCM_PATH", "OLCF_HIPFORT_ROOT"):
        if not os.environ.get(variable):
            raise SubmissionFailure("environment", f"module variable {variable} is not set")
    for command in (
        ["ftn", "--version"],
        ["cpp", "--version"],
        ["sbatch", "--version"],
        ["srun", "--version"],
    ):
        completed = _run(command, source_root)
        if completed.returncode != 0:
            raise SubmissionFailure("environment", f"{command[0]} is unavailable")
    return revision.stdout.strip()


def _stage_source(source_root: Path, artifact_root: Path) -> str:
    archive = artifact_root / "source.tar"
    with archive.open("wb") as stream:
        try:
            completed = subprocess.run(
                ["git", "archive", "--format=tar", "HEAD"],
                cwd=source_root,
                stdout=stream,
                stderr=subprocess.PIPE,
                check=False,
            )
        except OSError as error:
            raise SubmissionFailure("source", f"could not archive source: {error}") from error
    if completed.returncode != 0:
        raise SubmissionFailure(
            "source", f"git archive failed: {completed.stderr.decode(errors='replace').strip()}"
        )
    return _sha256(archive)


def _redacted_command(command: Sequence[str]) -> list[str]:
    redacted: list[str] = []
    sensitive_next = False
    for token in command:
        if sensitive_next:
            redacted.append("<supplied-at-submission>")
            sensitive_next = False
        elif token in ("--account", "--reservation"):
            redacted.append(token)
            sensitive_next = True
        else:
            redacted.append(token)
    return redacted


def classify_submission_failure(text: str) -> str:
    normalized = text.lower()
    if "source archive" in normalized or "source tree" in normalized:
        return "source"
    if any(
        marker in normalized
        for marker in ("invalid account", "invalid qos", "invalid partition")
    ):
        return "allocation"
    if any(
        marker in normalized
        for marker in ("socket timed out", "communication connection", "controller")
    ):
        return "facility"
    if any(
        marker in normalized
        for marker in ("cancelled", "deadline", "time limit", "node fail")
    ):
        return "queue"
    return "submission"


def collect_submission_diagnostics(
    completed: subprocess.CompletedProcess[str], artifact_root: Path
) -> str:
    """Combine sbatch client output with the Slurm job streams it redirects."""

    messages = [completed.stdout, completed.stderr]
    for name in ("slurm.stdout.txt", "slurm.stderr.txt"):
        path = artifact_root / name
        if path.is_file():
            try:
                messages.append(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeError):
                continue
    return "\n".join(messages)


def failed_submission_manifest(
    artifact_root: Path,
    *,
    source_sha: str,
    archive_sha256: str,
    category: str,
    message: str,
    job_id: str | None,
    partition: str,
    qos_supplied: bool,
    reservation_supplied: bool,
    cpus_per_task: int,
    time_limit: str,
) -> dict[str, object]:
    """Record a structured envelope when the job exits before its runner starts."""

    now = datetime.now(timezone.utc).isoformat()
    document: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "status": "failed",
        "failure": {
            "category": category,
            "phase": "pre-run",
            "message": message,
        },
        "source": {
            "sha": source_sha,
            "worktree_clean": True,
            "archive_sha256": archive_sha256,
            "archive_commit_sha": source_sha,
            "tree_sha256": None,
            "verified_before_build": False,
        },
        "environment": {},
        "slurm": {},
        "builds": {},
        "inputs": [],
        "checks": {},
        "artifact_inventory": [],
        "started_at_utc": now,
        "finished_at_utc": now,
        "runtime_seconds": 0.0,
    }
    return finalize_submission_manifest(
        document,
        artifact_root,
        job_id=job_id,
        partition=partition,
        qos_supplied=qos_supplied,
        reservation_supplied=reservation_supplied,
        cpus_per_task=cpus_per_task,
        time_limit=time_limit,
    )


def finalize_submission_manifest(
    document: dict[str, object],
    artifact_root: Path,
    *,
    job_id: str | None,
    partition: str,
    qos_supplied: bool,
    reservation_supplied: bool,
    cpus_per_task: int,
    time_limit: str,
) -> dict[str, object]:
    """Add redacted submission evidence after Slurm has finalized every artifact."""

    finalized = dict(document)
    finalized["slurm"] = {
        "job_id": job_id or "unknown",
        "account_supplied": True,
        "partition": partition,
        "qos_supplied": qos_supplied,
        "reservation_supplied": reservation_supplied,
        "nodes": 1,
        "tasks": 1,
        "cpus_per_task": cpus_per_task,
        "gpus_per_task": 1,
        "time_limit": time_limit,
    }
    finalized["artifact_inventory"] = inventory_regular_files(
        artifact_root,
        exclude=(Path("qualification_manifest.json"), Path("source.tar")),
        exclude_trees=(Path("source"),),
    )
    return finalized


def submit(arguments: argparse.Namespace) -> Path:
    source_root = arguments.source_root.resolve()
    artifact_root = arguments.artifact_root.resolve()
    if artifact_root == source_root or source_root in artifact_root.parents:
        raise SubmissionFailure("source", "artifact root must be outside the repository")
    if artifact_root.exists() and (not artifact_root.is_dir() or any(artifact_root.iterdir())):
        raise SubmissionFailure("source", "artifact root must be absent or empty")
    artifact_root.mkdir(parents=True, exist_ok=True)
    source_sha = _preflight(source_root)
    if arguments.expected_sha is not None and source_sha != arguments.expected_sha:
        raise SubmissionFailure(
            "source", f"source SHA {source_sha} differs from expected {arguments.expected_sha}"
        )
    archive_sha256 = _stage_source(source_root, artifact_root)
    job_script = source_root / "test" / "qualification" / "frontier" / "frontier_job.sh"
    if not job_script.is_file():
        raise SubmissionFailure("source", f"job script is missing from archive: {job_script}")

    command = [
        "sbatch",
        "--parsable",
        "--wait",
        "--job-name=xnet-frontier-qualification",
        "--nodes=1",
        "--ntasks=1",
        f"--cpus-per-task={arguments.cpus_per_task}",
        "--gpus-per-task=1",
        "--account",
        arguments.account,
        "--partition",
        arguments.partition,
        "--time",
        arguments.time,
        "--output",
        str(artifact_root / "slurm.stdout.txt"),
        "--error",
        str(artifact_root / "slurm.stderr.txt"),
    ]
    if arguments.qos:
        command.extend(("--qos", arguments.qos))
    if arguments.reservation:
        command.extend(("--reservation", arguments.reservation))
    command.extend(
        (
            str(job_script),
            str(artifact_root),
            source_sha,
            archive_sha256,
            str(arguments.build_jobs),
            arguments.time,
        )
    )
    (artifact_root / "sbatch.command.json").write_text(
        json.dumps(_redacted_command(command), indent=2) + "\n", encoding="utf-8"
    )
    completed = _run(command, source_root)
    (artifact_root / "sbatch.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (artifact_root / "sbatch.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    job_id = (
        completed.stdout.strip().split(";")[0].splitlines()[0]
        if completed.stdout.strip()
        else None
    )
    submission = {
        "submitted_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha": source_sha,
        "source_archive_sha256": archive_sha256,
        "job_id": job_id,
        "account_supplied": True,
        "partition": arguments.partition,
        "qos": arguments.qos,
        "reservation_supplied": arguments.reservation is not None,
        "nodes": 1,
        "tasks": 1,
        "cpus_per_task": arguments.cpus_per_task,
        "gpus_per_task": 1,
        "time_limit": arguments.time,
        "sbatch_return_code": completed.returncode,
    }
    (artifact_root / "submission.json").write_text(
        json.dumps(submission, indent=2) + "\n", encoding="utf-8"
    )
    manifest_path = artifact_root / "qualification_manifest.json"
    if manifest_path.is_file():
        try:
            document = json.loads(manifest_path.read_text(encoding="utf-8"))
            validate_manifest(document, require_pass=False)
            document = finalize_submission_manifest(
                document,
                artifact_root,
                job_id=job_id,
                partition=arguments.partition,
                qos_supplied=arguments.qos is not None,
                reservation_supplied=arguments.reservation is not None,
                cpus_per_task=arguments.cpus_per_task,
                time_limit=arguments.time,
            )
            manifest_path.write_text(
                json.dumps(document, indent=2) + "\n", encoding="utf-8"
            )
            validate_manifest(document, require_pass=False)
        except (OSError, UnicodeError, json.JSONDecodeError, FrontierFailure) as error:
            raise SubmissionFailure(
                "test", f"qualification manifest is invalid: {error}"
            ) from error
        if document["status"] != "passed":
            failure = document["failure"]
            raise SubmissionFailure(
                failure["category"], f"{failure['phase']}: {failure['message']}"
            )
        validate_manifest(document)
    elif completed.returncode != 0:
        combined = collect_submission_diagnostics(completed, artifact_root)
        category = classify_submission_failure(combined)
        document = failed_submission_manifest(
            artifact_root,
            source_sha=source_sha,
            archive_sha256=archive_sha256,
            category=category,
            message="Slurm returned before the qualification runner wrote a manifest",
            job_id=job_id,
            partition=arguments.partition,
            qos_supplied=arguments.qos is not None,
            reservation_supplied=arguments.reservation is not None,
            cpus_per_task=arguments.cpus_per_task,
            time_limit=arguments.time,
        )
        manifest_path.write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8"
        )
        validate_manifest(document, require_pass=False)
        raise SubmissionFailure(
            category,
            "Slurm returned without a qualification manifest; inspect sbatch/slurm logs",
        )
    else:
        raise SubmissionFailure("test", "Slurm returned zero without a qualification manifest")
    if completed.returncode != 0:
        raise SubmissionFailure("test", f"sbatch --wait returned {completed.returncode}")
    return manifest_path


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--account", required=True)
    parser.add_argument("--partition", default="batch")
    parser.add_argument("--qos")
    parser.add_argument("--reservation")
    parser.add_argument("--time", default="00:20:00")
    parser.add_argument("--cpus-per-task", type=int, default=7)
    parser.add_argument("--build-jobs", type=int, default=8)
    parser.add_argument("--expected-sha")
    arguments = parser.parse_args(argv)
    if not TIME_PATTERN.fullmatch(arguments.time):
        parser.error("--time must use Slurm D-HH:MM:SS or HH:MM:SS syntax")
    if arguments.cpus_per_task < 1 or arguments.build_jobs < 1:
        parser.error("CPU and build worker counts must be positive")
    return arguments


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        manifest = submit(arguments)
    except SubmissionFailure as error:
        print(f"Frontier submission failed [{error.category}]: {error}", file=sys.stderr)
        return 1
    print(f"Frontier qualification passed; retained artifact manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
