#!/bin/bash

set -u -o pipefail

artifact_root=$1
source_sha=$2
archive_sha256=$3
build_jobs=$4
time_limit=$5
source_archive="${artifact_root}/source.tar"
source_root="${artifact_root}/source"
launcher_status="${artifact_root}/srun.status.txt"

actual_archive_sha256=$(sha256sum "${source_archive}" | cut -d ' ' -f 1)
if [[ ${actual_archive_sha256} != "${archive_sha256}" ]]; then
  echo "Frontier qualification source archive SHA-256 differs before extraction" >&2
  exit 2
fi

archive_commit_sha=$(git get-tar-commit-id < "${source_archive}")
if [[ ${archive_commit_sha} != "${source_sha}" ]]; then
  echo "Frontier qualification source archive commit differs before extraction" >&2
  exit 2
fi

if [[ -e ${source_root} ]]; then
  echo "Frontier qualification source tree already exists before extraction" >&2
  exit 2
fi
if ! mkdir "${source_root}"; then
  echo "Frontier qualification could not create the source tree" >&2
  exit 2
fi
if ! tar -xf "${source_archive}" -C "${source_root}"; then
  echo "Frontier qualification could not extract the verified source archive" >&2
  exit 2
fi
runner="${source_root}/test/qualification/frontier/frontier_qualification.py"
export PYTHONDONTWRITEBYTECODE=1

srun \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${SLURM_CPUS_PER_TASK:-7}" \
  --gpus-per-task=1 \
  --gpu-bind=closest \
  python3 "${runner}" run \
    --source-root="${source_root}" \
    --artifact-root="${artifact_root}" \
    --source-sha="${source_sha}" \
    --archive-sha256="${archive_sha256}" \
    --build-jobs="${build_jobs}" \
    --time-limit="${time_limit}"
status=$?
printf '%s\n' "${status}" > "${launcher_status}"

if [[ ${status} -ne 0 && ! -f "${artifact_root}/qualification_manifest.json" ]]; then
  python3 "${runner}" failure \
    --artifact-root="${artifact_root}" \
    --source-sha="${source_sha}" \
    --archive-sha256="${archive_sha256}" \
    --time-limit="${time_limit}" \
    --category=allocation \
    --phase=slurm-step \
    --message="srun could not start or complete the qualification step"
fi

exit "${status}"
