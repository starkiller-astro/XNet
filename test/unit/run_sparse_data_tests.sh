#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 WORK_DIR SPARSE_DATA_EXE" >&2
  exit 1
fi

work_dir=$1
executable=$2
mutation_log="$work_dir/augmentation-mutation.log"

mkdir -p "$work_dir"
"$executable" "$work_dir"

if XNET_CRS_AUGMENTATION_MUTATION=missing-temperature-entry \
    "$executable" "$work_dir" >"$mutation_log" 2>&1; then
  echo "sparse_ind tests did not detect a missing CRS temperature entry" >&2
  exit 1
fi
if [[ $(grep -Fc '[FAILED]' "$mutation_log") -ne 1 ]] || \
    ! grep -Fq '... self-heating CRS augmentation and remapping [FAILED]' "$mutation_log" || \
    ! grep -Fq '1 test(s) failed' "$mutation_log"; then
  echo "sparse_ind mutation did not fail only the intended self-heating CRS test" >&2
  cat "$mutation_log" >&2
  exit 1
fi

echo "sparse_ind reader and self-heating CRS tests passed"
