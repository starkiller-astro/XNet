#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 WORK_DIR STARKILLER_EXE BAHCALL_EXE HELM_TABLE" >&2
  exit 2
fi

work_dir=$1
starkiller_exe=$2
bahcall_exe=$3
helm_table=$4

mkdir -p "$work_dir/starkiller" "$work_dir/bahcall"
ln -sf "$helm_table" "$work_dir/starkiller/helm_table.dat"

(
  cd "$work_dir/starkiller"
  "$starkiller_exe"
)
(
  cd "$work_dir/bahcall"
  "$bahcall_exe"
)

missing_table_dir=$(mktemp -d "$work_dir/missing-table.XXXXXX")
if (
  cd "$missing_table_dir"
  "$starkiller_exe"
) >"$missing_table_dir/output.log" 2>&1; then
  echo "STARKILLER initialization unexpectedly succeeded without helm_table.dat" >&2
  exit 1
fi
if ! grep -q "Failed to open helm_table.dat" "$missing_table_dir/output.log"; then
  echo "STARKILLER missing-table failure did not identify helm_table.dat" >&2
  exit 1
fi

echo "EOS provider contracts and missing-table handling passed"
