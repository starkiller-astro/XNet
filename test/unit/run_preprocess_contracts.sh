#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 FIXTURE_DIR WORK_DIR DIRECT_EXE NET_SETUP_EXE VERIFY_EXE" >&2
  exit 2
fi

fixture_dir=$1
work_dir=$2
direct_exe=$3
net_setup_exe=$4
verify_exe=$5

case "$work_dir" in
  */build/*/preprocess-work) ;;
  *) echo "refusing unexpected preprocessing work directory: $work_dir" >&2; exit 2 ;;
esac

rm -rf "$work_dir"
mkdir -p "$work_dir/direct" "$work_dir/standalone" "$work_dir/truncated" \
  "$work_dir/inconsistent" "$work_dir/unreadable"

copy_ascii_fixture() {
  local destination=$1
  cp "$fixture_dir/sunet" "$fixture_dir/netwinv" "$destination/"
}

copy_ascii_fixture "$work_dir/direct"
cp "$fixture_dir/netsu" "$work_dir/direct/netsu"
"$direct_exe" "$work_dir/direct" "preprocess contract fixture"

copy_ascii_fixture "$work_dir/standalone"
cp "$fixture_dir/netsu" "$work_dir/standalone/netsu"
printf '%s\n' "preprocess contract fixture" | (cd "$work_dir/standalone" && "$net_setup_exe")

"$verify_exe" "$work_dir/direct" "$work_dir/direct.summary"
"$verify_exe" "$work_dir/standalone" "$work_dir/standalone.summary"
cmp "$work_dir/direct.summary" "$work_dir/standalone.summary"

copy_ascii_fixture "$work_dir/truncated"
cp "$fixture_dir/netsu_truncated" "$work_dir/truncated/netsu"
if "$direct_exe" "$work_dir/truncated" "truncated reaction fixture"; then
  echo "truncated rate record unexpectedly succeeded" >&2
  exit 1
fi
if printf '%s\n' "truncated standalone fixture" | \
    (cd "$work_dir/truncated" && "$net_setup_exe"); then
  echo "standalone net_setup accepted a truncated rate record" >&2
  exit 1
fi

copy_ascii_fixture "$work_dir/inconsistent"
cp "$fixture_dir/netsu_inconsistent" "$work_dir/inconsistent/netsu"
if "$direct_exe" "$work_dir/inconsistent" "inconsistent participant fixture"; then
  echo "inconsistent participant record unexpectedly succeeded" >&2
  exit 1
fi

cp -R "$work_dir/standalone/." "$work_dir/unreadable/"
mv "$work_dir/unreadable/nets3" "$work_dir/unreadable/nets3.unreadable"
if "$verify_exe" "$work_dir/unreadable" "$work_dir/unreadable.summary"; then
  echo "unreadable generated reaction artifact unexpectedly loaded" >&2
  exit 1
fi

echo "preprocess process-integration contracts passed"
