#!/usr/bin/env bash
set -euo pipefail

build_dir=$1
reader_test=$2

reader_test=$(cd "$(dirname "$reader_test")" && pwd)/$(basename "$reader_test")

work_dir="$build_dir/work"
rm -rf "$work_dir"
mkdir -p "$work_dir"

(
  cd "$work_dir"
  "$reader_test"
)
