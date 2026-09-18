#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 PROVIDER WORK_DIR EXECUTABLE CONTROLS" >&2
  exit 1
fi

provider=$1
work_root=$2
executable=$3
tracked_controls=$4

case "$provider" in
  ma48)
    failure_diagnostic='Error during MA48'
    ;;
  pardiso-mkl)
    failure_diagnostic='PARDISO factorization failed'
    ;;
  *)
    echo "unsupported real sparse provider: $provider" >&2
    exit 1
    ;;
esac

run_case() {
  local mode=$1
  local work_dir="$work_root/$mode"
  local status

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    "$executable" "$mode" .
  ); then
    status=0
  else
    status=$?
  fi
  echo "$provider $mode process status=$status"
  return "$status"
}

run_tracked_controls() {
  local work_dir="$work_root/tracked-controls"
  local status

  mkdir -p "$work_dir"
  cp "$tracked_controls" "$work_dir/sparse_controls.nml"
  if (
    cd "$work_dir"
    XNET_USE_EXISTING_CONTROLS=1 "$executable" base .
  ); then
    status=0
  else
    status=$?
  fi
  echo "$provider tracked-controls process status=$status"
  return "$status"
}

expect_solver_failure() {
  local work_dir="$work_root/solver-failure"
  local log_file="$work_dir/output.log"
  local status

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    "$executable" failure .
  ) >"$log_file" 2>&1; then
    echo "$provider accepted controlled failure input" >&2
    exit 1
  else
    status=$?
  fi
  if ! grep -Fq "$failure_diagnostic" "$log_file"; then
    echo "$provider failure probe did not reach the expected production status path" >&2
    cat "$log_file" >&2
    exit 1
  fi
  echo "$provider controlled-failure process status=$status diagnostic='$failure_diagnostic'"
}

expect_residual_mutation_failure() {
  local work_dir="$work_root/residual-mutation"
  local log_file="$work_dir/output.log"
  local status

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    XNET_SPARSE_REAL_MUTATION=excessive_residual "$executable" base .
  ) >"$log_file" 2>&1; then
    echo "$provider residual check accepted a controlled result perturbation" >&2
    exit 1
  else
    status=$?
  fi
  echo "$provider residual-mutation process status=$status"
}

run_case base
run_case heat
run_tracked_controls
expect_solver_failure
expect_residual_mutation_failure

echo "$provider real-backend sparse contracts passed"
