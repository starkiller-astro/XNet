#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 6 ]]; then
  echo "usage: $0 WORK_DIR DENSE_EXE MA48_EXE PARDISO_EXE PARDISO_MKL_EXE CONTROLS" >&2
  exit 1
fi

work_root=$1
dense_exe=$2
ma48_exe=$3
pardiso_exe=$4
pardiso_mkl_exe=$5
tracked_controls=$6

run_case() {
  local name=$1
  local executable=$2
  local mode=$3
  local work_dir="$work_root/$name-$mode"

  mkdir -p "$work_dir"
  (
    cd "$work_dir"
    "$executable" "$mode" .
  )
}

expect_failure() {
  local name=$1
  local executable=$2
  local failure=$3
  local diagnostic=$4
  local work_dir="$work_root/failure-$name-$failure"
  local log_file="$work_dir/output.log"

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    XNET_SPARSE_FAIL="$failure" "$executable" base .
  ) >"$log_file" 2>&1; then
    echo "$name accepted injected failure $failure" >&2
    exit 1
  fi
  if ! grep -Fq "$diagnostic" "$log_file"; then
    echo "$name failure $failure did not produce an actionable solver diagnostic" >&2
    cat "$log_file" >&2
    exit 1
  fi
}

expect_input_failure() {
  local name=$1
  local executable=$2
  local mutation=$3
  local diagnostic=$4
  local work_dir="$work_root/input-failure-$name-$mutation"
  local log_file="$work_dir/output.log"

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    XNET_SPARSE_INPUT_MUTATION="$mutation" "$executable" base .
  ) >"$log_file" 2>&1; then
    echo "$name accepted malformed sparse_ind input $mutation" >&2
    exit 1
  fi
  if ! grep -Fq "$diagnostic" "$log_file"; then
    echo "$name sparse_ind input $mutation did not fail in the reader" >&2
    cat "$log_file" >&2
    exit 1
  fi
}

expect_mutation_failure() {
  local name=$1
  local executable=$2
  local mode=$3
  local mutation=$4
  local work_dir="$work_root/mutation-$name-$mutation"
  local log_file="$work_dir/output.log"

  mkdir -p "$work_dir"
  if (
    cd "$work_dir"
    XNET_SPARSE_MUTATION="$mutation" "$executable" "$mode" .
  ) >"$log_file" 2>&1; then
    echo "$name tests did not detect mutation $mutation" >&2
    exit 1
  fi
}

run_tracked_controls() {
  local name=$1
  local executable=$2
  local work_dir="$work_root/tracked-controls-$name"

  mkdir -p "$work_dir"
  cp "$tracked_controls" "$work_dir/sparse_controls.nml"
  (
    cd "$work_dir"
    XNET_USE_EXISTING_CONTROLS=1 "$executable" base .
  )
}

run_recovery_case() {
  local name=$1
  local executable=$2
  local recovery=$3
  local work_dir="$work_root/recovery-$name-$recovery"

  mkdir -p "$work_dir"
  (
    cd "$work_dir"
    XNET_SPARSE_RECOVERY="$recovery" "$executable" base .
  )
}

run_case dense "$dense_exe" base
run_case dense "$dense_exe" heat
run_case ma48 "$ma48_exe" base
run_case ma48 "$ma48_exe" heat
run_case pardiso "$pardiso_exe" base
run_case pardiso "$pardiso_exe" heat
run_case pardiso-mkl "$pardiso_mkl_exe" base
run_case pardiso-mkl "$pardiso_mkl_exe" heat

run_tracked_controls ma48 "$ma48_exe"
run_tracked_controls pardiso "$pardiso_exe"
run_tracked_controls pardiso-mkl "$pardiso_mkl_exe"

run_recovery_case ma48 "$ma48_exe" ma48_analysis_warning
run_recovery_case ma48 "$ma48_exe" ma48_storage_resize

for provider in ma48 pardiso pardiso-mkl; do
  case "$provider" in
    ma48) input_exe=$ma48_exe ;;
    pardiso) input_exe=$pardiso_exe ;;
    pardiso-mkl) input_exe=$pardiso_mkl_exe ;;
  esac
  expect_input_failure "$provider" "$input_exe" missing-file 'Failed to open sparse_ind file'
  expect_input_failure "$provider" "$input_exe" truncated-header 'Error reading sparse_ind header record'
done

expect_failure pardiso "$pardiso_exe" pardiso_init 'PARDISO initialization failed'
expect_failure pardiso "$pardiso_exe" pardiso_factor 'PARDISO factorization failed'
expect_failure pardiso "$pardiso_exe" pardiso_solve 'PARDISO solve failed'
expect_failure pardiso-mkl "$pardiso_mkl_exe" pardiso_factor 'PARDISO factorization failed'
expect_failure pardiso-mkl "$pardiso_mkl_exe" pardiso_solve 'PARDISO solve failed'
expect_failure ma48 "$ma48_exe" ma48_analysis 'Error during MA48AD'
expect_failure ma48 "$ma48_exe" ma48_storage 'Error during MA48AD'
expect_failure ma48 "$ma48_exe" ma48_factor 'Error during MA48BD'
expect_failure ma48 "$ma48_exe" ma48_singular 'Error during MA48BD'
expect_failure ma48 "$ma48_exe" ma48_solve 'Error during MA48CD'

expect_mutation_failure pardiso "$pardiso_exe" base shifted_row_pointer
expect_mutation_failure pardiso "$pardiso_exe" base wrong_index_base
expect_mutation_failure pardiso "$pardiso_exe" heat missing_self_heating_entry
expect_mutation_failure pardiso "$pardiso_exe" base wrong_reaction_map
expect_mutation_failure pardiso "$pardiso_exe" base incorrect_phase
expect_mutation_failure pardiso "$pardiso_exe" base result_offset
expect_mutation_failure pardiso "$pardiso_exe" base cross_zone_copy
expect_mutation_failure pardiso "$pardiso_exe" base excessive_residual

echo "sparse Jacobian and solver-adapter contracts passed"
