#!/bin/bash

# Preprocess XNet's variadic accelerator macros before invoking Cray Fortran.

set -u -o pipefail

compiler=${XNET_CRAY_FTN:-ftn}
preprocessor=${XNET_CPP:-cpp}
source_file=
source_form=
compiler_arguments=()
preprocessor_arguments=(-P -C -nostdinc)

for argument in "$@"; do
  case "${argument}" in
    -eZ)
      ;;
    -D*|-U*)
      preprocessor_arguments+=("${argument}")
      ;;
    -I*)
      preprocessor_arguments+=("${argument}")
      compiler_arguments+=("${argument}")
      ;;
    *.F90|*.F95|*.F03|*.F08|*.F18|*.FTN)
      if [[ -n ${source_file} ]]; then
        echo "crayftn_cpp.sh accepts one Fortran source per compilation" >&2
        exit 2
      fi
      source_file=${argument}
      source_form=free
      ;;
    *.F|*.FOR)
      if [[ -n ${source_file} ]]; then
        echo "crayftn_cpp.sh accepts one Fortran source per compilation" >&2
        exit 2
      fi
      source_file=${argument}
      source_form=fixed
      ;;
    *)
      compiler_arguments+=("${argument}")
      ;;
  esac
done

if [[ -z ${source_file} ]]; then
  exec "${compiler}" "$@"
fi

if [[ -n ${XNET_CPP_OUTPUT:-} ]]; then
  temporary_directory=${XNET_CPP_OUTPUT%/*}
  if [[ ${temporary_directory} == "${XNET_CPP_OUTPUT}" ]]; then
    temporary_directory=.
  fi
  mkdir -p -- "${temporary_directory}"
elif [[ -n ${XNET_CPP_OUTPUT_DIR:-} ]]; then
  temporary_directory=${XNET_CPP_OUTPUT_DIR}
  mkdir -p -- "${temporary_directory}"
else
  temporary_directory=$(mktemp -d "${TMPDIR:-/tmp}/xnet-crayftn-cpp.XXXXXX")
  trap 'rm -rf -- "${temporary_directory}"' EXIT HUP INT TERM
fi
source_basename=${source_file##*/}
source_stem=${source_basename%.*}
preprocessed_source="${temporary_directory}/${source_stem}.f90"
raw_preprocessed_source="${temporary_directory}/raw-${source_stem}.f90"
if [[ ${source_form} == fixed ]]; then
  preprocessed_source="${temporary_directory}/${source_stem}.f"
  raw_preprocessed_source="${temporary_directory}/raw-${source_stem}.f"
fi

if [[ -n ${XNET_CPP_OUTPUT:-} ]]; then
  preprocessed_source="${XNET_CPP_OUTPUT}.filtered.$$"
  raw_preprocessed_source="${XNET_CPP_OUTPUT}.raw.$$"
  trap 'rm -f -- "${preprocessed_source}" "${raw_preprocessed_source}"' EXIT HUP INT TERM
fi

if ! "${preprocessor}" "${preprocessor_arguments[@]}" \
    "${source_file}" > "${raw_preprocessed_source}"; then
  echo "crayftn_cpp.sh: preprocessing failed for ${source_file}" >&2
  exit 1
fi

# Drop OpenMP no-op sentinels and reconnect directives across ignored OpenACC present clauses.
if ! awk '
  /^[[:space:]]*!\$omp[[:space:]]+nothing([[:space:]]*&)?[[:space:]]*$/ {
    next
  }
  /^[[:space:]]*!\$omp[[:space:]]*$/ {
    next
  }
  /^[[:space:]]*!\$omp[[:space:]]*&[[:space:]]*!present([[:space:]]*&)?[[:space:]]*$/ {
    if ($0 !~ /&[[:space:]]*$/) {
      sub(/[[:space:]]*&[[:space:]]*$/, "", pending)
    }
    next
  }
  {
    if (have_pending) print pending
    pending = $0
    have_pending = 1
  }
  END {
    if (have_pending) print pending
  }
' "${raw_preprocessed_source}" > "${preprocessed_source}"; then
  echo "crayftn_cpp.sh: OpenMP continuation filtering failed for ${source_file}" >&2
  exit 1
fi

if [[ -n ${XNET_CPP_OUTPUT:-} ]]; then
  mv -- "${preprocessed_source}" "${XNET_CPP_OUTPUT}"
  rm -f -- "${raw_preprocessed_source}"
  exit 0
fi

"${compiler}" "${compiler_arguments[@]}" "${preprocessed_source}"
