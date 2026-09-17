#!/bin/sh

# Create config.txt for a new BUILD_DIR or reject incompatible reuse.

set -eu

config_file=$1
shift
config_dir=${config_file%/*}
temporary_file=${config_file}.tmp.$$

cleanup() {
    rm -f "${temporary_file}"
}
trap cleanup EXIT HUP INT TERM

if test -L "${config_dir}"; then
    echo "refusing symlink BUILD_DIR: ${config_dir}" >&2
    exit 2
fi

mkdir -p "${config_dir}"
{
    printf '%s\n' 'XNET_CONFIG_SCHEMA=1'
    for record do
        printf '%s\n' "${record}"
    done
} > "${temporary_file}"

if test -f "${config_file}"; then
    if ! cmp -s "${config_file}" "${temporary_file}"; then
        echo 'incompatible BUILD_DIR configuration; clean this BUILD_DIR or select another BUILD_DIR' >&2
        exit 2
    fi
else
    mv "${temporary_file}" "${config_file}"
fi
