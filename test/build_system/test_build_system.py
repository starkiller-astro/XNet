#!/usr/bin/env python3
"""Focused checks for the production GNU Make build."""
from __future__ import annotations

import os
import pathlib
import subprocess
import tempfile
from typing import Dict, Optional

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source"
DEPENDENCIES_MAKEFILE = SOURCE / "make" / "dependencies.mk"


def make(
    *arguments: str, environment: Optional[Dict[str, str]] = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["make", "-C", str(SOURCE), "--no-print-directory", *arguments],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
    )


def require_success(result: subprocess.CompletedProcess[str]) -> None:
    assert result.returncode == 0, result.stdout + result.stderr


def main() -> int:
    makefile = DEPENDENCIES_MAKEFILE.read_text(encoding="utf-8")
    sparse_jacobian_rule = next(
        line for line in makefile.splitlines() if "$(call solver_obj,$(JAC_SRC)):" in line
    )
    assert "$(call source_obj,$(MPI_SRC))" in sparse_jacobian_rule
    for name in ("xnet_jacobian_MA48.F90", "xnet_jacobian_PARDISO_MKL.F90"):
        assert "Use xnet_parallel" in (SOURCE / name).read_text(encoding="utf-8")

    with tempfile.TemporaryDirectory(prefix="xnet-build-system-") as temporary:
        build = pathlib.Path(temporary) / "gnu"
        tools = build.parent / "tools"
        tools.mkdir()
        uname = tools / "uname"
        uname.write_text(
            "#!/bin/sh\ncase $1 in -n) echo generic42.example.invalid;; -s) echo Linux;; esac\n",
            encoding="utf-8",
        )
        uname.chmod(0o755)
        hostname_environment = dict(os.environ)
        hostname_environment.pop("HOSTNAME", None)
        hostname_environment["PATH"] = f"{tools}{os.pathsep}{hostname_environment['PATH']}"
        architecture = make(
            f"BUILD_DIR={build.parent / 'hostname'}",
            "PE_ENV=GNU",
            "print-MACHINE",
            "print-ARCHOPT",
            environment=hostname_environment,
        )
        require_success(architecture)
        assert "MACHINE = generic" in architecture.stdout
        assert "ARCHOPT = -march=native" in architecture.stdout

        perlmutter_environment = dict(os.environ, LMOD_SYSTEM_NAME="perlmutter")
        perlmutter = make(
            f"BUILD_DIR={build.parent / 'perlmutter'}",
            "PE_ENV=GNU",
            "print-MACHINE",
            "print-FC",
            "print-LAPACK_VER",
            environment=perlmutter_environment,
        )
        require_success(perlmutter)
        assert "MACHINE = perlmutter" in perlmutter.stdout
        assert "FC = ftn" in perlmutter.stdout
        assert "LAPACK_VER = LIBSCI" in perlmutter.stdout

        summit_environment = dict(os.environ, LMOD_SYSTEM_NAME="summit")
        summit = make(
            f"BUILD_DIR={build.parent / 'summit'}",
            "PE_ENV=GNU",
            "print-MACHINE",
            "print-FC",
            "print-LAPACK_VER",
            "print-ARCHOPT",
            environment=summit_environment,
        )
        require_success(summit)
        assert "MACHINE = summit" in summit.stdout
        assert "FC = gfortran" in summit.stdout
        assert "LAPACK_VER = NETLIB" in summit.stdout
        assert "ARCHOPT = -mtune=native" in summit.stdout

        cori_environment = dict(os.environ, LMOD_SYSTEM_NAME="cori")
        cori = make(
            f"BUILD_DIR={build.parent / 'cori'}",
            "PE_ENV=INTEL",
            "print-MACHINE",
            "print-FC",
            "print-LAPACK_VER",
            "print-ARCHOPT",
            environment=cori_environment,
        )
        require_success(cori)
        assert "MACHINE = cori" in cori.stdout
        assert "FC = ftn" in cori.stdout
        assert "LAPACK_VER = LIBSCI" in cori.stdout
        assert "ARCHOPT = -align array64byte" in cori.stdout

        products = make(f"BUILD_DIR={build}", "-j4", "xnet", "xnse", "net_setup")
        require_success(products)
        assert "python" not in (products.stdout + products.stderr).lower()
        for executable in ("xnet", "xnse", "net_setup"):
            assert (build / "bin" / executable).is_file()

        mismatch = make(f"BUILD_DIR={build}", "CMODE=DEBUG", "xnet")
        assert mismatch.returncode != 0
        assert "incompatible BUILD_DIR configuration" in mismatch.stderr

        invalid_record = make(
            f"BUILD_DIR={build.parent / 'invalid-record'}",
            "FFLAGS=safe_flag\nINJECTED=1",
            "xnet",
        )
        assert invalid_record.returncode != 0
        assert "cannot represent single quotes or line breaks" in invalid_record.stderr

        cuda_selectors = (
            "GPU_MODE=ON",
            "GPU_BACKEND=CUDA",
            "GPU_LAPACK_VER=CUBLAS",
            "OPENACC_MODE=ON",
            "OPENMP_OL_MODE=OFF",
        )
        volta = make(
            f"BUILD_DIR={build.parent / 'cuda-volta'}",
            *cuda_selectors,
            "GPU_TARGET=Volta",
            "print-GPU_TARGET",
            "print-NVCCFLAGS",
        )
        require_success(volta)
        assert "GPU_TARGET = sm70" in volta.stdout
        assert "-gencode arch=compute_70,code=sm_70" in volta.stdout
        assert "compute_80" not in volta.stdout

        ampere = make(
            f"BUILD_DIR={build.parent / 'cuda-ampere'}",
            *cuda_selectors,
            "GPU_TARGET=Ampere",
            "print-GPU_TARGET",
            "print-NVCCFLAGS",
        )
        require_success(ampere)
        assert "GPU_TARGET = sm80" in ampere.stdout
        assert "-gencode arch=compute_80,code=sm_80" in ampere.stdout
        assert "compute_70" not in ampere.stdout

        direct_targets = make(
            f"BUILD_DIR={build.parent / 'cuda-direct-targets'}",
            *cuda_selectors,
            "GPU_TARGET=sm70 sm80",
            "print-GPU_TARGET",
            "print-NVCCFLAGS",
        )
        require_success(direct_targets)
        assert "GPU_TARGET = sm70 sm80" in direct_targets.stdout
        assert "-gencode arch=compute_70,code=sm_70" in direct_targets.stdout
        assert "-gencode arch=compute_80,code=sm_80" in direct_targets.stdout

        invalid_cuda_target = make(
            f"BUILD_DIR={build.parent / 'cuda-invalid-target'}",
            *cuda_selectors,
            "GPU_TARGET=sm90",
            "print-GPU_TARGET",
        )
        assert invalid_cuda_target.returncode != 0
        assert "unsupported CUDA GPU_TARGET 'sm90'" in invalid_cuda_target.stderr

        invalid_hip_target = make(
            f"BUILD_DIR={build.parent / 'hip-invalid-target'}",
            "GPU_MODE=ON",
            "GPU_BACKEND=HIP",
            "GPU_LAPACK_VER=ROCM",
            "OPENACC_MODE=OFF",
            "OPENMP_OL_MODE=ON",
            "GPU_TARGET=Ampere",
            "print-GPU_TARGET",
        )
        assert invalid_hip_target.returncode != 0
        assert "GPU_TARGET is valid only with GPU_BACKEND=CUDA" in invalid_hip_target.stderr

        clean = make(f"BUILD_DIR={build}", "clean")
        require_success(clean)
        assert not build.exists()
    print("production GNU Make checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
