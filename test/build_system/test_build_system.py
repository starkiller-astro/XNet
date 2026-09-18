#!/usr/bin/env python3
"""Focused behavior checks for the production GNU Make build."""
from __future__ import annotations

import os
import pathlib
import subprocess
import tempfile
import time
from typing import Mapping, Optional

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source"


def make(
    *arguments: str, environment: Optional[Mapping[str, str]] = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["make", "-C", str(ROOT), "--no-print-directory", *arguments],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=environment,
    )


def require_success(result: subprocess.CompletedProcess[str]) -> None:
    assert result.returncode == 0, result.stdout + result.stderr


def require_failure(result: subprocess.CompletedProcess[str], message: str) -> None:
    assert result.returncode != 0, result.stdout + result.stderr
    assert message in result.stdout + result.stderr, result.stdout + result.stderr


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="xnet-build-system-") as temporary:
        work = pathlib.Path(temporary)
        tools = work / "tools"
        tools.mkdir()
        uname = tools / "uname"
        uname.write_text(
            "#!/bin/sh\ncase $1 in -n) echo generic42.example.invalid;; -s) echo Linux;; esac\n",
            encoding="utf-8",
        )
        uname.chmod(0o755)
        for python_name in ("python", "python3"):
            python = tools / python_name
            python.write_text(
                "#!/bin/sh\necho 'production build invoked Python' >&2\nexit 97\n",
                encoding="utf-8",
            )
            python.chmod(0o755)
        generic_environment = dict(os.environ)
        generic_environment.pop("LMOD_SYSTEM_NAME", None)
        generic_environment.pop("HOSTNAME", None)
        generic_environment["PATH"] = f"{tools}{os.pathsep}{generic_environment['PATH']}"

        auto_base = work / "automatic-names"
        automatic_names = (
            (("PE_ENV=GNU",), "GNU-OPT"),
            (("PE_ENV=GNU", "CMODE=DEBUG"), "GNU-DEBUG"),
            (("PE_ENV=GNU", "MPI_MODE=ON"), "GNU-OPT-MPI"),
            (("PE_ENV=GNU", "OPENMP_MODE=ON"), "GNU-OPT-OPENMP"),
            (
                (
                    "PE_ENV=GNU",
                    "GPU_MODE=ON",
                    "GPU_BACKEND=CUDA",
                    "GPU_LAPACK_VER=CUBLAS",
                    "OPENACC_MODE=ON",
                ),
                "GNU-OPT-CUDA-CUBLAS-OPENACC",
            ),
            (("PE_ENV=GNU", "EOS=BAHCALL"), "GNU-OPT-BAHCALL"),
        )
        for options, expected_name in automatic_names:
            automatic = make(
                f"BUILD_BASE={auto_base}",
                *options,
                "print-BUILD_NAME",
                environment=generic_environment,
            )
            require_success(automatic)
            assert f"BUILD_NAME = {expected_name}" in automatic.stdout
            assert (auto_base / expected_name / "config.txt").is_file()

        ma48_dir = work / "ma48"
        ma48_dir.mkdir()
        (ma48_dir / "MA48.f").write_text("      end\n", encoding="utf-8")
        automatic_ma48 = make(
            f"BUILD_BASE={auto_base}",
            "PE_ENV=GNU",
            "MATRIX_SOLVER=MA48",
            f"MA48_DIR={ma48_dir}",
            "print-BUILD_NAME",
            environment=generic_environment,
        )
        require_success(automatic_ma48)
        assert "BUILD_NAME = GNU-OPT-MA48" in automatic_ma48.stdout

        explicit_name = make(
            f"BUILD_BASE={work / 'named-builds'}",
            "BUILD_NAME=my-debug-build",
            "CMODE=DEBUG",
            "print-BUILD_NAME",
            environment=generic_environment,
        )
        require_success(explicit_name)
        assert "BUILD_NAME = my-debug-build" in explicit_name.stdout

        explicit_directory_path = work / "explicit-directory"
        explicit_directory = make(
            f"BUILD_DIR={explicit_directory_path}",
            "CMODE=DEBUG",
            "print-BUILD_DIR",
            environment=generic_environment,
        )
        require_success(explicit_directory)
        assert f"BUILD_DIR = {explicit_directory_path}" in explicit_directory.stdout

        automatic_reuse_base = work / "automatic-reuse"
        automatic_reuse = make(
            f"BUILD_BASE={automatic_reuse_base}",
            "PE_ENV=GNU",
            "print-BUILD_NAME",
            environment=generic_environment,
        )
        require_success(automatic_reuse)
        changed_unencoded_setting = make(
            f"BUILD_BASE={automatic_reuse_base}",
            "PE_ENV=GNU",
            "EXTRA_FLAGS=-fno-inline",
            "print-BUILD_NAME",
            environment=generic_environment,
        )
        require_failure(
            changed_unencoded_setting, "incompatible BUILD_DIR configuration"
        )

        automatic_debug_clean = make(
            f"BUILD_BASE={auto_base}",
            "PE_ENV=GNU",
            "CMODE=DEBUG",
            "clean",
            environment=generic_environment,
        )
        require_success(automatic_debug_clean)
        assert not (auto_base / "GNU-DEBUG").exists()

        generic = make(
            f"BUILD_DIR={work / 'generic'}",
            "PE_ENV=GNU",
            "print-MACHINE",
            "print-FC",
            "print-ARCHOPT",
            environment=generic_environment,
        )
        require_success(generic)
        assert "MACHINE = generic42" in generic.stdout
        assert "FC = gfortran" in generic.stdout
        assert "ARCHOPT = -march=native" in generic.stdout

        explicit_machine = make(
            f"BUILD_DIR={work / 'manual-host'}",
            "PE_ENV=GNU",
            "MACHINE=manual-host",
            "print-MACHINE",
        )
        require_success(explicit_machine)
        assert "MACHINE = manual-host" in explicit_machine.stdout

        perlmutter_environment = dict(os.environ, LMOD_SYSTEM_NAME="perlmutter")
        perlmutter = make(
            f"BUILD_DIR={work / 'perlmutter'}",
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

        cuda_environment = dict(
            os.environ,
            LMOD_SYSTEM_NAME="perlmutter",
            CUDATOOLKIT_HOME="/facility/cudatoolkit",
            CRAY_CUDATOOLKIT_DIR="/facility/cray-cuda",
        )
        cuda_options = (
            "PE_ENV=CRAY",
            "GPU_MODE=ON",
            "GPU_BACKEND=CUDA",
            "GPU_LAPACK_VER=CUBLAS",
            "OPENACC_MODE=ON",
            "OPENMP_OL_MODE=OFF",
        )
        cuda = make(
            f"BUILD_DIR={work / 'perlmutter-cuda'}",
            *cuda_options,
            "print-CUDA_DIR",
            "print-GPU_INC",
            "print-GPU_LIBDIR",
            environment=cuda_environment,
        )
        require_success(cuda)
        assert "CUDA_DIR = /facility/cudatoolkit" in cuda.stdout
        assert "-I/facility/cudatoolkit/include" in cuda.stdout
        assert "-L/facility/cudatoolkit/lib64" in cuda.stdout

        cuda_environment.pop("CUDATOOLKIT_HOME")
        cray_cuda = make(
            f"BUILD_DIR={work / 'cray-cuda'}",
            *cuda_options,
            "print-CUDA_DIR",
            environment=cuda_environment,
        )
        require_success(cray_cuda)
        assert "CUDA_DIR = /facility/cray-cuda" in cray_cuda.stdout

        cray_openmp_offload = make(
            f"BUILD_DIR={work / 'cray-openmp-offload'}",
            "MACHINE=frontier",
            "PE_ENV=CRAY",
            "GPU_MODE=ON",
            "GPU_BACKEND=CUDA",
            "GPU_LAPACK_VER=CUBLAS",
            "OPENACC_MODE=OFF",
            "OPENMP_OL_MODE=ON",
            "CUDA_DIR=/facility/cuda",
            "print-CRAY_OMP_PREPROCESS",
            "print-COMPILE_FC",
        )
        require_success(cray_openmp_offload)
        assert "CRAY_OMP_PREPROCESS = yes" in cray_openmp_offload.stdout
        assert "COMPILE_FC = ftn" in cray_openmp_offload.stdout

        cray_openacc = make(
            f"BUILD_DIR={work / 'cray-openacc'}",
            "MACHINE=frontier",
            "PE_ENV=CRAY",
            "GPU_MODE=ON",
            "GPU_BACKEND=CUDA",
            "GPU_LAPACK_VER=CUBLAS",
            "OPENACC_MODE=ON",
            "OPENMP_OL_MODE=OFF",
            "CUDA_DIR=/facility/cuda",
            "print-CRAY_OMP_PREPROCESS",
        )
        require_success(cray_openacc)
        assert "CRAY_OMP_PREPROCESS = no" in cray_openacc.stdout

        optimized = work / "gnu-opt"
        products = make(
            f"BUILD_DIR={optimized}",
            "PE_ENV=GNU",
            "MACHINE=generic",
            "-j4",
            "xnet",
            "xnse",
            "net_setup",
            environment=generic_environment,
        )
        require_success(products)
        for executable in ("xnet", "xnse", "net_setup"):
            assert (optimized / "bin" / executable).is_file()

        config = (optimized / "config.txt").read_text(encoding="utf-8")
        assert config.startswith("XNET_CONFIG_SCHEMA=1\n")
        assert "MATRIX_SOLVER=dense\n" in config
        assert "JACOBIAN_SOURCE=" in config
        assert "GPU_TARGET=" not in config
        assert "NVCC" not in config
        assert "PROVIDER_ID" not in config

        executable = optimized / "bin" / "xnet"
        xnet_object = optimized / "obj" / "source" / "net.o"
        before = (executable.stat().st_mtime_ns, xnet_object.stat().st_mtime_ns)
        time.sleep(0.01)
        incremental = make(
            f"BUILD_DIR={optimized}",
            "PE_ENV=GNU",
            "MACHINE=generic",
            "-j4",
            "xnet",
        )
        require_success(incremental)
        after = (executable.stat().st_mtime_ns, xnet_object.stat().st_mtime_ns)
        assert after == before, incremental.stdout + incremental.stderr

        types_object = optimized / "obj" / "source" / "xnet_types.o"
        constants_object = optimized / "obj" / "source" / "xnet_constants.o"
        module_before = (
            types_object.stat().st_mtime_ns,
            constants_object.stat().st_mtime_ns,
            executable.stat().st_mtime_ns,
        )
        time.sleep(0.01)
        module_rebuild = make(
            f"BUILD_DIR={optimized}",
            "PE_ENV=GNU",
            "MACHINE=generic",
            "-j4",
            "-W",
            str(SOURCE / "xnet_types.F90"),
            "xnet",
        )
        require_success(module_rebuild)
        module_after = (
            types_object.stat().st_mtime_ns,
            constants_object.stat().st_mtime_ns,
            executable.stat().st_mtime_ns,
        )
        assert all(new > old for old, new in zip(module_before, module_after))

        macro_before = (xnet_object.stat().st_mtime_ns, executable.stat().st_mtime_ns)
        time.sleep(0.01)
        macro_rebuild = make(
            f"BUILD_DIR={optimized}",
            "PE_ENV=GNU",
            "MACHINE=generic",
            "-j4",
            "-W",
            str(SOURCE / "xnet_macros.fh"),
            "xnet",
        )
        require_success(macro_rebuild)
        macro_after = (xnet_object.stat().st_mtime_ns, executable.stat().st_mtime_ns)
        assert all(new > old for old, new in zip(macro_before, macro_after))

        require_failure(
            make(
                f"BUILD_DIR={optimized}",
                "PE_ENV=GNU",
                "MACHINE=generic",
                "CMODE=DEBUG",
                "xnet",
            ),
            "incompatible BUILD_DIR configuration",
        )
        require_failure(
            make(
                f"BUILD_DIR={work / 'invalid-record'}",
                "FFLAGS=safe_flag\nINJECTED=1",
                "xnet",
            ),
            "cannot represent single quotes or line breaks",
        )

        dense_sparse_selection = make(
            f"BUILD_DIR={work / 'dense-sparse-selection'}",
            "MATRIX_SOLVER=dense",
            "print-SPARSE_SRC",
        )
        require_success(dense_sparse_selection)
        assert "SPARSE_SRC = " in dense_sparse_selection.stdout
        assert "xnet_sparse.F90" not in dense_sparse_selection.stdout

        ma41_dir = work / "ma41"
        ma41_dir.mkdir()
        (ma41_dir / "MA41.f").write_text("      end\n", encoding="utf-8")
        # Isolate sparse-source selection from the separately supplied MA41
        # Jacobian implementation, which is not part of this staged PR.
        ma41_sparse_selection = make(
            f"BUILD_DIR={work / 'ma41-sparse-selection'}",
            "MATRIX_SOLVER=MA41",
            f"MA41_DIR={ma41_dir}",
            f"JAC_SRC={ROOT / 'source' / 'xnet_jacobian_dense.F90'}",
            "print-SPARSE_SRC",
        )
        require_success(ma41_sparse_selection)
        assert "SPARSE_SRC = " in ma41_sparse_selection.stdout
        assert "xnet_sparse.F90" not in ma41_sparse_selection.stdout

        ma48_sparse_selection = make(
            f"BUILD_DIR={work / 'ma48-sparse-selection'}",
            "MATRIX_SOLVER=MA48",
            f"MA48_DIR={ma48_dir}",
            "print-SPARSE_SRC",
        )
        require_success(ma48_sparse_selection)
        assert "xnet_sparse.F90" in ma48_sparse_selection.stdout

        pardiso = make(
            f"BUILD_DIR={work / 'pardiso-mkl'}",
            "MATRIX_SOLVER=PARDISO_MKL",
            "MKL_LIBS=-lmkl_rt",
            "print-JAC_SRC",
            "print-LAPACK_VER",
            "print-SPARSE_SRC",
        )
        require_success(pardiso)
        assert "xnet_jacobian_PARDISO_MKL.F90" in pardiso.stdout
        assert "xnet_sparse.F90" in pardiso.stdout
        assert "LAPACK_VER = MKL" in pardiso.stdout
        require_failure(
            make(
                f"BUILD_DIR={work / 'old-pardiso-name'}",
                "MATRIX_SOLVER=PARDISO",
                "print-JAC_SRC",
            ),
            "unsupported MATRIX_SOLVER 'PARDISO'",
        )
        require_failure(
            make(
                f"BUILD_DIR={work / 'pardiso-netlib'}",
                "MATRIX_SOLVER=PARDISO_MKL",
                "LAPACK_VER=NETLIB",
                "MKL_LIBS=-lmkl_rt",
                "print-JAC_SRC",
            ),
            "requires LAPACK_VER=MKL",
        )
        removed_alias = make(f"BUILD_DIR={work / 'solver-alias'}", "xnet_dense")
        assert removed_alias.returncode != 0, removed_alias.stdout + removed_alias.stderr

        require_failure(
            make(
                f"BUILD_DIR={work / 'inactive-cuda'}",
                "GPU_MODE=OFF",
                "GPU_BACKEND=CUDA",
                "print-GPU_BACKEND",
            ),
            "GPU_BACKEND is active while GPU_MODE=OFF",
        )
        require_failure(
            make(
                f"BUILD_DIR={work / 'hip-cublas'}",
                "GPU_MODE=ON",
                "GPU_BACKEND=HIP",
                "GPU_LAPACK_VER=CUBLAS",
                "OPENMP_OL_MODE=ON",
                "print-GPU_BACKEND",
            ),
            "HIP requires GPU_LAPACK_VER=ROCM",
        )
        require_failure(
            make(
                f"BUILD_DIR={work / 'two-directives'}",
                "GPU_MODE=ON",
                "GPU_BACKEND=CUDA",
                "GPU_LAPACK_VER=CUBLAS",
                "OPENACC_MODE=ON",
                "OPENMP_OL_MODE=ON",
                "print-GPU_BACKEND",
            ),
            "requires exactly one accelerator directive mode",
        )

        concurrent_opt = work / "concurrent-opt"
        concurrent_debug = work / "concurrent-debug"
        common_command = ["make", "-C", str(ROOT), "--no-print-directory"]
        opt_process = subprocess.Popen(
            common_command
            + [
                f"BUILD_DIR={concurrent_opt}",
                "PE_ENV=GNU",
                "MACHINE=generic",
                "-j4",
                "xnet",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        debug_process = subprocess.Popen(
            common_command
            + [
                f"BUILD_DIR={concurrent_debug}",
                "PE_ENV=GNU",
                "MACHINE=generic",
                "CMODE=DEBUG",
                "-j4",
                "xnet",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        opt_stdout, opt_stderr = opt_process.communicate()
        debug_stdout, debug_stderr = debug_process.communicate()
        assert opt_process.returncode == 0, opt_stdout + opt_stderr
        assert debug_process.returncode == 0, debug_stdout + debug_stderr
        assert (concurrent_opt / "bin" / "xnet").is_file()
        assert (concurrent_debug / "bin" / "xnet").is_file()

        clean = make(f"BUILD_DIR={optimized}", "clean")
        require_success(clean)
        assert not optimized.exists()

    print("production GNU Make checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
