# XNet Target Directory Structure

This document describes the target directory layout for the XNet modernization
effort. See `CLAUDE.md` Code Architecture [CURRENT STATE] for the current layout.

## Key differences from current state

- `source/` → `src/` with subdirectories reflecting the layer hierarchy — files
  are grouped by what they do, not piled flat
- `tools/` split into `extern/` + `tools/` — vendored build deps (LAPACK,
  starkiller-helmholtz) are not "tools"
- `cmake/` replaces `Makefile.internal` with composable compiler/site configs
- Out-of-source `build/` — source tree stays clean
- `test/` reorganized — `Data_*` → `test/networks/`, `Test_Problems/` →
  `test/problems/`, space for `test/unit/`
- `xnet_jacobian_common.F90` — shared physics extracted from the 4 Jacobian variants
- File names may be subject to change with further code restructuring.

## Target Layout

```
XNet/
  ├── src/                             # All Fortran source (library + programs)
  │   ├── core/                        # Foundation (Layers 0-2): no physics dependencies
  │   │   ├── xnet_types.F90
  │   │   ├── xnet_constants.F90
  │   │   ├── xnet_util.F90
  │   │   ├── xnet_timers.F90
  │   │   ├── xnet_fd.F90
  │   │   ├── xnet_controls.F90
  │   │   ├── xnet_conditions.F90
  │   │   └── xnet_abundances.F90
  │   │
  │   ├── data/                        # Nuclear/reaction data I/O and preprocessing
  │   │   ├── xnet_data.F90
  │   │   ├── xnet_preprocess.F90
  │   │   └── xnet_match.F90
  │   │
  │   ├── physics/                     # Rate evaluation, screening, EOS (Layer 3)
  │   │   ├── xnet_ffn.F90
  │   │   ├── xnet_nnu.F90
  │   │   ├── xnet_screening.F90
  │   │   ├── xnet_flux.F90
  │   │   ├── xnet_nse.F90
  │   │   └── eos/                     # EOS variants (compile-time selected)
  │   │       ├── xnet_eos_starkiller.F90
  │   │       ├── xnet_eos_helm.F90
  │   │       └── xnet_eos_bahcall.F90
  │   │
  │   ├── jacobian/                    # Jacobian build/solve (Layer 4)
  │   │   ├── xnet_jacobian_common.F90 # [NEW] Shared jacobian_build() physics
  │   │   ├── xnet_jacobian_dense.F90
  │   │   ├── xnet_jacobian_MA48.F90
  │   │   ├── xnet_jacobian_PARDISO.F90
  │   │   └── xnet_jacobian_PARDISO_MKL.F90
  │   │
  │   ├── integration/                 # Time integration (Layer 5)
  │   │   ├── xnet_integrate.F90
  │   │   ├── xnet_integrate_be.F90
  │   │   └── xnet_integrate_bdf.F90
  │   │
  │   ├── linalg/                      # Linear algebra abstractions + GPU dispatch
  │   │   └── xnet_linalg.F90
  │   │
  │   ├── parallel/                    # MPI layer (compile-time selected)
  │   │   ├── xnet_parallel.F90
  │   │   └── xnet_parallel_stubs.F90
  │   │
  │   ├── gpu/                         # GPU management + all vendor bindings
  │   │   ├── xnet_gpu.F90
  │   │   ├── xnet_macros.fh
  │   │   ├── cudaf.F90
  │   │   ├── cublasf.F90
  │   │   ├── cusolverf.F90
  │   │   ├── cusparsef.F90
  │   │   ├── hipf.F90
  │   │   ├── hipblasf.F90
  │   │   ├── rocblasf.F90
  │   │   ├── rocsolverf.F90
  │   │   ├── rocsparsef.F90
  │   │   ├── magmaf.F90
  │   │   ├── openaccf.F90
  │   │   └── openmpf.F90
  │   │
  │   ├── driver/                      # Main programs and driver-level I/O
  │   │   ├── xnet_evolve.F90
  │   │   ├── xnet_output.F90
  │   │   ├── model_input_ascii.F90
  │   │   ├── net.F90
  │   │   ├── net_setup.F90
  │   │   └── nse_slice.F90
  │   │
  │   └── CMakeLists.txt               # Source-level build definition
  │
  ├── build/                           # Out-of-source build directory (gitignored)
  │
  ├── cmake/                           # CMake modules
  │   ├── compilers/                   # Compiler-specific flag files (replaces Makefile.internal)
  │   │   ├── GNU.cmake
  │   │   ├── Intel.cmake
  │   │   ├── NVHPC.cmake
  │   │   └── Cray.cmake
  │   ├── sites/                       # Site-specific configuration (Frontier, Aurora, etc.)
  │   └── FindPARDISO.cmake, etc.      # Dependency finders
  │
  ├── test/                            # Test suite
  │   ├── CMakeLists.txt               # CTest integration
  │   ├── regression/                  # Current regression tests (test_xnet.sh + settings)
  │   │   ├── test_xnet.sh
  │   │   └── test_settings_*
  │   ├── networks/                    # Pre-built network data (currently Data_*)
  │   │   ├── Data_alpha/
  │   │   ├── Data_CNO/
  │   │   ├── Data_SN160/
  │   │   └── ...
  │   ├── problems/                    # Test problem definitions (currently Test_Problems/)
  │   │   ├── setup_*
  │   │   └── Results/
  │   └── unit/                        # [NEW] Unit tests (Fortran test files)
  │
  ├── extern/                          # Vendored external dependencies (currently tools/)
  │   ├── LAPACK/                      # NETLIB BLAS/LAPACK subset
  │   └── starkiller-helmholtz/        # Helmholtz EOS
  │
  ├── tools/                           # User-facing utilities (analysis, plotting, network gen)
  │   ├── analysis/
  │   ├── python/
  │   ├── thermo/
  │   ├── initial_abundance/
  │   └── build_net/                   # Network generation
  │
  ├── doc/
  │   ├── XNet_Formatting_Guidelines.md
  │   ├── xnet/
  │   ├── reaclib/
  │   ├── screening/
  │   └── solvers/
  │
  ├── CMakeLists.txt                   # Top-level CMake (calls into src/, test/, extern/)
  ├── Makefile                         # [TRANSITIONAL] Thin wrapper that calls CMake or legacy build
  ├── CLAUDE.md
  └── README.md
```
