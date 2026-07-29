# CLAUDE.md - XNet

This file provides guidance to Codex and other coding-agent workflows in this repository.

## Project Overview

**What this code does:**
XNet is a thermonuclear reaction network code for astrophysical nucleosynthesis simulations.
It is used primarily in one of two ways:

1. stand-alone post-processing nucleosynthesis calculations, typically from thermodynamic histories sampled by Lagrangian tracer particles;
2. integration with multi-physics simulation packages (e.g. Flash-X, CHIMERA).

The code evolves the abundances of an arbitrary set of nuclei that are linked by reaction rates that form a reaction network.
Time integration is performed using an implicit method --- Backward Euler (BE) or Backward Differentiation Formula (BDF)  --- as the coupled set of ODEs is quite numerically stiff.
Multiple linear algebra solver backends and compiler directives are used for portability and performance.

**Domain:**
Astrophysical reactive flow problems: core-collapse supernovae, thermonuclear supernovae, neutron star mergers, X-ray bursts, novae, etc.

**Scale / deployment context:**
From laptops to leadership-class HPC (Frontier, Aurora, Perlmutter); varying hardware architectures supported via CUDA/ROCm/oneMKL numerical libraries and OpenMP/OpenACC compiler directives.

**Primary language(s):**
Fortran 2003 (source code), GNU Make (build system), Python (analysis/plotting), bash (testing)

## Physics Domain Primer

This section provides enough physics context to write meaningful tests and
review code changes. It is not a textbook — it is a lookup reference.

### Key Physical Quantities

| Symbol | Variable(s) | Meaning | Units |
|--------|-------------|---------|-------|
| Y_i | `y(i)` | Number abundance of species i | mol/g |
| X_i = A_i * Y_i | — | Mass fraction of species i | dimensionless |
| A_i | `aa(i)` | Mass number (nucleons) of species i | — |
| Z_i | `zz(i)` | Proton number (charge) of species i | — |
| N_i = A_i - Z_i | `nn(i)` | Neutron number of species i | — |
| Y_e | `ye` | Electron fraction = sum(Z_i * Y_i) | dimensionless |
| T9 | `t9` | Temperature in units of 10^9 K (GK) | GK |
| rho | `rho` | Mass density | g/cm^3 |
| BE_i | `be(i)` | Nuclear binding energy of species i | MeV |
| dY_i/dt | `ydot(i)` | Time derivative of abundance | mol/g/s |

### Conservation Laws

These must hold at all times. Violations indicate bugs.

1. **Mass conservation**: `sum(A_i * Y_i) = 1` (sum of mass fractions = 1).
   Deviations indicate a bug in `yderiv()` or the Jacobian.
2. **Charge conservation**: `Y_e = sum(Z_i * Y_i)` is conserved in the absence
   of weak interactions (electron capture, beta decay, neutrino reactions).
   When weak rates are active, Y_e changes but must still equal `sum(Z_i * Y_i)`.
3. **Non-negativity**: `Y_i >= 0` always. Negative abundances indicate
   integrator failure or insufficient Newton-Raphson iterations.

### Reaction Types

XNet handles reactions with 1 to 4 reactants. The cross-section arrays
`csect1` through `csect4` correspond to:

- **1-body** (`csect1`): decays, photodisintegrations — (gamma,n), (gamma,p), (gamma,alpha)
- **2-body** (`csect2`): capture reactions — (p,gamma), (n,gamma), (alpha,gamma), (p,n), etc.
- **3-body** (`csect3`): e.g., triple-alpha (3 alpha -> 12C)
- **4-body** (`csect4`): rare; included for completeness

Forward and reverse rates are related by **detailed balance**: at thermodynamic
equilibrium, each reaction is exactly balanced by its reverse. The ratio of
forward to reverse rate is determined by Q-values and partition functions.

### Nuclear Statistical Equilibrium (NSE)

At sufficiently high temperature (T9 > ~5-7) and density, reactions are fast
enough that the composition reaches thermodynamic equilibrium (NSE). In NSE:

- Abundances are determined by T9, rho, Y_e, binding energies, and partition
  functions — not by individual reaction rates.
- The NSE solver (`xnet_nse.F90`) computes this equilibrium composition.
- XNet transitions between NSE and full network integration depending on
  thermodynamic conditions.
- At NSE, the Saha equation determines relative abundances.

### Key Physics Regimes

| Regime | Conditions | What to check |
|--------|-----------|---------------|
| NSE | T9 > 7, high rho | Abundances match Saha equation; Y_e consistent |
| Quasi-equilibrium | T9 ~ 3-6 | QSE clusters form (iron group, Si group) |
| Explosive burning | T9 ~ 2-5, short timescale | Alpha-chain reactions dominate; energy release |
| Freeze-out | Rapidly dropping T9/rho | Composition freezes; final Y_i converge |
| Hydrostatic burning | T9 ~ 0.1-2, long timescale | pp-chain, CNO, He burning |

### Key References

- Hix & Thielemann (1999), ApJ 511, 862 — foundational XNet methods paper
- Seitenzahl et al. (2009), ADNDT 95, 96 — NSE treatment (doi:10.1016/j.adt.2008.08.001)
- Timmes (1999), ApJS 124, 241 — integration benchmarks; also https://cococubed.com
- Longland, Martin & José (2014), A&A 563, A67 — modern treatment of implicit
  integration methods, timestep control, and convergence for reaction networks

See `doc/references.md` for the full reference list including community codes
and all references cited in source-code comments.

## Build System

All builds run from `source/`.
The `XNET_DIR` environment variable should point to the repo root.

```bash
# Default build (dense solver, optimized)
cd source && make

# Debug build
make CMODE=DEBUG

# Specific solver
make xnet MATRIX_SOLVER=dense   # Dense LAPACK solver
make xnet MATRIX_SOLVER=MA48    # HSL sparse solver
make xnet MATRIX_SOLVER=PARDISO # PARDISO sparse solver

# Parallel driver
make xnet MPI_MODE=ON # Use MPI to parallelize over zones/particles in driver

# Supporting tools
make net_setup        # Network setup pre-processing utility
make xnse             # Stand-alone NSE state calculator

# Clean
make clean
```

### Build Configuration

- Edit `source/Makefile.opt` or define with `make` command to change:
	- `CMODE` - OPT (default) or DEBUG
	- `PE_ENV` - Compiler: GNU (default), INTEL, PGI, CRAY, XL
	- `EOS` - STARKILLER (default), HELMHOLTZ, BAHCALL
	- `MATRIX_SOLVER` - dense (default), MA48, MA41, PARDISO, PARDISO_MKL
	- `MPI_MODE` / `OPENMP_MODE` / `GPU_MODE` - ON/OFF (default OFF)
	- `LAPACK_VER` - CPU BLAS/LAPACK: NETLIB (default), MKL, LIBSCI, ACCEL
	- `GPU_LAPACK_VER` - GPU BLAS/LAPACK: MAGMA, CUBLAS, ROCM
	- `GPU_BACKEND` - CUDA (default), HIP
- System/compiler-specific flags are inferred in `source/Makefile.internal`.

- `VPATH` used heavily to locate source files.

- This build system is quite dated and flawed and needs to be modernized (either with improved GNU Make framework or CMake)

## Testing

The existing test infrastructure is **effectively broken and needs modernization/replacement**.
Many of the individual tests are scientifically valid, but the test mechanism itself is flawed.

```bash
# From source/, these commands exist but have known issues:
make CMODE=DEBUG test_simple    # Self-heating alpha network
make CMODE=DEBUG test_setup     # Network setup
make CMODE=DEBUG test_nse       # NSE solver
make test                       # Runs test_serial, test_heat, test_setup, test_nse

# Individual test problems by ID (from test/ directory)
./test_xnet.sh <executable> <test_id>
```

### Why the existing tests are being scrapped

- `test_diff()` in `test_xnet.sh` does not exit nonzero on failure — CI passes
  even when results diverge. The tests are a no-op.
- Timer stripping uses hardcoded `sed +14` line count — breaks if timer output changes.
- Only 3 of 17+ test scenarios run in CI.
- Tests diff entire `net_diag` output byte-for-byte, including formatting and
  metadata irrelevant to correctness.
- Reference solutions are not reliably available for all test problems.
- No unit tests exist.

Do not invest time fixing `test_xnet.sh` or `makefile_tests.yml`. Build new
test infrastructure from scratch.

### Criteria for Retiring/Replacing Existing Tests

When replacing an existing regression test:
1. Identify what physics or code path it was covering
2. Ensure replacement coverage exists (unit test, physics-informed test, or
   improved regression test) before removing the old test
3. Document the retirement in the PR description

### Proposed New Test Ideas

Tests organized into three tiers:

#### Tier 1: Unit Tests (seconds, run on every push)

Test individual routines and modules in isolation.

| What to test | Module | Example |
|-------------|--------|---------|
| `safe_exp()` edge cases (overflow, underflow, zero) | `xnet_util` | `check(safe_exp(800.0_dp), ...)` |
| Physical constants sanity | `xnet_constants` | `check(amu, 1.66054e-24_dp, rel=.true., thr=1e-5_dp)` |
| Fermi-Dirac integrals vs tabulated values | `xnet_fd` | Compare F_k(eta) to published tables |
| Rate interpolation (REACLIB 7-parameter fit) | `xnet_data` | Evaluate known rate at known T9 |
| Partition function evaluation | `xnet_data` | Compare to published values at specific T9 |
| Screening corrections vs known limits | `xnet_screening` | Known Z1, Z2 at known T9/rho |
| EOS thermodynamic identity | `xnet_eos_*` | dP/dT at const rho, dP/drho at const T |
| Forward/reverse detailed balance | `xnet_match` | lambda_fwd / lambda_rev = expected from Q-value |
| Sparse indexing consistency | `xnet_jacobian_*` | Small network: sparse pattern matches dense |

#### Tier 2: Physics-Informed Tests (seconds to minutes, run on every PR)

Test physical invariants and known analytic results across integration.

| What to test | How |
|-------------|-----|
| **Mass conservation** | Run alpha network for N steps; verify `sum(A_i * Y_i) = 1` to tolerance |
| **Charge conservation** | Run without weak rates; verify `sum(Z_i * Y_i)` is constant |
| **Non-negative abundances** | Run any network; verify all `Y_i >= 0` at every output step |
| **NSE consistency** | At T9=10, rho=1e9: compare NSE solver to Saha equation for alpha network |
| **Detailed balance at equilibrium** | Initialize at NSE; run 1 step; verify ydot ~ 0 (net flux ~ 0) |
| **Jacobian accuracy** | Compare analytic Jacobian (`jacobian_build`) to finite-difference Jacobian |
| **Convergence order** | Run same problem at dt, dt/2, dt/4; verify BE is O(dt), BDF is higher order |
| **Solver equivalence** | Same problem with dense vs other solvers: final abundances agree to tolerance |
| **Energy consistency** | Verify edot = sum(dY_i/dt * BE_i * amu * c^2) matches reported energy generation |
| **Freeze-out** | Rapidly cool from T9=5 to 0.1: verify composition stops changing |
| **Integrator agreement** | Same initial conditions: BE and BDF give same answer (to tolerance) |

#### Tier 3: Integration/Regression Tests (minutes, run on merge to main)

Full-problem runs that verify end-to-end correctness. These replace the
current `test_xnet.sh` tests with properly tolerance-checked comparisons.

Representative problems (not exhaustive):
- Self-heating C/O detonation (alpha network) — fast, exercises self-heating
- Thermonuclear SN tracer (SN160 network) — exercises medium network
- Core-collapse SN tracer — exercises NSE transition, weak rates
- X-ray burst (CNO network) — exercises proton-rich burning
- Network setup round-trip — verifies preprocessing

### Proposed Tolerance Strategy

Tests do NOT require bitwise reproducibility. Tolerances depend on the quantity:

| Quantity | Method | Typical tolerance |
|----------|--------|-------------------|
| Conservation laws (mass, charge) | Absolute | 1e-12 (roundoff level) |
| Non-negativity | Absolute | Y_i >= -1e-30 (allow tiny roundoff) |
| Final abundances (regression) | Relative, per-species | 1e-6 dominant, 1e-3 trace |
| NSE vs Saha | Relative | 1e-4 |
| Jacobian vs finite difference | Relative, element-wise | 1e-4 |
| Convergence order | Slope of log-log | +/- 0.2 of expected order |
| Solver equivalence (dense vs sparse) | Relative | 1e-10 |
| Energy generation rate | Relative | 1e-6 |

For **monitoring changes over time** (performance, timestep count):
- Record current value as baseline
- Alert (not fail) if value changes by > 10%
- Fail if value changes by > 50%

## Code Architecture

### Directory Structure

```
XNet/
  ├── source/                          # Core Fortran source and build system
  │   ├── Makefile                     # Main build rules and dependency graph
  │   ├── Makefile.opt                 # User-selectable build options (compiler, solver, EOS, GPU)
  │   ├── Makefile.internal            # Compiler/system-specific flags (inferred from Makefile.opt)
  │   ├── Makefile.dev                 # Local developer overrides (optional, -include'd)
  │   ├── xnet_macros.fh               # Preprocessor macros abstracting OpenACC/OpenMP/serial directives
  │   │
  │   │ # --- Foundation (Layer 0-1) ---
  │   ├── xnet_types.F90               # Precision parameters (dp, sp, qp, i4, i8)
  │   ├── xnet_constants.F90           # Physical/numerical constants in CGS
  │   ├── xnet_util.F90                # Utilities: safe_exp(), file I/O helpers
  │   ├── xnet_timers.F90              # Wall-clock profiling
  │   ├── xnet_fd.F90                  # Fermi-Dirac integrals
  │   ├── xnet_parallel.F90            # MPI wrappers
  │   ├── xnet_parallel_stubs.F90      # Serial fallback (no MPI)
  │   │
  │   │ # --- Data (Layer 2) ---
  │   ├── xnet_controls.F90            # Input parameters, zone batching, integration control
  │   ├── xnet_data.F90                # Nuclear species (nuclear_data) + reaction networks (reaction_data)
  │   ├── xnet_conditions.F90          # T9/rho/Ye trajectories, timestep state
  │   ├── xnet_abundances.F90          # Y, Ydot arrays
  │   │
  │   │ # --- Physics (Layer 3) ---
  │   ├── xnet_ffn.F90                 # Weak interaction rates (FFN electron capture/beta decay)
  │   ├── xnet_nnu.F90                 # Neutrino-nucleus interactions
  │   ├── xnet_screening.F90           # Coulomb screening corrections
  │   ├── xnet_flux.F90                # Net reaction flux calculation
  │   ├── xnet_match.F90               # Forward/reverse reaction pairing
  │   ├── xnet_nse.F90                 # Nuclear statistical equilibrium
  │   ├── xnet_preprocess.F90          # ASCII→binary nuclear data conversion
  │   ├── xnet_eos_starkiller.F90      # Helmholtz EOS via starkiller-helmholtz (default)
  │   ├── xnet_eos_helm.F90            # Helmholtz EOS (standalone)
  │   ├── xnet_eos_bahcall.F90         # Simple Bahcall EOS
  │   │
  │   │ # --- Linear Algebra & Jacobian (Layer 4) ---
  │   ├── xnet_linalg.F90              # BLAS/LAPACK wrappers (CPU + GPU dispatch)
  │   ├── xnet_jacobian_dense.F90      # Dense Jacobian (LAPACK)
  │   ├── xnet_jacobian_MA48.F90       # HSL MA48 sparse Jacobian
  │   ├── xnet_jacobian_PARDISO.F90    # PARDISO sparse Jacobian
  │   ├── xnet_jacobian_PARDISO_MKL.F90 # PARDISO via Intel MKL
  │   │
  │   │ # --- Integration (Layer 5) ---
  │   ├── xnet_integrate.F90           # Shared: cross_sect(), yderiv(), timestep()
  │   ├── xnet_integrate_be.F90        # Backward Euler solver
  │   ├── xnet_integrate_bdf.F90       # BDF (Gear) multi-step solver
  │   ├── xnet_integrate_bd.F90        # (Unused/incomplete BD variant)
  │   │
  │   │ # --- Driver & I/O (Layer 6) ---
  │   ├── xnet_evolve.F90              # full_net() main time-stepping loop
  │   ├── xnet_output.F90              # Diagnostics and formatted output
  │   ├── model_input_ascii.F90        # ASCII thermodynamic profile reader
  │   ├── net.F90                      # Main program: init, zone loop, I/O
  │   ├── net_setup.F90                # Network setup preprocessing utility
  │   ├── nse_slice.F90                # Standalone NSE calculator
  │   │
  │   │ # --- GPU bindings ---
  │   ├── xnet_gpu.F90                 # Device init, stream management
  │   ├── cudaf.F90                    # CUDA Fortran bindings
  │   ├── cublasf.F90                  # cuBLAS Fortran bindings
  │   ├── cusolverf.F90                # cuSOLVER Fortran bindings
  │   ├── cusparsef.F90                # cuSPARSE Fortran bindings
  │   ├── hipf.F90                     # HIP Fortran bindings
  │   ├── hipblasf.F90                 # hipBLAS Fortran bindings
  │   ├── rocblasf.F90                 # rocBLAS Fortran bindings
  │   ├── rocsolverf.F90               # rocSOLVER Fortran bindings
  │   ├── rocsparsef.F90               # rocSPARSE Fortran bindings
  │   ├── magmaf.F90                   # MAGMA Fortran bindings
  │   ├── openaccf.F90                 # OpenACC helper routines
  │   └── openmpf.F90                  # OpenMP target helper routines
  │
  ├── test/                            # Regression test suite
  │   ├── test_xnet.sh                 # Main test driver (bash); defines all test IDs
  │   ├── test_xnet.csh                # Legacy csh test driver
  │   ├── test_settings*               # 9 settings variants (default, small, heat, bdf, logft, nse, xnse, batch, parallel)
  │   ├── build_net/                   # Network setup test data
  │   ├── Data_*/                      # Pre-built networks of various sizes
  │   └── Test_Problems/               # Test inputs (setup_*, ab_*, th_*) and Results/
  │
  ├── tools/                           # Supporting utilities and vendored dependencies
  │   ├── LAPACK/                      # NETLIB BLAS/LAPACK subset
  │   ├── starkiller-helmholtz/        # Helmholtz EOS tables + Fortran interface
  │   ├── build/                       # Network generation tools (Makefile_include)
  │   ├── analysis/                    # Post-processing: final abundances, HDF comparison
  │   ├── python/                      # Plotting & analysis scripts (plot_time_mf, draw_nz_flux, etc.)
  │   ├── thermo/                      # Thermodynamic trajectory utilities
  │   └── initial_abundance/           # Initial abundance file generation
  │
  ├── doc/                             # Documentation
  │   ├── XNet_Formatting_Guidelines.md # Fortran coding conventions
  │   ├── xnet/                        # Core XNet documentation (XNet.pdf)
  │   ├── reaclib/                     # Nuclear data format documentation
  │   ├── screening/                   # Screening physics references
  │   └── solvers/                     # Sparse solver documentation (HSL, PARDISO)
  │
  ├── .github/workflows/               # CI configuration
  │   └── makefile_tests.yml           # GitHub Actions: test_simple, test_setup, test_nse
  │
  ├── AGENTS.md                        # LLM coding agent guidance
  └── README.md                        # Project overview (minimal)
```

### Key Subsystems

- **Integration**: `xnet_evolve.F90` drives time evolution, calling `xnet_integrate.F90` (shared routines), `xnet_integrate_be.F90` (Backward Euler), `xnet_integrate_bdf.F90` (BDF)
- **Jacobian/Solvers**: `xnet_jacobian_dense.F90`, `xnet_jacobian_MA48.F90`, `xnet_jacobian_PARDISO.F90`, `xnet_jacobian_PARDISO_MKL.F90` — selected at build time via `MATRIX_SOLVER`
- **Physics**: `xnet_flux.F90` (reaction fluxes), `xnet_screening.F90` (screening corrections), `xnet_ffn.F90` (weak rates), `xnet_nnu.F90` (neutrino interactions), `xnet_nse.F90` (nuclear statistical equilibrium)
- **EOS**: `xnet_eos_starkiller.F90` (default Helmholtz), `xnet_eos_bahcall.F90`, `xnet_eos_helm.F90` — selected at build time via `EOS`
- **GPU**: `xnet_gpu.F90` manages GPU state; `xnet_macros.fh` provides preprocessor abstraction for OpenACC/OpenMP/serial GPU directives
- **Parallelization**: `xnet_parallel.F90` (MPI) or `xnet_parallel_stubs.F90` (serial fallback)

### Key Data Flow (per timestep)

```
full_net loop iteration:
  timestep() → cross_sect() → { update_eos(), screening(), partf(), ffn_rate(), nnu_rate() }
                             → computes csect1-4 (velocity-integrated cross sections)
             → yderiv()     → computes ydot from csect × abundances
  solve_be/bdf() → Newton-Raphson loop:
                     yderiv() → jacobian_build() → jacobian_solve() → update yt
```
