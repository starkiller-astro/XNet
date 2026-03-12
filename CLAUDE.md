# CLAUDE.md - XNet

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **How to read this file:** Sections marked **[CURRENT STATE]** describe the codebase as it actually exists today.
> Sections marked **[TARGET STATE]** describe where we are going.
> Sections marked **[MIGRATION POLICY]** tell you how to behave at the boundary.
> Never confuse current state with target state.
> When in doubt, read the code and ask.

---

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

* **[CURRENT STATE]:** Fortran 2003 (source code), GNU Make (build system), Python (analysis/plotting), bash (testing)
* **[TARGET STATE]:** Fortran 2018 (source code), CMake + GNU Make (build system), Python (analysis/plotting/testing)

---

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

---

## Build System [CURRENT STATE]

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

---

## Testing [CURRENT STATE]

The existing test infrastructure is **effectively broken and will be replaced**.

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
test infrastructure from scratch (see Testing [TARGET STATE]).

### Criteria for Retiring/Replacing Existing Tests

When replacing an existing regression test:
1. Identify what physics or code path it was covering
2. Ensure replacement coverage exists (unit test, physics-informed test, or
   improved regression test) before removing the old test
3. Document the retirement in the PR description

---

## Testing [TARGET STATE]

### Unit Testing Framework

XNet uses **test-drive** (fortran-lang/test-drive), a lightweight procedural
Fortran testing framework, pulled in as a git submodule. It is a single
Fortran source file with no external dependencies.

Key properties:
- Pure Fortran, no preprocessor — tests are standard `.F90` files
- `check()` assertions for logical, integer, real, complex, and string comparisons
- Real comparisons support absolute and relative tolerances
- Test discovery via `collect` subroutines; individual tests selectable by name
- Works with Make (no CMake required)
- Each test suite compiles to a standalone executable that exits nonzero on failure

### Test Tiers

Tests are organized into three tiers:

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

### Tolerance Strategy

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

---

## Code Architecture [CURRENT STATE]

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
  ├── CLAUDE.md                        # Claude Code guidance
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

---

## Fortran Coding Conventions [CURRENT STATE]

* See `doc/XNet_Formatting_Guidelines.md` for full details (somewhat dated). Key rules:

	- **Case**: lowercase except capitalized language keywords (`Call`, `Subroutine`, `EndDo`, `EndIf`, `Implicit None`)
	- **Indentation**: 2 spaces, no tabs, max 132 chars/line
	- **Line continuation**: `&` at both ends (`& var` on continuation line)
	- **Modules**: always `Use module, Only: needed_items`
	- **Arguments**: always declare `Intent(In/Out/InOut)`
	- **All routines**: include `Implicit None` even inside modules
	- **End routines with**: explicit `Return`
	- **Deallocate**: explicitly deallocate local allocatable variables
	- **No magic numbers**: use named `Parameter` constants
	- **Comparison operators**: use `>`, `<` not `.gt.`, `.lt.`
	- **Avoid**: post-Fortran2008 features, `While` loops
	- **Use**: `iso_fortran_env` variables (`output_unit` not unit 6)
	- **Use `Integer`** instead of `Logical` for execution flags

---

## Fortran Coding Conventions [TARGET STATE]

Extends [CURRENT STATE] conventions with modernization rules for new code:

- **Fortran 2008 is the floor**: new code may use any F2008 feature.
  F2018 features require explicit approval.
- **`error stop`**: prefer `error stop 'message'` over `stop` for abnormal
  termination. Use `stop 0` only for normal exit from a program unit.
- **No `block` constructs**: do not use `block`/`end block`.
- **`Integer` for flags**: retain existing convention (Integer, not Logical,
  for on/off execution flags).
- **Submodules**: consider for modules exceeding ~500 lines with multiple
  independent implementation sections, to reduce recompilation cascades.
- **`associate`**: use for long derived-type chains to improve readability,
  but do not use as a general-purpose rename mechanism.
- **`intent` on all arguments**: no exceptions, even for single-line wrappers.
- **Allocatable over pointer**: prefer `allocatable` arrays over `pointer`
  arrays unless pointer semantics (aliasing, target) are genuinely needed.
- **`iso_fortran_env`**: use `int32`, `int64`, `real64` etc. from
  `iso_fortran_env` in new code. Existing `dp`, `sp`, `i4` aliases in
  `xnet_types` remain authoritative within XNet source.

---

## Known Technical Debt (do not propagate, do not silently fix)

These are **known deficiencies**. Do not silently fix them unless the task
explicitly targets that debt. Do not propagate these patterns in new code.

### Code Structure
- Global state is extensively used via module-level allocatable arrays indexed
  by zone. Thread safety relies on zone-batch slicing, not encapsulation.
- ~3000 lines duplicated across 4 Jacobian variants (`xnet_jacobian_*.F90`).
  The `jacobian_build()` physics is nearly identical; only storage format and
  solver API calls differ.
- `xnet_data.F90` contains two modules (`nuclear_data`, `reaction_data`) in
  one file.
- Many subroutines are large (500+ lines) with side effects on module globals.

### I/O and Configuration
- Runtime parameters are read via fixed-format ASCII `control` files (not namelist).
- `xnet_output.F90` writes fixed-format `net_diag` files; format changes
  break all legacy regression tests.

### Build System
- `make clean` does not remove executables (xnetd, xnetm, xnetp, xnse, net_setup).
- `xnet_integrate_bd.F90` has Makefile dependency lines but is never compiled.

### Testing and CI
- The existing CI and test infrastructure (`test_xnet.sh`, `makefile_tests.yml`)
  is effectively non-functional and will be replaced entirely. See
  Testing [CURRENT STATE] for details.
- `safe_exp()` is used only in `xnet_nse.F90`; other modules use bare `exp()`.

### Numerics
- Hardcoded constants without explanation in several modules (e.g., `1e-42`
  in `xnet_nnu.F90`).

---

## Documentation Standards [TARGET STATE]

- Every module: header comment stating purpose and key references.
- Every public subroutine/function: brief description, intent of each argument,
  and (where applicable) the equation or algorithm being implemented with a
  literature citation.
- Physics routines should cite the relevant paper and equation number.
  Example:
  ```fortran
  ! Coulomb screening correction using the Graboske et al. (1973) prescription.
  ! See Eq. 4 of Hix & Thielemann (1999), ApJ 511, 862.
  ```
- Do NOT add boilerplate comments that restate what the code obviously does.
- When modifying a routine, bring its documentation up to this standard.
  Do not make documentation-only PRs for code you are not otherwise changing.
- **FORD** (FORtran Documenter) is the target for generated HTML documentation.
  FORD-compatible comment markup (`!>`, `!!`) should be used in new code and
  added to existing code when it is touched for other reasons.

---

## CI / GitHub Actions [TARGET STATE]

### Pipeline Structure

CI uses GitHub Actions with three tiers matching the test tiers:

| Tier | Trigger | Time budget | What runs |
|------|---------|-------------|-----------|
| **Fast** | Every push | < 2 min | Build (debug) + all Tier 1 unit tests |
| **Standard** | Every PR | < 10 min | Fast + Tier 2 physics tests + selected Tier 3 smoke tests |
| **Full** | Merge to main, weekly | < 30 min | Standard + full Tier 3 regression suite |

### Compiler Matrix

Starting with GFortran (latest) on ubuntu-latest. Will expand to include
Intel ifx (via oneAPI) and additional compilers/architectures over time.

### Requirements

- All test executables must exit nonzero on failure.
- CI must fail if any test fails — no silent warnings.
- Test output (diffs, logs) must be uploaded as GitHub Actions artifacts on failure.
- Performance baselines (timestep counts, wall time) tracked as artifacts
  for trend monitoring (alert, not fail).

---

## Code Architecture [TARGET STATE]

See `doc/target_architecture.md` for the full target directory tree.

---

## Modernization Roadmap [TARGET STATE]

See `doc/modernization_roadmap.md` for the phased checklist.

---

## Migration Policy [MIGRATION POLICY]

These rules govern how to behave when current state and target state conflict.
They apply to both human developers and Claude Code.

### General Rules

- **Do not move files speculatively.** File moves (`source/` → `src/physics/`) happen only
  as part of an explicit, planned migration step — never as a side effect of a bug fix
  or feature addition.
- **Current build system is authoritative.** Until the CMake migration is complete,
  `source/Makefile` is the source of truth. Do not add `CMakeLists.txt` files that
  contradict or duplicate it.
- **One structural change per commit.** A commit either moves/renames files OR changes
  logic, never both. This keeps git blame useful and reverts safe.

### File Moves and Renames

- When moving a file (e.g., `xnet_ffn.F90` → `src/physics/xnet_ffn.F90`), update:
  1. Makefile `VPATH` or explicit paths
  2. Any `Include` / `#include` directives referencing the old path
  3. `CLAUDE.md` directory structure (both `[CURRENT STATE]` and any file lists)
- Do NOT rename Fortran module names during a file move. Module renames are a
  separate step.

### Refactoring Boundaries

- **Jacobian deduplication**: Extract shared `jacobian_build()` physics into a common
  module (`xnet_jacobian_common.F90`) that each solver variant calls. Each variant
  retains only its storage-format-specific code and solver API calls. A later
  modernization step may introduce an abstract type with deferred procedures, but the
  common-module approach comes first.

- **EOS interface**: Keep compile-time selection (as today). The 3 EOS variants
  (starkiller, helm, bahcall) remain as separate modules with a common implicit
  interface. No abstract type unification is planned.

- **xnet_data.F90 split**: Leave as-is for now. The two modules (`nuclear_data`,
  `reaction_data`) stay in one file until a larger data-layer refactor is undertaken.

- **Global state reduction**: The long-term target is a state-object pattern — one or
  more derived types passed through the call chain instead of module-global allocatable
  arrays. This improves thread safety and testability. However, this is a large change
  that touches every routine signature, so it is **deferred until better infrastructure
  and test coverage are in place**. Approach will be incremental (module by module,
  starting from the leaves). New code should prefer explicit state passing where
  practical, but do not refactor existing global state without explicit instruction.

- **xnet_integrate_bd.F90**: This is a reference/template kept for historical purposes.
  It should be moved to an `archive/` or `deprecated/` directory to separate it from
  actively compiled source. Remove its dependency lines from the Makefile when moved.

### What Claude Code Should Do Autonomously

- Fix formatting violations (indentation, case, line length) in any file it touches
  for other reasons. Do not make formatting-only PRs.
- Update Makefile dependency lines when adding/removing `Use` statements.
- Add `Implicit None` to any routine found without it.
- Replace magic numbers with named parameters when the meaning is clear from context.
- Use `safe_exp()` instead of bare `exp()` where the argument could overflow.

### What Claude Code Should NOT Do Without Explicit Instruction

- Move or rename source files
- Change module names or public API signatures
- Add new module-global mutable state
- Modify the Makefile build targets or link rules
- Delete any file, even apparently dead code
- Change numerical algorithms or convergence criteria
- Introduce Fortran 2018 features (stick to F2008 unless told otherwise)

