# Current architecture overview

> Read this document when work spans modules, changes subsystem interaction,
> alters data flow, or requires architectural reasoning.

This document describes the current repository. Verify the relevant source
before making a change, and update this overview when verified implementation
details change. Proposed designs belong in their governing issue or a
dedicated design document.

## Repository areas

- `source/` contains production Fortran, the stand-alone driver, and utility
  programs. The production GNU Make entry point and configuration files are at
  the repository root, with included rules under `make/`.
- `tools/LAPACK/` contains the vendored NETLIB subset used by the tracked
  baseline configuration.
- `tools/starkiller-helmholtz/` supplies the tracked default EOS
  implementation and table.
- `test/` contains legacy problem inputs, network data, settings, and shell
  drivers.
- `tools/` also contains network-building, thermodynamic, analysis, and Python
  utilities.
- `doc/` contains formatting guidance and scientific or solver references.
- `docs/development/` contains task-oriented development guidance.

Production builds place objects, module files, retained preprocessing, and
executables under a readable build directory such as `build/GNU-OPT`, or one
selected explicitly with `BUILD_NAME` or `BUILD_DIR`.
The build compiles selected dependencies from `tools/` without generating
source-tree products.

## Production areas

| Area | Principal files | Current role |
| --- | --- | --- |
| Foundation | `xnet_types`, `xnet_constants`, `xnet_util`, `xnet_timers`, `xnet_fd` | Kinds, constants, utilities, timing, and Fermi-Dirac functions |
| Configuration and state | `xnet_controls`, `xnet_conditions`, `xnet_abundances`, `xnet_data` | Runtime controls, thermodynamic state, abundances, species data, and reaction data |
| Physics | `xnet_ffn`, `xnet_nnu`, `xnet_screening`, `xnet_flux`, `xnet_match`, `xnet_nse`, `xnet_eos_*` | Weak and neutrino rates, screening, fluxes, reaction matching, NSE, and EOS selection |
| Linear algebra | `xnet_linalg`, `xnet_jacobian_*` | BLAS/LAPACK dispatch and solver-specific Jacobian storage and solves |
| Integration | `xnet_integrate`, `xnet_integrate_be`, `xnet_integrate_bdf` | Rate updates, derivatives, timestep logic, and implicit integration |
| Driver and output | `net`, `model_input_ascii`, `xnet_evolve`, `xnet_output` | Initialization, zone input, evolution, diagnostics, and the `xnet` program |
| Utilities | `net_setup`, `nse_slice` | Network preprocessing and the `xnse` program |
| Portability | `xnet_parallel*`, `xnet_gpu`, `xnet_macros.fh`, vendor binding modules | Serial/MPI selection and accelerator interfaces |

This table is a navigation aid. Module dependencies and side effects cross
these areas, especially through shared module state.

## Program setup

`net.F90` is the stand-alone production program. Its setup path currently:

1. initializes serial or MPI execution, determines the OpenMP thread count,
   initializes GPU execution when selected, and starts setup timing;
2. reads `control` through `xnet_controls`;
3. preprocesses the requested network when needed, then reads nuclear,
   reaction, Jacobian, and match data;
4. initializes screening, flux evaluation, the selected EOS and integrator,
   NSE support, and shared run arrays;
5. reads thermodynamic histories and initial abundances through
   `model_input_ascii`;
6. evolves assigned zones and writes diagnostic or timestep output;
7. finalizes accelerator and parallel resources.

The exact order and conditional calls live in `source/net.F90`.

## Per-timestep flow

`full_net()` in `xnet_evolve.F90` coordinates a zone batch through the time
integration loop. The selected runtime integration method calls the shared
integration support and a build-selected Jacobian module.

The main numerical flow is:

```text
full_net()
  -> timestep()
       -> retain the prior timestep state and identify zones needing an initial estimate
       -> cross_sect() and yderiv() for zones needing an initial estimate
       -> timestep selection and trial-state preparation
  -> selected BE or BDF implicit solve
       -> Newton iterations
       -> cross_sect()
            -> EOS-dependent quantities
            -> screening and partition functions
            -> weak and neutrino rates
            -> velocity-integrated reaction-rate arrays
       -> yderiv()
            -> abundance derivatives
       -> jacobian_build()
       -> jacobian factorization/solve
       -> abundance and optional temperature update
  -> accept/retry logic and output
```

Inspect `xnet_evolve.F90`, `xnet_integrate.F90`, the selected integrator, and
the selected Jacobian file together when changing this path.

## Build-time selection

`Makefile.opt` exposes compiler environment, compile mode, MPI,
threading, accelerator mode, EOS, matrix solver, and CPU/GPU linear-algebra
choices. `Makefile.internal` maps those choices to compilers, flags,
sources, and libraries. The conventional fragments under `make/`
separate implementation selection, source lists, explicit module
prerequisites, and build rules. In particular, `providers.mk` selects the
source or library that supplies the requested MPI, EOS, Jacobian,
numerical-library, and accelerator implementation.
`machines.mk` retains hostname/LMOD detection and explicitly selects the
tracked generic or Cray Programming Environment defaults; compiler-family
flags remain in `Makefile.internal`.

Important selections change which file is intended to supply a common module
name:

- `xnet_parallel.F90` or `xnet_parallel_stubs.F90` supplies
  `xnet_parallel`.
- `xnet_eos_starkiller.F90`, `xnet_eos_helm.F90`, or
  `xnet_eos_bahcall.F90` supplies `xnet_eos`.
- one `xnet_jacobian_*.F90` file supplies `xnet_jacobian`.

Source inspection and tests must follow the configuration named by the task.

`xnet_linalg.F90` exposes one batched GPU LU/solve routine family:

```text
LinearSolveBatched_GPU
  -> LUDecompBatched_GPU
  -> LUBksubBatched_GPU
```

The vendor calling convention remains inside those routines. cuBLAS and MAGMA
use pointer-array batched calls; hipBLAS and oneMKL use strided-batched calls.
`xnet_jacobian_dense.F90` does not select between those representations.
cuBLAS retains XNet's implemented no-pivot path. XNet currently rejects a
no-pivot request for the hipBLAS, oneMKL, and MAGMA implementations rather
than silently changing the requested solve.

## Shared state and interfaces

Many modules own allocatable arrays and control values shared through `Use`
associations. Zone batching, thread ranges, and active-zone masks organize
work across these arrays. Changes to allocation, dimensions, initialization,
or batch indexing can affect distant routines.

`xnet_data.F90` currently defines both `nuclear_data` and `reaction_data`.
The selected EOS and Jacobian implementations use common module names, which
makes their interfaces compile-time requirements across implementations.

The stand-alone runtime interface includes:

- the ordered, block-oriented `control` file;
- nuclear data and pre-built network files under `Data_*` directories;
- thermodynamic histories and abundance inputs;
- fixed-format diagnostic and timestep outputs.

Treat ordering, units, dimensions, naming conventions, and output formats as
compatibility-sensitive behavior. Trace reads and writes in source when
changing any of them.

## Architectural investigation

For a change spanning modules:

1. identify the build-time implementations involved;
2. trace module variables, allocation, and initialization from the driver;
3. trace the active scalar, vector, batch, thread, MPI, and accelerator paths;
4. identify runtime and output interfaces affected;
5. characterize current behavior before structural edits;
6. state proposed structural changes separately from the current-state report.

Document uncertainty when the current implementation or supported
configurations have limited evidence.
