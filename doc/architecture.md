# Current architecture

This document describes the repository as it exists. It does not prescribe a
future directory layout.

## Repository areas

- `source/` contains the production Fortran and GNU Make build.
- `tools/LAPACK/` contains the vendored NETLIB subset used by the characterized
  serial build.
- `tools/starkiller-helmholtz/` supplies the STARKILLER EOS implementation and
  table.
- `test/` contains legacy problem inputs, network data, settings, and the
  historical shell driver.
- `tests/` contains the isolated serial unit-test harness.
- `tools/` also contains network-building, thermodynamic, analysis, and Python
  utilities.
- `doc/` contains current guidance and historical scientific or solver
  references.

Production builds are currently in-source: objects, module files, and
executables are written under `source/`. The Makefiles use `VPATH` to compile
some dependencies from `tools/`. A clean rebuild is required when changing
configuration.

## Production layers

| Area | Principal files | Role |
| --- | --- | --- |
| Foundation | `xnet_types`, `xnet_constants`, `xnet_util`, `xnet_timers`, `xnet_fd` | Kinds, constants, utilities, timing, Fermi-Dirac functions |
| Configuration and state | `xnet_controls`, `xnet_conditions`, `xnet_abundances`, `xnet_data` | Runtime controls, thermodynamic state, abundances, species and reaction data |
| Physics | `xnet_ffn`, `xnet_nnu`, `xnet_screening`, `xnet_flux`, `xnet_match`, `xnet_nse`, `xnet_eos_*` | Rates, screening, fluxes, reaction matching, NSE, and EOS selection |
| Linear algebra | `xnet_linalg`, `xnet_jacobian_*` | BLAS/LAPACK dispatch and solver-specific Jacobians |
| Integration | `xnet_integrate`, `xnet_integrate_be`, `xnet_integrate_bdf` | Cross sections, derivatives, timestep logic, and implicit solvers |
| Driver and output | `xnet_evolve`, `xnet_output`, `model_input_ascii`, `net` | Zone evolution, diagnostics, input, and the `xnet` program |
| Utilities | `net_setup`, `nse_slice` | Network preprocessing and the `xnse` program |
| Portability | `xnet_parallel*`, `xnet_gpu`, `xnet_macros.fh`, vendor binding modules | Serial/MPI selection and accelerator interfaces |

## Build-time selection

`source/Makefile.opt` exposes compiler environment, compile mode, EOS, matrix
solver, CPU/GPU linear algebra, and parallel modes.
`source/Makefile.internal` maps those choices to compiler flags, source
dependencies, and libraries. `source/Makefile` chooses the serial parallel
stubs or MPI wrappers, the EOS objects, and the solver-specific Jacobian.

The presence of a selection in these files means only that the path is
available. It does not establish that the path is characterized, scientifically
validated, or supported. Recorded evidence is listed in
[validation-plan.md](validation-plan.md).
