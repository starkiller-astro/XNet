# Build and test guidance

> Read this document when changing build logic, selecting a configuration,
> running legacy problems, adding tests, comparing numerical results, or
> measuring performance.

This document describes current repository behavior. Verify task-relevant
details in the Makefiles, test drivers, and source before relying on them.

## Build entry point

The production build uses GNU Make and writes its results below a
caller-selected build directory. From the repository root, build the tracked
default with:

```bash
make -j
```

With the tracked defaults, the executable is `build/GNU-OPT/bin/xnet`. The
major public selectors form a readable automatic directory name. Set
`BUILD_NAME` to choose a different name below `build/`, or set `BUILD_DIR`
directly. Relative paths are interpreted from the repository root:

```bash
make BUILD_NAME=gnu-debug CMODE=DEBUG -j xnet
make BUILD_DIR=/scratch/$USER/xnet-frontier -j xnet
```

The Makefile fragments have distinct roles:

- `Makefile` is the small public entry point. Its included `make/build.mk`
  defines the production build and its output layout.
- `make/configuration.mk` validates compiler/platform selectors.
- `make/providers.mk` selects the MPI, EOS, solver, accelerator, and
  numerical-library implementations.
- `make/sources.mk` lists sources and output files and records
  configuration reuse.
- `make/dependencies.mk` records explicit Fortran module prerequisites.
- `make/rules.mk` contains preprocessing, compilation, link, and public
  target rules.
- `Makefile.opt` defines tracked user-selectable defaults.
- `Makefile.internal` maps configuration choices to compilers, flags,
  libraries, source files, and solver objects.
- `make/machines.mk` detects the current host and explicitly selects
  the tracked generic or Cray Programming Environment defaults. Compiler
  defaults remain in `Makefile.internal`; add a machine fragment and one
  visible mapping entry only when a real repository-supported machine needs
  concrete overrides. Perlmutter and Frontier select the maintained HPE Cray
  Programming Environment settings; other hosts use the generic defaults.

Inspect the conditional path through these files for any configuration being
changed. Variable names and commented examples provide orientation; the
selected Make logic determines the build.

Each build directory contains a fixed-order `config.txt` record. Reusing it
with different selectors, selected libraries or sources, compiler commands,
or effective flags fails; clean that directory or choose another one. The
conventional POSIX record supports ordinary flag text but rejects effective
flags containing a single quote or line break. Do not run two top-level Make
processes in one build directory, and do not clean a directory while another
invocation uses it. Different build directories may build concurrently.

## Tracked defaults

The tracked selector defaults are:

| Setting | Value | Meaning |
| --- | --- | --- |
| `CMODE` | `OPT` | Optimized build |
| `PE_ENV` | `GNU` | GNU compiler configuration |
| `MPI_MODE` | `OFF` | Selects the serial parallel-interface stubs |
| `OPENMP_MODE` | `OFF` | OpenMP host threading disabled |
| `GPU_MODE` | `OFF` | Accelerator runtime and libraries disabled |
| `GPU_BACKEND` | `CUDA` | Accelerator vendor selection when GPU mode is enabled |
| `EOS` | `STARKILLER` | Starkiller Helmholtz EOS interface |
| `MATRIX_SOLVER` | `dense` | Dense Jacobian and linear solve |
| `LAPACK_VER` | `NETLIB` | Default on ordinary non-Cray systems |

These values describe configuration selection. Record actual validation in the
governing issue or PR with the compiler version, command, machine, relevant
environment, result, and date.

## Optional sparse-backend status

A named Make target records a build recipe, not a support claim. The current
serial CPU support status is:

| Implementation | Selection | Status and limit |
| --- | --- | --- |
| HSL MA48 2.2.0 | `MATRIX_SOLVER=MA48` with external `MA48.f` | Qualified on macOS arm64 with GNU Fortran 16.2.0 using the production build. The maintainer-supplied source is used under a maintainer-held non-redistributable HSL licence and must remain outside the repository. Other HSL versions, compilers, and platforms are unqualified. |
| Intel oneMKL PARDISO | `MATRIX_SOLVER=PARDISO_MKL LAPACK_VER=MKL` | Historical qualification used the `etacar` Linux x86_64 host with GNU Fortran 11.4.0 and oneMKL 2026.1. That run predates the current production Makefiles and does not qualify this tree. Other oneMKL versions, compilers, platforms, and parallel modes are also unqualified. |

MA48 source must not be copied, committed, archived, or attached to an issue or
pull request. HSL describes MA48 2.2.0 and its licensing restrictions in the
[official catalogue](https://www.hsl.rl.ac.uk/catalogue/ma48.html) and
[licensing overview](https://www.hsl.rl.ac.uk/). oneMKL is distributed under
Intel's [Simplified Software License](https://www.intel.com/content/www/us/en/developer/articles/tool/onemkl-license-faq.html).
The exact opt-in component and same-source dense comparison commands are in
[`test/qualification/sparse_backends/README.md`](../../test/qualification/sparse_backends/README.md).

## Platform status

| Configuration | Current status |
| --- | --- |
| GNU serial | The maintained serial regression suite passed in the optimized GNU configuration (194 pytest tests), and the component suite passed in both OPT and DEBUG configurations during the staged sequence on macOS arm64 with GNU Fortran 16.2.0. This is not a general cross-platform claim. |
| GNU MPI and OpenMP | Serial, two-rank MPI, and two-thread OpenMP results agreed on the ten-zone qualification problem on macOS arm64 with GNU Fortran 16.2.0 and Open MPI 5.0.10. This checks the selected functional and numerical behavior, not scaling, multi-node placement, binding performance, or hybrid MPI+OpenMP execution. |
| Frontier HIP/ROCm OpenMP offload | The accepted source `64951196032bf4622ee6c887323db33cf1de5beb`, included in this tree, passed Frontier job `5512553` on one MI250X with CPE 25.09, Cray Fortran 20.0.0, ROCm 6.4.2, hipfort 6.4.2, OpenMP target offload, Starkiller EOS, dense solver, and MPI off. The device-probe maximum residual was `3.552713678800501e-17` against a `1e-12` limit. The ten-zone partial batch and six-zone `heat_sn160` comparisons passed; the latter's maximum numerical-limit fraction was `0.16072834133922473`, with nonzero neutrino loss in every zone. CPE 26.03 with ROCm 7.0.2 encountered a hipfort/rocBLAS link incompatibility and is not qualified. |
| Perlmutter CUDA/OpenACC | Supplemental evidence for the same accepted source passed Perlmutter job `58588659` on one A100-SXM4-80GB with NVHPC 26.5, PrgEnv-nvidia 8.7.0, CUDA 13.2, OpenACC, and the cuBLAS pointer-array batched solve. Device-probe residuals were zero. The ten-zone partial batch and `heat_sn160` comparisons passed; the latter's maximum numerical-limit fraction was `0.07912200248018902`, with nonzero neutrino loss in every zone. This evidence applies only to that tested configuration. |
| Other host names | Use the generic machine defaults unless `MACHINE` or other build variables are supplied explicitly. |

The exact commands, source revision, and any launcher-specific options belong
in the issue or pull-request evidence for each run. A successful build alone
does not establish runtime or numerical agreement.

The Frontier and Perlmutter runs used source
`64951196032bf4622ee6c887323db33cf1de5beb`. Upstream staging merge
`1f41d9be38f8d3f171b1cc37f9977a3e815a53ad` has the same source tree.

The accelerator checks above do not establish performance, scaling,
multi-GPU or MPI+GPU behavior, BDF numerical behavior, other GPU generations,
newer ROCm stacks, arbitrary CUDA or NVHPC versions, or general accelerator
portability.

Make can display resolved values through the existing `print-%` target. For
example:

```bash
make --no-print-directory print-CMODE
make --no-print-directory print-MATRIX_SOLVER
```

## Configuration changes and clean builds

Each build directory contains predictable `obj/`, `mod/`, `pp/`, and
`bin/` subdirectories plus `config.txt`. The configuration record preserves
the exact effective selectors, commands, flags, selected library paths, and
external source paths. Reusing the directory with different material settings
fails before preprocessing or compilation and directs the user to clean or
choose another directory.

The line-oriented record accepts ordinary single-line Make values. Embedded
single quotes and line breaks in recorded commands, flags, libraries, or paths
are rejected with a clear diagnostic; select an equivalent spelling or another
build directory.

Use a different readable directory for a different configuration:

```bash
make BUILD_NAME=gnu-opt -j xnet
make BUILD_NAME=gnu-debug CMODE=DEBUG -j xnet
```

`clean` removes only the selected marked build directory and does not require
the old configuration to be restated. `clean-all` removes only marked direct
children of `BUILD_BASE` and requires `CONFIRM_CLEAN_ALL=yes`. Do not clean
while a build uses that directory; run cleaning and building as separate
commands. Mixed cleaning/product goals are rejected.

Pass local configuration choices on the Make command line and leave tracked
defaults unchanged. The main selection variables include:

- `CMODE` and `PE_ENV` for optimization/debug mode and compiler family;
- `MPI_MODE` and `OPENMP_MODE` for distributed and host-threaded execution;
- `GPU_MODE`, `GPU_BACKEND`, `OPENACC_MODE`, and `OPENMP_OL_MODE` for
  accelerator execution and directive model;
- `LAPACK_VER` and `GPU_LAPACK_VER` for CPU and accelerator numerical
  libraries;
- `MATRIX_SOLVER` for the Jacobian and linear solver implementation;
- `EOS` for the equation-of-state implementation.

Each selected path requires its compiler, headers, libraries, and runtime.
Validate support and numerical behavior for the exact combination used.

## Production and utility targets

The default target builds `build/GNU-OPT/bin/xnet` with the tracked defaults.
Common utility builds are:

```bash
make -j net_setup
make -j xnse
```

- `net_setup` preprocesses network data.
- `xnse` is the stand-alone NSE state calculator.

`all` builds the three canonical programs together. Select another solver with
`MATRIX_SOLVER`, for example `make MATRIX_SOLVER=MA48 MA48_DIR=... xnet`.
`xinab` and `xnet_gpu` are unsupported and fail early; accelerator builds use
`xnet` with explicit supported selectors. Target presence records a build
recipe, not a support claim.

## Focused component and executable tests

The Fortran component suite includes a complete serial XNet executable
smoke:

```bash
make -C test/unit
```

It uses the tracked GNU optimized configuration by default, compiles selected
production sources into the ignored `test/unit/build/` directory, and performs
no network access. Its build-net interoperability check always cleans and
builds the requested tracked configuration before resolving and running the
canonical XNet path, so an incompatible configuration cannot be reused. Run
the bounds-checking configuration with:

```bash
make -C test/unit clean test CMODE=DEBUG
```

See `test/unit/README.md` for the tested behavior, narrow test-only
state and stubs, vendored `test-drive` revision and license, update procedure,
and focused effectiveness checks.

## Maintained serial regression suite

The pytest suite under `test/regression/` is the normal serial CPU regression
path. It requires Python 3.11 or newer and the dependency recorded in
`test/regression/requirements.txt`. Build the production programs and pass
their paths explicitly:

```bash
python3 -m pip install -r test/regression/requirements.txt
make BUILD_NAME=regression-serial -j xnet xnse
python3 -m pytest test/regression \
    --xnet-executable="$PWD/build/regression-serial/bin/xnet" \
    --xnse-executable="$PWD/build/regression-serial/bin/xnse"
```

The suite runs `xnet` and `xnse` as external programs in isolated temporary
directories. It checks direct process status, required output, parsed
diagnostics, and numerical comparisons with stated limits. Its current 194
pytest tests include runner, parsing, and effectiveness checks as well as the
physical regression scenarios; they are not 194 separate scientific cases.
This evidence does not by itself establish scientific validity or portability.
See `test/regression/README.md` for the case definitions, requirements,
timeouts, reference provenance, and comparison policy.

## Runtime inputs

The stand-alone driver reads a file named `control` from its working directory.
`source/xnet_controls.F90` locates labeled blocks and reads the values within
each block in a specific order. Ordering and format changes can affect existing
inputs.

Legacy problems assemble a control file by joining a `test/test_settings*`
file with a matching `test/Test_Problems/setup_*` file. The setup file refers
to thermodynamic trajectories, initial abundances, and nuclear data under
`test/Data_*`. Source code remains authoritative for the values read and their
meaning.

`test/Data_*` directories contain pre-built nuclear networks. Network
preprocessing work should identify whether these tracked files are inputs,
generated results, or comparison data before changing them.

## Legacy shell test drivers

The older shell drivers remain useful for investigation and historical problem
runs, but they are secondary to the maintained pytest suite and have unreliable
pass/fail reporting.

`test/test_xnet.sh`:

- selects problems by numeric ID;
- combines settings and setup files into `test/control`;
- runs a supplied executable or `build/GNU-OPT/bin/xnet`;
- looks for an MPI build at `build/GNU-OPT-MPI/bin/xnet`; set `XNET_MPI` to use
  another predictable build directory;
- moves diagnostics into `test/Test_Results/`;
- removes timer sections before comparison;
- prints a warning and writes `diff_*` when results differ.

The wrapper invokes the selected executable without capturing or propagating
its exit status. Its final commands normally leave a zero wrapper status even
when the executable or intermediate file operations fail, and a request with
no recognized problem ID can run no problem. Treat the wrapper status as
uninformative. Confirm program invocation, direct program status, expected
diagnostics, output production, and numerical agreement separately.

A clean checkout currently supplies no tracked comparison files under
`test/Test_Problems/Results/`. The script creates that directory. Ordinary
problem paths then attempt comparison against absent results. The `xnse` path
copies the current result into the reference location when the expected file
is absent, and those diagnostic files are ignored by `.gitignore`. Establish
an independent comparison result with recorded provenance before claiming
numerical agreement.

Legacy problem drivers remain under `test/`; the production Makefile does not
duplicate them as recursive build targets. Some legacy drivers create a
Helmholtz-table symlink, write `control`, create result directories, move
diagnostic files, create comparison files, or run preprocessing inside a
tracked data directory.

`test/test_xnet.csh` is an older driver. Use it as historical information and
verify every command needed for a current task.

After any legacy run:

1. capture direct program status or state why the wrapper obscures it;
2. confirm that the requested problem ran and produced expected diagnostics;
3. inspect generated `diff_*` files and identify whether an independent
   comparison result exists;
4. record the comparison result's provenance;
5. state the quantities and tolerances used for the conclusion;
6. inspect ignored files as well as `git status` for generated or modified
   files;
7. preserve established reference data unless the issue explicitly changes
   it.

## Evidence for new work

For a defect fix, identify or add evidence that distinguishes the faulty
behavior from the corrected behavior. Prefer demonstrating that the check
fails before the fix when practical.

Choose the smallest check that exercises the requirement. Record:

- the exact command and working directory;
- compiler version and relevant build variables;
- input problem and data source;
- expected and observed behavior;
- comparison method and tolerances;
- generated files or reference results;
- checks that remain for other compilers, parallel modes, accelerators, or
  facilities.

Use `docs/development/scientific-validation.md` for changes that affect
physics, numerical behavior, tolerances, convergence, or performance.

## Generated and machine-specific files

Keep build objects, module files, executables, diagnostic outputs, comparison
files, temporary control files, local installation paths, and machine-specific
settings out of commits. Production objects, module files, preprocessing
output, `config.txt`, and executables remain below `BUILD_DIR`, but legacy
problem drivers can still write into `test/`; review ignored files as well as
`git status` after runs.
