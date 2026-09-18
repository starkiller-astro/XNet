# Optional sparse-backend qualification

These manual, dependency-bound checks qualify the real serial CPU backends
described below. They do not run as part of `make -C test/unit`, download
solver software, or turn an unavailable solver into a stubbed success.

The opt-in component targets compile the production `xnet_jacobian` implementation
against the real library and reuse the three-equation component fixture. Base
and self-heating solves must satisfy

```text
max(abs(matmul(A,x)-b)) <= 1e-12 * (1 + max(abs(b)))
```

for two distinct right-hand sides and two zones. The runner also requires the
tracked solver controls to work, a controlled real-solver error to reach the
solver implementation's production fatal-status path, and a `1e-11` controlled result
perturbation to fail the residual check. MA48 uses a singular matrix for the
status probe. oneMKL uses its matrix checker with a controlled invalid CRS
column index because its default pivot perturbation accepts a numerically zero
matrix. The result perturbation is applied only to the test's observed result
after the real solver returns.

The complete-executable check builds dense and sparse executables from one Git
revision, runs the existing six-zone `heat_sn160` self-heating case in separate
work directories, and reanchors the existing `xnet-comparison-v1` tolerances to
the newly generated dense endpoint. Non-exact difference bounds are multiplied
by two for this cross-direct-solver comparison; exact fields and the independent
mass-normalization requirement remain unchanged. This factor-of-two allowance
covers the additional Newton/timestep rounding path without defining
solver-specific canonical values. The runner writes an untracked JSON report with
process status, executable hashes, every compared field, complete-composition
norms, and the largest selected-species difference. It does not create or
update a canonical endpoint reference. It rejects dense and sparse executables
with identical SHA-256 hashes before staging or running the problem.

## HSL MA48

Use only MA48 source obtained by the maintainer under an applicable HSL
licence. The recorded qualification used the archived Fortran 77 MA48
2.2.0 source identified by the
[HSL catalogue](https://www.hsl.rl.ac.uk/catalogue/ma48.html). HSL source is
non-redistributable in this repository: keep it outside the worktree and never
commit, copy, archive, upload, or attach it.

From the repository root, with `HSL_MA48_DIR` naming the private source
directory:

```bash
make -C test/unit clean real-ma48-test \
  HSL_MA48_SOURCE="$HSL_MA48_DIR/MA48.f"

make BUILD_NAME=ma48-dense MATRIX_SOLVER=dense -j xnet
make BUILD_NAME=ma48-sparse MATRIX_SOLVER=MA48 MA48_DIR="$HSL_MA48_DIR" -j xnet

python3 test/qualification/sparse_backends/compare_heat_sn160.py \
  --backend=ma48 \
  --dense-executable="$PWD/build/ma48-dense/bin/xnet" \
  --sparse-executable="$PWD/build/ma48-sparse/bin/xnet" \
  --work-directory=/tmp/xnet-ma48
```

The recorded qualification applies only to the named source version, compiler,
host, and serial CPU configuration. Possession of MA48 source does not grant
redistribution rights.

## Intel oneMKL PARDISO

The recorded qualification used oneMKL PARDISO, not the distinct
standalone PARDISO ABI. Intel documents oneMKL licensing in its
[oneMKL License FAQ](https://www.intel.com/content/www/us/en/developer/articles/tool/onemkl-license-faq.html).
The oneMKL implementation keeps one solver handle per local batch slot. This preserves
each concurrently evolved zone's analysis/refactorization state independently;
sharing one handle across the two test zones corrupted the first
stored factorization under oneMKL 2026.1 even though both solver calls returned
success.
Initialize the locally installed oneAPI environment before each build. For
example, when using Intel's standard `setvars.sh` installation:

```bash
source /path/to/oneapi/setvars.sh

make -C test/unit clean real-pardiso-mkl-test LAPACK_VER=MKL

make BUILD_NAME=pardiso-mkl-dense MATRIX_SOLVER=dense LAPACK_VER=MKL -j xnet
make BUILD_NAME=pardiso-mkl-sparse MATRIX_SOLVER=PARDISO_MKL LAPACK_VER=MKL -j xnet

python3 test/qualification/sparse_backends/compare_heat_sn160.py \
  --backend=pardiso-mkl \
  --dense-executable="$PWD/build/pardiso-mkl-dense/bin/xnet" \
  --sparse-executable="$PWD/build/pardiso-mkl-sparse/bin/xnet" \
  --work-directory=/tmp/xnet-pardiso-mkl
```

`Makefile.internal` accepts the legacy
`$MKLROOT/tools/mkl_link_tool` location and the current
`$MKLROOT/bin/mkl_link_tool` location. Record the emitted compile/link lines so
the selected interface, sequential threading layer, and library version remain
reviewable.

## Standalone PARDISO support

Standalone PARDISO is unsupported and unqualified. The production build's
`PARDISO_MKL` selector names the maintained oneMKL implementation explicitly.
The maintainer does not have an approved compatible standalone installation,
and oneMKL success does not qualify the standalone ABI. Do not substitute the
component-test numerical stub or oneMKL library for a standalone qualification.

## Evidence limits

For each real qualification, record the exact Git revision, date, host,
compiler, library/source version, clean build commands and link line, component
process statuses, singular-path diagnostic, comparison report quantities, and
generated-file cleanup. These checks establish build, serial runtime,
solver setup and error behavior, known-system residual, and dense-versus-sparse endpoint
agreement. They are not independent scientific validation, performance or
scaling evidence, or qualification of MPI, OpenMP, accelerators, alternate
compilers, solver implementations, versions, or hosts.
