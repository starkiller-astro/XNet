# Validation plan

## Purpose and evidence maturity

XNet needs evidence that distinguishes “builds and runs” from “implements the
intended science.” This inventory records current evidence without treating
current output as a scientific reference answer or inventing tolerances.

The maturity terms used here are:

- **historical**: previously reported to work;
- **available**: the code path exists;
- **characterized**: reproduced in a recorded configuration;
- **validated**: compared with an accepted scientific reference answer or
  analytic result;
- **supported**: the project commits to maintaining it.

Evidence advances only when the exact configuration, command, result, and
expected-result source are recorded.

## Test categories

1. **Runner self-tests** verify that the test mechanism returns reliable
   process status.
2. **Unit tests** exercise dependency-light routines and implementation
   contracts.
3. **Build smoke** verifies that named artifacts compile and link in an exact
   configuration. It is not scientific validation.
4. **Characterization tests** record reproducible current behavior without
   claiming that behavior is correct.
5. **Scientific validation tests** compare a scientific quantity with an
   accepted analytic result, published benchmark, or independently verified
   result.
6. **Integration/regression tests** cover end-to-end behavior. A regression
   baseline is authoritative only when its provenance and scientific review
   are recorded.

## Reference-answer requirements

A scientific validation or regression proposal must identify the applicable
configurations; the reference answer and its source; input provenance; units;
the comparison quantity and norm; separately justified absolute and relative
tolerances where applicable; sensitivity near zero; expected floating-point
variation; runtime; and intended CI tier. Current output alone is not a source
of correctness.

## Evidence inventory

| Item | Category | Maturity | Evidence | Configuration | Automation | Limitation |
| --- | --- | --- | --- | --- | --- | --- |
| Test-drive pass/fail propagation | Runner self-test | characterized | Expected process status; nested Make target observes failure | GNU Fortran, serial | `serial-unit-tests` | Runner behavior only |
| `dp`, `sp`, `i4`, and `i8` kind/storage contracts | Unit | characterized | `iso_fortran_env` and `storage_size` | GNU Fortran, serial | `serial-unit-tests` | Exact kind/storage bug-fix characterization only |
| `real(dp)` evaluation of representative numerical and physical constants | Unit | characterized | Fortran kind semantics and the existing `third`, `ln_2`, and `clt` expressions and values | GNU Fortran, serial, ordinary default real kind; isolated `xnet_types` and `xnet_constants` compile | `serial-unit-tests` | Portability and implementation-correctness characterization only; not a scientific validation or update of physical values. |
| `fd0h` order-zero identity at `eta = [-4, -1, 0, 1, 4]` | Unit analytic validation | validated | [Fukushima (2015)](https://doi.org/10.1016/j.amc.2015.03.009), equation 17; section 3.7 limits use of the direct identity, and Table 27 reports order-zero errors in `2^-53` units. The test uses `32 * EPSILON(1.0_dp)`, with Fortran `EPSILON` equal to `2^-52` in the recorded binary64 configuration, to also cover rounding in the independently evaluated intrinsic reference. | GNU Fortran, serial, normal dependency-light unit configuration | `serial-unit-tests` | Limited to these five moderate inputs; not a full approximation error envelope, FFN/EOS/network validation, compiler/project support claim, or scientific baseline. |
| `reset_timers()` reset of public `timer_output` | Unit | characterized | Exact zero after a nonzero sentinel; reset contract for a reported summary accumulator. | GNU Fortran, serial, normal dependency-light unit configuration | `serial-unit-tests` | Implementation contract only; later summaries do not retain `Output` from an earlier reset epoch. |
| Scalar/vector `safe_exp()` ordinary and clamped behavior | Unit | characterized | Intrinsic `exp`, clamp parameters, IEEE finite check, epsilon-scaled comparisons | GNU Fortran, serial | `serial-unit-tests` | Dependency-light routine behavior only |
| Double-precision blank/tab numeric token parsing | Unit | characterized | Literal values after `replace_tabs()` and zero with `pos = 0` for exhausted input | GNU Fortran, serial | `serial-unit-tests` | Input-token behavior only |
| ASCII case folding before abundance-name lookup | Characterization | characterized | Exact `Fe56` to `fe56` result from `string_lc()` | GNU Fortran, serial | `serial-unit-tests` | ASCII example only; not a locale, Unicode, file, or lookup-path claim |
| Generated filename suffix formatting | Characterization | characterized | Literal one- and two-digit results from `name_ordered()` | GNU Fortran, serial | `serial-unit-tests` | Not arbitrary-width, path, filesystem, invalid-input, or compatibility validation |
| `xnet`, `net_setup`, and `xnse` compile/link | Build smoke | characterized | Clean command and executable presence | GNU Fortran; debug; serial; dense; NETLIB; STARKILLER | `serial-build-smoke` | Build evidence only |
| Legacy problem execution paths | Smoke/characterization | available | Inputs and scripts in `test/`; no tracked trusted results for exercised cases | Multiple historical configurations | Not a required gate | Historical use reported |
| Legacy numerical comparisons | Regression | available | Mismatches do not propagate failure; reference provenance unavailable | Historical | Not a required gate | Not trusted as a scientific gate |

Except for the `fd0h` row, these rows characterize only the exact configurations
exercised by recorded local or CI runs. The `fd0h` row validates only the stated
analytic identity at its five moderate inputs. No row is a broad/full-domain/
network scientific validation or support commitment.

## Commands

Run the characterized build from the repository root:

```bash
make -C source clean
make -C source -j2 \
  CMODE=DEBUG PE_ENV=GNU \
  MPI_MODE=OFF OPENMP_MODE=OFF GPU_MODE=OFF \
  MATRIX_SOLVER=dense LAPACK_VER=NETLIB EOS=STARKILLER \
  xnet net_setup xnse
```

Run the isolated unit and runner self-check targets with:

```bash
make -C source test_unit
make -C source test_unit_selfcheck
make -C tests strict_kind
```

Each unit command begins from a clean dedicated directory under
`tests/build/`, so the targets can run concurrently without sharing compiler
artifacts. `test_unit` runs both the promoted-default-real normal suite and the
isolated strict-kind constants check. The direct `strict_kind` command compiles
only vendored test-drive, `xnet_types.F90`, `xnet_constants.F90`, and the test
program with `-g -Og -cpp -fimplicit-none -fallow-argument-mismatch
-ffree-line-length-none`; it deliberately omits `-fdefault-real-8` and
`-fdefault-double-8`.

The self-check deliberately runs a failing fixture through a nested Make target
and succeeds only when that target returns nonzero.
