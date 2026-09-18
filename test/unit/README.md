# Component tests

Run the component suite from the repository root with:

```bash
make -C test/unit
```

The default uses the GNU optimized configuration. A bounds-checking run is:

```bash
make -C test/unit clean test CMODE=DEBUG
```

Generated files, test objects, and test executables are written below the
ignored `test/unit/build/` directory. The suite also builds a production
`xnet` executable at `build/unit-opt/bin/xnet` (or `unit-debug`) for the
generated-network smoke test.

## What is checked

The direct Fortran tests compile selected production routines with small test
modules that supply only the state needed by those routines. They check:

- safe exponentiation, abundance normalization, and output suffix formatting;
- scalar and vector trajectory interpolation, including inactive zones;
- scalar and vector abundance moments;
- neutrino-history interpolation at knots, endpoints, and outside the supplied
  time range;
- screened and unscreened NSE software behavior on a compact eight-species
  input, including mass, charge, finite composition, repeatability, and solver
  counters;
- STARKILLER and Bahcall EOS interfaces; and
- ASCII model-input success, missing-file status, and malformed input.

The compact NSE checks exercise software behavior and basic physical
invariants. A separate executable loads the retained 489-species `torch489`
network and three independently calculated complete compositions under
`test/nse_validation/`, then calls the production `nse_initialize` and
`nse_solve` routines. The scientific formulation, exact inputs, acceptance
limits, reproduction steps, and finite-network limitations are documented in
`test/nse_validation/README.md`.

The test programs use shared module state, so Test Drive runs them serially.
This restriction applies only to the test process; production sources are
compiled with the selected OpenMP options.

## Network preprocessing

The preprocessing checks use the synthetic public fixture under
`fixtures/preprocess/`. They compare direct `net_preprocess` results with the
stand-alone `net_setup` program and verify:

- nuclear and reaction translation;
- weak and reverse flags;
- repeated-participant multiplicities and recomputed Q values;
- reaction matching and sparse matrix indices;
- the generated files read by XNet; and
- nonzero status for truncated or inconsistent inputs.

The checks compare the contents interpreted by the readers rather than raw
compiler-dependent sequential-unformatted bytes. Temporary copies are mutated
to exercise failure paths; tracked network data are not changed.

The generated-file checks cover:

| File | Check |
| --- | --- |
| `nuc_data` | Species, nuclear values, temperature grid, partition functions, and spins loaded by the production nuclear-data reader. |
| `nets3`, `nets4` | Reactions, indices, participants, and multiplicities loaded by the production reaction reader. |
| `match_data` | Participants, signs, Q values, weak flags, and reverse associations loaded by the production match reader. |
| `sparse_ind` | Compressed-row indices and reaction-to-entry maps loaded by the PARDISO sparse-data reader. |
| `ab_blank` | Species order and zero abundances. |
| `match_read` | Participants, descriptors, endpoint coordinates, scales, and record order. |
| `matr_shape` | Trial right-hand side, coordinates, and values reconstructed from the sparse reaction maps. |
| `net_desc` | Description, dimensions, chapter ranges, non-REACLIB counts, and sparse widths. |
| `net_diag` | Q corrections, matched-reaction diagnostics, and sparse-row diagnostics. |

The PARDISO entry points are test stubs because these checks exercise only the
production sparse-file reader. They do not test a PARDISO solve.

## `build_net` and generated-network readers

`fixtures/build_net/README.md` describes the five-species network assembled
from small generated mass, partition-function, REACLIB, and weak-rate inputs.
The inputs also contain an unrequested species and a reaction with an
unavailable participant, allowing the checks to verify selection behavior.
No private neutrino data are used.

The runner builds identical positive cases twice and compares their parsed
contents. Its mass catalog includes a valid selected mass with a `#` in an
unrelated field, an unselected unavailable mass, a selected unavailable mass,
and a malformed row containing `#`; these distinguish an unavailable required
mass field from other uses of the character. It also checks a configuration
without weak rates and requires nonzero status for duplicate, blank,
unavailable, malformed, missing, or truncated inputs. Mutated generated files
exercise:

- reaction participants, species order, masses, and Q values;
- `sunet`/`netwinv` count and order mismatches, including the last species;
- `nets4` count and species-order mismatches; and
- malformed or inconsistent `match_data` reaction counts.

The accepted output is passed through the real `net_setup` program and the
production nuclear, reaction, match, and sparse-data readers. A short
one-zone production `xnet` run then checks that the generated network is
usable. The run requires normal target-time completion and final counters; it
is an interoperability smoke test, not a stored scientific endpoint.

## EOS checks

The STARKILLER and Bahcall implementations are compiled separately against the
same component tests. The checks cover initialization, finite thermodynamic
outputs, and expected failure for a missing Helmholtz table. The STARKILLER
fixtures include values from the retained Timmes reference driver; the fixture
README records how they were generated. These checks do not change the EOS
equations, interpolation, constants, or table interpretation.

## Test Drive dependency

The suite vendors the single-file Test Drive v0.5.0 release at commit
`fd66b4bca683c5fa5d92536075734f0792824d37`:

- `vendor/test-drive/testdrive.F90`, SHA-256
  `e8765129ba304f28c4bcfc20860cb49e0046e76527e4c240eeb54a5fea22837d`;
- `vendor/test-drive/LICENSE-MIT`, SHA-256
  `d34e0235cb56e251ea1c23f9c803857267d083459aeedcd06b538c0335d69e46`.

Test Drive permits redistribution under Apache-2.0 or MIT terms; XNet retains
the MIT license beside the source. Normal build and test execution makes no
network request. To update the dependency, review an upstream release, replace
the source and license, update the version, commit, and hashes above, then run
both optimized and debug component tests.

Test Drive v0.6.0 was evaluated with GNU Fortran 16.1.0 but test-suite
construction aborted with `SIGABRT`; no root cause is claimed here.

## Limits

These tests cover the GNU serial CPU configuration and small synthetic inputs.
They do not qualify MPI, accelerators, optional external sparse solvers, other
compilers, production nuclear-data quality, or long scientific calculations.
Passing the suite demonstrates the named component and interoperability
behavior only.
