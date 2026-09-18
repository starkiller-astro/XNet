# XNet pytest regression cases

This directory contains the maintained XNet regression cases. It
exercises the compiled XNet program as an external process; it is not a Python
binding, a scientific-validation suite, or a replacement for all legacy cases.
The migrated cases are the serial, CPU-only `tnsn_alpha` and `tnsn_torch47`
trajectory calculations, `heat_alpha`, `heat_sn160`, and `bdf_sn160`
self-heating calculations, the `batch_alpha` serial zone-batching case, and one
SN160 trajectory initialized through NSE. A separate process test exercises
the standalone `xnse` program.

## Prerequisites and command

Use Python 3.11 or newer and install the single test dependency:

```bash
python3 -m pip install -r test/regression/requirements.txt
```

Build XNet and `xnse` separately, then select both executables explicitly:

```bash
make BUILD_NAME=regression-serial -j xnet xnse
python3 -m pytest test/regression \
    --xnet-executable="$PWD/build/regression-serial/bin/xnet" \
    --xnse-executable="$PWD/build/regression-serial/bin/xnse"
```

There is no default executable name and no `XNET_EXECUTABLE` fallback. The
helper tests, which do not run XNet, can be selected independently:

```bash
python3 -m pytest test/regression/test_xnet_regression.py
```

The suite enforces a 30-second per-case process timeout. Use
`--xnet-timeout=SECONDS` to change it deliberately.

## Isolated execution and retained output

pytest supplies a new empty temporary directory for each case. The runner
copies the complete `control` input into it, creates a local writable network
directory containing the case's declared tracked source inputs, and stages
every trajectory, Helmholtz EOS table, and any explicitly declared nested
input at a safe relative destination. Duplicate destinations and unsafe paths
are rejected, so staging cannot silently replace an input. XNet runs with that
directory as its current directory.
Network-preprocessing products therefore stay inside the temporary directory
rather than mutating `test/Data_alpha`, `test/Data_torch47`, or
`test/Data_SN160`. The runner also reads each staged network's `sunet` and
requires its ordered, unique species list to match the case declaration before
execution. A nonempty work directory is a setup failure, so old diagnostics
cannot satisfy a new run.

Every invocation records `xnet.stdout.txt`, `xnet.stderr.txt`, and
`xnet.status.txt` beside the XNet outputs. Runs whose final diagnostic is
parsed successfully also record `composition_error_norms.json`. Failure
messages give the work directory path. pytest retains recent temporary
directories. To choose a stable diagnostic location for a run, use its
standard option, for example:

```bash
python3 -m pytest test/regression \
    --xnet-executable="$PWD/build/regression-serial/bin/xnet" \
    --xnse-executable="$PWD/build/regression-serial/bin/xnse" \
    --basetemp=/tmp/xnet-regression-output
```

The runner classifies invalid definitions, paths, and preparation as setup
failures; timeouts, signals, nonzero status, and missing/empty required output
as execution failures; missing or malformed diagnostic structure as parsing
failures; and out-of-policy values or invariants as comparison failures.

## Case inputs and legacy provenance

Legacy ID 1 in `test/test_xnet.sh` names `tnsn_alpha` and calls
`do_test_small`, which historically concatenates `test/test_settings_small`
and `test/Test_Problems/setup_tnsn_alpha` into `test/control`. It assumes the
script is run from `test/`, defaults silently to `../source/xnetp` if no
executable file argument is recognized, runs ten identical zones, reads
`Data_alpha/ab_co` and `Test_Problems/th_sn1aflame`, and expects
`net_diag01`, ten `ev_tnsn_alpha_*` ASCII histories, and ten
`ts_tnsn_alpha_*` binary histories. The legacy driver moves only
`net_diag01` into `Test_Results` for comparison.

The committed `cases/tnsn_alpha/control` is the one-time concatenated input.
Trailing whitespace was removed, and two deliberate path changes support
isolated execution: the `Test_Results/` prefixes were removed from the ASCII
and binary output roots, and each `Test_Problems/th_sn1aflame` path became
`th_sn1aflame`. All numerical controls, ten zone inputs, the `Data_alpha`
path, and the `Data_alpha/ab_co` paths are unchanged. The runner provides the
adjusted paths in the isolated directory. Normal execution does no fragment
concatenation.

The Starkiller EOS initializes even though `tnsn_alpha` disables self-heating,
so `tools/starkiller-helmholtz/helm_table.dat` is also an explicit required
input. Before the serial abort path returned a nonzero status, a missing table
could produce a fatal message, a partial `net_diag01`, and process status zero.
The runner still requires all expected output and complete final records
instead of relying on status or stderr keyword heuristics alone.

Legacy ID 51 in `test/test_xnet.sh` names `heat_alpha` and calls
`do_test_heat`, which historically concatenates `test/test_settings_heat` and
`test/Test_Problems/setup_heat_alpha` into `test/control`. It enables screening
and self-heating for six zones initialized from `Data_alpha/ab_co`. The zones
use `th_co_burn_1` through `th_co_burn_6`, covering densities from
`1.0e7` through `3.16227766e9 g/cm^3` and target times from 10 seconds through
`1.0e-5` seconds. The legacy run expects `net_diag01`, six
`ev_heat_alpha_1` through `ev_heat_alpha_6` ASCII histories, and six matching
`ts_heat_alpha_1` through `ts_heat_alpha_6` binary histories. It moves only
`net_diag01` for its warning-only comparison. No tracked historical
`net_diag_heat_alpha` reference is present in a clean checkout.

The committed `cases/heat_alpha/control` is the one-time concatenated input.
Trailing whitespace was removed, the `Test_Results/` prefixes were removed
from both output roots, and the six `Test_Problems/` prefixes were removed from
the trajectory paths. All numerical controls, zone ordering, `Data_alpha`
path, and abundance paths are unchanged. The runner stages the six adjusted
trajectory paths in the isolated directory. XNet uses the number of digits in
the largest zone number for output suffixes, so this six-zone case uses `_1`
through `_6`, while the ten-zone `tnsn_alpha` case uses `_01` through `_10`.

The paired SN160 study has separate cases for Backward Euler legacy ID 53,
`heat_sn160`, and BDF legacy ID 54, `bdf_sn160`. Both remain members of
aggregate self-heating ID 50. Each maintained integrator has its own
characterization reference and tolerance policy; the paired comparison is
diagnostic and does not make one integrator's execution history a requirement
for the other.

Legacy ID 53 calls `do_test_heat`, which concatenates
`test/test_settings_heat` and `test/Test_Problems/setup_heat_sn160`. It runs
six serial zones with weak reactions, screening, self-heating, runtime nuclear
data processing, no neutrino reactions, `Data_SN160/ab_co`, and ordered
trajectories `th_co_burn_1` through `th_co_burn_6`. The required fresh output
is `net_diag01`, six `ev_heat_sn160_1` through `_6` ASCII histories, and six
matching `ts_heat_sn160_1` through `_6` binary histories. The 14 species in the
ASCII-output control do not limit the diagnostic: `net_diag01` records all 160
species in the exact order of `test/Data_SN160/sunet`.

The committed `cases/heat_sn160/control` is that one-time concatenation with
trailing whitespace removed. Its only semantic-preserving path edits remove
the `Test_Results/` prefixes from both output roots and the `Test_Problems/`
prefixes from all six trajectories. Numerical controls, zone block size 1,
zone order, `Data_SN160`, abundance paths, and requested ASCII species are
unchanged. The case stages only `sunet`, `netsu`, `netweak`, `netwinv`, and
`ab_co` into a writable temporary `Data_SN160`; generated preprocessing files
never enter the tracked source directory. Screening and self-heating require
the tracked `tools/starkiller-helmholtz/helm_table.dat`.

The complete legacy control differences are:

| Control | Backward Euler ID 53 | BDF ID 54 | Classification |
| --- | ---: | ---: | --- |
| Integration choice | `1` | `3` | BDF-required numerical behavior |
| Maximum iterations per step | `5` | `10` | BDF-required numerical behavior |
| Convergence-condition flag | `0` | `3` | BDF-required numerical behavior |
| Lower abundance cutoff | `1e-30` | `1e-99` | BDF-required numerical behavior |
| Legacy zone block size | `1` | `4` | batching normalization to `1` below |
| ASCII/binary root | `heat_sn160` | `bdf_sn160` | output naming |

Legacy ID 54 calls `do_test_bdf`, which concatenates
`test/test_settings_bdf` and `test/Test_Problems/setup_bdf_sn160` once. The
committed `cases/bdf_sn160/control` removes trailing whitespace, removes
`Test_Results/` from both output roots, removes `Test_Problems/` from the six
trajectory paths, and intentionally changes the zone block size from 4 to 1.
The focused normalized-concatenation test permits exactly those path edits and
the batching change. All six zones, their order, the five SN160 source inputs,
`Data_SN160/ab_co`, the six trajectories, the EOS table, the 14 requested
ASCII species, and all other physical and numerical controls are unchanged.
The required fresh outputs are `net_diag01`, `ev_bdf_sn160_1` through `_6`,
and `ts_bdf_sn160_1` through `_6`.

The BDF iteration, convergence, and abundance-floor values are maintained
solver behavior. Holding block size at 1 keeps batching out of the integrator
comparison; no evidence from the isolated runs indicated a BDF semantic
dependency on the legacy block size of 4. Testing `szbatch > 1` remains
separate work. `source/xnet_evolve.F90` imports `solve_bdf` and dispatches
`isolv == 3` to it, `source/net.F90` calls `bdf_init` for the same choice, and
the normal root-level production build contains
`xnet_integrate_bdf.o`. `source/xnet_controls.F90` replaces both input
abundance- and temperature-change timestep limits with the effective value
`1e10` for `isolv == 3`; the committed reference records that effective state.

The maintained solver is Backward Differentiation Formula (BDF). Integration
choice 3 selects BDF; the obsolete Bader-Deuflhard (BD) choice 2 is not built.
The commented `Case (2)` in `source/xnet_evolve.F90` and the absent
`xnet_integrate_bd.o` production object provide additional confirmation.

Legacy ID 2 names `tnsn_torch47` and calls `do_test`, which concatenates
`test/test_settings` and `test/Test_Problems/setup_tnsn_torch47`. The settings
activate only zone 1, serial runtime network processing, Backward Euler, no
screening, and no self-heating. The setup nevertheless retains ten identical
`Data_torch47/ab_co` and `Test_Problems/th_sn1aflame` input pairs; the complete
standalone control preserves all ten records while XNet consumes the first
pair for the one active zone. Its ASCII history requests 14 selected species,
but `net_diag01` records the complete ordered 47-species network from
`test/Data_torch47/sunet`.

The committed `cases/tnsn_torch47/control` changes only the two output roots
by removing `Test_Results/` and the ten trajectory records by removing
`Test_Problems/`. Numerical and physical controls, `Data_torch47`, abundance
paths, repeated records, and the 14-species ASCII-history selection are
unchanged. The case stages `sunet`, `netsu`, `netweak`, `netwinv`, and `ab_co`
in a writable temporary `Data_torch47`. Runtime preprocessing generated
`ab_blank`, `match_data`, `match_read`, `matr_shape`, `net_desc`, `net_diag`,
`nets3`, `nets4`, `nuc_data`, and `sparse_ind` there without changing the
tracked inputs. The required fresh outputs are `net_diag01`,
`ev_tnsn_torch47_1`, and `ts_tnsn_torch47_1`.

No tracked `test/Test_Problems/Results/net_diag_tnsn_torch47` exists in a
clean checkout. The committed endpoint is therefore a characterization of the
recorded build and inputs, not historical or independent scientific truth.

### `batch_alpha`: serial zone-batching characterization

The maintained `batch_alpha` case corresponds to ID 61 within legacy aggregate
ID 60. The legacy driver maps ID 61 to `do_test_batch`, which concatenates
`test/test_settings_batch` with `test/Test_Problems/setup_batch_alpha`; ID 62
is the separate `batch_torch47` case. The former source-directory Makefile's
`test_batch` target invoked ID 62, so that target is historical context rather
than the definition of this regression.

The committed standalone control preserves ID 61's 16 serial zones,
`nzbatchmx = 4`, Backward Euler, self-heating, screening, weak reactions,
runtime nuclear-data processing, alpha network, and 0.1-second target per
zone. Its only edits remove `Test_Results/` from the history roots. It retains
the legacy prefix rules: XNet expands `Data_alpha/ab_batch/ab_batch_` and
`Test_Problems/th_batch/th_batch_` with two-digit global zone suffixes 01–16;
the output roots similarly produce `ev_batch_alpha_01`–`_16` and
`ts_batch_alpha_01`–`_16`.

The case declares and safely stages the four alpha network inputs, all 16
abundance files, all 16 trajectory files, and the Helmholtz table into an
empty work directory. Nested relative destinations are allowed; traversal,
absolute destinations, duplicates, missing sources, and replacement of an
already staged destination are rejected. Preprocessing writes only into the
isolated writable `Data_alpha` directory.

`nzbatchmx` is the configured maximum zones per batch; `nzevolve` is
`nzbatchmx * nthread`; `nzbatch` is the active count in the current batch;
and `szbatch` is that batch's global starting zone, not the configured block
size. `zb_lo` and `zb_hi` select a thread's local range in `nzevolve`, while
`lzactive` masks inactive local slots. `net.F90` uses ceiling division for the
batch count and maps a local index back to global zone number with
`szbatch`. The permanent case therefore requires precisely four ordered
diagnostic groups: zones 1–4, 5–8, 9–12, and 13–16. Each group has four
ordered End/composition records, one matching four-row Counters table, and
one delimited timer section. The parser rejects a sequence of one-zone groups
as evidence for this case.

End values are batch-level loop diagnostics, whereas the Counter values are
per-zone solver counters. In the canonical run every End value is 605; TS
counters repeat 553, 605, 1, and 573 by batch-local slot. They are recorded
in the characterization reference but are not numerical gates and are not
required to be equal. The four repeated input quartets (1/5/9/13,
2/6/10/14, 3/7/11/15, and 4/8/12/16) had identical parsed endpoints in fresh
canonical runs, so the runner requires exact endpoint equivalence as a
state-leakage and association assertion.

The committed reference is a serial-CPU characterization, not scientific
validation. It records all 16 complete ordered 14-species endpoints,
requested/achieved time, temperature, density, electron fraction, End and
solver-counter diagnostics, source provenance, and hashes. It uses the
existing `xnet-comparison-v1` exact policy for observed stable printed values,
selected species, and complete-vector L1/L-infinity gates; L2 remains
diagnostic-only. Normal execution only reads this reference.

Implementation evidence also ran the same inputs at `nzbatchmx = 1` and 6.
Block size 1 produced 16 singleton diagnostic groups and the same parsed
endpoints; block size 6 produced exactly groups 1–6, 7–12, and 13–16, with no
zones 17 or 18 and no inactive-slot records. These are temporary
software-equivalence checks, not registered portability or scientific cases.
Focused parser fixtures permanently cover the short final four-zone group.
The results are limited to the tested default serial CPU configuration.

## Comparison policy and reference status

The historical script deletes the `Timers Summary:` heading plus the next 14
lines, performs an exact whole-file diff, warns on differences, and normally
returns success. A clean checkout has no tracked
`test/Test_Problems/Results/net_diag_tnsn_alpha`, so no historical numerical
reference with compiler or platform provenance is available for either
trajectory case here.

Each case has a `reference/final_state.json` characterization baseline
generated from the tracked default serial build at its recorded revision.
These are not independently validated scientific results. Each JSON file
records the known compiler, platform, build selections, and input paths.
Normal tests only read these files and never create or replace them. Cases
that declare a reference schema bind the exact characterization status,
generating revision and date, toolchain and resolved build selection, legacy
assembly and solver provenance, documented input normalization, effective
controls, the complete input inventory, and SHA-256 hashes before
creating the work directory or executing XNet. Case identity and reference
schema association are also setup checks, so a swapped or incomplete reference
cannot be masked by a later process failure.

The parser requires ordered final records and zone-associated counters for
every case zone, the case-declared complete 14-, 47-, or 160-species
structure, one
delimited timer section per zone, and finite values. It has no default network
species list. Case setup independently requires the same ordered species in
`sunet`, so the network input, parsed diagnostic, and reference must agree. It
rejects negative mass fractions and applies a coarse
structural `|sum(X)-1|` check. The `tnsn_alpha` bound is `2.2e-8`.
`heat_alpha` records a zone-specific bound from `7.60556005e-9` through
`1.81105e-8`, and `tnsn_torch47` uses `1.7e-8`. Each bound is
the sum of the half-last-place rounding bounds for that zone's printed
baseline values. These bounds are not derived from
XNet's per-step Newton mass-convergence control and are not claimed as
scientific-validation thresholds.

XNet emits individual mass fractions, not a separate aggregate composition-sum
record. The regression parser sums those emitted values in output order and
compares the result with the canonical complete-vector sum through the normal
`mass_fraction_printed_sum` policy. It separately uses `math.fsum` for the
structural normalization-to-one check controlled by `mass_fraction_sum_atol`.
For `heat_sn160` and `bdf_sn160`, the printed-sum limits cover the documented
three-configuration envelope with a compact margin; the normalization bounds
retain the existing formatting-aware structural allowance. Neither is the BDF
`iconvc == 3` weighted RMS convergence norm or a check of `rtol`, `atol`, or
`ymin`.

Each reference records characterized final step counts for diagnosis:
`tnsn_alpha` records 2841 for every zone, `heat_alpha` records 654, 600, 553,
534, 532, and 540, `tnsn_torch47` records 2928, and `bdf_sn160` records 266,
274, 257, 231, 223, and 223. Step count is
diagnostic-only for every case and has
no cross-run pass/fail tolerance. Accepted-step count can vary when compiler,
library, architecture, optimization, or floating-point rounding changes which
tolerance-dependent convergence path is taken. Repeated results on one
configuration cannot justify either exact equality or a portable nonzero
bound. The parser retains the `End` step and all five source-labeled counters
separately. It requires the counter row to name the same zone, but does not
require TS to equal `End`: production writes `End` from `kstep`, while BDF TS
counts attempted timesteps and legitimately exceeds the accepted-step count.
Completion and selected endpoint values remain required.

The first two values in an `End` record are distinct: the first is the
requested target time (`tstop` in XNet), and the second is the achieved
integration time (`t`). The BDF reference records and compares target and
achieved time separately with their own printed-time policies, preserving a
genuinely distinct relationship if one occurs and avoiding compounded target
tolerance. Earlier references retain their established completion check
against the run's requested target.

The comparison also checks temperature, density, electron fraction, and an
independent composition selection for every zone. The established
silicon-burning anchors `si28`, `s32`, `ar36`, `ca40`, `ti44`, `cr48`, `fe52`,
and `ni56` are retained in each zone when they exist in the case network, even
when an anchor is below the general importance threshold. Every other species
gates comparison in a zone exactly when that zone's characterized endpoint
mass fraction satisfies `X >= 1e-4`. The inclusive boundary identifies species
that carry at least 0.01% of that zone's endpoint mass while leaving smaller
trace and intermediate abundances diagnostic-only there. Anchors are not
required: a CNO or nova network with none of them selects its own material
endpoint species by the same threshold.

Selection is deterministic. A single pass over the complete network vector
reported by XNet retains each anchor or above-threshold species in that vector's
order; anchor status changes membership, not position. The reference's
`mass_fraction_selection` object lists that exact ordered selection for every
expected zone. The loader rejects
missing or unknown zones, empty selections, invalid or duplicate species,
species absent from the complete composition vector, selected species without
tolerances, and tolerance entries unused by every zone. Case validation then
requires each committed list to equal the selection derived from that zone's
complete reference vector. Duplicate JSON object keys are rejected before
schema validation so a later zone, value, or tolerance cannot silently replace
an earlier one.

The per-zone change has the following effect. The table spells out ordering as
well as membership so reviewers can audit each list directly.

| Case | Zone | Previous case-wide selection | Per-zone selection |
| --- | --- | --- | --- |
| `tnsn_alpha` | 1-10 | `si28, s32, ar36, ca40, ti44, cr48, fe52, ni56` | unchanged |
| `tnsn_torch47` | 1 | `si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, s31, co55` | `si28, s31, s32, ar36, ca40, ti44, cr48, fe52, co55, ni56` |
| `heat_alpha` | 1 | `si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, he4, c12, o16, mg24, zn60` | `he4, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56` |
| `heat_alpha` | 2 | same case-wide list | `he4, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, zn60` |
| `heat_alpha` | 3 | same case-wide list | `he4, mg24, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, zn60` |
| `heat_alpha` | 4 | same case-wide list | `he4, o16, mg24, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, zn60` |
| `heat_alpha` | 5-6 | same case-wide list | `he4, c12, o16, mg24, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56, zn60` |

Thus the single-zone Torch47 case and all identical `tnsn_alpha` zones retain
the same pass/fail coverage. In `heat_alpha`, for example, `o16` gates zones
4-6 but remains diagnostic-only in zones 1-3. The focused tests apply the same
sum-preserving `o16`/`c12` perturbation to zones 4 and 3 to check both outcomes;
negative-fraction and composition-sum checks still cover the full vector.

Each `final_state.json` declares `xnet-comparison-v1`. Endpoint values remain
the one canonical `mac-gnu16` reference; comparison settings are properties of
that same case reference, never selected by host, platform, compiler, or
executable configuration. Scalar and selected-species settings explicitly
declare either `exact: true` with zero tolerances, or absolute and relative
limits. Values and limits can each be case-wide or a complete per-zone map.
The numerical-field criterion is

```text
base allowance = atol + rtol * abs(reference)
```

Exact policies still require exact equality and receive no numerical margin.
For non-exact scalar fields, the comparator adds only the binary ULPs involved
in representing the actual value, reference value, and calculated allowance.
This prevents a decimal value printed exactly at a declared boundary from
failing because the subtraction is performed in binary floating point. It is
not a stored decimal tolerance coefficient and has no practical effect away
from a boundary.

```text
effective allowance = base allowance
    + ulp(actual) + ulp(reference) + ulp(base allowance), when base allowance > 0
pass exactly when abs(actual - reference) <= effective allowance
```

An exact policy instead uses an effective allowance of exactly zero. Complete-
vector L1 and L-infinity gates compare directly to their stored limits; they do
not receive the scalar ULP guard. No selected-species, vector, printed-sum, or
normalization gate reuses another gate's slack. Every structural, association,
finite, nonnegative, ordering, process, timeout, required-output, provenance,
malformed-diagnostic, and independent-gate check is unchanged.

Each reference contains the final mass fraction for every case-network species
and every zone. The separate `mass_fraction_selection` and
`mass_fraction_tolerances` mappings identify field-aware species checks without
repeating endpoint values. A compact `all_selected` setting can cover every
selected species, while future cases may use per-species settings when evidence
requires them. `composition_norm_limits` optionally adds complete
vector `L1` and `L-infinity` gates; `L2` remains diagnostic-only. Selected
species and vector gates are independent, so either can reject a result. An
`L-infinity` failure reports the responsible species.

`composition_error_norms.json` records complete-vector `L1`, `L2`, and
`L-infinity` norms, the `L-infinity` species, and any reference-owned vector
limits. It is test output, not a reference-update mechanism.

A three-configuration study (`mac-gnu16`, `mac-llvm`, and `etacar-gnu16`)
provides characterization evidence for the five current cases under the
documented serial optimized configuration. This is characterization, not
scientific validation, and it does not establish portability outside that
matrix. Future cases should keep one canonical endpoint, choose the coarsest
quantitatively justified settings, preserve exact invariants, and document the
supporting evidence rather than adding study archives or platform branches to
this tree.

A later hosted GNU 13.3.0 Ubuntu 24.04 observation extended that empirical
envelope after exposing optimized-run variation. For each affected case, zone,
and comparison category, the stored coefficient is derived from

```text
1.5 * max(abs(observation - canonical mac-gnu16 value))
    + half the final printed decimal unit
```

The maximum is taken over the union of the three-configuration observations
and the preserved `NUM-OBS-001` observations. The resulting value is rounded upward to
two significant decimal digits. The hosted observation controls each revised
coefficient below; the earlier maxima were smaller. The structural
normalization row retains the canonical sum's offset from one because that gate
compares directly with one, then adds the same cross-canonical envelope.

| Policy category | Controlling maximum deviation | Half printed unit | Unrounded result | Stored coefficient |
| --- | ---: | ---: | ---: | ---: |
| `tnsn_torch47` z1 selected species (`co55`) | `2.800e-9` | `5e-12` | `4.2050e-9` | `4.3e-9` |
| `tnsn_torch47` z1 L1 | `6.786e-9` | `5e-13` | `1.01795e-8` | `1.1e-8` |
| `tnsn_torch47` z1 L-infinity (`co55`) | `2.800e-9` | `5e-13` | `4.2005e-9` | `4.3e-9` |
| `tnsn_torch47` z1 printed sum | `3.131e-9` | `5e-13` | `4.6970e-9` | `4.7e-9` |
| `bdf_sn160` z1 printed sum | `4.866186044e-7` | `5e-14` | `7.299279566e-7` | `7.3e-7` |
| `bdf_sn160` z3 printed sum | `4.696526643e-7` | `5e-13` | `7.0447949645e-7` | `7.1e-7` |
| `bdf_sn160` z3 normalization to one | canonical offset `3.01014263357e-5` plus `4.696526643e-7` | `5e-13` | `3.080590583215e-5` | `3.1e-5` |

The `tnsn_torch47` selected setting remains one zone-wide `all_selected`
coefficient because that is the existing smallest supported category. All
unaffected BDF zones retain their prior printed-sum and normalization values.
Canonical endpoint values, exact policies, and unrelated coefficients are
unchanged. Zero-valued non-exact references receive only their absolute
allowance; relative scale contributes zero. Trace values remain structurally
required but gate as selected species only under the existing per-zone policy.
Near-one printed sums and normalization use their independently stored
coefficients. This remains empirical regression characterization for the
recorded configurations, not scientific validation or a general
compiler/platform support claim.

Timer exclusion is structural rather than line-count based: only a
`Timers Summary:` heading and immediately following timer-name/numeric-value
rows are removed. The first non-timer row ends the exclusion, so unrelated
diagnostic content is not hidden.

## SN160 Backward Euler characterization evidence

No tracked `test/Test_Problems/Results/net_diag_heat_sn160` exists, and a
repository-history search found no historical endpoint with usable compiler,
platform, input, or scientific provenance. The committed result is therefore
a new characterization, not historical truth. `net_diag01` retains enough
printed precision for the selected endpoint policy, while each `ts_*` file
remains a required fresh output and is neither decoded nor compared.

The reference was generated on 2026-08-05 from production and input revision
`01dd4963e9b9677f64711c90e08f50d468bc99a4` on macOS 26.6 arm64 with GNU
Fortran 16.1.0. The following commands and in-source executable path record
that revision's historical build; they are not current build instructions:

```bash
make -C source clean
make -C source -j
```

Resolved selections were `EXE=xnet`, `CMODE=OPT`, `PE_ENV=GNU`,
`MPI_MODE=OFF`, `OPENMP_MODE=OFF`, `GPU_MODE=OFF`, `EOS=STARKILLER`,
`MATRIX_SOLVER=dense`, and `LAPACK_VER=NETLIB`. The known nominal-serial Make
dependency also compiled `xnet_parallel.F90` with `mpifort` and emitted its
existing argument and rank mismatch warnings; the linked executable used
`xnet_parallel_stubs.o`. MPI, OpenMP, accelerators, other matrix solvers,
other libraries, other compilers, and other platforms were not validated.

Three isolated optimized runs returned direct process status 0. The historical
focused command below was run three times with distinct external pytest base
directories:

```bash
.venv/bin/python -m pytest -q test/regression/test_regression.py::test_heat_sn160 --xnet-executable="$PWD/source/xnet" --basetemp=<external-directory> --durations=1
```

The three retained `net_diag01` files were parsed independently with
`parse_diagnostic`, and equality of the resulting `FinalState` sequences
confirmed identical endpoints, all 160 mass fractions, `End` steps, and all
five solver counters. Pytest call times were 2.04, 2.02, and 2.02 seconds,
well inside the unchanged 30-second timeout. Step is diagnostic-only and is
not compared. Requested and achieved times were identical at printed
precision in every zone:

| Zone | `End` step | Time (s) | Final T (GK) | Density (g/cm3) | Ye | Printed sum(X) | Sum bound | Compared species |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 668 | 10 | 4.5362097 | 1.0000000e7 | 0.49886745 | 0.999999989648 | 2.7904e-8 | 26 |
| 2 | 602 | 1 | 5.4842926 | 3.1622777e7 | 0.49912204 | 1.000000021899 | 4.2621e-8 | 38 |
| 3 | 548 | 0.1 | 6.3425061 | 1.0000000e8 | 0.49928188 | 1.000000113104 | 1.3592e-7 | 52 |
| 4 | 531 | 0.001 | 7.1755363 | 3.1622777e8 | 0.49953200 | 1.000000033670 | 5.2818e-8 | 65 |
| 5 | 548 | 0.0001 | 8.1086726 | 1.0000000e9 | 0.49953803 | 1.000000031608 | 5.0945e-8 | 76 |
| 6 | 564 | 0.00001 | 9.2354937 | 3.1622777e9 | 0.49953848 | 1.000000007620 | 2.7537e-8 | 84 |

The parser preserves the source-labeled counter values separately from the
`End` step. `xnet_output.F90` writes `End` from `kstep` and writes the five
counter columns from `ktot(1:5)`. For Backward Euler, `TS` accumulates trial
timestep attempts, `NR` accumulates Newton-Raphson iterations, and `Jacobian`,
`Deriv`, and `CrossSect` count the corresponding builds or evaluations. These
are solver diagnostics rather than pass/fail values. All three optimized runs
produced the same records:

| Zone | `End` | TS | NR | Jacobian | Deriv | CrossSect |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 668 | 668 | 1238 | 1238 | 1239 | 1239 |
| 2 | 602 | 602 | 1103 | 1103 | 1104 | 1104 |
| 3 | 548 | 548 | 1050 | 1050 | 1051 | 1051 |
| 4 | 531 | 531 | 1013 | 1013 | 1014 | 1014 |
| 5 | 548 | 548 | 1048 | 1048 | 1049 | 1049 |
| 6 | 564 | 564 | 1084 | 1084 | 1085 | 1085 |

Every zone stores all 160 values, including every trace abundance printed by
this baseline. This particular reference contains no printed zeros and does
not exercise values near the Backward Euler `1e-30` molar-abundance cutoff;
the smallest printed mass fraction is zone 1 `ca48 = 1.9843634e-25`.
Pass/fail species are the established silicon-burning anchors when present
plus every species whose zone-specific baseline has `X >= 1e-4`; the reference
records the exact ordered lists summarized by the final column. Each selected
value and scalar field uses half its baseline's last printed place as `atol`
and `5e-8` as `rtol`, matching the diagnostic's eight-significant-digit
representation. Complete vectors still receive exact identity/order,
uniqueness, finite, nonnegative, normalization, and diagnostic norm checks.
Cross-integrator agreement is not evaluated by this regression.

The required outputs totaled 30,663,926 bytes: 47,898 bytes for `net_diag01`,
694,600 bytes for six ASCII histories, and 29,921,428 bytes for six binary
histories. Isolated preprocessing created `ab_blank`, `match_data`,
`match_read`, `matr_shape`, `net_desc`, `net_diag`, `nets3`, `nets4`,
`nuc_data`, and `sparse_ind`, totaling 847,933 bytes. The committed complete
JSON reference is 52,447 bytes. Hashes for the control, five network sources,
six trajectories, and EOS table are recorded in the reference; the tracked
inputs remained unchanged after the runs. The post-run ignored-file check
was scoped to every repository input and case directory that the execution
could mutate:

```bash
git status --short --ignored -- test/Data_SN160 test/Test_Problems tools/starkiller-helmholtz test/regression/cases/heat_sn160
# no entries
```

Thus neither tracked nor ignored runtime or preprocessing outputs appeared in
those source directories. At that revision, expected build products remained
confined to the ignored `source/` directory, and pytest outputs remained under
the explicitly external `/private/tmp` base directories.

An alternate clean GNU `CMODE=DEBUG` build was also exercised with these
historical commands before restoring that revision's tracked-default optimized
build:

```bash
make -C source clean
make -C source -j CMODE=DEBUG
.venv/bin/python -m pytest -q test/regression/test_regression.py::test_heat_sn160 --xnet-executable="$PWD/source/xnet" --xnet-timeout=120 --durations=1
make -C source clean
make -C source -j
```

The debug case had first exceeded the normal 30-second timeout, then completed
directly with status 0 in the explicit 120-second pytest run, whose case call
took 37.19 seconds. It retained the optimized run's step counts but failed the
narrow endpoint characterization; the largest complete-composition difference
was `6.33e-7` for zone 3 `ni56`. Its zone 3 L1, L2, and L-infinity differences
were `3.733e-6`, `1.058e-6`, and `6.33e-7`. These observations do not establish
a scientifically acceptable cross-mode tolerance, so the reference was not
widened merely to make the debug build pass. Portability beyond the optimized
configuration is not established.

As a controlled end-to-end effectiveness check, a temporary reference changed
zone 3 `ni56` from `0.065404983` to `0.07`. The real pytest case returned
status 1 and reported a `4.595e-3` difference against an allowed `4.000e-9`.
The reference value was then restored. Normal execution has no reference
creation or update path. After restoration and review fixes, the focused
helpers passed 70 tests and the complete optimized suite passed 74 tests in
4.94 seconds.

## SN160 Backward Differentiation Formula characterization evidence

No tracked `test/Test_Problems/Results/net_diag_bdf_sn160` exists. A
path-scoped search of all repository history found no prior BDF SN160 endpoint
reference with usable revision, compiler, platform, inputs, or scientific
provenance. This reference is therefore a new characterization of current
maintained BDF behavior, not a historical truth reference and not independent
scientific validation.

The reference was generated on 2026-08-05 from production and input revision
`a8b64764a6d614f406da6c897e6b051fb3e1972d` on macOS 26.6 arm64 with GNU
Fortran 16.1.0, Python 3.13.0, and pytest 9.1.1. The complete control is bound
by its SHA-256 hash in the reference. The following commands and executable
path record that revision's historical build; they are not current build
instructions:

```bash
make -C source clean
make -C source -j
```

Resolved selections were `EXE=xnet`, `CMODE=OPT`, `PE_ENV=GNU`,
`MPI_MODE=OFF`, `OPENMP_MODE=OFF`, `GPU_MODE=OFF`, `EOS=STARKILLER`,
`MATRIX_SOLVER=dense`, and `LAPACK_VER=NETLIB`. The exact executable was
`source/xnet`. The known nominal-serial dependency also compiled
`xnet_parallel.F90` with GNU `mpifort` 16.1.0 and emitted the existing MPI
argument type and rank mismatch warnings; the executable linked
`xnet_parallel_stubs.o`. MPI execution, OpenMP, accelerators, other matrix
solvers and libraries, other compilers, other platforms, and DEBUG mode were
not checked for this BDF characterization.

The legacy provenance and normalized comparison are executable checks, not
only prose. `test_bdf_control_is_the_normalized_legacy_id_54_concatenation`
forms the same `test_settings_bdf + setup_bdf_sn160` concatenation as
`do_test_bdf`, strips trailing whitespace, applies the two isolated path
normalizations, changes only the block-size line from 4 to 1, and requires an
exact match with the committed control. The reference records SHA-256 hashes
for that control, the five network sources, abundance input, six trajectories,
and EOS table. It also records integration choice 3, the maintained solver
identity, and the effective `changemx = changemxt = 1e10` state imposed by
XNet.

Three isolated optimized pytest runs used the unchanged default 30-second
timeout and distinct external base directories. This historical command was
run once for each directory:

```bash
.venv/bin/python -m pytest -q test/regression/test_regression.py::test_bdf_sn160 --xnet-executable="$PWD/source/xnet" --basetemp=<external-directory> --durations=1
```

Each XNet subprocess recorded direct return status 0. Pytest call times were
0.76, 0.75, and 0.76 seconds, well inside 30 seconds. Independent parsing of
the three retained `net_diag01` files produced identical targets, achieved
times, scalar endpoints, all 160 mass fractions, `End` steps, and all five
solver counters. Requested and achieved times were identical at printed
precision in every zone:

| Zone | `End` | Time (s) | Final T (GK) | Density (g/cm3) | Ye | Printed sum(X) | Sum bound | Compared species |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 266 | 10 | 4.5376972 | 1.0000000e7 | 0.49885684 | 1.000000279473 | 2.9702e-7 | 26 |
| 2 | 274 | 1 | 5.4873128 | 3.1622777e7 | 0.49912074 | 1.000000414592 | 4.3527e-7 | 38 |
| 3 | 257 | 0.1 | 6.3448519 | 1.0000000e8 | 0.49928049 | 1.000030101426 | 3.0124e-5 | 52 |
| 4 | 231 | 0.001 | 7.1773063 | 3.1622777e8 | 0.49952960 | 1.000000560949 | 5.8010e-7 | 65 |
| 5 | 223 | 0.0001 | 8.1101200 | 1.0000000e9 | 0.49953348 | 1.000001564188 | 1.5835e-6 | 76 |
| 6 | 223 | 0.00001 | 9.2368042 | 3.1622777e9 | 0.49953675 | 1.000001415571 | 1.4355e-6 | 84 |

Every zone stores all 160 values in exact `sunet` order. The pass/fail
selection is made independently per zone: all available silicon-burning
anchors plus each species with characterized `X >= 1e-4`, retained in the
complete-vector order. Each scalar and selected composition value uses half
its baseline's last printed place as `atol` and `5e-8` as `rtol`. The printed
sum bound is the absolute baseline printed-sum offset shown above plus one summed
half-last-place bound for the complete 160-value vector. In particular, the
zone 3 printed-sum offset is recorded directly. The maintainer assessed the
reported value as acceptable. This endpoint check is separate from the BDF
`iconvc == 3` weighted RMS convergence norm.

The reference records `End` and all solver counters for diagnosis. No counter
has a pass/fail threshold:

| Zone | `End` | TS | NR | Jacobian | Deriv | CrossSect |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 266 | 280 | 416 | 280 | 417 | 417 |
| 2 | 274 | 289 | 423 | 289 | 424 | 424 |
| 3 | 257 | 271 | 416 | 271 | 417 | 417 |
| 4 | 231 | 245 | 359 | 245 | 360 | 360 |
| 5 | 223 | 237 | 357 | 237 | 358 | 358 |
| 6 | 223 | 237 | 351 | 237 | 352 | 352 |

All six TS values exceed `End`, directly exercising the shared parser's
solver-independent representation. The parser still requires the exact
`Counters: Zone TS NR Jacobian Deriv CrossSect` heading, six nonnegative
integer fields, correct zone association, completeness, and order.

The required BDF outputs totaled 13,154,474 bytes: 47,890 bytes for
`net_diag01`, 297,200 bytes for the six ASCII histories, and 12,809,384 bytes
for the six binary histories. Isolated preprocessing created `ab_blank`,
`match_data`, `match_read`, `matr_shape`, `net_desc`, `net_diag`, `nets3`,
`nets4`, `nuc_data`, and `sparse_ind`, totaling 847,933 bytes. The committed
JSON reference is 54,190 bytes. Binary histories are required fresh and
nonempty but are not decoded, compared, or committed.

A controlled end-to-end failure temporarily changed selected zone 3 `ni56`
from `0.064921366` to `0.07`. Both the focused BDF command and the complete
pytest command returned status 1. They reported the zone, species, actual
value `6.492136600e-02`, reference value `7.000000000e-02`, absolute
difference `5.079e-03`, and allowance `4.000e-09`; the complete run otherwise
passed all other tests. The reference was restored, the focused case passed
again, and normal test execution still has no reference creation or update
path.

### Backward Euler/BDF diagnostic comparison

The following same-platform comparison uses the accepted Backward Euler
characterization and the BDF characterization above. Target and achieved times
and densities agree at printed precision. The other differences are diagnostic
only; no shared endpoint, step, counter, or history criterion has been
established.

| Zone | T BE / BDF (GK) | Ye BE / BDF | `End` BE / BDF | TS BE / BDF | NR BE / BDF | Jacobian BE / BDF | Deriv BE / BDF | CrossSect BE / BDF |
| ---: | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.5362097 / 4.5376972 | 0.49886745 / 0.49885684 | 668 / 266 | 668 / 280 | 1238 / 416 | 1238 / 280 | 1239 / 417 | 1239 / 417 |
| 2 | 5.4842926 / 5.4873128 | 0.49912204 / 0.49912074 | 602 / 274 | 602 / 289 | 1103 / 423 | 1103 / 289 | 1104 / 424 | 1104 / 424 |
| 3 | 6.3425061 / 6.3448519 | 0.49928188 / 0.49928049 | 548 / 257 | 548 / 271 | 1050 / 416 | 1050 / 271 | 1051 / 417 | 1051 / 417 |
| 4 | 7.1755363 / 7.1773063 | 0.49953200 / 0.49952960 | 531 / 231 | 531 / 245 | 1013 / 359 | 1013 / 245 | 1014 / 360 | 1014 / 360 |
| 5 | 8.1086726 / 8.1101200 | 0.49953803 / 0.49953348 | 548 / 223 | 548 / 237 | 1048 / 357 | 1048 / 237 | 1049 / 358 | 1049 / 358 |
| 6 | 9.2354937 / 9.2368042 | 0.49953848 / 0.49953675 | 564 / 223 | 564 / 237 | 1084 / 351 | 1084 / 237 | 1085 / 352 | 1085 / 352 |

Representative selected abundant species and complete-vector differences are:

| Zone | he4 BE / BDF | si28 BE / BDF | ni56 BE / BDF | L1 | L2 | Linf (species) |
| ---: | --- | --- | --- | ---: | ---: | --- |
| 1 | 6.5052911e-3 / 6.6237558e-3 | 2.3790347e-3 / 2.1938516e-3 | 5.8064083e-1 / 5.8052829e-1 | 3.2274e-3 | 1.0058e-3 | 5.6476e-4 (`ni58`) |
| 2 | 4.6308135e-2 / 4.6664867e-2 | 5.1601149e-3 / 5.1783213e-3 | 2.0222165e-1 / 2.0068812e-1 | 3.9905e-3 | 1.7702e-3 | 1.5335e-3 (`ni56`) |
| 3 | 1.0357065e-1 / 1.0402255e-1 | 8.2611060e-3 / 8.2685346e-3 | 6.5404983e-2 / 6.4921366e-2 | 2.6063e-3 | 8.2952e-4 | 4.8362e-4 (`ni56`) |
| 4 | 1.4310443e-1 / 1.4347333e-1 | 9.3280462e-3 / 9.3293077e-3 | 3.1695134e-2 / 3.1519896e-2 | 1.8395e-3 | 5.4246e-4 | 3.6890e-4 (`he4`) |
| 5 | 1.6456740e-1 / 1.6484535e-1 | 9.2537712e-3 / 9.2534273e-3 | 1.8892735e-2 / 1.8814958e-2 | 1.3423e-3 | 3.8126e-4 | 2.7795e-4 (`he4`) |
| 6 | 1.7204490e-1 / 1.7226093e-1 | 8.2250493e-3 / 8.2247830e-3 | 1.3048592e-2 / 1.3005155e-2 | 1.0613e-3 | 2.9026e-4 | 2.1603e-4 (`he4`) |

| Quantity | Backward Euler | BDF |
| --- | ---: | ---: |
| Repeated pytest call time | 2.04, 2.02, 2.02 s | 0.76, 0.75, 0.76 s |
| Required output | 30,663,926 bytes | 13,154,474 bytes |
| Generated preprocessing | 847,933 bytes | 847,933 bytes |
| Committed reference | 52,447 bytes | 54,190 bytes |

These runtime and size observations describe one optimized serial CPU
configuration and are not performance benchmarks. Differing timestep,
convergence, and abundance-floor policies make the different endpoints and
execution histories expected characterization evidence, not evidence that
either integrator is scientifically superior. Scientific interpretation of
the endpoint differences remains a human-maintainer decision.

## Torch47 characterization evidence

The Torch47 reference was generated on 2026-08-05 from revision
`5e7e1543d432f3c2792e40e271816ecaf8184fad` on macOS 26.6 arm64 with GNU
Fortran 16.1.0. The following commands and in-source executable path record
that revision's historical build; they are not current build instructions:

```bash
make -C source clean
make -C source -j
```

Resolved selections were `CMODE=OPT`, `PE_ENV=GNU`, `MPI_MODE=OFF`,
`OPENMP_MODE=OFF`, `GPU_MODE=OFF`, `EOS=STARKILLER`,
`MATRIX_SOLVER=dense`, and `LAPACK_VER=NETLIB`. The known Make dependency
limitation also compiled `xnet_parallel.F90` with `mpifort` and emitted its
existing argument-mismatch warnings; the linked executable used
`xnet_parallel_stubs.o`.

An isolated direct Torch47 process returned status 0 in 0.37 seconds. Three
subsequent pytest runs each completed in 0.40 seconds including test overhead
and produced identical parsed 47-species endpoints. This repeated result only
characterizes one compiler, build, and machine. The required outputs totaled
5,877,984 bytes: 5,176 for `net_diag01`, 586,000 for the ASCII history, and
5,286,808 for the binary history. The committed JSON reference is 4,098 bytes.
For comparison on the same configuration, the pytest case times were 0.52
seconds for `tnsn_alpha` and 0.37 seconds for `heat_alpha`; their required
outputs totaled 14,470,756 and 1,757,702 bytes respectively. Torch47 remains
well inside the unchanged 30-second timeout and is suitable for the fast local
suite on this configuration.

At that revision, the focused helper command passed 39 tests and the complete
suite passed 42:

```bash
python -m pytest -q test/regression/test_xnet_regression.py
python -m pytest -q test/regression \
    --xnet-executable="$PWD/source/xnet"
```

A temporary copy of the Torch47 reference changed the network-specific
selected product `co55` from `0.00076324081` to `0.001`. The real end-to-end
runner then failed with pytest status 1 and reported an absolute difference of
`2.368e-4` against an allowed `5.500e-11`. Normal execution has no
reference-writing path.

The third explicit Python registration remains a short `RegressionCase`
declaration and shares the existing loader-free validation path. A TOML
manifest would duplicate these values while adding schema and loading code, so
all three cases remain registered directly in Python.

## NSE runtime-behavior coverage

The suite includes exactly one complete ordinary NSE-initialized evolution. The
`nse_sn160` case uses the tracked `th_frohlich2006_nse` trajectory, whose
initial temperature exceeds the 8 GK threshold and whose subsequent evolution
cools below it. The supplied `ab_ye49` file differs from the trajectory's
initial `Ye = 0.55`. The runner therefore requires the production diagnostic
marker `Initial abundances from NSE` in addition to the normal complete-state
comparison; silently selecting the file state cannot satisfy the case.

The standalone `xnse` test supplies three ordered density/temperature/Ye rows
and checks their association with complete state and counter output. Two
failure probes require nonzero status for malformed input and for a
deliberately out-of-domain electron fraction that drives NSE nonconvergence.
These checks cover program behavior and provide characterization evidence; they
do not scientifically validate NSE physics or screening models.

## Current limits and next cases

These cases establish runtime and software-behavior checks plus narrow
numerical characterization. They do not establish broad scientific validity,
portability, performance benchmarking, CI suitability, or support for MPI,
threading, accelerators, log-ft rates, or networks larger than SN160. NSE
coverage is limited to direct software invariants, one ordinary unscreened
trajectory characterization, and deterministic process behavior; it does not
validate screened physics. BDF coverage is limited to the characterized
`bdf_sn160` endpoint on the recorded optimized configuration.

The self-heating comparisons cover final `net_diag01` endpoints only. They do
not inspect the evolution history in `ev_*` or binary `ts_*` output. Binary
time-series data would require a separate decoding and comparison study before
they could affect regression pass/fail. The Torch47 and SN160 cases change
endpoint coverage only; their larger `ev_*` and `ts_*` outputs remain required
for freshness but are not parsed or compared.

Per-zone importance selection changes only regression classification. It does
not change XNet calculations, establish scientific validity, or extend the
single-compiler and single-platform characterization evidence described above.
