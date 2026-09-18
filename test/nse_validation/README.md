# Independent unscreened NSE validation

This directory compares XNet's existing, unscreened ideal-gas NSE calculation with
three compositions computed independently from published equilibrium
equations, using the same finite species set and explicitly reconciled nuclear
inputs.  It is not a general NSE package and it does not change XNet physics.

## Verification record (before reference generation)

The statements below were checked at historical `jaharris87/XNet` fork commit
`65271bcbeea430534c1adc92748ef13bea10c228`.  This record was written before
freezing expected compositions or running XNet against them; that commit is
not expected to exist in a normal `starkiller-astro/XNet` clone.

### XNet software facts

- The component test directly calls the public `nse_initialize` and
  `nse_solve` procedures from the production `source/xnet_nse.F90`, then reads
  public `xnse`, `ynse`, `unse`, and `knrtot`.  The test build copies that
  production source into an isolated test-drive build and supplies only the
  surrounding XNet module interfaces.  The scientific comparison uses the
  same production entry points.
- The direct component fixture has eight deliberately synthetic species with
  unit spins and partition functions and rounded binding energies.  It is a
  software-test fixture, not an independent scientific reference, and is not
  used to define the retained expected values.
- XNet represents the NSE result as mass fractions `X_i`; it also sets
  `Y_i = X_i / (m_i N_A)`.  With the active translational-mass convention
  `m_i = A_i/N_A`, this is `Y_i = X_i/A_i`.
- `rho` is in g cm^-3, `T9` is temperature in GK, and
  `Ye = sum_i (Z_i/A_i) X_i`.  Free neutrons and protons are ordinary entries
  in the finite network and their stored mass excesses define every binding
  energy through `B_i = N_i MEn + Z_i MEp - ME_i`.
- The active NSE translational mass is `A_i/N_A`; commented alternatives based
  on actual nuclear masses are not active.  The selected `netwinv` file
  supplies ground-state spin, mass excess, and 24 normalized partition-factor
  values.  XNet multiplies the normalized factor by `2J_i+1` and interpolates
  linearly in `ln(G)` versus `T9`.
- The constants that enter the NSE abundance equation are the binary64 values
  compiled from `pi`, `hbar` (MeV s), `avn` (mol^-1), `bok` (MeV GK^-1), and
  `epmev` (erg MeV^-1) in `source/xnet_constants.F90`.  Tracked GNU, Intel,
  NVHPC, Cray, and AMD compiler configurations use 64-bit default real values;
  therefore unsuffixed decimal constants and the encoded partition-temperature
  grid are binary64 in supported builds.
- For `iscrn = 0`, XNet sets every Coulomb correction `h_i` to zero.  The final
  abundance equation then contains no EOS call or other nonideal correction.
- `nse_solve` exposes the two chemical-potential-like roots, iteration and
  evaluation counters, and the final mass fractions.  Its local `info` status
  and residual vector are written only through optional diagnostics and are
  not public module state.  The residuals are `sum(X)-1` and
  `sum((Z/A-Ye)X)`.  The configured function tolerance is `1e-8`; a separate
  `1e-14` step-size exit can return `info = 2`, so the scientific test must
  calculate and gate the physical residuals itself instead of treating either
  positive status as sufficient.
- Exponentials are protected by `safe_exp`, and each resulting mass fraction
  is then clipped to `[0,1]`.  The validation therefore compares the complete
  vector, including stored tiny values, and preflights the reference states so
  no scientifically material species relies on either clipping boundary.
- The current test-drive dependency is version 0.5.0.  Its Makefile builds
  separate test executables, copies production Fortran into the isolated build
  to select the stub module files, and treats `build/` plus generated module,
  object, and executable files as untracked generated files.

### Network candidates recorded before the freeze

No network was selected merely because the planning report proposed SN160.
All tracked `test/Data_*` directories with both `sunet` and `netwinv` were
surveyed.  The networks considered are:

| Network | Species | Useful property | Material concern |
| --- | ---: | --- | --- |
| `Data_SN160` | 160 | Smallest tracked dense iron-group network with n and p | Must show negligible boundary sensitivity against a larger network |
| `Data_SN231` | 231 | Broader isotopic coverage with the same input-data format | A larger committed vector and slower test |
| `Data_torch47` | 47 | Small and includes n and p | Sparse off the alpha chain; may truncate neutron- or proton-rich equilibrium |
| `Data_Nova` | 169 | Broad `Z/A` range | Light-network emphasis and maximum A = 54 |
| `test/build_net/sunet.torch489` | 489 | Strictly contains SN160 and SN231; broad light and iron-group coverage | Requires a reproducible `netwinv` build |
| `Data_Reaclib20180621` | 7852 | Full tracked comparison set | Too large for the smallest routine validation fixture |

`Data_alpha` cannot represent `Ye != 0.5` and omits free nucleons;
`Data_CNO` is too light; and the 304--7852-species networks are not the
smallest credible first choice.  The choice was made from high-precision,
XNet-free preflight calculations.  In particular, SN160 was acceptable only
if SN231 showed that mass in additional species and changes to shared species
were immaterial at all three states.  The disposition appears below.

### Candidate state design recorded before the freeze

Exactly three states were to be retained.  Before selection, the independent
calculator examined the report's `(T9,rho,Ye)` values
`(7,1e9,0.50)`, `(7,1e9,0.45)`, and `(9,1e7,0.50)` together with bounded
stress candidates at lower NSE-domain temperature/density, proton-rich
`Ye`, and more neutron-rich `Ye`.  Selection criteria are:

1. the state remains in a regime where NSE is a physically meaningful target,
2. the chosen finite network represents the result without material boundary
   mass,
3. analytic- and finite-difference-Jacobian variants and multiple starts converge,
4. no material abundance is controlled by XNet's exponential/clipping guards,
5. the three states are scientifically nonredundant and numerically durable.

A `Ye` below the minimum bound-nucleus `Z/A` is not intrinsically
unrepresentable because free neutrons have `Z/A = 0`; it is nevertheless a
network-boundary stress test and will only be retained if expanded-network
preflight supports it.  Similarly, a proton-rich state may be more demanding
but will not be chosen solely because it makes the solver fail.

### Scientific authority and implementation independence

The scientific authority is the Maxwell-Boltzmann chemical-equilibrium and
mass/charge formulation in Seitenzahl et al. (2009), equations (2)--(10),
with the constraint notation corroborated by Hix & Meyer (2006).  Lippuner &
Roberts (2017), Appendix B, supplies an independent published NSE derivation
and numerical-method precedent.  The generator is an independent numerical
realization of those sources, not itself scientific authority.  Its equations
are implemented in `reference_solver.py` with Python's `Decimal` arithmetic.
That code imports no XNet routine, does not execute XNet, and is organized
around log-sum-exp constrained equilibrium rather than the implementation
structure of `source/xnet_nse.F90`.

The supplied paper files used in this verification have these SHA-256 hashes:

| Paper | SHA-256 |
| --- | --- |
| Seitenzahl et al. (2009) | `5d1c684abf36463dc88e0f3e6aedf0234f95af2b0e85acfe01621b4e7beb86f4` |
| Hix & Meyer (2006) | `59da109818ec8e053af07f9eae7a4d0bfa326c0a969613e5e4889cf80bd9bd95` |
| Lippuner & Roberts (2017) | `246c62865d292df60046572084c9f543d7d575efbc4af65099ec2d50a926140e` |
| Reichert et al. (2023) | `8e20bf9bf43e38c26dfe7416129d6f22f0ef576061866e68663960dcff01510d` |

External-code roles are intentionally qualified:

The source survey was rechecked on 2026-08-10 at pynucastro commit
`a7268f86f42c556172578ca53293cf00b6c539ab`, Microphysics commit
`6fb41b5f7475b42a06eb5b09ff0520c9f08aa7f0`, WinNet commit
`de0c9c852907a8714444287f3ebe0bdedf328397`, `wneq` commit
`d0bede932936d5b30a0d73e7a435f1c3d3f33ca4`, and `wnnet` commit
`4aef8352373f144bb8c58516a7067fc656991315`.  SkyNet is pinned through its
published v1.0 archive, DOI `10.5281/zenodo.1008754`.  These versions record
what was evaluated; only the sources identified below as scientific authority
or numerical corroboration support the retained result.

- Frank Timmes's public `nse` code is useful historical and methodological
  evidence, but XNet documents historical inspiration from that implementation;
  it is not a fully independent numerical authority for this test.  The
  official page describes a 47-isotope illustrative solver; its linked archive
  returned HTTP 403 during this verification, so no unaudited copy was used.
- No XNet--pynucastro developer overlap was found.  Austin Harris made one
  2017 Microphysics commit (`9c08513b`) that declared an EOS name public in
  seven EOS files; it did not touch NSE, nuclear data, or a solver.  No
  XNet--Microphysics NSE source lineage was found.  Pynucastro and Microphysics
  do share developers, network-generation machinery, and nuclear-data paths
  with each other, so agreement from both would be one family of secondary
  corroboration, not two independent results.
- SkyNet and WinNet provide relevant published verification methods.  Their
  default nuclear data, species sets, mass conventions, partition functions,
  and nonideal corrections are not automatically compatible with this finite
  XNet problem.  They may corroborate a reconciled state, but they do not
  replace the primary, explicit-input calculation.

After the expected data and tolerances were frozen, pynucastro commit
`a7268f86f42c556172578ca53293cf00b6c539ab` supplied a secondary numerical
check.  Its default problem was not compared directly: `ar45` lacks a default
spin, and its mass, translational-mass, partition, and constant choices differ.
The optional `corroborate_pynucastro.py` adapter instead supplies the exact
retained inputs to pynucastro's independently maintained NSE equation and
SciPy solver.  It obtains complete-vector L1 differences of `5.29e-14`,
`8.74e-14`, and `3.63e-13` for the three states.  This result corroborates the
primary calculation but is neither the expected-data generator nor a gating
test.  Microphysics was not counted as another numerical corroboration because
it belongs to the same developer and nuclear-data family as pynucastro.

### Input-data reconciliation

`extract_inputs.py` is data tooling, not an NSE solver.  It checks exact
`sunet`/`netwinv` identity and order, records both source decimals and the
binary64 values consumed by tracked builds, and recomputes the active XNet
binding and translational-mass inputs.  Every selected spin and normalized
partition factor is checked against
`test/build_net/partf_data/winvne_JINAv22`.  Every mass excess is checked
against the mass source selected by that file after the network builder's
eight-decimal MeV storage rounding.  This reconciles the exact retained input
bytes; it does not claim that those nuclear inputs are exact measurements or
that every mass table belongs to REACLIB v2.2.

The retained scientific-input identity covers the network bytes and species
order, effective binary64 constants and temperature grid, effective
per-species nuclear inputs, and stated conventions.  It deliberately excludes
equivalent source-literal spellings, Git commits, paths, raw production-source
hashes, reconciliation annotations, the top-level raw-data provenance record,
and generator provenance.  The older `canonical_input_sha256` and those source
hashes remain frozen historical metadata describing how the reference was
generated; they are not compared with the current implementation during an
ordinary test.  Production changes remain covered by compiling and running the
current `xnet_nse` calculation against the independently retained compositions
and unchanged numerical gates.

The public `jaharris87/build_net` archive fixes both central raw inputs at its
initial database commit `77141ca2a3dfc9fa9fd52ef0fcf39a49d74c08e1`
(2017-01-24).  At that commit, `mass_reac1.dat` has Git blob
`acd42b416e52edea020990d68631fb2b8063275d` and retained SHA-256
`0d068a92c6694485e117df2da8b081fe1a7cd7b80757303d1f72d4bfc16fff50`;
`winvne_JINAv22` has Git blob
`5bac4be2a2095bfd6a6fe3b1282bea4bcd6040ea` and retained SHA-256
`5fc838645cac2f5e2eb9ddfd5a60bb41df8a6378b8f28748b1b47ed6148555aa`.
XNet imported those exact blobs at subtree commit
`90e9363d5f9443a8ad2d5e986232c1f60bb5b96a`.  JINA's public snapshot archive
dates REACLIB v2.2 to 2016-11-14.  The `mass_reac1.dat` header identifies the
JINA Nuclide Database evaluation label `reac1`, but the original JINA
per-record publication or snapshot identifier is not recoverable from the
retained metadata.  `reference.json` records that limitation instead of
inventing a stronger source attribution.

The proton record uses the neutral-hydrogen atomic mass excess.  Because the
same atomic-mass convention appears with the same `Z` in both sides of the
binding-energy difference, electron rest masses cancel in `B_i`.  The
published equation is reconciled to XNet's explicitly approximate
translational mass `A_i/N_A`, rather than silently substituting actual nuclear
masses.

### Preflight decisions

- SN160 is a candidate, not a requirement.  It must pass an expanded-network
  boundary comparison; SN231 will be used if it does not.
- The three initially proposed states were not treated as fixed answers.  The
  final set will consider proton-rich and harder thermodynamic conditions as
  alternatives while retaining exactly three robust states.
- Pynucastro/Microphysics overlap is primarily within that software family.
  The one identified XNet-author contribution to historical Microphysics EOS
  declarations is unrelated to NSE.  Scientific implementation lineage and a
  literal contributor-list intersection are recorded separately.
No expected composition, validation tolerance, or observed XNet discrepancy
was used in the network/state decision recorded above.

## Final design after independent preflight

### Network selection and boundary evidence

The validation uses the tracked `test/build_net/sunet.torch489` selection,
built to the retained `network/sunet` and `network/netwinv` files.  It contains
489 unique ordered species: n, p, d, t, and bound nuclei through `tc91`
(`Z <= 43`, `A <= 91`).  It strictly contains every SN160 and SN231 species.
The retained file hashes are:

- `sunet`: `930fcd3e43b9fbbac1ce2441b69368b341baf488e29edff53935d9f3cbf05d1f`;
- `netwinv`: `efd4551e7597c732a49931e699453cf469e9fd9f9296bc1b3e841bd21ed42379`;
- ordered species: `8620804c59002cad785f539d1beeddcc13e4ec5aa3b77fb7fd0d223120795038`.

This choice follows quantitative preflight, not network size alone:

| Comparison using identical common nuclear inputs | State | Mass outside smaller set | Shared L1 | Shared L-infinity |
| --- | --- | ---: | ---: | ---: |
| SN160 versus SN231 | `(7,1e9,0.50)` | `1.20e-3` in SN231-only species | `1.44e-3` | `4.25e-4` |
| SN160 versus SN231 | `(7,1e9,0.45)` | `1.41e-1` in SN231-only species | `1.76e-1` | `4.68e-2` |
| torch489 versus full 7852 | `(7,1e9,0.50)` | `2.254e-6` | `2.276e-6` | `5.504e-7` |
| torch489 versus full 7852 | `(7,1e9,0.45)` | `1.911e-6` | `5.856e-6` | `8.657e-7` |
| torch489 versus full 7852 | `(6.5,1e7,0.55)` | `7.241e-6` | `7.241e-6` | `5.636e-6` |

These aggregate figures are non-gating network/state-selection evidence, not
part of the ordinary three-state XNet test.  The retained preflight program now
refuses a pairwise comparison unless all shared species have identical masses,
bindings, spins, partition data, translational masses, constants, temperature
grid, and conventions.  The verified shared-input fingerprints are
`685793a11a2b68ce46b4a5fd2ad5bd026aeb735dced82078b8cd32f1fd83aff6`
for SN160/SN231 and
`cf42edfae59d56c1b3d25abc114152e0811c9eba9993ef6e3dbf0c8e743fac00`
for torch489/full.  The exact rerun reconstructs the exploratory manifests and
diagnostic JSON; the full network intentionally skips raw-source reconciliation
because it is not the retained reference network:

```bash
python3 test/nse_validation/extract_inputs.py test/Data_SN160 /tmp/sn160.json
python3 test/nse_validation/extract_inputs.py test/Data_SN231 /tmp/sn231.json
python3 test/nse_validation/extract_inputs.py \
  test/nse_validation/network /tmp/torch489.json
python3 test/nse_validation/extract_inputs.py \
  test/Data_Reaclib20180621 /tmp/full7852.json --skip-raw-reconciliation
python3 test/nse_validation/preflight_states.py \
  /tmp/sn160.json /tmp/sn231.json /tmp/torch489.json /tmp/full7852.json \
  > /tmp/nse-preflight.json
```

SN160 therefore fails the bounded-network requirement, especially at the
report's neutron-rich point.  SN231 improves the result but still leaves about
one percent of the full-network mass outside its set at that point.  Torch489
reduces every retained-state boundary mass below `1e-5`; the symmetric and
proton-rich differences are dominated by the transient `be8` entry in the
full REACLIB set.  A `Ye=0.40` dense trial was rejected because about 70% of
the full-network mass lies outside torch489.  The Seitenzahl 443-species list
was not available as a retained list in the supplied article, and its actual
masses and Rauscher--Thielemann partition inputs would not define XNet's
finite-data problem.  It remains equation and methodology evidence, not a
substitute network.

The build exposed one necessary tooling defect: the mass reader stopped on
the documented `#` unavailable-value rows even when those rows were not
selected.  `test/build_net/partf_module.f90` now skips a row only when its
required mass field is exactly `#`, while malformed rows remain errors and the
existing later error still rejects a requested species with no mass.  Focused
checks cover all three cases.  `network/build_input.namelist` retains the
exact torch489, REACLIB, partition, mass, and disabled weak/neutrino settings
used for the snapshot.

### Exact states

| ID | `T9` [GK] | `rho` [g cm^-3] | `Ye` | Role | Jacobian infinity condition estimate |
| --- | ---: | ---: | ---: | --- | ---: |
| `symmetric` | 7.0 | `1e9` | 0.50 | dense, approximately symmetric iron-group equilibrium | `3.159e3` |
| `neutron_rich` | 7.0 | `1e9` | 0.45 | dense neutron-rich iron-group equilibrium | `8.321e3` |
| `proton_rich_low_density` | 6.5 | `1e7` | 0.55 | proton-rich, 100-times lower density, lower temperature, and partition-function interpolation | `9.353e1` |

The third state replaces the report's hot symmetric point.  It supplies a
proton-rich regime, a material thermodynamic change, and a non-grid-node
partition-function interpolation.  All three are static equilibrium benchmark
points in a high-temperature range commonly relevant to NSE; whether a dynamic
system attains NSE also depends on density, reaction rates, and its available
timescale.  All three states converged from four
starting offsets using both analytic- and numerical-Jacobian Newton routes.
XNet is also run at every retained state from its default guess and from a
supplied root offset by `(+2,-2)`; both results must independently pass the
reference gates.  This is numerical-robustness evidence, not scientific
authority.
The more extreme neutron-rich trials were useful preflight but were rejected
for finite-network boundary sensitivity, not because XNet was consulted.

### Reference quality, stored data, and tolerances

Generator `xnet-independent-nse-reference-v3` uses the standard library only;
`generate_reference.py` calls
`reference_solver.py`; neither imports or executes XNet.  The solver uses
CPython 3.13.0 with `decimal` 1.70/libmpdec 4.0.1 and 50 requested decimal
digits with 20--30 working guard digits.  Checks at 35
and 65 requested digits produce
identical stored binary64 mass fractions, so the report's proposed
80/120/180-digit ladder was unnecessary.  The largest retained reference
mass or charge residual is `2.25e-28`; the largest analytic- versus
finite-difference-Jacobian L1 disagreement is `2.85e-49`.  The smallest
retained mass fraction is
`5.25e-44`, far above XNet's protected-exponential floor, and the largest is
`8.80e-1`, below its upper clip.

The scientific dataset SHA-256 is
`c1651b543e821b05ba9c9d4449011df2d65db06e47a75c0a670e60083f5de4ac`.
The Fortran-facing `reference.dat` SHA-256 is
`df5b50465152e99aa97b4cb275cb938ea25f019bb8154fe8f7d24b7ea151dffd`.
`reference.json` records the generator/source hashes, exact binary64 state
inputs and constants, residuals, starts, precision checks, complete vectors,
and derivation details.  The scientific hash covers that complete durable
record except its two derived hash fields and the separately checked,
provenance-only scientific-input identity, so changing a residual or other
numerical-quality diagnostic changes the fingerprint.

The frozen JSON retains the original fork issue URL in its historical
generation metadata. That field is part of the recorded dataset identity; it
is provenance, not a current issue or build instruction.

The residual gate starts from XNet's configured `1e-8` function tolerance.
An inspected upper bound of 32 correctly-rounded binary64 operation
equivalents per species plus complete serial accumulation gives
`gamma_16136 = 1.792e-12`, rounded upward to a `2e-12` arithmetic budget.
This is an explicit portability assumption rather than a rigorous bound on
every compiler's optimized transcendental library.  For each state the
generator solves the four corners and four edge midpoints sampled on the
resulting mass/XNet-charge residual box.  The composition gate is the maximum
sampled displacement plus analytic/numerical-Jacobian, stored-binary64, and
arithmetic budgets, rounded upward to four significant digits.  These are
sampled, GNU-serial-configuration-qualified acceptance envelopes, not rigorous
global sensitivity bounds over the full residual box.  That distinction does
not alter the frozen gate values.

| Gating quantity (dimensionless absolute error) | symmetric | neutron rich | proton-rich low density |
| --- | ---: | ---: | ---: |
| `abs(sum(X)-1)` | `1.0002e-8` | `1.0002e-8` | `1.0002e-8` |
| `abs(sum((Z/A-Ye)X))` | `1.0002e-8` | `1.0002e-8` | `1.0002e-8` |
| reconstructed `abs(sum((Z/A)X)-Ye)` | `1.5005e-8` | `1.4505e-8` | `1.5506e-8` |
| complete-vector L1 | `4.192e-7` | `9.608e-7` | `4.689e-8` |
| complete-vector L-infinity | `9.529e-8` | `1.050e-7` | `2.057e-8` |

Finiteness, nonnegativity, unique name and `A/Z/N` identity, completeness, and
exact order have no fallback tolerance.  The dominant species are reported by
name and their errors are checked through the complete-vector L-infinity gate;
no looser dominant-only gate exists.  Candidate XNet output is never
renormalized.

### Reproduction and ordinary test separation

The retained JSON's generator/source payload was frozen at historical
`jaharris87/XNet` fork commit
`66e3ee7399e522011aae17fc714a942511f47041`; the earlier verification record
began at fork commit `65271bcbeea430534c1adc92748ef13bea10c228`.
Neither commit is expected to exist in a normal `starkiller-astro/XNet` clone.
The later `scientific_input_sha256` field is the separately checked
preservation identity documented above.  Reproducing the historical payload
and `reference.dat` therefore requires fetching the fork's development history,
then using the recorded CPython and libmpdec versions:

```bash
git fetch https://github.com/jaharris87/XNet.git development
git cat-file -e 65271bcbeea430534c1adc92748ef13bea10c228^{commit}
git cat-file -e 66e3ee7399e522011aae17fc714a942511f47041^{commit}
mkdir /tmp/xnet-nse-reference-snapshot
git archive 66e3ee7399e522011aae17fc714a942511f47041 | \
  tar -x -C /tmp/xnet-nse-reference-snapshot
cd /tmp/xnet-nse-reference-snapshot
PYTHONDONTWRITEBYTECODE=1 python3 test/nse_validation/generate_reference.py \
  test/nse_validation/network /tmp/reference.json /tmp/reference.dat
shasum -a 256 /tmp/reference.json /tmp/reference.dat
cmp /tmp/reference.dat /path/to/current/XNet/test/nse_validation/reference.dat
python3 - /path/to/current/XNet/test/nse_validation/reference.json \
  /tmp/reference.json <<'PY'
import json
import sys

retained = json.load(open(sys.argv[1], encoding="utf-8"))
generated = json.load(open(sys.argv[2], encoding="utf-8"))
retained["network"].pop("scientific_input_sha256")
if retained != generated:
    raise SystemExit("retained historical payload differs from generation snapshot")
PY
```

Running the generator from a later checkout is a manual scientific-result
recalculation, not a byte reproduction of historical provenance.  Its
`reference.dat` may be compared with the retained data, but its JSON records
the later generator and production-source hashes and must not replace the
retained JSON merely to make those provenance fields current.

The ordinary test target only verifies the retained hashes/data and runs the
XNet comparator.  It has no rule that invokes `generate_reference.py`, no
update mode, and no candidate-output input to the generator.

The optional secondary check requires pynucastro installed from the pinned
checkout and verifies both that checkout's commit and the imported material
source files before using it:

```bash
python3 test/nse_validation/corroborate_pynucastro.py \
  /path/to/pynucastro-checkout-at-a7268f86
```

It is deliberately absent from the ordinary test target and never writes the
retained expected data.

### Effectiveness coverage and supported claim

Focused tests reject controlled changes to a dominant expected species, the
complete vector, reconstructed `Ye`, XNet's charge residual, normalization,
species order, duplicate identity, missing identity, physical `A/Z/N`
identity, NaN, negativity, and a 10 keV `co55` binding-energy input.  Separate
metric-comparator boundary tests set each already computed numeric metric to
the next representable binary64 value above its gate; they are not constructed
489-species physical-vector tests.  The independent Python test changes the
same binding input in the generator path and requires both composition norms
to move beyond their gates.  These mutations operate on parsed, identity-valid
data or in-memory scientific inputs; none depends on a metadata parse failure.

The strongest supported claim is limited to the complete, unscreened static
ideal-NSE composition for this exact 489-species set and these three states,
using the retained XNet mass, partition, constant, and translational-mass
conventions.  It does not validate screening, reaction or weak rates, the
timescale for reaching NSE, other networks/states, or the absolute accuracy of
the underlying nuclear data.  The several-parts-per-million full-network
boundary differences are larger than the XNet/reference discrepancies and
remain outside this finite-network claim.  The recorded provenance and test
results do not replace independent human scientific review.
