# EOS scientific-reference fixture

This test validates XNet's default STARKILLER Helmholtz EOS against five
independently calculated states from the direct Timmes EOS. These are distinct
implementations. The direct Timmes EOS evaluates the electron-positron
integrals and is described by
[Timmes & Arnett (1999)](https://scixplorer.org/abs/1999ApJS..125..277T/abstract).
The Helmholtz EOS interpolates a tabulated Helmholtz free energy and is
described by
[Timmes & Swesty (2000)](https://scixplorer.org/abs/2000ApJS..126..501T/abstract).

The authoritative source is
[Cococubed's stellar-EOS page](https://cococubed.com/code_pages/eos.shtml),
which distributes the two implementations separately as `timmes.tbz` and
`helmholtz.tbz`. Cococubed states that the Timmes EOS is its comparison
reference and that the table used by the Helmholtz EOS is calculated from the
Timmes EOS. The reference calculation here uses `eosfxt` from the Cococubed
`timmes.tbz` archive and does not use XNet's tracked table or adapted
STARKILLER implementation.

The archives downloaded on 2026-08-08 have these SHA-256 values:

| Archive | Authoritative source | SHA-256 |
| --- | --- | --- |
| `timmes.tbz` | [Cococubed](https://cococubed.com/codes/eos/timmes.tbz) | `e20c7d27e66c240486a3397a649c49673e33100284f0905cd6fa9893dbad30a9` |
| `helmholtz.tbz` | [Cococubed](https://cococubed.com/codes/eos/helmholtz.tbz) | `fc55ca3b188598ed19f9dbf63bacf1033676a22d2065f420b67375a493b91eeb` |

The SHA-256 values below are for the immutable raw files, before any harness
edits:

| File in `timmes.tbz` | SHA-256 |
| --- | --- |
| `eosfxt.f90` | `4c3e45924c00cff751885377c936cdd44793838c3b38db6ba4c5fec1c58710d3` |
| `const.dek` | `bb6089c86b183d119a8e03ec5d893626e8a4a2f168bb201a09060a4996fc8ab1` |
| `vector_eos.dek` | `a76e72e174aa1440a2c8321968313f98f8d2cb8f488deb4283e3a2785058d427` |
| `implno.dek` | `e62d95098f636eb3d90871f3e1f07f462f184af4557a23cfdf3d7228a723b2ea` |

For the isolated reference run, the bundled example `program teos` was renamed
to an unused subroutine so that `timmes_reference_driver.F90` could provide the
main program. This harness-only edit does not modify `eosfxt` calculations.
GNU Fortran 16.1.0 compiled and ran the driver on 2026-08-08.

The reference can be reproduced from an XNet checkout with:

```bash
reference_dir=$(mktemp -d)
cd "$reference_dir"
curl -fL https://cococubed.com/codes/eos/timmes.tbz -o timmes.tbz
shasum -a 256 timmes.tbz
tar -xjf timmes.tbz
cd timmes
shasum -a 256 eosfxt.f90 const.dek vector_eos.dek implno.dek
perl -0pi -e 's/\A      program teos/      subroutine teos_unused/' eosfxt.f90
cp <XNET_CHECKOUT>/test/unit/fixtures/eos/timmes_reference_driver.F90 .
gfortran -O0 -I. eosfxt.f90 timmes_reference_driver.F90 -o timmes_reference_driver
./timmes_reference_driver
```

After the documented harness edit, the four SHA-256 values are
`8604501a0eeba64d71f2fe2ff9d99c6e707c4c9d3e41f1e12dc4ffb44470cc00`,
`bb6089c86b183d119a8e03ec5d893626e8a4a2f168bb201a09060a4996fc8ab1`,
`a76e72e174aa1440a2c8321968313f98f8d2cb8f488deb4283e3a2785058d427`,
and `e62d95098f636eb3d90871f3e1f07f462f184af4557a23cfdf3d7228a723b2ea`
in the table's file order. A clean reproduction with those inputs returned the
values below exactly. The tracked `timmes_reference_driver.F90` used for that
run has SHA-256
`571332fdb0ce65424cebdad6d7039c142db9b9d0a553021a44b02e9c36b558c9`.

The raw independent results are:

| State | T (K) | rho (g cm^-3) | abar | zbar | cv (erg g^-1 K^-1) | eta | d eta/dT (K^-1) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.2e7 | 2.5e2 | 1.2958963282937364 | 1.1015118790496758 | 2.0653230156968951e8 | -1.8807503609806382 | -4.9592562974862513e-8 |
| 2 | 1.7e8 | 3.0e6 | 4 | 2 | 4.6199576060239784e7 | 18.319585249344883 | -1.0876372950003569e-7 |
| 3 | 8.0e8 | 2.0e9 | 12 | 6 | 1.8932479108972270e7 | 67.68924196485338 | -8.4721540555149227e-8 |
| 4 | 4.0e9 | 7.0e7 | 56 | 26 | 9.8736503214137957e7 | 2.8245988707608238 | -1.0485292044799550e-9 |
| 5 | 9.0e9 | 1.0e10 | 19.764705882352942 | 9.529411764705882 | 4.3195927677688472e7 | 10.296830235239677 | -1.2091574746119109e-9 |

These states cover negative-eta nondegenerate material, degenerate and strongly
degenerate material, a hot moderately degenerate state, and a hot dense state.
They exercise mixtures as well as single-species compositions.

The test maps Timmes `cv` to XNet's exposed MeV nucleon^-1 GK^-1 field with
`cv * xnet_constants::amu * 1e9`, maps `d eta/dT` to GK^-1 with a factor of
`1e9`, and compares eta directly. The relative tolerances are `1e-5` for
`cv`, `2e-5` plus `5e-6` absolute for eta, and `3e-4` plus `1e-6` absolute for
the temperature derivative. The state-by-state differences from the tracked
STARKILLER table are:

| State | cv relative | eta absolute | eta relative | d eta/dT9 relative |
| --- | ---: | ---: | ---: | ---: |
| 1 | 1.0961e-6 | 2.0118e-6 | 1.0697e-6 | 6.6475e-5 |
| 2 | 2.9728e-8 | 1.7977e-4 | 9.8127e-6 | 1.5113e-4 |
| 3 | 4.5224e-8 | 3.6200e-5 | 5.3480e-7 | 1.4077e-4 |
| 4 | 5.6173e-6 | 6.0867e-7 | 2.1549e-7 | 8.5062e-5 |
| 5 | 9.4704e-7 | 5.7354e-6 | 5.5701e-7 | 1.7546e-4 |

The limits cover the observed tabulation/interpolation differences with modest
compiler margin and remain well below the controlled one-percent reference
perturbation.

The tracked table used by the STARKILLER Helmholtz implementation has SHA-256
`c9a57c26c6fd2b2b378b9d5295ca1214022f6fec6289d038b47bf8c8938881a1`.
It is byte-for-byte identical to `helmholtz/helm_table.dat` in the Cococubed
`helmholtz.tbz` archive above. This establishes the table's authoritative
Helmholtz provenance; the direct Timmes run remains the independent numerical
reference. The identity check is reproducible with:

```bash
curl -fL https://cococubed.com/codes/eos/helmholtz.tbz -o helmholtz.tbz
shasum -a 256 helmholtz.tbz
tar -xjf helmholtz.tbz
shasum -a 256 helmholtz/helm_table.dat
cmp -s <XNET_CHECKOUT>/tools/starkiller-helmholtz/helm_table.dat \
  helmholtz/helm_table.dat
```
