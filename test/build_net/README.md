REACLIB Reaction Network Generator for XNet
============

Tools for creating customized reaction networks for XNet from JINA REACLIB rates and FFN weak rates.

## Getting started

The file `input.namelist` contains controls for incorporating non-REACLIB rates and the path for the new `netsu`, `netwinv`, `netweak`, and `netneutr` files to be used by XNet.

The `./*_data` directories contain the rate databases:
* `./reaclib_data`: ReaclibVX.X rates in slightly modified Reaclib2 format w/o ch. 9/10/11
  * Data Files:
    * `reaclib_JINAvXX` - Parameterized Rates
    * `winvne_JINAvXX`  - Tabulated Partition Functions
  * References:
    * Cyburt et al., ApJS 189 (2010) 240
    * https://groups.nscl.msu.edu/jina/reaclib/db/library.php?action=viewsnapshots
* `./partf_data`: Tabulated partition function data
  * Data Files:
    * `winvne_JINAvXX`  - Tabulated Partition Functions
  * References:
    * Cyburt et al., ApJS 189 (2010) 240
    * https://groups.nscl.msu.edu/jina/reaclib/db/library.php?action=viewsnapshots
* `./mass_data`: Atomic mass evaluations used from JINA Nuclide Database
  * Data Files:
    * `mass_X.dat`      - Tabulated Atomic Mass Tables
  * References:
    * https://groups.nscl.msu.edu/jina/nucdatalib/evaluations/3
* `./ffn_data`: Tabulated EC/PC rates in FFN-style formatting
  * Files:
    * `lmpffnoda.data` - Tabulated Rates
  * References:
    * Fuller et al., ApJ 293 (1985) 1
    * Oda et al., ADNDT 56 (1994) 231
    * Langanke & Martinez-Pinedo, ADNDT 79 (2001) 1
* `./neutrio_data`: Tabulated neutrino capture rates (not public and not
  included in this repository)
  * Files:
    * `neutrino.data` - Tabulated Rates
  * References:
    * private communication with Carla Froehlich (2015)


## To run:

1. Compile the code:
  * `make`
2. Create a list of nuclei in the file `sunet` in the run directory that is a subset of nuclei in the REACLIB database
  * see `sunet.example` for correct format
3. Configure `input.namelist`
4. Run the code:
  * `./build_net`


## Example(s)

See `sunet.example` for a simple example of an alpha-network with neutrons and protons

Some other useful sunet files are also included

## Retained data provenance

The files retained in this tree have exact-byte provenance independent of
claims about their original scientific evaluations. XNet imported the
`build_net` subtree at commit
`90e9363d5f9443a8ad2d5e986232c1f60bb5b96a` on 2019-03-05. The subtree
content came from `build_net` split commit
`79a2a46654aea68161de0dee482e79880ca9a25d` and the public
`jaharris87/build_net` database commit
`77141ca2a3dfc9fa9fd52ef0fcf39a49d74c08e1` (2017-01-24).

The primary retained inputs are bound to the following upstream Git blobs
and local SHA-256 values:

| Retained file | Upstream blob | SHA-256 |
| --- | --- | --- |
| `mass_data/mass_reac1.dat` | `acd42b416e52edea020990d68631fb2b8063275d` | `0d068a92c6694485e117df2da8b081fe1a7cd7b80757303d1f72d4bfc16fff50` |
| `partf_data/winvne_JINAv22` | `5bac4be2a2095bfd6a6fe3b1282bea4bcd6040ea` | `5fc838645cac2f5e2eb9ddfd5a60bb41df8a6378b8f28748b1b47ed6148555aa` |
| `reaclib_data/reaclib_JINAv22` | `0a492936abae2808d3f27b4ed71d04294a7bb70a` | `15e670a5b39a0ffce6a3f06db0d76c355dff7c8b3ff2280dd819004d9756214b` |

The official JINA REACLIB site identifies ReaclibV2.2 as a public snapshot
dated 2016-11-14 and cites Cyburt et al., *ApJS* 189 (2010) 240. This is
contextual scientific provenance for the `JINAv22` filenames:
<https://reaclib.jinaweb.org/library.php?action=details&libindex=ReaclibV2.2>.
The retained `reaclib_JINAv22` and `winvne_JINAv22` filenames associate them
with that snapshot at the library level. The Git and SHA-256 values above
prove their retained bytes only against `build_net@77141ca`; no immutable
official release artifact or checksum was retained to establish byte identity
with JINA or a complete derivative chain.

The original JINA download timestamp and per-record release or publication
identifier for the `reac1` mass evaluation are not preserved. The AME03,
AME11, extrapolated, and FRDM files retain only their evaluation/model labels;
those labels do not establish an exact upstream release artifact. Later
`ffngeff.dat` and `updated_rate_table.txt.gz` additions likewise have no
equivalent authoritative release binding in this repository. These gaps are
provenance limitations, not reasons to refresh the retained data in place.
