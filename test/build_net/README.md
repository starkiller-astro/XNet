REACLIB Reaction Network Generator for XNet
============

Tools for creating customized reaction networks for XNet from JINA REACLIB rates and FFN weak rates.

## Getting started

The file `input` contains controls for incorporating non-REACLIB rates and the path for the new `netsu`, `netwinv`, `netweak`, and `netneutr` files to be used by XNet.

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
* `./neutrio_data`: Tabulated neutrino capture rates (not public)
  * Files:
    * `neutrino.data` - Tabulated Rates
  * References:
    * private communication with Carla Froehlich (2015)


## To run:

1. Compile the code:
  * `make`
2. Create a list of nuclei in the file `sunet` in the run directory that is a subset of nuclei in the REACLIB database
  * see `sunet.example` for correct format
3. Configure `input` file
4. Run the code:
  * `./build_net`


## Example(s)

See `sunet.example` for a simple example of an alpha-network with neutrons and protons

Some other useful sunet files are also included