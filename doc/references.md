# XNet References

## XNet Core References

These papers describe methods and physics directly implemented in XNet.

- **Hix & Thielemann (1999)**, ApJ 511, 862 — Foundational XNet methods paper:
  implicit integration, Jacobian construction, sparse solver strategy.
- **Hix & Meyer (2006)**, Nuclear Physics A 777, 188 — Background on screening
  and NSE methods in reaction networks (specific implementations in XNet may
  differ; see comments in `xnet_screening.F90` and `xnet_nse.F90`).
  Referenced in `xnet_data.F90`, `xnet_integrate.F90`.
- **Seitenzahl et al. (2009)**, ADNDT 95, 96 (doi:10.1016/j.adt.2008.08.001) —
  In-depth treatment of NSE in the context of thermonuclear supernovae.
  Referenced in `xnet_nse.F90`.
- **Fuller, Fowler & Newman (1982, 1985)** — Electron and positron capture rates
  (FFN format). Referenced in `xnet_ffn.F90`, `xnet_data.F90`.
  See also Fuller et al. (1985), ApJ 293, 1.
- **Arnett (1996)**, *Supernovae and Nucleosynthesis: An Investigation of the
  History of Matter, from the Big Bang to the Present*, Princeton Univ. Press —
  General reference for nucleosynthesis theory.

## Numerical Methods

### Integration

- **Longland, Martin & José (2014)**, A&A 563, A67
  (doi:10.1051/0004-6361/201321958) — Modern treatment of implicit integration
  methods, timestep control, and convergence testing for thermonuclear reaction
  networks. Referenced in `xnet_integrate_bdf.F90`.
- **Timmes (1999)**, ApJS 124, 241 — Benchmark of integration methods and
  solvers for thermonuclear reaction networks. Referenced in
  `xnet_integrate_bd.F90`. Also see:
  - https://cococubed.com/research_pages/research.shtml
  - https://cococubed.com/code_pages/codes.shtml
- **Brown, Byrne & Hindmarsh (1989)**, SIAM J. Sci. Stat. Comput. 10, 1038
  (doi:10.1137/0910062) — VODE: A variable-coefficient ODE solver.
  Referenced in `xnet_integrate_bdf.F90`.
- **Jackson & Sacks-Davis (1980)**, ACM Trans. Math. Softw. 6, 295
  (doi:10.1145/355900.355903) — Alternative implementation of variable
  step-size multistep formulas for stiff ODEs.
  Referenced in `xnet_integrate_bdf.F90`.
- **Byrne & Hindmarsh (1975)**, ACM Trans. Math. Soft. 1, 71
  (doi:10.1145/355626.355636) — Polyalgorithm for numerical solution of ODEs.
  Referenced in `xnet_integrate_bdf.F90`.
- **Gear (1971)**, Commun. ACM 14, 176 (doi:10.1145/362566.362571) — Automatic
  integration of ordinary differential equations.
  Referenced in `xnet_integrate_bdf.F90`.
- **Nordsieck (1962)**, Mathematics of Computation 16, 22
  (doi:10.1090/S0025-5718-1962-0136519-5) — Numerical integration of ODEs.
  Referenced in `xnet_integrate_bdf.F90`.

### Nonlinear Solvers

- **Dennis & Schnabel (1996)**, *Numerical Methods for Unconstrained
  Optimization and Nonlinear Equations* (doi:10.1137/1.9781611971200) —
  Newton-Raphson methodology. Referenced in `xnet_nse.F90`.
- **Press, Teukolsky, Vetterling & Flannery (1992, 1996)**, *Numerical
  Recipes* — General numerical methods. Referenced in `xnet_integrate_bd.F90`.

## Screening and Coulomb Physics

- **Salpeter (1954)**, Aust. J. Phys. 7, 353 — Electron screening ratio.
  Referenced in `xnet_screening.F90`, `xnet_eos_helm.F90`,
  `xnet_eos_starkiller.F90`, `xnet_eos_bahcall.F90`.
- **DeWitt, Graboske & Cooper (1973)**, ApJ 181, 439 — Intermediate screening
  factors. Referenced in `xnet_screening.F90`, `xnet_eos_helm.F90`,
  `xnet_eos_starkiller.F90`.
- **Graboske, DeWitt, Grossman & Cooper (1973)**, ApJ 181, 457 — Intermediate
  screening factors (Table 4). Referenced in `xnet_screening.F90`,
  `xnet_constants.F90`.
- **DeWitt & Slattery (2003)**, Contrib. Plasma Phys. 43, 279 — Strong
  screening fitting coefficients (Eq. 4). Referenced in `xnet_screening.F90`,
  `xnet_constants.F90`, `xnet_nse.F90`.

## Equation of State

- **Timmes & Swesty (1999)**, ApJS — Helmholtz free energy EOS for stellar
  matter. Referenced in `xnet_screening.F90`.
- **Fukushima (2015)**, Appl. Math. Comput. 259, 708 — Rational function
  expansions for Fermi-Dirac integrals. Referenced extensively in `xnet_fd.F90`,
  `xnet_eos_helm.F90`, `xnet_eos_starkiller.F90`, `xnet_eos_bahcall.F90`.

## Physical Constants

- **Mohr, Newell & Taylor (2016)**, Rev. Mod. Phys. 88, 035009 — CODATA 2014
  recommended values of fundamental constants. Referenced in
  `xnet_constants.F90`. See also http://physics.nist.gov/cuu/Constants/

## Weak Interaction and Neutrino Physics

- **Fuller, Fowler & Newman (1985)**, ApJ 293, 1 — FFN electron/positron
  capture rate tables. Referenced in `xnet_ffn.F90`, `xnet_data.F90`,
  `test/build_net/README.md`.
- **Oda et al. (1994)**, ADNDT 56, 231 — Weak interaction rates for
  sd-shell nuclei. Referenced in `test/build_net/README.md`.
- **Langanke & Martínez-Pinedo (2001)**, ADNDT 79, 1 — Log(ft) electron
  capture rates for pf-shell nuclei. Referenced in `test/build_net/README.md`,
  `test/Data_SN231_logft/README`.
- **Fröhlich (2007)**, PhD thesis — Neutrino-nucleus interaction rate
  implementation. Referenced in `xnet_nnu.F90`.
- **Zinner & Langanke** — Neutrino-induced reaction rates.
  Referenced in `xnet_nnu.F90`.

## Nuclear Data

- **Cyburt et al. (2010)**, ApJS 189, 240 — JINA REACLIB: reaction rates and
  partition functions. Referenced in `test/build_net/README.md`.
  Database: https://groups.nscl.msu.edu/jina/reaclib/db/
- **Anders & Grevesse (1989)**, Geochim. Cosmochim. Acta 53, 197 — Solar
  isotopic abundances. Referenced in `tools/initial_abundance/anders_grevesse_89.txt`.

## Community Codes and Methods

These are similar codes or complementary tools that may be useful as references
for verification tests, alternative implementations, or community best practices.

- **SkyNet** — Lippuner & Roberts (2017), ApJS 233, 18
  (doi:10.3847/1538-4365/aa94cb). Code comparison benchmarks, verification tests.
- **WinNet** — Reichert et al. (2023), ApJS 268, 66
  (doi:10.3847/1538-4365/acf033). Comprehensive verification test suite.
- **AMReX-Astro Microphysics** — https://github.com/AMReX-Astro/Microphysics.
  Helmholtz EOS, reaction networks integrated with AMReX hydro codes.
  XNet's BDF integrator is based on the VBDF code from this repository
  (original by Matthew Emmett; modifications by Adam Jacobs and Mike Zingale).
- **pynucastro** — https://github.com/pynucastro/pynucastro.
  Python library for nuclear reaction rate evaluation and network generation.

## Educational and Review Material

These are referenced in `README.md` as background reading.

- **Hix, W. R.** — "Thermonuclear Reaction Networks for Astrophysics"
  (lecture notes): http://astro.phys.utk.edu/_media/gallery:networks.pdf
- **Timmes, F. X.** — "Tools and Toys for Nuclear Astrophysics: Nuclear
  Reaction Network Techniques" (JINA lecture):
  http://cococubed.asu.edu/talk_pages/jina.shtml
- **Hix & Meyer (2006)**, Nuclear Physics A 777, 188
  (doi:10.1016/j.nuclphysa.2004.10.009) — "Thermonuclear Kinetics in
  Astrophysics" (review).

## External Libraries

- **HSL MA48** — http://www.hsl.rl.ac.uk. Sparse direct solver.
  Used in `xnet_jacobian_MA48.F90`.
- **MAGMA** — Matrix Algebra on GPU and Multicore Architectures (2017).
  Used in `magmaf.F90`.
- **LAPACK** — https://github.com/reference-lapack/lapack.
  NETLIB reference BLAS/LAPACK (subset vendored in `tools/LAPACK/`).
