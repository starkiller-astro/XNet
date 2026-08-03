Welcome to XNet 
======

XNet is a thermonuclear reaction network for astrophysical applications,
using a variety of temporal integration techniques. It is written in portable modern Fortran and makes use of a variety of matrix solving packages to obtain excellent speed (as much as 50% of peak) on as many platforms as possible.

The version contained in the repository is targeted at standalone
execution (post-processing nucleosynthesis). However, with the appropriate
interfaces, XNet has been used in a variety of large scale simulation
packages, including FLASH, CHIMERA and several other codes.

Before you dig deep into the source code, we suggest you learn the
basics of how thermonuclear reaction networks function, from one of the following sources.

"Thermonuclear Reaction Networks for Astrophysics" by Raph Hix, http://astro.phys.utk.edu/_media/gallery:networks.pdf

Tools and Toys for Nuclear Astrophysics: Nuclear Reaction Network Techniques by Frank Timmes, http://cococubed.asu.edu/talk_pages/jina.shtml

Thermonuclear Kinetics in Astrophysics by Raph Hix and Brad Meyer http://dx.doi.org/10.1016/j.nuclphysa.2004.10.009

Supernovae and Nucleosynthesis: An Investigation of the History of Matter, from the Big Bang to the Present by Dave Arnett http://press.princeton.edu/titles/5859.html

Someday, we'll make a list of publications using XNet, but that will wait...

## Characterized GNU serial build

The repository records a build-smoke configuration using GNU Fortran in debug
mode with serial execution, the dense solver, vendored NETLIB, and the
STARKILLER EOS. Run it from the repository root:

```bash
make -C source clean
make -C source -j2 \
  CMODE=DEBUG PE_ENV=GNU \
  MPI_MODE=OFF OPENMP_MODE=OFF GPU_MODE=OFF \
  MATRIX_SOLVER=dense LAPACK_VER=NETLIB EOS=STARKILLER \
  xnet net_setup xnse
```

The production build is in-source, so clean before changing configuration.
A successful build is compilation and link evidence for this exact
configuration; it is not scientific validation or a support claim.

## Dependency-light unit tests

The isolated GNU serial unit-test targets are:

```bash
make -C source test_unit
make -C source test_unit_selfcheck
```

The self-check deliberately exercises a failing fixture and verifies that the
failure propagates through Make. See [the current architecture](doc/architecture.md)
and [the validation plan](doc/validation-plan.md) for scope and limitations.


