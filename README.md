Welcome to XNet 
======

XNet is a thermonuclear reaction network for astrophysical applications,
using a variety of temporal integration techniques. It is written in portable modern Fortran and makes use of a variety of matrix solving packages to obtain excellent speed (as much as 50% of peak) on as many platforms as possible.

The version contained in the repository is targeted at standalone
execution (post-processing nucleosynthesis). However, with the appropriate
interfaces, XNet has been used in a variety of large scale simulation
packages, including FLASH, CHIMERA and several other codes.

Building XNet
-------------

Run GNU Make from the repository root. The default GNU optimized build places
programs under `build/GNU-OPT/bin/`:

```bash
make -j
```

Major build choices automatically select a different readable directory. For
example, this command uses `build/GNU-DEBUG/bin/`:

```bash
make CMODE=DEBUG -j xnet
```

The automatic name reflects the compiler family, build mode, MPI/OpenMP,
accelerator, non-default EOS, and non-default matrix-solver choices. The
directory's `config.txt` still rejects reuse when other material settings
change.

Explicit `BUILD_NAME` and `BUILD_DIR` remain available when a developer wants
to choose the location. The build produces object files, module files,
preprocessed sources, and programs under the selected directory. See
`Makefile.opt` for the principal compiler, parallel, solver, EOS, and
accelerator choices.

Before you dig deep into the source code, we suggest you learn the
basics of how thermonuclear reaction networks function, from one of the following sources.

"Thermonuclear Reaction Networks for Astrophysics" by Raph Hix, http://astro.phys.utk.edu/_media/gallery:networks.pdf

Tools and Toys for Nuclear Astrophysics: Nuclear Reaction Network Techniques by Frank Timmes, http://cococubed.asu.edu/talk_pages/jina.shtml

Thermonuclear Kinetics in Astrophysics by Raph Hix and Brad Meyer http://dx.doi.org/10.1016/j.nuclphysa.2004.10.009

Supernovae and Nucleosynthesis: An Investigation of the History of Matter, from the Big Bang to the Present by Dave Arnett http://press.princeton.edu/titles/5859.html

Someday, we'll make a list of publications using XNet, but that will wait...
