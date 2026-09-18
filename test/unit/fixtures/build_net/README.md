# Synthetic `build_net` sources

`run_build_net_contracts.py` creates the fixture below each configuration's
ignored `build-net-work/` directory. The source catalog is intentionally
defined in code so the fixed-width REACLIB, partition, mass, and weak records
remain reviewable without tracking hundreds of repeated table values.

The master partition/mass catalog contains `n`, `p`, `he4`, `c12`, `o16`, and
the deliberately unrequested `ne20`. The requested file uses `N01`, `H1`, and
upper-case isotope names to exercise normalization. The REACLIB source has one
retained forward/reverse pair and one `ne20` rate that must be excluded. The
public synthetic weak source has one `p`/`n` pair and 143 constant table
points. All six required mass filenames contain the same six synthetic
records; `ame11` is selected by the partition records.
The raw partition records deliberately use different per-species mass offsets
from the selected AME values so the weak-Q checks distinguish those contracts.

No production database or private neutrino source is copied or opened. The
runner also constructs duplicate, blank, unavailable, malformed-namelist,
malformed-REACLIB, and missing-mass negative cases from this definition.
