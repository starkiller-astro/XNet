# Network preprocessing fixture

This fixture is a synthetic six-species network for the `net_preprocess` and
`net_setup` contract tests. It contains no private data. Nuclear masses are
small public-format test values and every partition function is deliberately
one, because the checks cover preprocessing behavior rather than rate-library
quality or partition-function physics.

The positive `netsu` contains:

- a REACLIB-style electron-capture `n -> p` record;
- matched forward and reverse `he4 + c12 <-> o16` records;
- matched forward and reverse triple-alpha records, which exercise repeated
  participant multiplicities and the one- and three-reactant layouts; and
- a syntactically valid `he4 + ne20 -> mg24` record that must be dropped
  because `mg24` is not in `sunet`.

The positive rate headers deliberately contain incorrect Q values. The
semantic verifier calculates the expected values from the fixture mass
excesses, so copying rather than recomputing a header Q value is detectable.

`netsu_truncated` ends after a reaction header without its coefficient record.
`netsu_inconsistent` leaves a required chapter-4 reactant blank. Both must
produce a detectable nonzero process result.
