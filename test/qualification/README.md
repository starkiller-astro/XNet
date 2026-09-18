# MPI and OpenMP zone-semantics qualification

This directory contains a configuration qualification for XNet's current MPI
and coarse-batch OpenMP zone calculations. It is not part of the ordinary
physical regression list and stores no MPI- or OpenMP-specific reference result.

The fixture has ten zones and `nzbatchmx=4`, so execution has three batches and
the final batch has two inactive lanes. Every zone has a distinct abundance
file and a distinct constant thermodynamic history. The abundance inputs vary
the active species and reaction paths, while temperature, density, and target
time vary the outer-zone state. The alpha network retains the representative
marked outer-zone, inner-species/reaction, mask, and reduction paths without
adding another physical scenario to the regression catalog.

The runner executes one serial reference, one two-rank MPI run, a two-rank
non-root missing-abundance check, and three two-thread OpenMP runs. It requires:

- exactly one nonempty ASCII and binary history for each global zone 1-10;
- no output for inactive zones 11-12;
- runtime topology records proving ranks 0-1 or threads 1-2 actually ran;
- complete diagnostic records for each zone exactly once across all workers;
- each ASCII filename's final row to agree with that zone's diagnostic state;
- exact serial/MPI/OpenMP equality for the normalized diagnostic endpoints and
  final ASCII energy-generation, neutrino-loss, and timestep fields, using the
  established `batch_alpha` endpoint policy and ignoring worker order;
- exact equality across the three OpenMP repetitions; and
- nonzero status within the requested timeout when zone 5's required abundance
  file is removed. With the current three rank-strided batches and two ranks,
  zones 5-8 belong to rank 1.

The present `OPENMP_MODE=ON` implementation distributes the coarse `ibatch`
loop. It is characterized here, not prescribed. When CPU OpenMP moves to the
computational loops marked through `xnet_macros.fh` with `!XDIR XLOOP*`, rerun
this same qualification; worker assignment and diagnostic order are
not assertions.

## Focused helper checks

From the repository root:

```bash
python3 -m pytest -q test/qualification/test_parallel_zones.py
```

These checks prove that the qualification helpers reject duplicate/missing
zones, an off-by-one zone 11 output, filename/result misassociation, copied
endpoint or ASCII energy values, missing MPI/OpenMP workers, a zero-status
failure probe, and a process-group timeout.

## Configuration run

Use separate build directories for the serial, MPI, and OpenMP configurations:

```bash
make BUILD_NAME=parallel-serial -j \
  CMODE=OPT PE_ENV=GNU MPI_MODE=OFF OPENMP_MODE=OFF xnet

make BUILD_NAME=parallel-mpi -j \
  CMODE=OPT PE_ENV=GNU MPI_MODE=ON OPENMP_MODE=OFF xnet

make BUILD_NAME=parallel-openmp -j \
  CMODE=OPT PE_ENV=GNU MPI_MODE=OFF OPENMP_MODE=ON xnet

python3 test/qualification/parallel_zones.py \
  --serial-executable="$PWD/build/parallel-serial/bin/xnet" \
  --mpi-executable="$PWD/build/parallel-mpi/bin/xnet" \
  --openmp-executable="$PWD/build/parallel-openmp/bin/xnet" \
  --mpi-launcher="$(command -v mpiexec)" \
  --work-root=/tmp/xnet-parallel-zone-work \
  --timeout=60
```

Use fresh empty output paths for each run. The runner retains stdout, stderr,
status, command, diagnostics, histories, and a final JSON summary below the
chosen work root. Record the compiler, MPI implementation, OpenMP settings,
commands, host, date, and any untested environments in the PR; a successful
run qualifies only that exact configuration.

The manual Frontier HIP/ROCm qualification that reuses this fixture is
documented in [`frontier/README.md`](frontier/README.md). It is a separate GPU
check and does not change this CPU/MPI/OpenMP comparison policy.

The example intentionally uses the launcher's ordinary slot and binding
policy. If a scheduler allocation exposes fewer than two slots, request an
allocation with enough slots rather than allowing oversubscription. Allowing
oversubscription can change Open MPI's binding policy even when two ranks fit
on the host. `--mpi-launcher-argument` remains available for site-required
launcher options.
