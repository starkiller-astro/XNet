# Manual Frontier GPU correctness qualification

This manual qualification checks one exact, clean XNet commit on one Frontier
MI250X GPU. It is not CI, a performance benchmark, a scaling claim, or a
permanent GPU endpoint baseline.

A human completes Frontier's interactive RSA authentication and starts the
submission. No credential, token, account name, reservation name, or private
path belongs in the repository. The submitter supplies allocation values on
the command line, and the retained manifest records only whether an account,
QoS, or reservation was supplied. The full run directory remains outside
the checkout.

The package exercises three independent requirements:

- `gpu_linalg_probe.F90` maps two small nonsingular systems to the device and
  calls the public production `LinearSolveBatched` path. It requires a real
  OpenMP target device, successful HIP/rocBLAS handles and factor/solve status,
  and a relative residual no larger than `1e-12` for both systems.
- the ten-zone, `nzbatchmx=4` fixture runs once with the CPU executable and
  once with the GPU executable. The existing parallel-zone runner
  requires zones 1-10 exactly once, correct filename/content association, and
  no output for inactive lanes 11-12.
- `heat_sn160` runs with both executables to reach the larger network,
  self-heating, screening, default Starkiller Helmholtz EOS, dense Jacobian,
  integrator, runtime preprocessing, and accelerator data paths. Its final
  ASCII histories must also exercise and compare at least one nonzero
  neutrino-loss diagnostic.

Both XNet comparisons normalize diagnostic endpoints by global zone and use
the CPU result from the same source archive as the reference. The tracked
`comparison_policy.json` supplies scalar, selected material-species,
complete-vector, normalization, neutrino-loss, and timestep limits. The final
ASCII energy-generation rate is required to be finite and is retained with
its CPU/GPU difference, but it is report-only because its difference of two
accepted-step energies can be ill-conditioned when the true rate is near zero.
This package does not qualify the correctness of that final diagnostic; the
claimed CPU/GPU support is limited to the checked endpoint and ASCII fields.
The policy never creates or updates a canonical GPU result.

## Package contents

- `submit_frontier.py` verifies a clean commit, creates and hashes a `git
  archive`, stages it outside the checkout, submits one Slurm job, and
  waits for the returned report. After Slurm closes the job logs, it finalizes
  the redacted resource evidence and complete output-file inventory.
- `frontier_job.sh` verifies the staged archive hash and embedded Git commit,
  extracts it only after the queued job starts, and then starts one GPU-bound
  Slurm step. Before building, the runner independently verifies that every
  extracted file, symlink, executable bit, and content hash matches the
  archive.
- `frontier_qualification.py` captures the environment, builds CPU and GPU
  configurations from the archive, runs the three checks, compares results,
  inventories every regular output file, and writes
  `qualification_manifest.json`.
- `manifest.schema.json` defines the versioned evidence envelope. The runner's
  stricter standard-library validator requires the successful build, device,
  zone, comparison, hash, and inventory fields.
- `comparison_policy.json` is the reviewed numerical policy, not a reference
  result.

Generated run files, including the validated manifest, remain outside the
repository. Record the source commit, tested software configuration, job,
numerical results, review outcome, and limitations concisely with the pull
request or release record. Do not commit executables, object/module files,
runtime histories, raw environment paths, or account data.

## Human login and submission

Use an explicitly recorded Frontier module set and consult the current OLCF
documentation before a future rerun. OLCF's Frontier guide requires the
`craype-accel-amd-gfx90a` module for HPE Cray Programming Environment OpenMP
offload and documents hipfort as an OLCF module. The last successfully
recorded qualification used CPE 25.09, Cray Fortran 20.0.0, ROCm 6.4.2, and
hipfort 6.4.2.

After authenticating interactively on Frontier:

```bash
module purge
module load PrgEnv-cray
module load cpe/25.09
module load rocm/6.4.2
module load craype-accel-amd-gfx90a
module load hipfort/6.4.2
module load cray-python/3.12.12
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

python3 test/qualification/frontier/submit_frontier.py \
  --account=<allocation> \
  --partition=batch \
  --artifact-root=<fresh-external-directory> \
  --expected-sha="$(git rev-parse HEAD)"
```

This block records the last stack that passed; it is not a requirement to
retain those versions indefinitely. A newer CPE 26.03 / ROCm 7.0.2 attempt
currently encounters a hipfort/rocBLAS link incompatibility and is not
qualified. For a later current-stack run, update the module versions, keep the
same explicit recording and validation, and report the exact tested stack with
the result.

Use `--qos` or `--reservation` only when the facility requires it. The default
request is one node, one task, seven CPUs, one GPU, and 20 minutes. The script
uses `sbatch --wait`; queue wait is not part of the runtime limit. The artifact
root must be absent or empty and outside the repository.

The source worktree must be clean because `git archive HEAD` is the executable
source of record. The archive SHA-256, embedded archive commit, and verified
extracted-tree hash bind every build and test to that commit. Extraction occurs
inside the allocated job, after queue wait. CPU and GPU builds use separate
caller-selected build directories. Different directories may run concurrently;
do not run two top-level Make invocations or a clean concurrently in one
affected directory.

The qualification tools require Python 3.10 or newer. The current Frontier
procedure loads `cray-python/3.12.12`; the job checks the interpreter version
before invoking the runner and retains `python.version.txt` with the other run
files.

The explicit build selections are:

| Selection | CPU | GPU |
| --- | --- | --- |
| Build directory | `build/frontier-cpu` | `build/frontier-gpu` |
| Compiler/mode | `PE_ENV=CRAY`, `CMODE=OPT` | same |
| Parallel modes | MPI/OpenMP off | MPI/host OpenMP off |
| Accelerator | off | `GPU_MODE=ON`, `GPU_BACKEND=HIP` |
| Directives | off | `OPENMP_OL_MODE=ON`, OpenACC off |
| GPU linear algebra | none | `GPU_LAPACK_VER=ROCM` |
| EOS/solver | Starkiller, dense, LibSci | same |

The loaded accelerator target module makes `ftn -fopenmp` target the MI250X
`gfx90a` device. GPU runs set `OMP_TARGET_OFFLOAD=MANDATORY`, so host fallback
is a failure rather than a misleading pass.

The runner passes `BUILD_NAME=frontier-cpu` and `BUILD_NAME=frontier-gpu` to
the root production build through this directory's small qualification
Makefile. That file adds only the GPU linear-algebra probe and reuses the
production objects and build settings. The runner collects the selected
executables from the corresponding `build/<name>/bin/` directories in the
extracted source tree.

OpenMP-offload builds currently treat the shared `XASYNC` and `XWAIT` markers
as no-ops and execute synchronously. Mapping XNet's queue-oriented OpenACC
behavior to OpenMP tasks and dependencies is outside this qualification.

Cray's native Fortran preprocessor requires a fixed argument count for
function-like macros, while XNet's shared accelerator-directive layer uses
variadic macros. For Cray GPU builds, `make/crayftn_cpp.sh` therefore runs
the system `cpp -P -C -nostdinc` first and passes the resulting Fortran
source to the selected Cray compiler wrapper. The preprocessing wrapper is
selected after any MPI compiler override so an MPI-enabled Cray GPU build does
not silently return to native preprocessing. This procedure qualifies only the
explicit MPI-off configuration above. Comment preservation retains
Fortran `//` concatenation; disabling standard include directories avoids
injecting C system-header text. The qualification records both compiler and
preprocessor versions.

## Evidence and result review

On return, inspect and validate the manifest:

```bash
python3 test/qualification/frontier/frontier_qualification.py validate \
  <artifact-root>/qualification_manifest.json
```

The manifest records:

- full source SHA, clean-worktree assertion, source-archive and extracted-tree
  hashes, embedded archive commit, and pre-build verification status;
- exact loaded modules, compiler and preprocessor summaries, ROCm and hipfort
  versions, GPU model, and hashes of the raw environment reports;
- Slurm job ID and account-neutral resource parameters;
- every requested and resolved build variable, path-neutral compiler/flag
  summary, build status/runtime, executable size/hash, and dynamic-library
  link summary;
- hashes of every control, trajectory, abundance, network, and EOS input;
- direct status, timeout behavior, observed residual, zone list, numerical
  differences, policy fractions, and runtime output inventories; and
- a complete relative-path, size, and SHA-256 inventory of qualification
  output files except the self-referential manifest and the separately hashed
  source archive/build tree.

Review the observed scalar differences, selected-species fraction of allowed,
complete-vector `L1`/`L-infinity` values, and the report-only final
energy-generation-rate difference before accepting the policy. A policy
change requires a numerical explanation, a controlled perturbation that the
new limit still rejects, and a final rerun from the exact source commit.
Do not derive limits automatically from the current output.

Keep the validated manifest and full run directory outside the repository.
Record concise, path-neutral source, configuration, job, numerical, and review
results with the pull request or release record. Do not publish raw histories,
binaries, build products, absolute link paths, source archives, or private
allocation information.

## Focused local checks

These checks need no Frontier access:

```bash
python3 -m pytest -q \
  test/qualification/test_frontier_qualification.py \
  test/qualification/test_parallel_zones.py
```

They prove rejection of missing/duplicate/off-by-one zones, filename/result
misassociation, endpoint state leakage, material CPU/GPU perturbation, host
fallback, absent device evidence, handle/factor status failure, a large solve
residual, incomplete manifests, and queue/allocation/facility
misclassification.

## Failure categories and troubleshooting

Every job-side failure manifest names a category and phase:

- `source`: dirty checkout, wrong SHA, missing input, invalid policy, or archive
  mismatch;
- `environment`: missing/wrong module, compiler wrapper, ROCm, or hipfort;
- `submission`: malformed `sbatch` request or unavailable command;
- `queue`: cancellation, deadline, time-limit, or node failure before useful
  test evidence;
- `allocation`: invalid account/partition/QOS, failed `srun`, no MI250X, no HIP
  device, or OpenMP host fallback;
- `facility`: scheduler/controller communication failure;
- `build`: clean, compile, link, resolved-variable, or `ldd` failure;
- `test`: nonzero/timeout, device/handle/factor failure, residual failure,
  missing/wrong output association, or incomplete diagnostic; and
- `comparison`: a CPU/GPU endpoint exceeds the reviewed numerical limits.

Inspect `sbatch.stderr.txt` and `slurm.stderr.txt` first for submission,
queue, allocation, and facility failures. For builds, inspect the relevant
`build/*/build.stderr.txt` and resolved-variable evidence. For runtime checks,
use the per-run command, stdout, stderr, status, and complete output inventory.
Do not widen a comparison limit to hide a device, initialization, convergence,
or output-association defect.

## Required reruns

Repeat this manual qualification after changes to accelerator memory mapping or
lifetimes, `xnet_macros.fh` directive expansion, OpenMP target/compiler/module
selection, HIP/ROCm bindings or handles, batched GPU linear algebra, zone
batching or active masks, the BE/BDF integrators, GPU EOS/self-heating paths,
or the CPU/GPU endpoint policy. Also rerun before a release claims the
HIP/ROCm Frontier path is qualified. A new Frontier compiler, ROCm, hipfort,
GPU runtime, or materially changed facility module stack is qualification
drift and requires a new recorded run; it does not silently inherit an older
result.

OLCF platform guidance: [Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html).
