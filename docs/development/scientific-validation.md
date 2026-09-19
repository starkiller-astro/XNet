# Scientific validation guidance

> Read this document when changing physics, reaction rates, EOS behavior,
> network preprocessing, integrators, Jacobians, linear solvers, convergence
> criteria, tolerances, parallel numerical behavior, or performance-sensitive
> code.

This document begins with general evidence requirements. Add XNet-specific
checks as real development work establishes their definitions, inputs,
tolerances, and trusted comparisons.

## Separate the questions being answered

Evaluate each relevant category explicitly:

- **Build:** Does the requested configuration compile and link?
- **Runtime execution:** Does the program start, complete the intended path,
  and produce the expected diagnostics or outputs?
- **Software behavior:** Does the implementation meet the specified interface
  and error-handling requirements?
- **Numerical agreement:** Do computed quantities agree with a justified
  comparison within stated tolerances?
- **Scientific validity:** Does the calculation represent the intended
  physical model and regime?
- **Portability:** Which compilers, libraries, parallel modes, and hardware
  configurations produce acceptable behavior?
- **Performance:** Does the change preserve or improve the relevant runtime,
  memory use, transfer cost, or scaling behavior?

Report evidence by category. A build result answers the build question.
Runtime completion requires direct program status and expected runtime
evidence. A numerical comparison requires its own quantities, reference, and
tolerances.

## Evidence for defect fixes

For a defect fix, identify or add evidence that distinguishes the faulty
behavior from the corrected behavior. Prefer demonstrating that the check
fails before the fix when practical.

Useful demonstrations include:

- running the check on the pre-fix revision;
- temporarily reverting the critical implementation change;
- constructing a focused input that reaches the faulty path;
- introducing a controlled mutation during local investigation;
- showing that the assertion measures the scientific or interface
  requirement directly.

Choose the method that produces clear evidence with reasonable effort. Legacy
or exploratory repairs may require characterization before a focused
regression check is possible.

## Validation methods

Select methods from the scientific and numerical requirement. Useful options
for XNet work may include:

- small deterministic calculations;
- conservation or balance checks derived for the selected network and
  physical setup;
- analytic solutions or independently derived limiting cases;
- convergence studies across timestep or solver tolerances;
- comparisons among BE and BDF integration where both apply;
- comparisons among dense and sparse Jacobian or linear-solver paths;
- serial, MPI, threaded, batch, or accelerator comparisons;
- NSE limiting behavior and comparisons with independently established NSE
  results;
- network-preprocessing round trips and consistency checks;
- restart or decomposition consistency where the relevant workflow supports
  them;
- comparisons with trusted data or historical outputs whose provenance and
  applicability are known.

These are possible methods. The governing issue defines the applicable
physical invariant, comparison quantity, regime, and expected result. Record
uncertainty when current behavior has limited characterization.

## Characterization before refactoring

Capture current behavior before restructuring code whose behavior lacks clear
coverage. Characterization can include output quantities, conserved values,
timestep counts, iteration counts, selected rates, or focused intermediate
results.

Describe which captured behavior represents a scientific requirement, which
represents an interface expectation, and which simply records the current
implementation. This distinction keeps historical output from being treated
as evidence of scientific validity without justification.

## Tolerances and floating-point behavior

Choose tolerances from:

- the mathematical scale of the quantity;
- the physical importance of the difference;
- expected discretization or iteration error;
- floating-point ordering and reduction behavior;
- the precision and conditioning of the algorithm;
- variation across supported compilers, libraries, and hardware.

Document absolute and relative comparisons precisely, including treatment of
zero, small trace abundances, underflow, and inactive species when relevant.
Any tolerance change requires a numerical explanation and examples of the
results it accepts and rejects.

Use bitwise comparison when exact reproducibility is an established
requirement. Use justified numerical tolerances when valid execution orders or
implementations can change floating-point rounding.

## Recording a numerical change

Record the following in the issue or PR:

- the physical or numerical requirement;
- the affected code path and configuration;
- compiler, version, libraries, and build variables;
- network, input history, initial abundances, and runtime controls;
- quantities compared and their units;
- origin and provenance of the comparison result;
- absolute or relative comparison method and tolerance;
- observed differences, including trace and dominant quantities as relevant;
- changes to timestep selection, nonlinear iterations, or solver behavior;
- configurations and regimes that remain to be checked.

Preserve scripts or small inputs that make the evidence reproducible when they
fit the repository and issue scope.

## Performance evidence

Measure performance-sensitive changes with a representative problem. Record:

- hardware and relevant system libraries;
- compiler, flags, parallel mode, thread/process placement, and accelerator
  selection;
- network size, zone count, trajectory length, and other workload dimensions;
- warm-up and repetition method;
- runtime variability and comparison baseline;
- memory use, transfers, or scaling data when those motivate the change.

Compare fixed revisions and equivalent workloads. Explain the expected source
of improvement or regression. Keep numerical validation alongside performance
measurements so the benchmark exercises acceptable results.

## Scientific and architecture decisions

Developers should prepare evidence, identify assumptions, and expose
uncertainty. The maintainer decides whether a physical model, algorithm,
tolerance, reference result, or maintenance tradeoff is scientifically
appropriate.

Escalate questions involving physical interpretation, accepted regimes,
scientific provenance, or consequential algorithm selection. Record the
maintainer's decision in the governing issue, PR, or an appropriate design
document.
