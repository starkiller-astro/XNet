"""High-precision ideal-NSE solver derived from published equilibrium equations.

The implementation uses dimensionless neutron and proton kinetic chemical
potentials and solves normalized log-sum-exp constraints.  It does not import,
call, or translate XNet's NSE routines.  See Seitenzahl et al. (2009), eqs.
(2)-(10), and Lippuner & Roberts (2017), Appendix B.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from typing import Any, Iterable, Sequence


class ReferenceSolveError(RuntimeError):
    """Raised when a high-precision equilibrium solve does not converge."""


def decimal_from_binary64(entry: dict[str, str]) -> Decimal:
    return Decimal.from_float(float.fromhex(entry["hex"]))


@dataclass(frozen=True)
class State:
    state_id: str
    rho: Decimal
    t9: Decimal
    ye: Decimal

    @classmethod
    def from_strings(cls, state_id: str, rho: str, t9: str, ye: str) -> "State":
        # XNet receives binary64 values.  Lift those exact values into Decimal.
        return cls(
            state_id,
            Decimal.from_float(float(rho)),
            Decimal.from_float(float(t9)),
            Decimal.from_float(float(ye)),
        )


@dataclass(frozen=True)
class Species:
    name: str
    a: int
    z: int
    n: int
    log_coefficient: Decimal

    @property
    def charge_per_baryon(self) -> Decimal:
        return Decimal(self.z) / Decimal(self.a)


@dataclass(frozen=True)
class Evaluation:
    residual: tuple[Decimal, Decimal]
    jacobian: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]
    log_mass_sum: Decimal
    normalized_composition: tuple[Decimal, ...]
    log_composition: tuple[Decimal, ...]


@dataclass(frozen=True)
class Solution:
    state: State
    eta_n: Decimal
    eta_p: Decimal
    composition: tuple[Decimal, ...]
    mass_residual: Decimal
    charge_residual: Decimal
    xnet_charge_residual: Decimal
    iterations: int
    route: str


def _partition_factor(
    manifest: dict[str, Any], species: dict[str, Any], t9: Decimal
) -> Decimal:
    grid = [decimal_from_binary64(item) for item in manifest["temperature_grid_gk"]]
    values = [decimal_from_binary64(item) for item in species["partition_factors"]]
    upper = next((index for index, value in enumerate(grid) if value >= t9), len(grid))
    if upper == 0:
        return values[0]
    if upper == len(grid):
        return values[-1]
    fraction = (t9 - grid[upper - 1]) / (grid[upper] - grid[upper - 1])
    return (
        fraction * values[upper].ln()
        + (Decimal(1) - fraction) * values[upper - 1].ln()
    ).exp()


def build_species(
    manifest: dict[str, Any], state: State, precision: int
) -> tuple[Species, ...]:
    with localcontext() as context:
        context.prec = precision + 20
        constants = manifest["constants"]
        pi = decimal_from_binary64(constants["pi"])
        hbar = decimal_from_binary64(constants["hbar"])
        bok = decimal_from_binary64(constants["bok"])
        epmev = decimal_from_binary64(constants["epmev"])
        bkt = state.t9 * bok * epmev
        quantum = bkt / (Decimal(2) * pi * hbar * hbar * epmev * epmev)
        common = Decimal("1.5") * quantum.ln() - state.rho.ln()
        result = []
        for item in manifest["species"]:
            mass = decimal_from_binary64(item["translational_mass_g"])
            degeneracy = decimal_from_binary64(item["ground_state_degeneracy"])
            partition = _partition_factor(manifest, item, state.t9)
            binding = decimal_from_binary64(item["binding_energy_mev"])
            log_coefficient = (
                common
                + degeneracy.ln()
                + partition.ln()
                + Decimal("2.5") * mass.ln()
                + binding * epmev / bkt
            )
            result.append(
                Species(
                    name=item["name"],
                    a=item["a"],
                    z=item["z"],
                    n=item["n"],
                    log_coefficient=+log_coefficient,
                )
            )
        return tuple(result)


def evaluate(
    species: Sequence[Species],
    eta_n: Decimal,
    eta_p: Decimal,
    ye: Decimal,
    *,
    mass_sum_target: Decimal = Decimal(1),
    charge_residual_target: Decimal = Decimal(0),
) -> Evaluation:
    logs = tuple(
        item.log_coefficient + Decimal(item.n) * eta_n + Decimal(item.z) * eta_p
        for item in species
    )
    maximum = max(logs)
    scaled = tuple((value - maximum).exp() for value in logs)
    scaled_sum = sum(scaled, Decimal(0))
    normalized = tuple(value / scaled_sum for value in scaled)
    log_mass_sum = maximum + scaled_sum.ln()

    mean_n = sum(
        (Decimal(item.n) * value for item, value in zip(species, normalized, strict=True)),
        Decimal(0),
    )
    mean_z = sum(
        (Decimal(item.z) * value for item, value in zip(species, normalized, strict=True)),
        Decimal(0),
    )
    mean_q = sum(
        (item.charge_per_baryon * value for item, value in zip(species, normalized, strict=True)),
        Decimal(0),
    )
    mean_qn = sum(
        (
            item.charge_per_baryon * Decimal(item.n) * value
            for item, value in zip(species, normalized, strict=True)
        ),
        Decimal(0),
    )
    mean_qz = sum(
        (
            item.charge_per_baryon * Decimal(item.z) * value
            for item, value in zip(species, normalized, strict=True)
        ),
        Decimal(0),
    )
    return Evaluation(
        residual=(
            log_mass_sum - mass_sum_target.ln(),
            mean_q - ye - charge_residual_target / mass_sum_target,
        ),
        jacobian=(
            (mean_n, mean_z),
            (mean_qn - mean_q * mean_n, mean_qz - mean_q * mean_z),
        ),
        log_mass_sum=log_mass_sum,
        normalized_composition=normalized,
        log_composition=logs,
    )


def _solve_step(
    jacobian: tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]],
    residual: tuple[Decimal, Decimal],
) -> tuple[Decimal, Decimal]:
    a, b = jacobian[0]
    c, d = jacobian[1]
    determinant = a * d - b * c
    if determinant == 0:
        raise ReferenceSolveError("singular two-variable reference Jacobian")
    f, g = residual
    return ((-d * f + b * g) / determinant, (c * f - a * g) / determinant)


def initial_guess(species: Sequence[Species], ye: Decimal) -> tuple[Decimal, Decimal]:
    neutron = next(item for item in species if item.n == 1 and item.z == 0)
    proton = next(item for item in species if item.n == 0 and item.z == 1)
    return (
        (Decimal(1) - ye).ln() - neutron.log_coefficient,
        ye.ln() - proton.log_coefficient,
    )


def _numeric_jacobian(
    species: Sequence[Species],
    eta_n: Decimal,
    eta_p: Decimal,
    ye: Decimal,
    step: Decimal,
    mass_sum_target: Decimal,
    charge_residual_target: Decimal,
) -> tuple[tuple[Decimal, Decimal], tuple[Decimal, Decimal]]:
    evaluation_args = {
        "mass_sum_target": mass_sum_target,
        "charge_residual_target": charge_residual_target,
    }
    n_plus = evaluate(species, eta_n + step, eta_p, ye, **evaluation_args).residual
    n_minus = evaluate(species, eta_n - step, eta_p, ye, **evaluation_args).residual
    p_plus = evaluate(species, eta_n, eta_p + step, ye, **evaluation_args).residual
    p_minus = evaluate(species, eta_n, eta_p - step, ye, **evaluation_args).residual
    denominator = Decimal(2) * step
    return (
        (
            (n_plus[0] - n_minus[0]) / denominator,
            (p_plus[0] - p_minus[0]) / denominator,
        ),
        (
            (n_plus[1] - n_minus[1]) / denominator,
            (p_plus[1] - p_minus[1]) / denominator,
        ),
    )


def solve(
    manifest: dict[str, Any],
    state: State,
    precision: int,
    *,
    route: str = "analytic-newton",
    offset: tuple[str, str] = ("0", "0"),
    constraint_target: tuple[str, str] = ("0", "0"),
    max_iterations: int = 240,
) -> Solution:
    if route not in {"analytic-newton", "numeric-newton"}:
        raise ValueError(f"unknown reference route {route}")
    with localcontext() as context:
        context.prec = precision + 24
        species = build_species(manifest, state, precision)
        eta_n, eta_p = initial_guess(species, state.ye)
        eta_n += Decimal(offset[0])
        eta_p += Decimal(offset[1])
        mass_sum_target = Decimal(1) + Decimal(constraint_target[0])
        charge_residual_target = Decimal(constraint_target[1])
        if mass_sum_target <= 0:
            raise ValueError("shifted mass-sum target must be positive")
        target = Decimal(10) ** -min(precision // 2, 70)
        finite_difference_step = Decimal(10) ** -max(18, precision // 3)
        for iteration in range(max_iterations + 1):
            current = evaluate(
                species,
                eta_n,
                eta_p,
                state.ye,
                mass_sum_target=mass_sum_target,
                charge_residual_target=charge_residual_target,
            )
            norm = max(abs(current.residual[0]), abs(current.residual[1]))
            if norm <= target:
                raw = tuple(value.exp() for value in current.log_composition)
                mass = sum(raw, Decimal(0))
                charge = sum(
                    (
                        item.charge_per_baryon * value
                        for item, value in zip(species, raw, strict=True)
                    ),
                    Decimal(0),
                )
                return Solution(
                    state=state,
                    eta_n=+eta_n,
                    eta_p=+eta_p,
                    composition=tuple(+value for value in raw),
                    mass_residual=+(mass - Decimal(1)),
                    charge_residual=+(charge - state.ye),
                    xnet_charge_residual=+(
                        charge - state.ye * mass
                    ),
                    iterations=iteration,
                    route=route,
                )
            jacobian = (
                current.jacobian
                if route == "analytic-newton"
                else _numeric_jacobian(
                    species,
                    eta_n,
                    eta_p,
                    state.ye,
                    finite_difference_step,
                    mass_sum_target,
                    charge_residual_target,
                )
            )
            delta_n, delta_p = _solve_step(jacobian, current.residual)
            damping = Decimal(1)
            for _ in range(100):
                trial_n = eta_n + damping * delta_n
                trial_p = eta_p + damping * delta_p
                trial = evaluate(
                    species,
                    trial_n,
                    trial_p,
                    state.ye,
                    mass_sum_target=mass_sum_target,
                    charge_residual_target=charge_residual_target,
                )
                trial_norm = max(abs(trial.residual[0]), abs(trial.residual[1]))
                if trial_norm < norm:
                    eta_n, eta_p = trial_n, trial_p
                    break
                damping /= 2
            else:
                raise ReferenceSolveError(
                    f"{route} line search failed for {state.state_id} at iteration {iteration}"
                )
        raise ReferenceSolveError(
            f"{route} exceeded {max_iterations} iterations for {state.state_id}"
        )


def composition_norms(left: Iterable[Decimal], right: Iterable[Decimal]) -> tuple[Decimal, Decimal]:
    differences = tuple(abs(a - b) for a, b in zip(left, right, strict=True))
    return sum(differences, Decimal(0)), max(differences, default=Decimal(0))


def dominant_species(
    manifest: dict[str, Any], solution: Solution, threshold: Decimal = Decimal("1e-3")
) -> list[tuple[str, Decimal]]:
    return [
        (item["name"], value)
        for item, value in zip(manifest["species"], solution.composition, strict=True)
        if value >= threshold
    ]
