"""Exact repeated-daughter closure check; no material imports or computations.

Assumed Markov Fock jumps: (-1,+2) at c*x*(y+1)*(y+2),
(+1,-2) at c*(x+1)*y*(y-1). This model is an additional kinetic
assumption, not the dynamics of an isolated finite cubic Hamiltonian.
Requires Python >=3.10 and SymPy. Run python -B C2_closure_check.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import sympy as sp

x, y = sp.symbols("x y", integer=True, nonnegative=True)
a, b, c = sp.symbols("a b c", positive=True)
z1, z2 = sp.symbols("delta_p delta_d", real=True)
rate_forward = c * x * (y + 1) * (y + 2)
rate_reverse = c * (x + 1) * y * (y - 1)


def generator(observable):
    forward = observable.subs({x: x - 1, y: y + 2}, simultaneous=True)
    reverse = observable.subs({x: x + 1, y: y - 2}, simultaneous=True)
    return sp.expand(
        rate_forward * (forward - observable)
        + rate_reverse * (reverse - observable)
    )


@lru_cache(maxsize=None)
def geometric_moment(power, mean):
    # Independently obtain raw moments from (r d/dr)^k [1/(1-r)].
    r = sp.Symbol("r")
    result = 1 / (1 - r)
    for _ in range(power):
        result = r * sp.diff(result, r)
    return sp.factor(((1 - r) * result).subs(r, mean / (1 + mean)))


def geometric_expectation(polynomial):
    total = 0
    for (px, py), coefficient in sp.Poly(polynomial, x, y).terms():
        total += coefficient * geometric_moment(px, a) * geometric_moment(py, b)
    return sp.factor(total)


def exact_zero(expression):
    if isinstance(expression, sp.MatrixBase):
        return all(sp.simplify(entry) == 0 for entry in expression)
    return sp.simplify(expression) == 0


def matrix_strings(matrix):
    return [[str(sp.simplify(entry)) for entry in row] for row in matrix.tolist()]


def symbolic_checks():
    flux = geometric_expectation(rate_forward - rate_reverse)
    expected_flux = 2 * c * (a * (1 + 2 * b) - b**2)
    assert exact_zero(flux - expected_flux)
    incidence = sp.Matrix([-1, 2])
    energy = sp.Matrix([2, 1])
    flow = incidence * flux
    assert exact_zero(generator(2 * x + y))
    assert exact_zero(energy.T * flow)
    jacobian = flow.jacobian([a, b])
    equilibrium_parent = b**2 / (1 + 2 * b)
    loss = sp.simplify(-jacobian.subs(a, equilibrium_parent))
    susceptibility = sp.diag(
        equilibrium_parent * (1 + equilibrium_parent), b * (1 + b)
    )
    bose_factor = equilibrium_parent * (1 + b)**2
    r = -incidence
    kappa = 2 * c
    rank_one_loss = kappa * bose_factor * r * r.T * susceptibility.inv()
    assert exact_zero(loss - rank_one_loss)
    assert exact_zero(energy.T * loss)
    assert exact_zero(loss * susceptibility * energy)
    entropy_hessian = susceptibility.inv()
    entropy_matrix = sp.simplify(entropy_hessian * loss)
    assert exact_zero(entropy_matrix - entropy_matrix.T)
    displacement = sp.Matrix([z1, z2])
    entropy_quadratic = sp.factor((displacement.T * entropy_matrix * displacement)[0])
    square_form = kappa * bose_factor * (r.T * entropy_hessian * displacement)[0]**2
    assert exact_zero(entropy_quadratic - square_form)
    positive_eigenvalue = sp.factor(sp.trace(loss))
    assert exact_zero(sp.det(loss))
    assert exact_zero(
        positive_eigenvalue - 2 * c * (1 + 8 * b + 8 * b**2) / (1 + 2 * b)
    )

    # Independently derived Jacobian compared with D's stated source comparator.
    # No source import/re-audit is claimed for the comparator in this script.
    d_comparator = sp.Matrix([
        kappa * (1 + 2 * b),
        (2 * kappa * (b - a)).subs(a, equilibrium_parent),
    ])
    d_ratios = [
        sp.factor(loss[index, index] / d_comparator[index]) for index in range(2)
    ]
    assert d_ratios == [1, 2]

    # One exact positive equilibrium gives a symmetric susceptibility-normalized matrix.
    sample_subs = {b: sp.Integer(1), c: sp.Integer(1)}
    sample_loss = loss.subs(sample_subs)
    sample_w = susceptibility.subs(sample_subs)
    root_w = sp.diag(*[sp.sqrt(sample_w[j, j]) for j in range(2)])
    sample_symmetric = sp.simplify(root_w.inv() * sample_loss * root_w)
    sample_energy_kernel = root_w * energy
    assert exact_zero(sample_symmetric - sample_symmetric.T)
    assert exact_zero(sample_symmetric * sample_energy_kernel)
    assert sample_symmetric.eigenvals() == {sp.Integer(0): 1, sp.Rational(34, 3): 1}

    moment_derivatives = {
        "parent_factorial2": geometric_expectation(generator(x * (x - 1))),
        "daughter_factorial2": geometric_expectation(generator(y * (y - 1))),
        "parent_daughter": geometric_expectation(generator(x * y)),
    }
    parent_defect = sp.factor(moment_derivatives["parent_factorial2"] - 4 * a * flow[0])
    daughter_defect = sp.factor(moment_derivatives["daughter_factorial2"] - 4 * b * flow[1])
    covariance_derivative = sp.factor(
        moment_derivatives["parent_daughter"] - b * flow[0] - a * flow[1]
    )
    assert exact_zero(parent_defect)
    assert exact_zero(daughter_defect - 2 * (2 * b + 1) * flux)
    assert exact_zero(covariance_derivative - 2 * (a - b) * flux)
    reciprocal_probability_ratio = b**2 * (1 + a) / (a * (1 + b)**2)
    assert exact_zero(reciprocal_probability_ratio.subs(a, equilibrium_parent) - 1)
    # The edge rates equal their destination's reverse edge rates exactly.
    assert exact_zero(
        rate_forward - rate_reverse.subs({x: x - 1, y: y + 2}, simultaneous=True)
    )

    sample_nonequilibrium = {a: sp.Rational(1, 3), b: sp.Rational(1, 2), c: 1}
    return {
        "raw_geometric_moments": {str(k): str(geometric_moment(k, a)) for k in range(4)},
        "flux": str(flux),
        "flow_jacobian": matrix_strings(jacobian),
        "equilibrium_parent_mean": str(equilibrium_parent),
        "equilibrium_loss": matrix_strings(loss),
        "susceptibility": matrix_strings(susceptibility),
        "equilibrium_bose_factor": str(sp.factor(bose_factor)),
        "entropy_quadratic": str(entropy_quadratic),
        "nonzero_equilibrium_eigenvalue": str(positive_eigenvalue),
        "exact_energy_left_and_right_kernel_checks": True,
        "D_report_comparator": {
            "formula": matrix_strings(d_comparator),
            "event_diagonal_over_comparator": [str(value) for value in d_ratios],
            "scope": "D report formula only; source not independently imported or audited",
        },
        "sample_equilibrium": {
            "means": ["1/3", "1"], "c": "1",
            "loss": matrix_strings(sample_loss),
            "symmetric_loss": matrix_strings(sample_symmetric),
            "symmetric_energy_kernel": matrix_strings(sample_energy_kernel),
            "eigenvalues": ["0", "34/3"],
        },
        "geometric_family_tangency": {
            "moment_derivatives": {key: str(value) for key, value in moment_derivatives.items()},
            "parent_factorial2_defect": str(parent_defect),
            "daughter_factorial2_defect": str(daughter_defect),
            "covariance_derivative": str(covariance_derivative),
            "destination_over_initial_probability": str(reciprocal_probability_ratio),
            "equilibrium_edge_detailed_balance_exact": True,
            "sample_means": ["1/3", "1/2"],
            "sample_c": "1",
            "sample_flux": str(flux.subs(sample_nonequilibrium)),
            "sample_daughter_factorial2_derivative": str(
                moment_derivatives["daughter_factorial2"].subs(sample_nonequilibrium)),
            "sample_geometric_tangent_derivative": str((4 * b * flow[1]).subs(sample_nonequilibrium)),
            "sample_daughter_tangency_defect": str(daughter_defect.subs(sample_nonequilibrium)),
            "sample_covariance_derivative": str(covariance_derivative.subs(sample_nonequilibrium)),
        },
        "all_symbolic_identities_exact": True,
    }


def integer_rates(px, dy):
    return px * (dy + 1) * (dy + 2), (px + 1) * dy * (dy - 1)


def generator_value(observable, px, dy):
    forward, reverse = integer_rates(px, dy)
    value = 0
    current = observable(px, dy)
    if forward:
        value += forward * (observable(px - 1, dy + 2) - current)
    if reverse:
        value += reverse * (observable(px + 1, dy - 2) - current)
    return value


def finite_same_mean_checks():
    distributions = {
        "daughter_fixed_one": {(1, 1): Fraction(1)},
        "daughter_zero_or_two": {(1, 0): Fraction(1, 2), (1, 2): Fraction(1, 2)},
    }
    result = {}
    for label, distribution in distributions.items():
        assert sum(distribution.values()) == 1
        stats = {key: Fraction(0) for key in [
            "mean_parent", "mean_daughter", "daughter_factorial2", "forward", "reverse", "flux"
        ]}
        for (px, dy), probability in distribution.items():
            forward, reverse = integer_rates(px, dy)
            contributions = {
                "mean_parent": px, "mean_daughter": dy,
                "daughter_factorial2": dy * (dy - 1),
                "forward": forward, "reverse": reverse, "flux": forward - reverse,
            }
            for key, value in contributions.items():
                stats[key] += probability * value
            assert generator_value(lambda u, v: 2 * u + v, px, dy) == 0
        assert stats["mean_parent"] == stats["mean_daughter"] == 1
        result[label] = {key: str(value) for key, value in stats.items()}
    assert result["daughter_fixed_one"]["flux"] == "6"
    assert result["daughter_zero_or_two"]["flux"] == "5"
    result["geometric_daughter_mean_one_reference"] = {
        "mean_parent": "1", "mean_daughter": "1", "daughter_factorial2": "2", "flux": "4",
        "note": "independent parent with mean one; parent need not be geometric for this mean flux",
    }
    return result


def rational_partial_sum_checks():
    # Independent direct state enumeration, without clipping destination states.
    parent_ratio, daughter_ratio = Fraction(1, 4), Fraction(1, 3)
    targets = {
        "flux": Fraction(5, 6),
        "parent_derivative": Fraction(-5, 6),
        "daughter_derivative": Fraction(5, 3),
        "parent_factorial2_derivative": Fraction(-10, 9),
        "daughter_factorial2_derivative": Fraction(20, 3),
        "parent_daughter_derivative": Fraction(-5, 36),
    }
    observables = {
        "parent_derivative": lambda u, v: u,
        "daughter_derivative": lambda u, v: v,
        "parent_factorial2_derivative": lambda u, v: u * (u - 1),
        "daughter_factorial2_derivative": lambda u, v: v * (v - 1),
        "parent_daughter_derivative": lambda u, v: u * v,
    }
    results = []
    for max_parent, max_daughter in [(8, 12), (16, 24), (32, 48)]:
        sums = {key: Fraction(0) for key in targets}
        mass = Fraction(0)
        for px in range(max_parent + 1):
            p_parent = (1 - parent_ratio) * parent_ratio**px
            for dy in range(max_daughter + 1):
                probability = p_parent * (1 - daughter_ratio) * daughter_ratio**dy
                mass += probability
                forward, reverse = integer_rates(px, dy)
                sums["flux"] += probability * (forward - reverse)
                for key, observable in observables.items():
                    sums[key] += probability * generator_value(observable, px, dy)
        omitted = 1 - mass
        assert omitted == 1 - (
            (1 - parent_ratio**(max_parent + 1)) * (1 - daughter_ratio**(max_daughter + 1))
        )
        errors = {key: abs(sums[key] - targets[key]) for key in targets}
        results.append({
            "max_parent": max_parent, "max_daughter": max_daughter,
            "omitted_probability_exact": str(omitted),
            "omitted_probability": float(omitted),
            "partial_generator_expectations": {key: str(value) for key, value in sums.items()},
            "absolute_errors": {key: float(value) for key, value in errors.items()},
            "maximum_absolute_error_exact": str(max(errors.values())),
            "maximum_absolute_error": float(max(errors.values())),
        })
    assert max(errors.values()) < Fraction(1, 10**14)
    return {
        "ratios_exact": ["1/4", "1/3"], "means_exact": ["1/3", "1/2"], "c": "1",
        "targets_exact": {key: str(value) for key, value in targets.items()},
        "convergence": results,
        "boundary_treatment": "unnormalized partial initial distribution; destination states not clipped",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpreter": sys.executable, "python": platform.python_version(), "sympy": sp.__version__,
        "symbolic": symbolic_checks(),
        "mean_only_counterexamples": finite_same_mean_checks(),
        "rational_partial_sums": rational_partial_sum_checks(),
        "all_assertions_passed": True,
        "scope": "assumed kinetic Fock jump model; no material, source, or finite Hamiltonian dynamics claim",
    }
    destination = args.output_dir / "C2_closure_check.json"
    destination.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(destination), "script_sha256": results["script_sha256"],
        "all_assertions_passed": True,
        "D_comparator_ratios": results["symbolic"]["D_report_comparator"]["event_diagonal_over_comparator"],
        "daughter_tangency_defect": results["symbolic"]["geometric_family_tangency"]["daughter_factorial2_defect"],
        "sample_daughter_tangency_defect": results["symbolic"]["geometric_family_tangency"]["sample_daughter_tangency_defect"],
        "final_partial_sum_max_absolute_error": results["rational_partial_sums"]["convergence"][-1]["maximum_absolute_error"],
    }, indent=2))


if __name__ == "__main__":
    main()

