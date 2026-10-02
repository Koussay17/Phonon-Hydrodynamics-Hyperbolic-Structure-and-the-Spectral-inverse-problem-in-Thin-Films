"""Standalone finite entropy weak-form checks, not production code.

Python >=3.10 and NumPy only. Positive volume and event quadratures,
common entropy-variable reconstruction, complex-step derivative checks,
exact two-node references, and a small periodic regular-surface example.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np


if not __debug__:
    raise SystemExit("C2_weak_form: checks use assert; refusing to run under python -O")


def bose(xi):
    return 1.0 / np.expm1(xi)


def log_bose(xi):
    """log n = -xi - log(1-exp(-xi)) for Re(xi)>0, without overflow."""
    return -xi - np.log(-np.expm1(-xi))


def log_one_plus_bose(xi):
    """log(1+n) = -log(1-exp(-xi)) for Re(xi)>0."""
    return -np.log(-np.expm1(-xi))


def exprel_stable(z):
    """(e^z-1)/z. Series below |z|<1e-2 (2 Oct 2026; was 1e-5 with six terms).

    The wider series removes the complex-step derivative cancellation of
    expm1(z)/z, which lost up to five digits for 1e-5<=|z|<=1e-2. Truncation
    after z^10/11! is below 3e-28 relative at |z|=1e-2.
    """
    value = np.asarray(z)
    result = np.empty_like(value, dtype=np.result_type(value, float))
    small = np.abs(value) < 1e-2
    v = value[small]
    coefficients = [1/math.factorial(k+1) for k in range(11)]
    series = np.zeros_like(v)
    for c in reversed(coefficients):
        series = series*v + c
    result[small] = series
    result[~small] = np.expm1(value[~small]) / value[~small]
    return result


@dataclass
class EntropyWeakForm:
    Phi: np.ndarray
    masses: np.ndarray
    parent: np.ndarray
    daughter_a: np.ndarray
    daughter_b: np.ndarray
    weights: np.ndarray

    def __post_init__(self):
        for name in ("Phi", "masses", "parent", "daughter_a", "daughter_b", "weights"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        assert np.all(self.masses > 0) and np.all(self.weights > 0)
        assert np.linalg.matrix_rank(self.Phi) == self.Phi.shape[1]
        assert self.parent.shape == self.daughter_a.shape == self.daughter_b.shape
        self.b = self.parent - self.daughter_a - self.daughter_b

    def minimum_xi(self, alpha):
        return float(min(np.min(evaluation @ alpha) for evaluation in
                         (self.Phi, self.parent, self.daughter_a, self.daughter_b)))

    def moments(self, alpha):
        return self.Phi.T @ (self.masses * bose(self.Phi @ alpha))

    def capacity(self, alpha):
        n = bose(self.Phi @ alpha)
        return self.Phi.T @ ((self.masses*n*(1+n))[:, None] * self.Phi)

    def reaction_data(self, alpha):
        # Log-space evaluation (2 Oct 2026 fix). The frozen candidate formed
        # B=(1+n_p)n_a n_b and B*exprel(-d); for widely separated entropy
        # variables B underflows to 0 while exprel(-d) overflows (0*inf=NaN).
        # Here d=log B-log A=b@alpha, F=A-B and Lambda=(A-B)/(log A-log B).
        # The branch on Re(d) keeps each expression analytic for complex steps:
        # d>=0: F=B*expm1(-d), Lambda=B*exprel(-d); d<0: F=-A*expm1(d),
        # Lambda=A*exprel(d). Both branches are equal because A=B*exp(-d).
        affinity = self.b @ alpha
        with np.errstate(over="ignore", invalid="ignore"):
            affinity, flux, mobility = self._log_space_reaction(alpha, affinity)
        # errstate only silences the discarded np.where branch; the selected
        # outputs must still be finite (2 Oct 2026 numerical review).
        if not (np.all(np.isfinite(flux)) and np.all(np.isfinite(mobility))):
            raise FloatingPointError("non-finite selected flux or mobility; check xi>0 and range")
        return affinity, flux, mobility

    def _log_space_reaction(self, alpha, affinity):
        log_a = (log_bose(self.parent @ alpha) + log_one_plus_bose(self.daughter_a @ alpha)
                 + log_one_plus_bose(self.daughter_b @ alpha))
        log_b = (log_one_plus_bose(self.parent @ alpha) + log_bose(self.daughter_a @ alpha)
                 + log_bose(self.daughter_b @ alpha))
        reverse_wins = np.real(affinity) >= 0
        flux = np.where(reverse_wins, np.exp(log_b)*np.expm1(-affinity),
                        -np.exp(log_a)*np.expm1(affinity))
        mobility = np.where(reverse_wins, np.exp(log_b)*exprel_stable(-affinity),
                            np.exp(log_a)*exprel_stable(affinity))
        return affinity, flux, mobility

    def rhs(self, alpha, direct_products=False):
        if direct_products:
            np_ = bose(self.parent @ alpha)
            na = bose(self.daughter_a @ alpha)
            nb = bose(self.daughter_b @ alpha)
            flux = np_*(1+na)*(1+nb) - (1+np_)*na*nb
        else:
            _, flux, _ = self.reaction_data(alpha)
        return -self.b.T @ (self.weights*flux)

    def collision_matrix(self, alpha):
        _, _, mobility = self.reaction_data(alpha)
        return self.b.T @ ((self.weights*mobility)[:, None]*self.b)

    def coefficient_rhs(self, alpha, direct_products=False):
        return -np.linalg.solve(self.capacity(alpha), self.rhs(alpha, direct_products))

    def entropy(self, alpha):
        xi = self.Phi @ alpha
        n = bose(xi)
        return np.dot(self.masses, xi*n - np.log(-np.expm1(-xi)))

    def energy(self, alpha, energy_coefficients):
        return energy_coefficients @ self.moments(alpha)


def complex_jacobian(function, alpha, step=1e-28):
    columns = []
    for j in range(len(alpha)):
        z = np.asarray(alpha, dtype=complex).copy()
        z[j] += 1j*step
        columns.append(np.atleast_1d(np.imag(function(z))/step))
    return np.column_stack(columns)


def directional_derivative(function, alpha, velocity, step=1e-28):
    return float(np.imag(function(np.asarray(alpha, complex) + 1j*step*velocity))/step)


def relative_error(actual, expected):
    return float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected), 1e-30))


def structural_checks(model, energy, beta, perturbation):
    alpha_eq = beta*energy
    alpha = alpha_eq + perturbation
    assert model.minimum_xi(alpha_eq) > 0 and model.minimum_xi(alpha) > 0
    M = model.capacity(alpha_eq)
    K = model.collision_matrix(alpha_eq)
    dU = complex_jacobian(model.moments, alpha_eq)
    dG = complex_jacobian(lambda z: model.rhs(z, direct_products=True), alpha_eq)
    d_alpha_rhs = complex_jacobian(lambda z: model.coefficient_rhs(z, True), alpha_eq)
    mvalues, mvecs = np.linalg.eigh(M)
    invroot = (mvecs*(1/np.sqrt(mvalues))) @ mvecs.T
    root = (mvecs*np.sqrt(mvalues)) @ mvecs.T
    C = invroot @ K @ invroot
    ceig = np.linalg.eigvalsh(C)
    normalized_energy_kernel = relative_error(C @ (root @ energy), np.zeros(len(energy)))
    # Use operator/vector scales for the zero-kernel residual.
    normalized_energy_kernel = float(
        np.linalg.norm(C @ (root @ energy))
        / max(np.linalg.norm(C)*np.linalg.norm(root @ energy), 1e-30)
    )
    G = model.rhs(alpha)
    Koff = model.collision_matrix(alpha)
    velocity = model.coefficient_rhs(alpha)
    affinity, flux, mobility = model.reaction_data(alpha)
    event_entropy = float(np.dot(model.weights*mobility, affinity**2))
    event_energy = float(-np.dot(model.weights*flux, model.b @ energy))
    energy_square = float(np.dot(model.weights*mobility, (model.b @ energy)**2))
    entropy_derivative = directional_derivative(model.entropy, alpha, velocity)
    energy_derivative = directional_derivative(lambda z: model.energy(z, energy), alpha, velocity)
    off_jac = complex_jacobian(lambda z: model.rhs(z, True), alpha)
    off_M_derivative = complex_jacobian(model.moments, alpha)
    finite_difference_errors = []
    for step in [1e-3, 1e-4, 1e-5, 1e-6]:
        columns = []
        for j in range(len(alpha_eq)):
            direction = np.zeros(len(alpha_eq))
            direction[j] = step
            columns.append((model.moments(alpha_eq+direction)-model.moments(alpha_eq-direction))/(2*step))
        finite_difference_errors.append({
            "step": step, "relative_error": relative_error(np.column_stack(columns), -M)
        })
    result = {
        "number_coefficients": len(energy), "volume_nodes": len(model.masses),
        "events": len(model.weights), "volume_mass": float(sum(model.masses)),
        "reaction_mass": float(sum(model.weights)),
        "minimum_xi_equilibrium": model.minimum_xi(alpha_eq),
        "minimum_xi_perturbed": model.minimum_xi(alpha),
        "capacity_eigenvalues": mvalues.tolist(), "capacity_condition_number": float(mvalues[-1]/mvalues[0]),
        "normalized_collision_eigenvalues": ceig.tolist(),
        "collision_rank": int(np.linalg.matrix_rank(K, tol=max(np.linalg.norm(K),1e-30)*1e-10)),
        "relative_capacity_derivative_error_equilibrium": relative_error(dU, -M),
        "relative_capacity_derivative_error_perturbed": relative_error(off_M_derivative, -model.capacity(alpha)),
        "relative_collision_jacobian_error_equilibrium": relative_error(dG, K),
        "relative_coefficient_jacobian_error_equilibrium": relative_error(d_alpha_rhs, -np.linalg.solve(M,K)),
        "equilibrium_rhs_norm": float(np.linalg.norm(model.rhs(alpha_eq))),
        "normalized_energy_kernel_residual": normalized_energy_kernel,
        "maximum_reaction_energy_mismatch": float(max(abs(model.b @ energy))),
        "off_equilibrium": {
            "rhs_vs_K_alpha_relative_error": relative_error(G, Koff @ alpha),
            "direct_vs_stable_rhs_relative_error": relative_error(model.rhs(alpha,True),G),
            "jacobian_vs_K_relative_difference": relative_error(off_jac,Koff),
            "entropy_production_event": event_entropy,
            "entropy_directional_derivative": entropy_derivative,
            "entropy_relative_identity_error": abs(entropy_derivative-event_entropy)/event_entropy,
            "energy_drift_event": event_energy,
            "energy_directional_derivative": energy_derivative,
            "energy_drift_from_moments": float(energy @ G),
            "energy_square": energy_square,
            "energy_cauchy_schwarz_bound": math.sqrt(max(0,energy_square*event_entropy)),
        },
        "capacity_centered_difference_checks": finite_difference_errors,
    }
    assert result["relative_capacity_derivative_error_equilibrium"] < 1e-12
    assert result["relative_capacity_derivative_error_perturbed"] < 1e-12
    assert result["relative_collision_jacobian_error_equilibrium"] < 1e-11
    assert result["relative_coefficient_jacobian_error_equilibrium"] < 1e-11
    assert normalized_energy_kernel < 1e-12
    assert min(ceig) >= -1e-12*max(1,max(ceig))
    assert event_entropy > 0
    assert result["off_equilibrium"]["entropy_relative_identity_error"] < 1e-11
    assert abs(energy_derivative-event_energy) < 1e-12
    assert result["off_equilibrium"]["rhs_vs_K_alpha_relative_error"] < 1e-12
    return result


def two_node_model():
    return EntropyWeakForm(np.eye(2), np.ones(2), [[0.5,0.5]], [[1,0]], [[1,0]], [1])


def occupation_interpolation_rhs(model, alpha):
    nodes = bose(alpha)
    np_ = model.parent @ nodes
    na = model.daughter_a @ nodes
    nb = model.daughter_b @ nodes
    flux = np_*(1+na)*(1+nb)-(1+np_)*na*nb
    return -model.b.T @ (model.weights*flux), flux


def two_node_checks():
    model = two_node_model()
    energy = np.array([1.,3.])
    beta = math.log(2)
    eq = beta*energy
    desired_off_eq = np.array([0.8,1.7])
    result = structural_checks(model,energy,beta,desired_off_eq-eq)
    exact_M = np.diag([2.,8/49])
    exact_K = np.array([[3.,-1.],[-1.,1/3]])
    assert relative_error(model.capacity(eq),exact_M) < 1e-14
    assert relative_error(model.collision_matrix(eq),exact_K) < 1e-14
    assert abs(max(result["normalized_collision_eigenvalues"])-85/24) < 1e-13
    bad_rhs, bad_flux = occupation_interpolation_rhs(model,eq)
    bad_alpha = eq + 1e-3*model.b[0]
    near_bad_rhs, near_bad_flux = occupation_interpolation_rhs(model,bad_alpha)
    bad_velocity = -np.linalg.solve(model.capacity(bad_alpha), near_bad_rhs)
    bad_entropy_derivative = directional_derivative(model.entropy,bad_alpha,bad_velocity)
    assert abs(float(bad_flux[0])-5/7) < 1e-14
    assert bad_entropy_derivative < 0
    result["exact_references"] = {
        "M": [["2","0"],["0","8/49"]], "K": [["3","-1"],["-1","1/3"]],
        "Lambda": "4/3", "normalized_eigenvalues": ["0","85/24"]
    }
    result["occupation_interpolation_failure"] = {
        "false_equilibrium_flux": float(bad_flux[0]),
        "false_equilibrium_rhs": bad_rhs.tolist(),
        "false_equilibrium_energy_drift": float(energy @ bad_rhs),
        "perturbation": "alpha=beta*(1,3)+1e-3*(-3/2,1/2)",
        "nearby_flux": float(near_bad_flux[0]),
        "nearby_affinity": float(model.b[0] @ bad_alpha),
        "nearby_entropy_derivative": bad_entropy_derivative,
        "nearby_energy_drift": float(energy @ near_bad_rhs),
    }
    surrogate_cases = []
    for change in [-0.2,0.2]:
        physical_energy = energy + np.array([0.,change])
        alpha_physical = beta*physical_energy
        G = model.rhs(alpha_physical)
        _,_,mobility = model.reaction_data(alpha_physical)
        mismatch = model.b @ physical_energy
        expected_heating = float(beta*np.dot(model.weights*mobility,mismatch**2))
        velocity_off = model.coefficient_rhs(desired_off_eq)
        Goff = model.rhs(desired_off_eq)
        physical_off_derivative = directional_derivative(
            lambda z:model.energy(z,physical_energy),desired_off_eq,velocity_off)
        surrogate_cases.append({
            "physical_energy_coefficients": physical_energy.tolist(),
            "surrogate_energy_coefficients": energy.tolist(),
            "physical_mismatch": mismatch.tolist(),
            "surrogate_drift_at_physical_Bose": float(energy @ G),
            "physical_drift_at_physical_Bose": float(physical_energy @ G),
            "positive_square_prediction": expected_heating,
            "surrogate_Bose_rhs_norm": float(np.linalg.norm(model.rhs(eq))),
            "physical_occupation_error_at_surrogate_Bose":
                float(max(abs(bose(eq)-bose(alpha_physical)))),
            "physical_drift_off_equilibrium": float(physical_energy @ Goff),
            "physical_directional_derivative_off_equilibrium": physical_off_derivative,
        })
        assert expected_heating > 0
        assert abs(float(physical_energy @ G)-expected_heating) < 1e-14
        assert abs(float(energy @ G)) < 1e-14
        assert abs(physical_off_derivative-float(physical_energy @ Goff)) < 1e-13
    result["physical_vs_surrogate_energy"] = surrogate_cases
    return result


def evaluate_basis(branch, q):
    q = np.atleast_1d(q)
    result = np.zeros((len(q),9))
    result[:,3*branch] = 1
    result[:,3*branch+1] = np.cos(q)
    result[:,3*branch+2] = np.sin(q)
    return result


def periodic_model(root_residual=0.0, volume_nodes=64, parent_nodes=16):
    u = 0.6
    grid = -math.pi + 2*math.pi*np.arange(volume_nodes)/volume_nodes
    Phi = np.vstack([evaluate_basis(branch,grid) for branch in range(3)])
    masses = np.full(3*volume_nodes,1/volume_nodes)
    pgrid = -math.pi + 2*math.pi*(np.arange(parent_nodes)+0.17)/parent_nodes
    roots = [-math.acos(u-root_residual),math.acos(u-root_residual)]
    pvals,a_vals,bvals = [],[],[]
    for p in pgrid:
        for k in roots:
            pvals.append(evaluate_basis(0,[p])[0])
            a_vals.append(evaluate_basis(1,[k])[0])
            bvals.append(evaluate_basis(2,[p-k])[0])
    weights = np.full(2*parent_nodes,1/(parent_nodes*2*math.pi*math.sqrt(1-u*u)))
    energy = np.array([5+u,0,0,3,1,0,2,0,0.])
    model = EntropyWeakForm(Phi,masses,pvals,a_vals,bvals,weights)
    return model,energy


def periodic_checks():
    model,energy = periodic_model()
    perturbation = np.array([.04,-.03,.02,.02,.01,-.015,-.03,.02,.01])
    result = structural_checks(model,energy,.7,perturbation)
    target_mass = 1/(math.pi*.8)
    result["analytic_reaction_mass"] = target_mass
    result["reaction_mass_relative_error"] = sum(model.weights)/target_mass-1
    assert abs(result["reaction_mass_relative_error"]) < 1e-14
    base,base_energy = periodic_model(volume_nodes=32,parent_nodes=8)
    beta = .7
    _,_,base_mobility = base.reaction_data(beta*base_energy)
    base_total_mobility = float(np.dot(base.weights,base_mobility))
    residual_results = []
    for magnitude in [1e-2,1e-3,1e-4,1e-5,1e-6]:
        for sign in [-1,1]:
            residual = sign*magnitude
            perturbed,e = periodic_model(residual,volume_nodes=32,parent_nodes=8)
            alpha = beta*e
            G = perturbed.rhs(alpha)
            affinity,flux,mobility = perturbed.reaction_data(alpha)
            delta = perturbed.b @ e
            energy_event = float(-np.dot(perturbed.weights*flux,delta))
            positive_square = float(beta*np.dot(perturbed.weights*mobility,delta**2))
            entropy = float(np.dot(perturbed.weights*mobility,affinity**2))
            bound = math.sqrt(float(np.dot(perturbed.weights*mobility,delta**2))*entropy)
            thermal_roundoff = float(np.dot(perturbed.weights*mobility,
                np.abs(delta)*np.abs(affinity-beta*delta)))
            thermal_identity_bound = thermal_roundoff + 32*np.finfo(float).eps*(abs(energy_event)+abs(positive_square))
            residual_results.append({
                "prescribed_residual": residual,
                "maximum_actual_residual": float(max(abs(delta))),
                "energy_drift_event": energy_event,
                "energy_drift_moment_dot": float(e @ G),
                "energy_positive_square": positive_square,
                "thermal_identity_absolute_discrepancy": abs(energy_event-positive_square),
                "thermal_identity_roundoff_bound": thermal_identity_bound,
                "max_affinity_distributivity_residual": float(max(abs(affinity-beta*delta))),
                "energy_drift_over_quadratic_limit":
                    energy_event/(beta*base_total_mobility*magnitude*magnitude),
                "rhs_norm": float(np.linalg.norm(G)),
                "rhs_norm_over_absolute_residual": float(np.linalg.norm(G)/magnitude),
                "entropy_production": entropy, "energy_cauchy_schwarz_bound": bound,
                "bound_gap": bound-abs(energy_event),
            })
            # Discriminating replacement (2 Oct 2026 numerical review): the earlier
            # bound added the residual to its own tolerance and accepted states that
            # violate alpha=beta*e by 30-1580%. Since F=-Lambda*d exactly,
            # E-PS = sum w Lambda Delta (d-beta Delta) =: T. Check (i) the a-priori
            # rounding bound on d-beta*Delta and (ii) E-PS-T at a few ulps of E.
            unit = np.finfo(float).eps/2
            n_terms = perturbed.b.shape[1]+1
            gamma = n_terms*unit/(1-n_terms*unit)
            distributivity_bound = gamma*(np.abs(perturbed.b) @ np.abs(alpha))
            exact_T = math.fsum(perturbed.weights*mobility*delta*(affinity-beta*delta))
            identity_residual = abs(energy_event-positive_square-exact_T)
            residual_results[-1]["identity_residual_after_exact_T"] = identity_residual
            assert energy_event > 0
            assert np.all(np.abs(affinity-beta*delta) <= distributivity_bound)
            assert identity_residual <= 16*unit*abs(energy_event)
            assert abs(energy_event-positive_square) <= thermal_identity_bound  # diagnostic, not discriminating
            assert abs(energy_event) <= bound*(1+1e-10)  # consistency only: holds by Cauchy-Schwarz
    result["root_residual_tests"] = {
        "weights": "held fixed at nominal exact coarea weights to isolate support error",
        "equilibrium_total_mobility": base_total_mobility,
        "cases": residual_results,
    }
    return result


def stable_near_zero_checks():
    result = []
    for raw in ["0","1e-16","-1e-16","1e-12","-1e-12","1e-8","-1e-8","0.1","-0.1"]:
        with localcontext() as context:
            context.prec = 80
            affinity = Decimal(raw)
            reverse = Decimal("1.5")
            flux_reference = reverse*((-affinity).exp()-1)
            lambda_reference = -flux_reference/affinity if affinity else reverse
        af = float(raw)
        stable_flux = 1.5*np.expm1(-af)
        stable_lambda = 1.5*float(exprel_stable(np.array([-af]))[0])
        direct_flux = 1.5*np.exp(-af)-1.5
        flux_relative = abs(stable_flux-float(flux_reference))/abs(float(flux_reference)) if af else 0.
        lambda_relative = abs(stable_lambda-float(lambda_reference))/abs(float(lambda_reference))
        result.append({
            "affinity": raw,"reference_flux_80_digits": str(flux_reference),
            "stable_flux": float(stable_flux),"direct_subtraction_flux": float(direct_flux),
            "stable_flux_relative_error": flux_relative,
            "stable_lambda_relative_error": lambda_relative,
            "direct_flux_relative_error": abs(direct_flux-float(flux_reference))/abs(float(flux_reference)) if af else 0.,
        })
        assert flux_relative < 5e-15 and lambda_relative < 5e-15
    return result


def decimal_reaction_reference(xi_parent, xi_a, xi_b, prec=120):
    """Independent Decimal reference at the EXACT binary values of the float inputs."""
    with localcontext() as context:
        context.prec = prec
        def n(x):
            return 1/(Decimal(float(x)).exp()-1)
        npar, na, nb = n(xi_parent), n(xi_a), n(xi_b)
        A = npar*(1+na)*(1+nb)
        B = (1+npar)*na*nb
        affinity = Decimal(float(xi_parent))-Decimal(float(xi_a))-Decimal(float(xi_b))
        if affinity == 0:
            # Exact resonance: A=B mathematically; rounded A-B would be noise.
            return 0.0, float(A), affinity
        flux = A-B
        mobility = flux/(A.ln()-B.ln())
        return float(flux), float(mobility), affinity


def separated_variable_checks():
    """Regression for the NaN found by the 28 Sep numerical audit.

    The frozen candidate returned NaN mobility at xi=(1,400,400): its reverse
    product underflowed to zero while exprel(-affinity) overflowed.
    Exact resonance must give F=0 exactly. Elsewhere the flux error is judged
    against the conditioning of forming d=b@alpha, kappa=sum|xi|/|d| (2 Oct
    review): near resonance the loss is inherent, not an instability.
    """
    identity = np.eye(3)
    model = EntropyWeakForm(Phi=identity, masses=np.ones(3), parent=identity[[0]],
                            daughter_a=identity[[1]], daughter_b=identity[[2]],
                            weights=np.ones(1))
    unit = np.finfo(float).eps/2
    rows = []
    cases = [(1.0, 400.0, 400.0), (1.0, 30.0, 30.0), (0.5, 2.0, 3.0), (5.0, 2.0, 3.0),
             (800.0, 1.0, 1.0), (1e-3, 700.0, 705.0), (705.0, 700.0, 4.0), (3.0, 1.0, 2.0),
             (3.0+1e-9, 1.0, 2.0), (2.0, 1.0, 1.0-1e-12), (0.25, 0.125, 0.125+3e-6)]
    for xi in cases:
        affinity, flux, mobility = model.reaction_data(np.array(xi))
        flux_ref, mobility_ref, d_exact = decimal_reaction_reference(*xi)
        mobility_error = abs(float(mobility[0])-mobility_ref)/abs(mobility_ref)
        if d_exact == 0:
            flux_error, flux_tolerance = abs(float(flux[0])), 0.0
        else:
            kappa = sum(abs(x) for x in xi)/abs(float(d_exact))
            flux_error = abs(float(flux[0])-flux_ref)/abs(flux_ref)
            flux_tolerance = 1e-12 + 8*unit*kappa
        rows.append({"xi": list(xi), "exact_affinity": str(d_exact), "flux": float(flux[0]),
                     "flux_reference": flux_ref, "mobility": float(mobility[0]),
                     "mobility_reference": mobility_ref, "flux_relative_error": flux_error,
                     "flux_tolerance": flux_tolerance, "mobility_relative_error": mobility_error})
        assert np.isfinite(flux[0]) and np.isfinite(mobility[0]) and mobility[0] > 0
        assert mobility_error < 1e-12
        assert flux_error <= flux_tolerance
    return rows


def thermal_identity_negative_control():
    """The discriminating check must reject a state violating alpha=beta*e (2 Oct 2026).

    The superseded triangle-inequality bound accepted alpha=1.3*beta*e, whose
    thermal identity is violated by about 30%.
    """
    beta = .7
    model, e = periodic_model(1e-3, volume_nodes=32, parent_nodes=8)
    unit = np.finfo(float).eps/2
    n_terms = model.b.shape[1]+1
    gamma = n_terms*unit/(1-n_terms*unit)
    out = {}
    for scale in (1.0, 1.3):
        alpha = scale*beta*e
        affinity, _, _ = model.reaction_data(alpha)
        delta = model.b @ e
        out[str(scale)] = bool(np.all(np.abs(affinity-beta*delta) <= gamma*(np.abs(model.b) @ np.abs(alpha))))
    assert out["1.0"] and not out["1.3"]
    return {"distributivity_check_passes": out, "expected": {"1.0": True, "1.3": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    results = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpreter": sys.executable,"python":platform.python_version(),"numpy":np.__version__,
        "two_node":two_node_checks(),"periodic":periodic_checks(),
        "near_zero_affinity":stable_near_zero_checks(),
        "separated_entropy_variables":separated_variable_checks(),
        "thermal_identity_negative_control":thermal_identity_negative_control(),
        "all_assertions_passed":True,
        "scope":"Finite quadrature weak-form identities only; no time integrator, global invariance, or material validation.",
    }
    destination=args.output_dir/"C2_weak_form.json"
    destination.write_text(json.dumps(results,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "output":str(destination),"script_sha256":results["script_sha256"],
        "all_assertions_passed":True,
        "two_node_capacity_derivative_error":results["two_node"]["relative_capacity_derivative_error_equilibrium"],
        "periodic_collision_jacobian_error":results["periodic"]["relative_collision_jacobian_error_equilibrium"],
        "bad_occupation_interpolation_entropy_derivative":
            results["two_node"]["occupation_interpolation_failure"]["nearby_entropy_derivative"],
        "periodic_normalized_energy_kernel_residual":results["periodic"]["normalized_energy_kernel_residual"],
    },indent=2))


if __name__ == "__main__":
    main()


