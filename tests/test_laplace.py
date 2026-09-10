"""Unit tests for src/laplace.py.

The decisive test is `test_semi_infinite_pulse_matches_closed_form`: it
validates the inversion together with the quadrupole chain, against a solution
known in closed form.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import laplace as lp  # noqa: E402
import quadrupoles as q  # noqa: E402


LAM = 30.0
RHO_C = 3.0e6
E = 5.0e-6

TIMES = np.logspace(-6, -1, 30)


# --------------------------------------------------------------------------
# The weights themselves
# --------------------------------------------------------------------------

def test_odd_n_rejected():
    with pytest.raises(ValueError):
        lp.stehfest_coefficients(11)


def test_weights_sum_to_zero():
    """A consequence of the construction; it is what makes the method exact
    on constants."""
    for n in (6, 8, 10, 12, 14):
        assert lp.stehfest_coefficients(n).sum() == pytest.approx(0.0, abs=1e-6)


def test_weights_alternate_and_grow():
    """The alternating signs and rapid growth are the source of the
    cancellation that limits accuracy in double precision."""
    v = lp.stehfest_coefficients(12)
    assert np.all(np.diff(np.sign(v)) != 0)
    assert np.abs(v).max() > 1e4


# --------------------------------------------------------------------------
# Analytical transform pairs
# --------------------------------------------------------------------------

def test_step_function():
    """F(p) = 1/p  ->  f(t) = 1."""
    got = lp.stehfest_inverse(lambda p: 1.0 / p, TIMES)
    assert np.allclose(got, 1.0, rtol=1e-8)


def test_ramp():
    """F(p) = 1/p^2  ->  f(t) = t.

    The relative error settles at about 1e-6, which is the intrinsic accuracy
    of the method at n = 12 in double precision.
    """
    got = lp.stehfest_inverse(lambda p: 1.0 / p ** 2, TIMES)
    assert np.allclose(got, TIMES, rtol=1e-5)


def test_decaying_exponential():
    """F(p) = 1/(p+a)  ->  f(t) = exp(-a t), while the signal stays within one
    decade of its maximum.
    """
    a = 1.0e3
    t = np.logspace(-6, -1, 60)
    t = t[np.exp(-a * t) > 1e-1]
    got = lp.stehfest_inverse(lambda p: 1.0 / (p + a), t)
    assert np.allclose(got, np.exp(-a * t), rtol=1e-3)


def test_accuracy_degrades_with_signal_decay():
    """Quantifies the usable range of the method.

    Measured relative error on exp(-a t) at n = 12, as a function of how far
    the signal has decayed from its maximum:

        above 1e-1 of the peak   ->  ~5e-4
        above 1e-2               ->  ~7e-3
        above 1e-3               ->  ~1e-1

    Practical rule: the inversion holds over roughly two decades of decay.
    Beyond that the result is dominated by cancellation between the
    alternating weights, and a fit performed on such a tail is driven by
    numerical noise rather than by the physics.
    """
    a = 1.0e3
    errors = {}
    for floor in (1e-1, 1e-2, 1e-3):
        t = np.logspace(-6, -1, 60)
        t = t[np.exp(-a * t) > floor]
        got = lp.stehfest_inverse(lambda p: 1.0 / (p + a), t)
        errors[floor] = np.max(np.abs(got - np.exp(-a * t)) / np.exp(-a * t))

    assert errors[1e-1] < 1e-3
    assert errors[1e-2] < 1e-2
    assert errors[1e-3] > 1e-2          # the degradation is real
    assert errors[1e-1] < errors[1e-2] < errors[1e-3]


def test_rejects_non_positive_time():
    with pytest.raises(ValueError):
        lp.stehfest_inverse(lambda p: 1.0 / p, 0.0)


# --------------------------------------------------------------------------
# Coupling with the quadrupoles
# --------------------------------------------------------------------------

def test_semi_infinite_pulse_matches_closed_form():
    """Front face of a semi-infinite medium after a Dirac pulse.

    The Laplace-domain expression comes from the quadrupole of a slab closed
    by an impedance of the same effusivity; the time-domain counterpart is
    known exactly. Validates the whole chain.
    """
    b = q.effusivity(LAM, RHO_C)

    def theta_p(p):
        m = q.homogeneous_wall(p, LAM, RHO_C, E)
        z = q.semi_infinite_impedance(p, b)
        return q.front_face_temperature(m, z, power=1.0, h=0.0).real

    got = lp.stehfest_inverse(theta_p, TIMES)
    expected = lp.semi_infinite_pulse_response(TIMES, b)
    assert np.allclose(got, expected, rtol=1e-5)


def test_adiabatic_slab_long_time_limit():
    """A slab insulated on both faces heated by a pulse of energy density Q
    ends up uniformly heated:  theta -> Q / (rho c e).

    Checks the low-frequency behaviour of the quadrupole through the inversion.
    """
    b = q.effusivity(LAM, RHO_C)
    a = q.diffusivity(LAM, RHO_C)
    xi1 = q.xi_from_thickness(E, a)

    def theta_p(p):
        # rear face adiabatic: theta_0 / phi_0 = A / C = coth(sqrt(p) xi1) / (b sqrt(p))
        sp = np.sqrt(p)
        return (np.cosh(sp * xi1) / np.sinh(sp * xi1) / (b * sp)).real

    t_long = np.array([50.0 * E ** 2 / a])
    got = lp.stehfest_inverse(theta_p, t_long)
    expected = 1.0 / (RHO_C * E)
    assert got[0] == pytest.approx(expected, rel=1e-4)


# --------------------------------------------------------------------------
# Accuracy versus the number of terms
# --------------------------------------------------------------------------

def test_accuracy_is_not_monotonic_in_n():
    """Accuracy improves with n up to about 12, then degrades as cancellation
    takes over. Knowing this prevents raising n in the false belief that
    precision follows.
    """
    t = np.array([1.0e-3])
    exact = 1.0 / np.sqrt(np.pi * t[0])

    errors = {}
    for n in (6, 8, 10, 12, 14, 16, 18, 20):
        got = lp.stehfest_inverse(lambda p: 1.0 / np.sqrt(p), t, n=n)[0]
        errors[n] = abs(got - exact) / exact

    assert errors[12] < errors[6]
    assert errors[20] > errors[12]
