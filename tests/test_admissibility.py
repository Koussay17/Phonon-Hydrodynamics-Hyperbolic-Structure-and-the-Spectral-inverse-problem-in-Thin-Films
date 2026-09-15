"""Admissibility: entropy on one side, propagation speed on the other.

The point these tests record is that the two criteria are independent. All
three laws admit a convex entropy with non-negative production; only Cattaneo
propagates at finite speed. The nonlocal term of Guyer-Krumhansl is diffusive
in the flux, and above the frequency where it dominates it restores the
infinite speed that Cattaneo had removed.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import admissibility as adm  # noqa: E402

A = 60.0 / 2.41e6          # diffusivity of the reference film
TAU_R = 1.0e-11
TAU_L = 1.0e-10
LAM = 60.0
T0 = 300.0


def _exponent(omega, k):
    """Local slope of log |Re k| against log omega."""
    return float(np.polyfit(np.log(omega), np.log(np.abs(np.real(k))), 1)[0])


# --------------------------------------------------------------------------
# Dispersion
# --------------------------------------------------------------------------

def test_fourier_wavenumber_grows_as_the_square_root():
    omega = np.logspace(8, 14, 20)
    k = adm.wavenumber(omega, A)
    assert _exponent(omega, k) == pytest.approx(0.5, abs=1e-6)


def test_cattaneo_wavenumber_grows_linearly():
    """A linear growth means a constant phase velocity: a genuine wave."""
    omega = np.logspace(13, 17, 20)
    k = adm.wavenumber(omega, A, TAU_R)
    assert _exponent(omega, k) == pytest.approx(1.0, abs=1e-3)


def test_guyer_krumhansl_returns_to_the_square_root():
    """The main negative result of this file.

    The exponent falls back to one half, the signature of diffusion, because
    the nonlocal term dominates the term in tau_R at high frequency.
    """
    omega = np.logspace(13, 17, 20)
    k = adm.wavenumber(omega, A, TAU_R, TAU_L)
    assert _exponent(omega, k) == pytest.approx(0.5, abs=1e-2)


def test_cattaneo_phase_velocity_saturates():
    omega = np.logspace(13, 17, 20)
    v = adm.phase_velocity(omega, A, TAU_R)
    expected = adm.cattaneo_speed(A, TAU_R)
    assert np.allclose(v, expected, rtol=1e-3)
    assert expected == pytest.approx(np.sqrt(A / TAU_R))


def test_guyer_krumhansl_phase_velocity_diverges():
    """The velocity grows as the square root of the frequency, without bound,
    and overtakes the Cattaneo limit by orders of magnitude."""
    omega = np.logspace(13, 17, 20)
    v = adm.phase_velocity(omega, A, TAU_R, TAU_L)
    slope = float(np.polyfit(np.log(omega), np.log(v), 1)[0])
    assert slope == pytest.approx(0.5, abs=1e-2)
    assert v[-1] > 1000.0 * adm.cattaneo_speed(A, TAU_R)


def test_fourier_phase_velocity_diverges():
    omega = np.logspace(8, 16, 20)
    v = adm.phase_velocity(omega, A)
    slope = float(np.polyfit(np.log(omega), np.log(v), 1)[0])
    assert slope == pytest.approx(0.5, abs=1e-6)


def test_dispersion_reduces_to_fourier():
    omega = np.logspace(6, 12, 20)
    assert np.allclose(adm.wavenumber(omega, A, 0.0, 0.0),
                       adm.wavenumber(omega, A), rtol=1e-14)


def test_equal_times_reproduce_fourier_dispersion():
    """The blind diagonal again, seen from the dispersion relation."""
    omega = np.logspace(6, 16, 30)
    for tau in (1e-12, 1e-10, 1e-8):
        assert np.allclose(adm.wavenumber(omega, A, tau, tau),
                           adm.wavenumber(omega, A), rtol=1e-12), tau


# --------------------------------------------------------------------------
# Characteristic analysis
# --------------------------------------------------------------------------

def test_characteristic_speeds_are_real_and_opposite():
    s = np.sort(adm.characteristic_speeds(A, TAU_R))
    assert np.all(np.isreal(s))
    expected = np.sqrt(A / TAU_R)
    assert s[0] == pytest.approx(-expected)
    assert s[1] == pytest.approx(+expected)


def test_strict_hyperbolicity_requires_a_positive_relaxation_time():
    with pytest.raises(ValueError):
        adm.characteristic_matrix(A, 0.0)
    assert np.isinf(adm.cattaneo_speed(A, 0.0))


def test_speed_matches_the_debye_relation():
    """With a = v^2 tau_R / 3 the limiting speed is v / sqrt(3)."""
    v = 6000.0
    tau = 1.1e-11
    a = v ** 2 * tau / 3.0
    assert adm.cattaneo_speed(a, tau) == pytest.approx(v / np.sqrt(3.0))


# --------------------------------------------------------------------------
# Entropy
# --------------------------------------------------------------------------

def test_entropy_production_is_non_negative():
    """A sum of squares with positive coefficients, whatever the fields."""
    rng = np.random.default_rng(0)
    q = rng.standard_normal(500) * 1e4
    dq = rng.standard_normal(500) * 1e10
    for ell_sq in (0.0, 1e-16, 1e-14):
        assert np.all(adm.entropy_production(q, dq, LAM, T0, ell_sq) >= 0.0)


def test_cattaneo_production_has_the_classical_form():
    q = np.array([1.0e4])
    got = adm.entropy_production(q, np.array([1e12]), LAM, T0, 0.0)
    assert got[0] == pytest.approx(q[0] ** 2 / (LAM * T0 ** 2))


def test_nonlocal_term_adds_a_square():
    q, dq, ell_sq = 1.0e4, 1.0e12, 1.0e-15
    base = adm.entropy_production(q, dq, LAM, T0, 0.0)
    full = adm.entropy_production(q, dq, LAM, T0, ell_sq)
    assert full - base == pytest.approx(3.0 * ell_sq * dq ** 2 / (LAM * T0 ** 2))


def test_entropy_flux_coefficient():
    """Without this extra flux the production of Guyer-Krumhansl keeps a cross
    term of indefinite sign and is not a sum of squares."""
    ell_sq = 1.0e-15
    assert adm.entropy_flux_coefficient(LAM, T0, ell_sq) == pytest.approx(
        3.0 * ell_sq / (LAM * T0 ** 2))
    assert adm.entropy_flux_coefficient(LAM, T0, 0.0) == 0.0


# --------------------------------------------------------------------------
# Symbolic derivation, when sympy is available
# --------------------------------------------------------------------------

def test_entropy_production_derivation():
    """Keep gradients while expanding to quadratic order about equilibrium."""
    import sympy as sp
    z, eps = sp.symbols("z eps")
    theta = sp.Function("theta")(z)
    j = sp.Function("j")(z)
    lam, tau, t0, l2, beta = sp.symbols("lam tau T0 l2 beta", positive=True)
    T, flux = t0 + eps * theta, eps * j
    def production(extra, wrong_flux=False):
        qt = (-flux-lam*sp.diff(T,z)+3*l2*sp.diff(flux,z,2))/tau
        sdot = -sp.diff(flux,z)/T - tau*flux*qt/(lam*t0**2)
        js = flux/(t0 if wrong_flux else T) + extra*flux*sp.diff(flux,z)
        return sp.expand(sp.series(sdot+sp.diff(js,z),eps,0,3).removeO()).coeff(eps,2)
    raw = production(beta)
    solution = sp.solve(sp.Eq(raw.coeff(sp.diff(j,z,2)),0),beta)[0]
    assert sp.simplify(solution-3*l2/(lam*t0**2)) == 0
    target = (j*j+3*l2*sp.diff(j,z)**2)/(lam*t0**2)
    assert sp.simplify(raw.subs(beta,solution)-target) == 0
    assert sp.simplify(production(0).subs(l2,0)-j*j/(lam*t0**2)) == 0
    # This mutation must fail; replacing T by T0 too early used to hide it.
    assert sp.simplify(production(solution,True)-target) != 0


def test_table_is_printable():
    text = adm.admissibility_table()
    assert "Cattaneo" in text and "independent" in text
