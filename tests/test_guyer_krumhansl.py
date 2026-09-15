"""Guyer--Krumhansl: external validation, invariance, and the blind diagonal.

The first test of this file is the one that matters most, and its absence
allowed an error to survive a hundred others. Every earlier test compared the
model to itself: unimodularity, homogeneous limit, composition, splitting. All
of them passed on a model whose non-Fourier phase ran the wrong way, because
the error was internally consistent.

What caught it was a comparison with a closed form taken from outside the
codebase. That comparison is now the first thing this file does.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import quadrupoles as q  # noqa: E402

B = 12025.0
XI1 = 500e-9 / np.sqrt(60.0 / 2.41e6)


# --------------------------------------------------------------------------
# External validation
# --------------------------------------------------------------------------

def test_semi_infinite_phase_matches_published_closed_form():
    """Phase of a semi-infinite Cattaneo medium under modulated heating.

    Camacho de la Rosa, Esquivel-Sirvent and Becerril (2025) give, in the
    one-dimensional limit of their equation 15,

        phi(omega) = -45 degrees + (1/2) arctan(omega tau).

    It runs from -45 degrees in the diffusive limit to 0 in the wave limit,
    passing through -22.5 at omega tau = 1.

    Substituting p by p(1 + tau p) throughout the Fourier quadrupole produces
    the opposite behaviour, from -45 towards -90 degrees. The distinction is
    the whole point of this test.
    """
    tau = 1.0e-10
    for omega_tau in (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3):
        omega = omega_tau / tau
        p = 1j * omega
        _, coeff = fm.relaxation_pair(p, B, 1.0, tau, 0.0)
        phase = np.degrees(np.angle(1.0 / coeff))
        expected = -45.0 + 0.5 * np.degrees(np.arctan(omega_tau))
        assert phase == pytest.approx(expected, abs=1e-6), omega_tau


def test_amplitude_matches_published_closed_form():
    """The same comparison on the modulus."""
    tau = 1.0e-10
    omega = np.logspace(6, 13, 20)
    p = 1j * omega
    _, coeff = fm.relaxation_pair(p, B, 1.0, tau, 0.0)
    got = np.abs(1.0 / coeff)
    expected = np.abs(np.sqrt(1.0 + 1j * omega * tau) / (B * np.sqrt(1j * omega)))
    assert np.allclose(got, expected, rtol=1e-12)


# --------------------------------------------------------------------------
# Reduction to the simpler laws
# --------------------------------------------------------------------------

def test_both_times_zero_gives_fourier():
    p = 2j * np.pi * np.logspace(4, 11, 30)
    arg, coeff = fm.relaxation_pair(p, B, XI1, 0.0, 0.0)
    assert np.allclose(arg, XI1 * np.sqrt(p), rtol=1e-14)
    assert np.allclose(coeff, B * np.sqrt(p), rtol=1e-14)


def test_nonlocal_time_zero_gives_cattaneo():
    p = 2j * np.pi * np.logspace(4, 11, 30)
    tau = 1e-11
    arg, coeff = fm.relaxation_pair(p, B, XI1, tau, 0.0)
    assert np.allclose(arg, XI1 * np.sqrt(p * (1.0 + tau * p)), rtol=1e-12)
    assert np.allclose(coeff, B * np.sqrt(p / (1.0 + tau * p)), rtol=1e-12)


def test_sample_without_relaxation_is_unchanged():
    """The correction must leave every Fourier result untouched."""
    p = 2j * np.pi * np.logspace(4, 10, 30)
    s = fm.Sample()
    m = fm.stack_matrix(p, s)
    ref = q.homogeneous_wall(p, s.film_lam, s.film_rho_c, s.thickness)
    assert np.allclose(m, ref, rtol=1e-13)


def test_graded_layer_with_relaxation_is_refused():
    """Not treated, and refused rather than silently wrong."""
    s = fm.Sample(film_lam_back=120.0, film_rho_c_back=2.41e6, graded_xi1=500e-9/np.sqrt(60/2.41e6),
                  relaxation_time=1e-11)
    with pytest.raises(NotImplementedError):
        fm.stack_matrix(2j * np.pi * 1e8, s)


# --------------------------------------------------------------------------
# The blind diagonal
# --------------------------------------------------------------------------

def test_equal_times_reproduce_fourier_exactly():
    """When the two times coincide the square roots cancel.

    Such a medium is thermally indistinguishable from a Fourier medium,
    whatever the common value of the times. Checked to machine precision over
    four decades of relaxation time and eight of frequency.
    """
    p = 2j * np.pi * np.logspace(4, 12, 40)
    for tau in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8):
        arg, coeff = fm.relaxation_pair(p, B, XI1, tau, tau)
        assert np.allclose(arg, XI1 * np.sqrt(p), rtol=1e-13), tau
        assert np.allclose(coeff, B * np.sqrt(p), rtol=1e-13), tau


def test_equal_times_response_is_fourier():
    """The same statement on the complete front-face response."""
    p = 2j * np.pi * np.logspace(4, 11, 30)
    tau = 1e-10
    blind = fm.Sample(relaxation_time=tau, nonlocal_time=tau)
    fourier = fm.Sample()
    assert np.allclose(fm.response(p, blind), fm.response(p, fourier), rtol=1e-12)


def test_times_must_differ_to_be_separable():
    """A quantified criterion: below roughly a sixth of relative difference the
    two times cannot be told apart at one per cent noise."""
    import inversion as inv

    def sigma(ratio):
        tau_r = 1e-10
        s = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9,
                      relaxation_time=tau_r, nonlocal_time=tau_r * ratio)
        fc = np.log10(s.characteristic_frequency)
        f = np.logspace(fc, fc + 4, 80)
        _, cov = inv.fisher_analysis(
            f, s, ["relaxation_time", "nonlocal_time"], 0.01, 0.1)
        return np.sqrt(np.diag(cov)).max()

    assert sigma(0.99) > 10.0 * sigma(0.5)      # divergence at the diagonal
    assert sigma(0.5) < 0.05                    # comfortably separable
    assert sigma(3.0) < 0.05


# --------------------------------------------------------------------------
# The invariance survives the nonlocal term
# --------------------------------------------------------------------------

def _euler_residual(tau_r, tau_l, with_length_term=True, h=1e-6):
    """Residual of  e dG/de + 2a dG/da + lam dG/dlam ( + ell dG/dell ).

    The parameters are the physical triplet plus the nonlocal length, the
    latter entering through tau_l = 3 ell^2 / a.
    """
    e0, kappa0 = 500e-9, 60.0
    alpha0 = kappa0 / 2.41e6
    ell0 = np.sqrt(tau_l * alpha0 / 3.0) if tau_l > 0 else 0.0
    # The band is capped: at high frequency the argument tends to
    # xi_1 sqrt(tau_R / tau_l) sqrt(p), which overflows when the two times are
    # far apart.
    f_top = min(1e11, ((400.0 / (e0 / np.sqrt(alpha0)
                                 * np.sqrt(max(tau_r / max(tau_l, 1e-30), 1.0))))
                       ** 2) / (2.0 * np.pi))
    f = np.logspace(4, np.log10(max(f_top, 1e6)), 30)
    p = 2j * np.pi * f
    b_sub = 10298.0

    def g(e, alpha, kappa, ell):
        b = kappa / np.sqrt(alpha)
        xi1 = e / np.sqrt(alpha)
        tl = 3.0 * ell ** 2 / alpha
        m = fm.relaxation_wall(p, b, xi1, tau_r, tl)
        z = 1.0 / (b_sub * np.sqrt(p))
        return q.front_face_temperature(m, z)

    base = [e0, alpha0, kappa0, ell0]
    weights = [1.0, 2.0, 1.0, 1.0 if with_length_term else 0.0]

    terms = []
    for i, w in enumerate(weights):
        if w == 0.0 or base[i] == 0.0:
            continue
        up, dn = list(base), list(base)
        up[i] *= 1.0 + h
        dn[i] *= 1.0 - h
        d = (g(*up) - g(*dn)) / (2.0 * h * base[i])
        terms.append(w * base[i] * d)

    scale = sum(np.abs(t) for t in terms)
    return float(np.max(np.abs(sum(terms)) / scale))


def test_invariance_survives_the_nonlocal_term():
    """Main negative result, extended.

    The scale invariance established under Fourier, and shown to survive
    Cattaneo, survives Guyer-Krumhansl as well. The reason is structural: the
    nonlocal length never enters alone, only through tau_l = 3 ell^2 / a,
    which is a time. No length appears separately in the response, and that is
    precisely what leaves the depth free to be relabelled.
    """
    for tau_r, tau_l in ((1e-11, 1e-10), (1e-10, 1e-11), (1e-9, 1e-12)):
        assert _euler_residual(tau_r, tau_l) < 1e-5, (tau_r, tau_l)


def test_the_length_term_is_required():
    """Omitting the derivative with respect to the nonlocal length breaks the
    identity, which shows that the length does transform under the group."""
    tau_r, tau_l = 1e-11, 1e-10
    with_it = _euler_residual(tau_r, tau_l, with_length_term=True)
    without = _euler_residual(tau_r, tau_l, with_length_term=False)
    assert with_it < 1e-5
    assert without > 0.05
    assert without > 1000.0 * with_it


# --------------------------------------------------------------------------
# The microscopic link
# --------------------------------------------------------------------------

def test_nonlocal_time_equals_nine_fifths_of_the_normal_time():
    """With ell^2 = v^2 tau_N tau_R / 5 and a = v^2 tau_R / 3,

        tau_l = 3 ell^2 / a = (9/5) tau_N ,

    so the blind condition tau_R = tau_l reads tau_R = 1.8 tau_N: a ratio of
    collision times, which varies with temperature.
    """
    v, tau_n, tau_r = 6000.0, 3.0e-12, 1.1e-11
    ell_sq = v ** 2 * tau_n * tau_r / 5.0
    a = v ** 2 * tau_r / 3.0
    assert 3.0 * ell_sq / a == pytest.approx(1.8 * tau_n, rel=1e-12)
