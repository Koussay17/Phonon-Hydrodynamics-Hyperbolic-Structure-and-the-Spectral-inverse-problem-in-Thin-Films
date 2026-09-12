"""Identifiability under the Cattaneo constitutive law.

This file records three results, each verified independently of the
analytical derivation that produced it.

1. Negative. The relaxation term does not lift the degeneracy established by
   Krapez and Rigollet under Fourier. Their scaling invariance survives
   intact, because the transit time, the effusivity and the relaxation time
   are all invariant under the group.

2. Positive. The relaxation time enters as a genuinely independent parameter.
   Estimating it degrades neither the effusivity nor the transit time.

3. Quantitative. Its relative uncertainty follows an inverse scaling law in
   the product of the highest measured angular frequency and the relaxation
   time, with a threshold at one. That threshold coincides with the spectral
   transition where the energy of the Schrodinger analogy leaves the
   imaginary axis.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import inversion as inv  # noqa: E402

SIGMA_REL = 0.01
SIGMA_PHASE = 0.1

# Film properties expressed in the triplet of the 2017 note.
L0 = 500e-9                      # thickness
KAPPA0 = 60.0                    # conductivity
ALPHA0 = KAPPA0 / 2.41e6         # diffusivity


def sample_from_triplet(L, alpha, kappa, tau=0.0, **kw):
    """Build a Sample from (thickness, diffusivity, conductivity)."""
    return fm.Sample(film_lam=kappa, film_rho_c=kappa / alpha,
                     thickness=L, relaxation_time=tau, **kw)


def band_above_diffusion(s, decades=2.0, n=80):
    fc = np.log10(s.characteristic_frequency)
    return np.logspace(fc, fc + decades, n)


# --------------------------------------------------------------------------
# The extension reduces to Fourier
# --------------------------------------------------------------------------

def test_zero_relaxation_recovers_fourier():
    s_f = fm.Sample()
    s_c = fm.Sample(relaxation_time=0.0, sub_relaxation_time=0.0)
    p = 2j * np.pi * np.logspace(4, 9, 30)
    assert np.allclose(fm.response(p, s_f), fm.response(p, s_c), rtol=1e-14)


def test_effective_parameter_definition():
    p = 2j * np.pi * 1e9
    tau = 1e-10
    assert fm.effective_p(p, 0.0) == pytest.approx(p)
    assert fm.effective_p(p, tau) == pytest.approx(p * (1.0 + tau * p))


def test_relaxation_frequency():
    tau = 1e-10
    s = fm.Sample(relaxation_time=tau)
    assert s.relaxation_frequency == pytest.approx(1.0 / (2.0 * np.pi * tau))
    assert np.isinf(fm.Sample().relaxation_frequency)


# --------------------------------------------------------------------------
# Result 1, negative: the degeneracy survives
# --------------------------------------------------------------------------

def _euler_residual(tau, freq, h=1e-7):
    """Relative residual of  L dG/dL + 2 alpha dG/dalpha + kappa dG/dkappa.

    The combination vanishes identically under Fourier, which is the analytic
    content of the 2017 note. The question is whether a relaxation time makes
    it depart from zero.
    """
    def g(L, alpha, kappa):
        s = sample_from_triplet(L, alpha, kappa, tau)
        return fm.response(2j * np.pi * freq, s)

    def d(i):
        args = [L0, ALPHA0, KAPPA0]
        up, dn = list(args), list(args)
        up[i] *= 1.0 + h
        dn[i] *= 1.0 - h
        return (g(*up) - g(*dn)) / (2.0 * h * args[i])

    terms = [L0 * d(0), 2.0 * ALPHA0 * d(1), KAPPA0 * d(2)]
    total = sum(terms)
    scale = sum(np.abs(t) for t in terms)
    return float(np.max(np.abs(total) / scale))


def test_euler_relation_holds_under_fourier():
    """Reproduces the central relation of Krapez and Rigollet (2017)."""
    freq = np.logspace(4, 10, 40)
    assert _euler_residual(0.0, freq) < 1e-6


def test_relaxation_does_not_lift_the_degeneracy():
    """The main negative result.

    The residual stays at the level of the finite-difference truncation error
    whatever the relaxation time, up to omega tau of several hundred. The
    scaling invariance is therefore untouched: the relaxation time is a time,
    not a length, and a depth rescaling leaves it alone, so the invariants of
    the group remain invariant.
    """
    freq = np.logspace(4, 10, 40)
    reference = _euler_residual(0.0, freq)
    for tau in (1e-12, 1e-11, 1e-10, 1e-9, 1e-8):
        omega_tau = 2.0 * np.pi * freq.max() * tau
        residual = _euler_residual(tau, freq)
        assert residual < 10.0 * reference + 1e-6, (tau, omega_tau, residual)


# --------------------------------------------------------------------------
# Result 2, positive: the relaxation time is an independent parameter
# --------------------------------------------------------------------------

def test_estimating_relaxation_does_not_degrade_the_others():
    """Adding the relaxation time to the fit leaves the uncertainties on the
    conductivity and the heat capacity essentially unchanged.

    The thickness is held fixed, as it must be: it is known independently, and
    without it the triplet is not identifiable at all.
    """
    tau = 1e-10
    s = sample_from_triplet(L0, ALPHA0, KAPPA0, tau)
    f = band_above_diffusion(s)

    _, cov2 = inv.fisher_analysis(f, s, ["film_lam", "film_rho_c"],
                                  SIGMA_REL, SIGMA_PHASE)
    _, cov3 = inv.fisher_analysis(f, s, ["film_lam", "film_rho_c",
                                         "relaxation_time"],
                                  SIGMA_REL, SIGMA_PHASE)

    s2 = np.sqrt(np.diag(cov2))
    s3 = np.sqrt(np.diag(cov3))[:2]
    assert np.all(s3 < 2.0 * s2)


def test_relaxation_is_identifiable_inside_the_band():
    tau = 1e-9
    s = sample_from_triplet(L0, ALPHA0, KAPPA0, tau)
    f = band_above_diffusion(s)
    _, cov = inv.fisher_analysis(f, s, ["film_lam", "film_rho_c",
                                        "relaxation_time"],
                                 SIGMA_REL, SIGMA_PHASE)
    assert np.sqrt(np.diag(cov))[2] < 0.02          # better than 2 per cent


# --------------------------------------------------------------------------
# Result 3, quantitative: the scaling law and its threshold
# --------------------------------------------------------------------------

def _sigma_tau(tau, n=80, decades=2.0):
    s = sample_from_triplet(L0, ALPHA0, KAPPA0, tau)
    f = band_above_diffusion(s, decades, n)
    _, cov = inv.fisher_analysis(f, s, ["film_lam", "film_rho_c",
                                        "relaxation_time"],
                                 SIGMA_REL, SIGMA_PHASE)
    return np.sqrt(np.diag(cov))[2], 2.0 * np.pi * f.max() * tau


def test_uncertainty_scales_as_the_inverse_of_omega_tau():
    """Below the threshold the product of the relative uncertainty and of
    omega tau is constant over decades.

    The constant depends on the noise level and on the number of points; the
    inverse scaling and the position of the threshold do not.
    """
    products = []
    for tau in (1e-14, 1e-13, 1e-12, 1e-11):
        sigma, omega_tau = _sigma_tau(tau)
        assert omega_tau < 0.2
        products.append(sigma * omega_tau)

    products = np.array(products)
    assert products.std() / products.mean() < 0.05


def test_scaling_exponent_is_exactly_minus_one_below_the_threshold():
    """The local exponent of the uncertainty against omega tau is -1 to three
    decimal places over three decades, then leaves that value at the knee."""
    def exponent(t1, t2):
        s1, w1 = _sigma_tau(t1)
        s2, w2 = _sigma_tau(t2)
        return np.log(s2 / s1) / np.log(w2 / w1)

    assert exponent(1e-15, 1e-14) == pytest.approx(-1.0, abs=1e-3)
    assert exponent(1e-14, 1e-13) == pytest.approx(-1.0, abs=1e-3)
    assert exponent(1e-13, 1e-12) == pytest.approx(-1.0, abs=1e-3)

    # the knee is already felt one decade later
    assert exponent(1e-12, 1e-11) > -1.0 + 1e-3

    # past the knee the gain per decade collapses
    assert exponent(3e-10, 3e-9) > -0.5


def test_knee_is_located_at_omega_tau_of_order_one():
    """Scanning the exponent locates the knee: it is the last point where the
    inverse law still holds, and it sits at omega tau of order one."""
    taus = np.logspace(-14, -8, 13)
    data = [_sigma_tau(t) for t in taus]
    wt = np.array([d[1] for d in data])
    sg = np.array([d[0] for d in data])
    exponents = np.diff(np.log(sg)) / np.diff(np.log(wt))

    last_inverse = np.where(exponents < -0.9)[0]
    knee = wt[last_inverse[-1] + 1]
    assert 0.1 < knee < 10.0


def test_an_optimum_exists_beyond_the_threshold():
    """Pushing the band ever higher does not help indefinitely: past roughly
    two decades above the threshold the uncertainty starts growing again, so
    there is a finite optimum rather than a monotone gain.
    """
    taus = np.logspace(-12, -7.5, 13)
    data = [_sigma_tau(t) for t in taus]
    wt = np.array([d[1] for d in data])
    sg = np.array([d[0] for d in data])

    best = int(np.argmin(sg))
    assert 0 < best < len(sg) - 1                # interior minimum
    assert wt[best] > 1.0                        # beyond the threshold
    assert sg[-1] > sg[best]                     # degradation afterwards


def test_relaxation_unmeasurable_far_below_the_threshold():
    sigma, omega_tau = _sigma_tau(1e-15)
    assert omega_tau < 1e-4
    assert sigma > 1.0                  # worse than one hundred per cent


# --------------------------------------------------------------------------
# The spectral counterpart of the same threshold
# --------------------------------------------------------------------------

def test_analogy_energy_argument():
    """The energy of the Schrodinger analogy reaches exactly -45 degrees at
    omega tau = 1, which is the threshold found above from estimation theory.

    Two independent routes, one spectral and one statistical, single out the
    same value of the product omega tau.
    """
    tau = 1e-10
    for omega_tau, expected in ((0.001, -89.94), (1.0, -45.0), (1000.0, -0.057)):
        f = omega_tau / (2.0 * np.pi * tau)
        e = fm.analogy_energy(f, tau)
        assert np.degrees(np.angle(e)) == pytest.approx(expected, abs=0.01)


def test_analogy_energy_depends_only_on_the_product():
    """Different relaxation times at matched omega tau give the same argument."""
    args = []
    for tau in (1e-12, 1e-10, 1e-8):
        f = 1.0 / (2.0 * np.pi * tau)          # omega tau = 1
        args.append(np.degrees(np.angle(fm.analogy_energy(f, tau))))
    assert np.allclose(args, -45.0, atol=1e-9)


# --------------------------------------------------------------------------
# A constraint the relaxation term introduces on the formalism itself
# --------------------------------------------------------------------------

def test_relaxation_shrinks_the_usable_band():
    """The transfer matrix overflows sooner when a relaxation time is present.

    The effective parameter behaves as tau omega squared once omega tau
    exceeds one, so the argument of the hyperbolic functions grows linearly
    with frequency instead of as its square root. The band that double
    precision can carry therefore shrinks as the relaxation time grows.

    This is a limitation of the transfer-matrix representation, not of the
    physics, and it bounds the range over which the optimum of the previous
    test can be sought.
    """
    s_f = sample_from_triplet(L0, ALPHA0, KAPPA0, 0.0)
    f = band_above_diffusion(s_f)
    p = 2j * np.pi * f

    # Fourier: the band is carried without difficulty
    fm.response(p, s_f)

    # a large relaxation time overflows on the very same band
    with pytest.raises(OverflowError):
        fm.response(p, sample_from_triplet(L0, ALPHA0, KAPPA0, 1e-5))


def test_overflow_threshold_follows_the_predicted_law():
    """The argument scales as omega sqrt(tau) xi1, so the largest usable
    relaxation time falls as the inverse square of the highest frequency."""
    s = sample_from_triplet(L0, ALPHA0, KAPPA0, 0.0)
    f = band_above_diffusion(s)
    omega_max = 2.0 * np.pi * f.max()
    predicted = (700.0 / (omega_max * s.xi1)) ** 2

    ok = sample_from_triplet(L0, ALPHA0, KAPPA0, 0.3 * predicted)
    fm.response(2j * np.pi * f, ok)

    with pytest.raises(OverflowError):
        fm.response(2j * np.pi * f,
                    sample_from_triplet(L0, ALPHA0, KAPPA0, 3.0 * predicted))
