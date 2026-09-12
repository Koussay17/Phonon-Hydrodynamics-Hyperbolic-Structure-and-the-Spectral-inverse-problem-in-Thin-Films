"""Unit tests for src/bayesian.py.

The chains are deliberately short: these tests check correctness, not
convergence quality. A production run would use far more steps.

The test that justifies the module is
`test_posterior_is_asymmetric_when_badly_conditioned`: it shows the posterior
departing from the ellipse that Levenberg-Marquardt assumes.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

pytest.importorskip("emcee")

import bayesian as bay  # noqa: E402
import forward_model as fm  # noqa: E402
import inversion as inv  # noqa: E402

SIGMA_REL = 0.01
SIGMA_PHASE = 0.1


def truth():
    return fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9)


def band(sample, lo, hi, n=40):
    fc = np.log10(sample.characteristic_frequency)
    return np.logspace(fc + lo, fc + hi, n)


def synthetic(sample, freq, rng=None):
    amp, ph = fm.modulated_response(freq, sample)
    if rng is None:
        return amp, ph
    return fm.add_noise(amp, ph, SIGMA_REL, SIGMA_PHASE, rng=rng)


# --------------------------------------------------------------------------
# Basic behaviour
# --------------------------------------------------------------------------

def test_posterior_recovers_the_truth():
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=11)

    res = bay.sample_posterior(f, amp, ph, s, ["film_lam", "film_rho_c"],
                               sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE,
                               n_walkers=24, n_steps=1500, seed=1)

    ci = res.credible_interval(0.95)
    assert ci[0, 0] < 60.0 < ci[0, 1]
    assert ci[1, 0] < 2.41e6 < ci[1, 1]


def test_acceptance_is_reasonable():
    """A healthy ensemble sampler accepts between roughly 15 and 60 per cent."""
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=12)
    res = bay.sample_posterior(f, amp, ph, s, ["film_lam", "film_rho_c"],
                               n_walkers=24, n_steps=1200, seed=2)
    assert 0.10 < res.acceptance < 0.75


def test_sampling_is_reproducible():
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=13)
    kw = dict(n_walkers=16, n_steps=600, seed=7)
    a = bay.sample_posterior(f, amp, ph, s, ["film_lam"], **kw)
    b = bay.sample_posterior(f, amp, ph, s, ["film_lam"], **kw)
    assert np.allclose(a.median, b.median)


def test_rejects_non_positive_parameter():
    s = fm.Sample(contact_resistance=0.0)
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f)
    with pytest.raises(ValueError):
        bay.sample_posterior(f, amp, ph, s, ["contact_resistance"],
                             n_walkers=8, n_steps=50)


def test_prior_bounds_are_enforced():
    """No sample may fall outside the prior support."""
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=14)
    res = bay.sample_posterior(f, amp, ph, s, ["film_lam"],
                               bounds=[(50.0, 70.0)],
                               n_walkers=16, n_steps=800, seed=3)
    assert res.chain.min() >= 50.0 and res.chain.max() <= 70.0


# --------------------------------------------------------------------------
# Agreement with, and departure from, the linearised estimate
# --------------------------------------------------------------------------

def test_posterior_matches_fisher_when_well_conditioned():
    """In the favourable band the model is nearly linear over the region the
    data allow, so both approaches must agree."""
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=21)
    names = ["film_lam", "film_rho_c"]

    lm = inv.fit_modulated(f, amp, ph, s, names=names,
                           sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    post = bay.sample_posterior(f, amp, ph, s, names,
                                sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE,
                                n_walkers=32, n_steps=4000, seed=4)

    spread = post.chain.std(axis=0, ddof=1)
    for got, expected in zip(spread, lm.std_errors):
        assert got == pytest.approx(expected, rel=0.35)


def test_posterior_is_asymmetric_when_badly_conditioned():
    """Far above the transition only the effusivity is measured, and the
    posterior becomes a long curved ridge rather than an ellipse.

    The asymmetry of the credible interval is what the Gaussian error bar of
    the linearised estimate cannot express.
    """
    s = truth()
    f_good = band(s, 0, 1)
    f_bad = band(s, 2, 3)
    names = ["film_lam", "film_rho_c"]

    amp_g, ph_g = synthetic(s, f_good, rng=31)
    amp_b, ph_b = synthetic(s, f_bad, rng=31)

    kw = dict(sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE,
              n_walkers=32, n_steps=4000, seed=5)
    good = bay.sample_posterior(f_good, amp_g, ph_g, s, names, **kw)
    bad = bay.sample_posterior(f_bad, amp_b, ph_b, s, names, **kw)

    # The badly conditioned case is far wider
    assert bad.chain.std(axis=0).max() > 20.0 * good.chain.std(axis=0).max()
    # and markedly more asymmetric
    assert np.abs(bad.asymmetry).max() > 3.0 * np.abs(good.asymmetry).max()


def test_correlation_is_strong_and_negative():
    """Conductivity and heat capacity trade off against each other, since the
    measurement constrains their product far better than their ratio."""
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=41)
    res = bay.sample_posterior(f, amp, ph, s, ["film_lam", "film_rho_c"],
                               n_walkers=32, n_steps=3000, seed=6)
    assert res.correlation[0, 1] < -0.85


def test_summary_is_printable():
    s = truth()
    f = band(s, 0, 1)
    amp, ph = synthetic(s, f, rng=51)
    res = bay.sample_posterior(f, amp, ph, s, ["film_lam", "film_rho_c"],
                               n_walkers=16, n_steps=600, seed=8)
    text = res.summary()
    assert "film_lam" in text and "acceptance" in text
