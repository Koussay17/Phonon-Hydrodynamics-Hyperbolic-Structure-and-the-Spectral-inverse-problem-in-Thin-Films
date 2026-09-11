"""Unit tests for src/inversion.py.

Two tests carry the weight.

`test_noiseless_recovery_is_exact` checks that the estimator finds the truth
when there is nothing to obscure it. If this fails, nothing else matters.

`test_fisher_uncertainties_match_monte_carlo` checks that the reported
uncertainties are real. A fit that returns the right value with a wrong error
bar is worse than useless, because it looks trustworthy. The only honest way
to verify an error bar is to repeat the experiment many times and compare the
observed spread with the predicted one.
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


def truth():
    return fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9)


def wide_band(sample, decades_below=2.0, decades_above=2.0, n=60):
    """Frequency range centred on the characteristic frequency of the film."""
    fc = sample.characteristic_frequency
    return np.logspace(np.log10(fc) - decades_below,
                       np.log10(fc) + decades_above, n)


def band(sample, lo, hi, n=60):
    """Frequency band, in decades relative to the characteristic frequency.

    band(s, 0, 1) covers fc to 10 fc; band(s, -3, -2) covers fc/1000 to fc/100.
    """
    fc = np.log10(sample.characteristic_frequency)
    return np.logspace(fc + lo, fc + hi, n)


def synthetic(sample, freq, rng=None):
    amp, ph = fm.modulated_response(freq, sample)
    if rng is None:
        return amp, ph
    return fm.add_noise(amp, ph, relative_amplitude=SIGMA_REL,
                        absolute_phase_deg=SIGMA_PHASE, rng=rng)


# --------------------------------------------------------------------------
# Recovery
# --------------------------------------------------------------------------

def test_noiseless_recovery_is_exact():
    """With clean data the estimator must return the truth to high precision."""
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f)

    start = fm.Sample(film_lam=25.0, film_rho_c=1.0e6, thickness=500e-9)
    res = inv.fit_modulated(f, amp, ph, start,
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)

    assert res.success
    assert res.values[0] == pytest.approx(60.0, rel=1e-6)
    assert res.values[1] == pytest.approx(2.41e6, rel=1e-6)
    assert res.chi2 < 1e-12


def test_recovery_from_a_distant_start():
    """Levenberg-Marquardt in log space must tolerate a starting point off by
    a factor of several."""
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f)

    start = fm.Sample(film_lam=6.0, film_rho_c=1.0e7, thickness=500e-9)
    res = inv.fit_modulated(f, amp, ph, start,
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    assert res.success
    assert res.values[0] == pytest.approx(60.0, rel=1e-5)
    assert res.values[1] == pytest.approx(2.41e6, rel=1e-5)


def test_noisy_recovery_stays_within_error_bars():
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f, rng=7)

    res = inv.fit_modulated(f, amp, ph, truth(),
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    assert res.success
    for got, true, err in zip(res.values, (60.0, 2.41e6), res.std_errors):
        assert abs(got - true) < 4.0 * err


def test_reduced_chi2_is_about_one_on_noisy_data():
    """A reduced chi-square far from one means the noise model is wrong, not
    that the fit is good or bad."""
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f, rng=3)

    res = inv.fit_modulated(f, amp, ph, truth(),
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    assert 0.5 < res.reduced_chi2 < 2.0


# --------------------------------------------------------------------------
# Are the error bars real
# --------------------------------------------------------------------------

def test_fisher_uncertainties_match_monte_carlo():
    """Repeat the experiment many times and compare the observed spread of the
    estimates with the uncertainty predicted from the Fisher matrix.

    This is the only honest validation of an error bar.
    """
    s = truth()
    f = wide_band(s)

    n_runs = 120
    estimates = np.empty((n_runs, 2))
    predicted = None

    for i in range(n_runs):
        amp, ph = synthetic(s, f, rng=1000 + i)
        res = inv.fit_modulated(f, amp, ph, truth(),
                                names=["film_lam", "film_rho_c"],
                                sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
        estimates[i] = res.values
        if predicted is None:
            predicted = res.std_errors

    observed = estimates.std(axis=0, ddof=1)
    for obs, pred in zip(observed, predicted):
        assert obs == pytest.approx(pred, rel=0.25)


def test_estimator_is_unbiased():
    """The mean of many estimates must sit on the truth, within the error of
    the mean."""
    s = truth()
    f = wide_band(s)

    n_runs = 120
    estimates = np.empty((n_runs, 2))
    for i in range(n_runs):
        amp, ph = synthetic(s, f, rng=2000 + i)
        res = inv.fit_modulated(f, amp, ph, truth(),
                                names=["film_lam", "film_rho_c"],
                                sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
        estimates[i] = res.values

    mean = estimates.mean(axis=0)
    sem = estimates.std(axis=0, ddof=1) / np.sqrt(n_runs)
    for m, true, e in zip(mean, (60.0, 2.41e6), sem):
        assert abs(m - true) < 4.0 * e


# --------------------------------------------------------------------------
# Identifiability
# --------------------------------------------------------------------------

def test_band_far_above_transition_loses_the_diffusivity():
    """Once the thermal wave is confined to the film, only the effusivity is
    measured, and conductivity and heat capacity become degenerate.

    Measured on the default sample: the condition number of the Fisher matrix
    reaches 2e11 two decades above the characteristic frequency, against 35 in
    the optimal band. The correlation between the two parameters reaches -1.
    """
    s = truth()
    names = ["film_lam", "film_rho_c"]

    best = band(s, 0, 1)                        # fc to 10 fc
    far = band(s, 2, 3)                         # 100 fc to 1000 fc

    f_best, cov_best = inv.fisher_analysis(best, s, names, SIGMA_REL, SIGMA_PHASE)
    f_far, cov_far = inv.fisher_analysis(far, s, names, SIGMA_REL, SIGMA_PHASE)

    assert np.linalg.cond(f_far) > 1e6 * np.linalg.cond(f_best)
    assert np.sqrt(np.diag(cov_far)).max() > 100.0 * np.sqrt(np.diag(cov_best)).max()

    r = cov_far[0, 1] / np.sqrt(cov_far[0, 0] * cov_far[1, 1])
    assert r < -0.9999


def test_degenerate_direction_is_the_diffusivity():
    """The eigenvectors of the Fisher matrix name what the measurement
    determines.

    Far above the transition the well-determined direction is (1, 1) in
    logarithmic coordinates, that is log(lambda) + log(rho c) = 2 log(b): the
    effusivity. The orthogonal direction, the diffusivity, carries almost no
    information.

    This is why fitting (effusivity, diffusion time) is far better conditioned
    than fitting (conductivity, heat capacity): the former are the eigen-
    directions of the problem, the latter are not.
    """
    s = truth()
    far = band(s, 2, 3)
    fisher, _ = inv.fisher_analysis(far, s, ["film_lam", "film_rho_c"],
                                    SIGMA_REL, SIGMA_PHASE)

    w, v = np.linalg.eigh(fisher)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]

    strong = v[:, 0]
    assert abs(abs(strong[0]) - abs(strong[1])) < 0.02      # direction (1, 1)
    assert strong[0] * strong[1] > 0                        # same sign
    assert w[0] / w[1] > 1e6                                # one direction only


def test_optimum_sits_just_above_the_characteristic_frequency():
    """Scanning one-decade bands, the smallest uncertainties are obtained just
    above the characteristic frequency, not far from it in either direction.
    """
    s = truth()
    names = ["film_lam", "film_rho_c"]

    worst_low = band(s, -3, -2)                 # far below
    best = band(s, 0, 1)                        # just above
    worst_high = band(s, 2, 3)                  # far above

    def worst_sigma(band):
        _, cov = inv.fisher_analysis(band, s, names, SIGMA_REL, SIGMA_PHASE)
        return np.sqrt(np.diag(cov)).max()

    assert worst_sigma(best) < worst_sigma(worst_low)
    assert worst_sigma(best) < worst_sigma(worst_high)


def test_fisher_prediction_matches_the_fit():
    """The covariance predicted before the experiment must agree with the one
    obtained after it."""
    s = truth()
    f = wide_band(s)
    names = ["film_lam", "film_rho_c"]

    _, cov_log = inv.fisher_analysis(f, s, names, SIGMA_REL, SIGMA_PHASE)
    values = np.array([getattr(s, n) for n in names])
    predicted = np.sqrt(np.diag(cov_log)) * values

    amp, ph = synthetic(s, f)
    res = inv.fit_modulated(f, amp, ph, s, names=names,
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)

    assert np.allclose(predicted, res.std_errors, rtol=1e-3)


def test_correlation_is_symmetric_and_unit_diagonal():
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f)
    res = inv.fit_modulated(f, amp, ph, s,
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    c = res.correlation
    assert np.allclose(np.diag(c), 1.0)
    assert np.allclose(c, c.T)
    assert np.all(np.abs(c) <= 1.0 + 1e-12)


def test_rejects_non_positive_parameter():
    s = fm.Sample(contact_resistance=0.0)
    f = wide_band(s)
    amp, ph = synthetic(s, f)
    with pytest.raises(ValueError):
        inv.fit_modulated(f, amp, ph, s, names=["contact_resistance"])


def test_summary_is_printable():
    s = truth()
    f = wide_band(s)
    amp, ph = synthetic(s, f)
    res = inv.fit_modulated(f, amp, ph, s,
                            names=["film_lam", "film_rho_c"],
                            sigma_rel=SIGMA_REL, sigma_phase=SIGMA_PHASE)
    text = res.summary()
    assert "film_lam" in text and "correlation" in text
