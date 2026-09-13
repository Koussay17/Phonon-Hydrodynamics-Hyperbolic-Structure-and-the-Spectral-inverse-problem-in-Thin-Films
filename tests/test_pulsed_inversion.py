"""Parameter estimation in the pulsed regime.

The test that matters most here is `test_euler_relation_holds_in_time_domain`.
Krapez and Rigollet (2017), section 4.4, state a transposition principle: a
combination of parameters that is unidentifiable in the frequency domain for
every frequency is equally unidentifiable in the time domain for every
instant, and conversely. The reason is that both regimes are representations
of one forward model, so a rank degeneracy independent of the spectral
parameter survives the transform.

That principle is stated, not proved, in the note. It is checked here.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import inversion as inv  # noqa: E402

SIGMA_REL = 0.01

L0 = 500e-9
KAPPA0 = 60.0
ALPHA0 = KAPPA0 / 2.41e6


def truth():
    return fm.Sample(film_lam=KAPPA0, film_rho_c=KAPPA0 / ALPHA0, thickness=L0)


def times(n=40):
    t0, t1 = inv.usable_time_window(truth())
    return np.logspace(np.log10(t0), np.log10(t1), n)


def synthetic(s, t, rng=None):
    temperature = fm.pulsed_response(t, s)
    if rng is None:
        return temperature
    g = np.random.default_rng(rng)
    return temperature * (1.0 + SIGMA_REL * g.standard_normal(temperature.shape))


# --------------------------------------------------------------------------
# The window
# --------------------------------------------------------------------------

def test_window_brackets_the_diffusion_time():
    s = truth()
    t0, t1 = inv.usable_time_window(s)
    assert t0 < s.diffusion_time < t1
    assert t1 / t0 == pytest.approx(1e4)


def test_response_decreases_over_the_window():
    assert np.all(np.diff(synthetic(truth(), times())) < 0.0)


# --------------------------------------------------------------------------
# Recovery
# --------------------------------------------------------------------------

def test_noiseless_recovery_is_exact():
    t = times()
    data = synthetic(truth(), t)
    start = fm.Sample(film_lam=25.0, film_rho_c=1.0e6, thickness=L0)

    res = inv.fit_pulsed(t, data, start, names=["film_lam", "film_rho_c"],
                         sigma_rel=SIGMA_REL)
    assert res.success
    assert res.values[0] == pytest.approx(KAPPA0, rel=1e-5)
    assert res.values[1] == pytest.approx(KAPPA0 / ALPHA0, rel=1e-5)


def test_noisy_recovery_stays_within_error_bars():
    t = times()
    data = synthetic(truth(), t, rng=5)
    res = inv.fit_pulsed(t, data, truth(), names=["film_lam", "film_rho_c"],
                         sigma_rel=SIGMA_REL)
    for got, true, err in zip(res.values, (KAPPA0, KAPPA0 / ALPHA0),
                              res.std_errors):
        assert abs(got - true) < 4.0 * err


def test_rejects_non_positive_parameter():
    t = times()
    data = synthetic(truth(), t)
    with pytest.raises(ValueError):
        inv.fit_pulsed(t, data, fm.Sample(contact_resistance=0.0),
                       names=["contact_resistance"])


def test_degenerate_case_does_not_crash():
    """Same guard as in the modulated regime: an excursion towards a vanishing
    parameter must be reported, not raised."""
    t = times()
    data = synthetic(truth(), t, rng=6)
    res = inv.fit_pulsed(t, data, truth(),
                         names=["film_lam", "film_rho_c", "thickness"],
                         sigma_rel=SIGMA_REL)
    assert np.all(np.isfinite(res.values))


# --------------------------------------------------------------------------
# The transposition principle
# --------------------------------------------------------------------------

def _euler_residual_time(t, h=1e-4):
    """Relative residual of  L dG/dL + 2 alpha dG/dalpha + kappa dG/dkappa,
    evaluated on the pulsed response."""
    def g(L, alpha, kappa):
        s = fm.Sample(film_lam=kappa, film_rho_c=kappa / alpha, thickness=L)
        return np.asarray(fm.pulsed_response(t, s), dtype=float)

    def d(i):
        args = [L0, ALPHA0, KAPPA0]
        up, dn = list(args), list(args)
        up[i] *= 1.0 + h
        dn[i] *= 1.0 - h
        return (g(*up) - g(*dn)) / (2.0 * h * args[i])

    terms = [L0 * d(0), 2.0 * ALPHA0 * d(1), KAPPA0 * d(2)]
    total = sum(terms)
    scale = sum(np.abs(x) for x in terms)
    return float(np.max(np.abs(total) / scale))


def test_euler_relation_holds_in_time_domain():
    """The degeneracy transposes to the time domain.

    Checks the transposition principle of Krapez and Rigollet, section 4.4.
    The residual is far larger than in the frequency domain, not because the
    principle holds less well but because the numerical inversion limits how
    accurately any derivative can be formed. See the next test.
    """
    assert _euler_residual_time(times(20)) < 1e-3


def test_inversion_noise_floors_the_derivatives():
    """Refining the differentiation step degrades the derivative rather than
    improving it, which is the signature of a noise floor.

    The inversion carries a relative noise near 1e-6; a central difference of
    step h amplifies it by 1/(2h). The optimum therefore sits at a step far
    coarser than usual, and below it the derivative is pure noise.

    This bounds the accuracy of every Jacobian, hence of every covariance,
    obtained in the pulsed regime.
    """
    t = times(20)
    coarse = _euler_residual_time(t, h=1e-4)
    fine = _euler_residual_time(t, h=1e-7)

    assert coarse < 1e-3            # usable
    assert fine > 0.1               # pure noise
    assert fine > 100.0 * coarse


def test_thickness_alone_is_not_identifiable_either():
    """A corollary: with the three parameters free, the covariance blows up,
    exactly as it does in the frequency domain."""
    t = times(20)
    s = truth()
    _, cov = inv.fisher_analysis_pulsed(t, s,
                                        ["thickness", "film_lam", "film_rho_c"],
                                        SIGMA_REL)
    assert np.sqrt(np.abs(np.diag(cov))).max() > 1.0     # worse than 100 per cent


def test_two_parameters_are_identifiable_once_thickness_is_known():
    t = times(40)
    _, cov = inv.fisher_analysis_pulsed(t, truth(),
                                        ["film_lam", "film_rho_c"], SIGMA_REL)
    assert np.sqrt(np.diag(cov)).max() < 0.1            # better than ten per cent


# --------------------------------------------------------------------------
# Are the error bars real
# --------------------------------------------------------------------------

def test_fisher_uncertainties_match_monte_carlo():
    """Fewer repetitions than in the modulated case, the numerical inversion
    being far more costly, but the same check."""
    t = times(30)
    n_runs = 40
    estimates = np.empty((n_runs, 2))
    predicted = None

    for i in range(n_runs):
        data = synthetic(truth(), t, rng=3000 + i)
        res = inv.fit_pulsed(t, data, truth(),
                             names=["film_lam", "film_rho_c"],
                             sigma_rel=SIGMA_REL)
        estimates[i] = res.values
        if predicted is None:
            predicted = res.std_errors

    observed = estimates.std(axis=0, ddof=1)
    for obs, pred in zip(observed, predicted):
        assert obs == pytest.approx(pred, rel=0.5)


def test_fisher_prediction_matches_the_fit():
    t = times()
    names = ["film_lam", "film_rho_c"]
    _, cov_log = inv.fisher_analysis_pulsed(t, truth(), names, SIGMA_REL)
    values = np.array([getattr(truth(), n) for n in names])
    predicted = np.sqrt(np.diag(cov_log)) * values

    res = inv.fit_pulsed(t, synthetic(truth(), t), truth(), names=names,
                         sigma_rel=SIGMA_REL)
    # Looser than in the modulated regime, and for a reason: both covariances
    # rest on Jacobians formed through the numerical inversion, whose noise
    # floor limits their accuracy.
    assert np.allclose(predicted, res.std_errors, rtol=0.25)
