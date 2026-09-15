"""Posterior sampling of thermal parameters by Markov chain Monte Carlo.

Why this exists alongside `inversion.py`
----------------------------------------
Levenberg-Marquardt returns an estimate and a covariance obtained by
linearising the model around the optimum. That covariance describes an
ellipse: symmetric, Gaussian, fully specified by a width and a correlation.

The actual posterior is only an ellipse when the model is close to linear
over the region the data allow. For a film on a substrate it is not. The
measurement constrains the effusivity tightly and the diffusivity loosely,
which produces a long curved ridge in parameter space. An ellipse fitted at
the summit of that ridge understates the uncertainty on the poorly determined
direction and hides its asymmetry.

Sampling the posterior shows the ridge as it is. The price is computation
    time; useful intervals require convergence and prior-sensitivity checks.

Parameterisation and priors
---------------------------
Parameters are sampled in logarithmic space, as in `inversion.py`. The prior
is uniform in the logarithm within the given bounds. This is an explicit
prior choice; its bounds can control results along a weakly constrained ridge.

The likelihood uses the same residuals as the least-squares estimation, so
the two approaches differ only in how they explore parameter space, never in
what they assume about the data.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

import forward_model as fm

__all__ = ["PosteriorResult", "sample_posterior"]


# --------------------------------------------------------------------------
# Result container
# --------------------------------------------------------------------------

@dataclass
class PosteriorResult:
    """Outcome of a posterior sampling."""

    names: list[str]
    chain: np.ndarray            # shape (n_samples, n_params), linear units
    acceptance: float
    autocorr_time: np.ndarray | None
    autocorr_reliable: bool = False

    @property
    def median(self) -> np.ndarray:
        return np.median(self.chain, axis=0)

    def credible_interval(self, level: float = 0.68) -> np.ndarray:
        """Percentile interval containing `level` of the posterior mass.

        Returns an array of shape (n_params, 2). The interval is not forced to
        be symmetric about the median, which is the whole point.
        """
        if not 0 < level < 1:
            raise ValueError("level must lie strictly between zero and one")
        lo = 50.0 * (1.0 - level)
        hi = 100.0 - lo
        return np.percentile(self.chain, [lo, hi], axis=0).T

    @property
    def asymmetry(self) -> np.ndarray:
        """Relative difference between the upper and lower halves of the 68 %
        interval.

        Zero for a symmetric posterior. A large value means the Gaussian error
        bar of the linearised estimate is misleading.
        """
        ci = self.credible_interval(0.68)
        med = self.median
        upper = ci[:, 1] - med
        lower = med - ci[:, 0]
        return (upper - lower) / (upper + lower)

    @property
    def correlation(self) -> np.ndarray:
        return np.corrcoef(np.log(self.chain).T)

    def summary(self) -> str:
        ci = self.credible_interval(0.68)
        med = self.median
        lines = [
            f"acceptance     : {self.acceptance:.3f}",
            f"autocorrelation estimate reliable: {self.autocorr_reliable}",
            "",
        ]
        for i, n in enumerate(self.names):
            lines.append(
                f"  {n:<18s} {med[i]:11.5e}   "
                f"+{ci[i, 1] - med[i]:.3e} / -{med[i] - ci[i, 0]:.3e}   "
                f"(asym {100 * self.asymmetry[i]:+.1f} %)"
            )
        return "\n".join(lines)


# --------------------------------------------------------------------------
# Sampling
# --------------------------------------------------------------------------

def _log_probability(log_values, names, sample, freq, amp, phase,
                     sigma_rel, sigma_phase, log_bounds):
    if np.any(log_values < log_bounds[:, 0]) or np.any(log_values > log_bounds[:, 1]):
        return -np.inf

    try:
        s = replace(sample, **{n: float(np.exp(v)) for n, v in zip(names, log_values)})
        a_mod, p_mod = fm.modulated_response(freq, s)
    except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
        return -np.inf
    if not np.all(np.isfinite(a_mod)) or np.any(a_mod <= 0.0):
        return -np.inf

    r_amp = (np.log(a_mod) - np.log(amp)) / sigma_rel
    r_phase = ((p_mod - phase + 180) % 360 - 180) / sigma_phase
    return -0.5 * float(np.sum(r_amp ** 2) + np.sum(r_phase ** 2))


def sample_posterior(freq, amplitude, phase_deg, initial: fm.Sample, names,
                     bounds=None, sigma_rel=0.01, sigma_phase=0.1,
                     n_walkers=32, n_steps=3000, n_burn=None, seed=None,
                     progress=False) -> PosteriorResult:
    """Sample the posterior of `names` given modulated data.

    Parameters
    ----------
    bounds : sequence of (low, high), optional
        Prior support in linear units. Defaults to a factor of 100 either side
        of the starting value. These bounds are part of the prior and may
        materially affect a weakly identified direction; check sensitivity.
    n_burn : int, optional
        Steps discarded at the start. Defaults to a fifth of `n_steps`.
    """
    import emcee

    names = list(names)
    freq = np.asarray(freq, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    phase_deg = np.asarray(phase_deg, dtype=float)

    start = np.array([getattr(initial, n) for n in names], dtype=float)
    if not names or len(set(names)) != len(names) or np.any(~np.isfinite(start)) or np.any(start <= 0.0):
        raise ValueError("all sampled parameters must be strictly positive")

    if bounds is None:
        bounds = [(v / 100.0, v * 100.0) for v in start]
    bounds_array = np.asarray(bounds, dtype=float)
    if (bounds_array.shape != (len(names), 2) or
        np.any(~np.isfinite(bounds_array)) or np.any(bounds_array <= 0) or
        np.any(bounds_array[:, 0] >= bounds_array[:, 1])):
        raise ValueError("bounds must be finite, positive and ordered")
    if np.any(start <= bounds_array[:, 0]) or np.any(start >= bounds_array[:, 1]):
        raise ValueError("starting values must be strictly inside the prior bounds")
    if sigma_rel <= 0 or sigma_phase <= 0:
        raise ValueError("noise levels must be positive")
    log_bounds = np.log(bounds_array)

    n_dim = len(names)
    if n_burn is None:
        n_burn = n_steps // 5

    if not 0 <= n_burn < n_steps:
        raise ValueError("require 0 <= n_burn < n_steps")
    rng = np.random.default_rng(seed)
    p0 = np.log(start) + 1e-3 * rng.standard_normal((n_walkers, n_dim))

    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, _log_probability,
        args=(names, initial, freq, amplitude, phase_deg,
              sigma_rel, sigma_phase, log_bounds),
    )
    if seed is not None:
        sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, n_steps, progress=progress)

    flat = sampler.get_chain(discard=n_burn, flat=True)

    try:
        tau = sampler.get_autocorr_time(discard=n_burn, quiet=False)
        reliable = bool(np.all(np.isfinite(tau)))
    except emcee.autocorr.AutocorrError as exc:
        tau = exc.tau
        reliable = False

    return PosteriorResult(
        names=names,
        chain=np.exp(flat),
        acceptance=float(np.mean(sampler.acceptance_fraction)),
        autocorr_time=tau,
        autocorr_reliable=reliable,
    )
