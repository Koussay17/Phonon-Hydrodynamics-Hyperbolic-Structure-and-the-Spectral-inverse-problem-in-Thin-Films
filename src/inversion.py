"""Parameter estimation from photothermal data, with uncertainty analysis.

Given a measurement and a forward model, recover the sample properties and
quantify how well they are determined.

Parameterisation
----------------
Parameters are fitted in logarithmic space. This is not cosmetic. Thermal
properties span orders of magnitude, and the logarithm makes a relative
perturbation of one parameter comparable to a relative perturbation of
another, which is what a Levenberg-Marquardt step assumes. It also enforces
positivity without a constraint.

More importantly, the choice of parameters carries the identifiability
structure. The modulated response determines the film effusivity at high
frequency and the diffusion time through the transition; conductivity and
heat capacity are recovered only through their combination. Fitting
(effusivity, diffusion time) is therefore far better conditioned than fitting
(conductivity, heat capacity), even though the two pairs describe the same
sample. The condition number reported by `FitResult` measures exactly that.

Residuals
---------
Amplitude noise from a lock-in measurement is multiplicative, phase noise
additive. The residuals are built accordingly: logarithm of the amplitude
divided by the relative noise level, and phase difference divided by the
absolute noise level. Both are then dimensionless and of unit variance, so
the covariance of the estimate is the inverse of J transpose J with no
further scaling.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

import forward_model as fm

__all__ = ["FitResult", "fisher_analysis", "fit_modulated",
           "usable_time_window", "fisher_analysis_pulsed", "fit_pulsed"]


# --------------------------------------------------------------------------
# Result container
# --------------------------------------------------------------------------

@dataclass
class FitResult:
    """Outcome of an estimation."""

    names: list[str]
    values: np.ndarray
    covariance: np.ndarray
    sample: fm.Sample
    chi2: float
    n_data: int
    success: bool
    message: str

    @property
    def std_errors(self) -> np.ndarray:
        """One-sigma uncertainties, in the same units as the parameters."""
        return np.sqrt(np.diag(self.covariance))

    @property
    def relative_errors(self) -> np.ndarray:
        return self.std_errors / np.abs(self.values)

    @property
    def correlation(self) -> np.ndarray:
        """Correlation matrix of the estimate.

        Off-diagonal terms close to plus or minus one mean the corresponding
        parameters are not separately determined: only a combination of them
        is constrained by the data.
        """
        s = self.std_errors
        return self.covariance / np.outer(s, s)

    @property
    def condition_number(self) -> float:
        """Condition number of the Fisher information matrix.

        A large value means the data constrain some directions of parameter
        space far more tightly than others. It is the numerical counterpart
        of the identifiability question.
        """
        fisher = np.linalg.inv(self.covariance)
        return float(np.linalg.cond(fisher))

    @property
    def reduced_chi2(self) -> float:
        dof = self.n_data - len(self.values)
        return self.chi2 / dof if dof > 0 else np.nan

    def summary(self) -> str:
        lines = [
            f"success        : {self.success}",
            f"reduced chi2   : {self.reduced_chi2:.4f}",
            f"condition num. : {self.condition_number:.3e}",
            "",
        ]
        for i, n in enumerate(self.names):
            lines.append(
                f"  {n:<20s} {self.values[i]:12.5e}  "
                f"+/- {self.std_errors[i]:.3e}  "
                f"({100 * self.relative_errors[i]:.2f} %)"
            )
        if len(self.names) > 1:
            lines.append("")
            lines.append("  correlation")
            for i, n in enumerate(self.names):
                row = "  ".join(f"{c:+.3f}" for c in self.correlation[i])
                lines.append(f"  {n:<20s} {row}")
        return "\n".join(lines)


# --------------------------------------------------------------------------
# Residuals
# --------------------------------------------------------------------------

def _apply(sample: fm.Sample, names, log_values) -> fm.Sample:
    """Rebuild a Sample from logarithmic parameter values."""
    return replace(sample, **{n: float(np.exp(v)) for n, v in zip(names, log_values)})


_INVALID = 1e6


def _residuals(log_values, names, sample, freq, amp, phase,
               sigma_rel, sigma_phase):
    """Weighted residuals, guarded against non-physical excursions.

    On a badly conditioned problem the optimiser drifts along the
    unconstrained direction until a parameter underflows to exactly zero, at
    which point the forward model divides by zero. Returning a large finite
    residual instead of raising keeps the optimiser inside the physical region
    and lets it report failure rather than crash.
    """
    s = _apply(sample, names, log_values)
    try:
        a_mod, p_mod = fm.modulated_response(freq, s)
    except (ZeroDivisionError, OverflowError, FloatingPointError):
        return np.full(2 * freq.size, _INVALID)

    if not np.all(np.isfinite(a_mod)) or np.any(a_mod <= 0.0):
        return np.full(2 * freq.size, _INVALID)

    r_amp = (np.log(a_mod) - np.log(amp)) / sigma_rel
    r_phase = (p_mod - phase) / sigma_phase
    out = np.concatenate([r_amp, r_phase])
    return np.where(np.isfinite(out), out, _INVALID)


# --------------------------------------------------------------------------
# Identifiability without fitting
# --------------------------------------------------------------------------

def fisher_analysis(freq, sample: fm.Sample, names,
                    sigma_rel=0.01, sigma_phase=0.1, step=1e-5):
    """Fisher information of a planned measurement.

    Requires no data: it answers what a given frequency range and noise level
    would allow, before the experiment is run. Returns the information matrix
    and the covariance it implies.
    """
    log0 = np.array([np.log(getattr(sample, n)) for n in names])

    def model(lv):
        s = _apply(sample, names, lv)
        a, p = fm.modulated_response(freq, s)
        return np.concatenate([np.log(a) / sigma_rel, p / sigma_phase])

    base = model(log0)
    jac = np.empty((base.size, log0.size))
    for i in range(log0.size):
        up = log0.copy()
        up[i] += step
        jac[:, i] = (model(up) - base) / step

    fisher = jac.T @ jac
    return fisher, np.linalg.inv(fisher)


# --------------------------------------------------------------------------
# Estimation
# --------------------------------------------------------------------------

def fit_modulated(freq, amplitude, phase_deg, initial: fm.Sample, names,
                  sigma_rel=0.01, sigma_phase=0.1, bounds=None,
                  **kwargs) -> FitResult:
    """Estimate parameters from modulated data by Levenberg-Marquardt.

    Parameters
    ----------
    freq : array-like
        Modulation frequencies, in hertz.
    amplitude, phase_deg : array-like
        Measured amplitude and phase.
    initial : Sample
        Starting point. Only the fields named in `names` are adjusted.
    names : sequence of str
        Fields of `Sample` to estimate. They must be strictly positive.
    sigma_rel, sigma_phase : float
        Noise levels: relative on amplitude, absolute in degrees on phase.
    bounds : sequence of (low, high), optional
        Bounds in linear units. Recommended whenever the problem may be badly
        conditioned: without them the optimiser can wander arbitrarily far
        along a direction the data do not constrain.
    """
    names = list(names)
    freq = np.asarray(freq, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    phase_deg = np.asarray(phase_deg, dtype=float)

    for n in names:
        if getattr(initial, n) <= 0.0:
            raise ValueError(f"{n} must be strictly positive to be fitted in log space")

    log0 = np.array([np.log(getattr(initial, n)) for n in names])

    args = (names, initial, freq, amplitude, phase_deg, sigma_rel, sigma_phase)
    if bounds is None:
        out = least_squares(_residuals, log0, args=args, method="lm", **kwargs)
    else:
        lb = np.log(np.asarray(bounds, dtype=float))
        out = least_squares(_residuals, log0, args=args, method="trf",
                            bounds=(lb[:, 0], lb[:, 1]), **kwargs)

    jac = out.jac
    fisher = jac.T @ jac
    try:
        cov_log = np.linalg.inv(fisher)
    except np.linalg.LinAlgError:
        cov_log = np.full_like(fisher, np.inf)

    values = np.exp(out.x)
    # Delta method: a relative uncertainty in log space becomes an absolute
    # uncertainty proportional to the value itself.
    scale = np.outer(values, values)
    cov = cov_log * scale

    return FitResult(
        names=names,
        values=values,
        covariance=cov,
        sample=_apply(initial, names, out.x),
        chi2=float(2.0 * out.cost),
        n_data=int(2 * freq.size),
        success=bool(out.success),
        message=str(out.message),
    )


# --------------------------------------------------------------------------
# Pulsed regime
# --------------------------------------------------------------------------
#
# The pulsed and modulated regimes carry the same information about the
# sample: they are two representations of one forward model, related by an
# integral transform. A rank degeneracy that is independent of the spectral
# parameter is therefore invariant under that transform, which is the
# transposition principle stated by Krapez and Rigollet (2017), section 4.4.
#
# What differs is practical. The pulsed response must be reconstructed by
# numerical inversion, whose accuracy holds over roughly two decades of decay
# (see laplace.py). Choosing the time window is therefore part of the
# estimation, not a detail of it.


def usable_time_window(s: fm.Sample, decades_before=2.0, decades_after=2.0):
    """Time window centred on the diffusion time of the film.

    Outside it the numerical inversion is the limiting factor rather than the
    physics: too early and the response is that of the film alone, too late
    and the signal has decayed past what the inversion can carry.
    """
    t0 = s.diffusion_time
    return t0 * 10.0 ** (-decades_before), t0 * 10.0 ** decades_after


def _pulsed_model(log_values, names, sample, times, n_stehfest):
    s = _apply(sample, names, log_values)
    try:
        out = fm.pulsed_response(times, s, energy=1.0, n_stehfest=n_stehfest)
    except (ZeroDivisionError, OverflowError, FloatingPointError):
        return None
    out = np.asarray(out, dtype=float)
    if not np.all(np.isfinite(out)) or np.any(out <= 0.0):
        return None
    return out


def _residuals_pulsed(log_values, names, sample, times, temperature,
                      sigma_rel, n_stehfest):
    """Weighted residuals in the time domain.

    The noise is taken multiplicative, as it is on an amplitude, so the
    residual is built on the logarithm. A flash measurement spans decades of
    signal level, and an absolute weighting would let the earliest points
    dominate entirely.
    """
    model = _pulsed_model(log_values, names, sample, times, n_stehfest)
    if model is None:
        return np.full(times.size, _INVALID)
    out = (np.log(model) - np.log(temperature)) / sigma_rel
    return np.where(np.isfinite(out), out, _INVALID)


def fisher_analysis_pulsed(times, sample: fm.Sample, names,
                           sigma_rel=0.01, n_stehfest: int = 12, step=1e-4):
    """Fisher information of a planned pulsed measurement.

    Same role as `fisher_analysis` for the modulated regime: it answers what a
    given time window and noise level would allow, before the experiment runs.

    The default differentiation step is larger here than in the modulated
    case, and deliberately so. The numerical inversion carries a relative
    noise floor near 1e-6, which a central difference of step h amplifies by
    1/(2h). Refining the step therefore degrades the derivative instead of
    improving it: measured on the Euler combination, the residual falls to
    1.5e-4 at h = 1e-4 and rises to 1 at h = 1e-7, where the derivative is
    pure noise. The same calculation in the frequency domain reaches 4e-10.

    The practical consequence is that a Jacobian obtained through the
    inversion is noisy, and that this propagates into the covariance. Pulsed
    error bars deserve more caution than modulated ones, and the agreement
    between the predicted and the fitted covariance is correspondingly
    looser.
    """
    times = np.asarray(times, dtype=float)
    log0 = np.array([np.log(getattr(sample, n)) for n in names])

    def model(lv):
        out = _pulsed_model(lv, names, sample, times, n_stehfest)
        if out is None:
            raise ValueError("the forward model is not defined at this point")
        return np.log(out) / sigma_rel

    base = model(log0)
    jac = np.empty((base.size, log0.size))
    for i in range(log0.size):
        up = log0.copy()
        up[i] += step
        jac[:, i] = (model(up) - base) / step

    fisher = jac.T @ jac
    return fisher, np.linalg.inv(fisher)


def fit_pulsed(times, temperature, initial: fm.Sample, names,
               sigma_rel=0.01, n_stehfest: int = 12, bounds=None,
               diff_step=1e-4, **kwargs) -> FitResult:
    """Estimate parameters from a pulsed measurement.

    Parameters
    ----------
    times : array-like
        Measurement times, strictly positive. `usable_time_window` gives a
        sensible range for a given sample.
    temperature : array-like
        Measured front-face temperature, strictly positive.
    n_stehfest : int
        Number of terms of the numerical inversion. Twelve is the practical
        optimum in double precision; raising it degrades rather than improves.
    diff_step : float
        Relative step of the numerical Jacobian. Chosen larger than the
        default for the reason given in `fisher_analysis_pulsed`: below about
        1e-5 the derivative is dominated by the noise of the inversion.
    """
    names = list(names)
    times = np.asarray(times, dtype=float)
    temperature = np.asarray(temperature, dtype=float)

    for n in names:
        if getattr(initial, n) <= 0.0:
            raise ValueError(f"{n} must be strictly positive to be fitted in log space")

    log0 = np.array([np.log(getattr(initial, n)) for n in names])
    args = (names, initial, times, temperature, sigma_rel, n_stehfest)

    if bounds is None:
        out = least_squares(_residuals_pulsed, log0, args=args,
                            method="lm", diff_step=diff_step, **kwargs)
    else:
        lb = np.log(np.asarray(bounds, dtype=float))
        out = least_squares(_residuals_pulsed, log0, args=args, method="trf",
                            bounds=(lb[:, 0], lb[:, 1]),
                            diff_step=diff_step, **kwargs)

    jac = out.jac
    fisher = jac.T @ jac
    try:
        cov_log = np.linalg.inv(fisher)
    except np.linalg.LinAlgError:
        cov_log = np.full_like(fisher, np.inf)

    values = np.exp(out.x)
    cov = cov_log * np.outer(values, values)

    return FitResult(
        names=names,
        values=values,
        covariance=cov,
        sample=_apply(initial, names, out.x),
        chi2=float(2.0 * out.cost),
        n_data=int(times.size),
        success=bool(out.success),
        message=str(out.message),
    )
