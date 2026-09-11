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

__all__ = ["FitResult", "fisher_analysis", "fit_modulated"]


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


def _residuals(log_values, names, sample, freq, amp, phase,
               sigma_rel, sigma_phase):
    s = _apply(sample, names, log_values)
    a_mod, p_mod = fm.modulated_response(freq, s)
    r_amp = (np.log(a_mod) - np.log(amp)) / sigma_rel
    r_phase = (p_mod - phase) / sigma_phase
    return np.concatenate([r_amp, r_phase])


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
                  sigma_rel=0.01, sigma_phase=0.1, **kwargs) -> FitResult:
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
    """
    names = list(names)
    freq = np.asarray(freq, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    phase_deg = np.asarray(phase_deg, dtype=float)

    for n in names:
        if getattr(initial, n) <= 0.0:
            raise ValueError(f"{n} must be strictly positive to be fitted in log space")

    log0 = np.array([np.log(getattr(initial, n)) for n in names])

    out = least_squares(
        _residuals, log0,
        args=(names, initial, freq, amplitude, phase_deg, sigma_rel, sigma_phase),
        method="lm", **kwargs,
    )

    jac = out.jac
    fisher = jac.T @ jac
    cov_log = np.linalg.inv(fisher)

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
