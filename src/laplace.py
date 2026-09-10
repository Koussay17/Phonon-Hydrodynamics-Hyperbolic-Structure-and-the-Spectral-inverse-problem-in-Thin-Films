"""Numerical inversion of the Laplace transform by the Gaver-Stehfest method.

The transformed quantities produced by `quadrupoles.py` live in the Laplace
domain. This module brings them back to the time domain, which is what a
pulsed photothermal experiment measures.

Method
------
    f(t)  ~=  (ln 2 / t) * sum_{k=1..N} V_k * F(k ln2 / t)

The weights V_k depend only on N, which must be even. The transform F is
evaluated at *real* positive arguments only.

Limitations
-----------
The weights alternate in sign and grow rapidly with N, so the sum suffers
from catastrophic cancellation in double precision. Accuracy improves with N
up to roughly N = 12, then degrades.

Measured relative error on exp(-a t) at N = 12, as a function of how far the
signal has decayed from its maximum:

    above 1e-1 of the peak   ->  ~5e-4
    above 1e-2               ->  ~7e-3
    above 1e-3               ->  ~1e-1

Practical rule: the inversion holds over roughly two decades of decay. A fit
performed on a tail beyond that is driven by numerical noise rather than by
the physics. This bears directly on experiment design: acquiring a long
decaying tail does not add usable information if the model is evaluated
through this inversion.

The method also assumes a smooth, non-oscillatory f(t): step discontinuities
away from the origin and oscillations are not recovered. For those cases the
De Hoog method, based on Fourier series with epsilon acceleration, is the
usual replacement. Krapez (2018) uses De Hoog.

Reference: Maillet, Andre, Batsale, Degiovanni, Moyne, *Thermal Quadrupoles*,
Wiley (2000), chapter 9, section 9.3.2 and appendix 1.1.
"""

from __future__ import annotations

from math import factorial, log

import numpy as np

__all__ = [
    "stehfest_coefficients",
    "stehfest_inverse",
    "semi_infinite_pulse_response",
]

_LN2 = log(2.0)
_CACHE: dict[int, np.ndarray] = {}


def stehfest_coefficients(n: int = 12) -> np.ndarray:
    """Stehfest weights V_k for k = 1..n. `n` must be even.

    The weights sum to zero for n >= 2, which is a useful sanity check: it is
    what makes the method exact on constants.
    """
    if n % 2 or n < 2:
        raise ValueError("n must be an even integer >= 2")
    if n in _CACHE:
        return _CACHE[n]

    half = n // 2
    v = np.zeros(n, dtype=float)
    for k in range(1, n + 1):
        total = 0.0
        for j in range((k + 1) // 2, min(k, half) + 1):
            total += (
                j ** half
                * factorial(2 * j)
                / (
                    factorial(half - j)
                    * factorial(j)
                    * factorial(j - 1)
                    * factorial(k - j)
                    * factorial(2 * j - k)
                )
            )
        v[k - 1] = (-1.0) ** (k + half) * total

    _CACHE[n] = v
    return v


def stehfest_inverse(f_laplace, t, n: int = 12):
    """Invert `f_laplace` at time(s) `t`.

    Parameters
    ----------
    f_laplace : callable
        Function of a single real positive argument p, returning a real value
        or an array of real values.
    t : float or array-like
        Times at which the original function is wanted. Must be strictly
        positive.
    n : int
        Number of terms, even. 12 is a reasonable default in double precision.
    """
    t = np.asarray(t, dtype=float)
    if np.any(t <= 0.0):
        raise ValueError("t must be strictly positive")

    v = stehfest_coefficients(n)
    scalar = t.ndim == 0
    t = np.atleast_1d(t)

    out = np.zeros_like(t)
    for i, ti in enumerate(t):
        a = _LN2 / ti
        s = 0.0
        for k in range(1, n + 1):
            s += v[k - 1] * np.real(f_laplace(k * a))
        out[i] = a * s

    return float(out[0]) if scalar else out


def semi_infinite_pulse_response(t, b: float, energy: float = 1.0):
    """Front-face temperature of a semi-infinite medium after a Dirac pulse.

        theta(t) = energy / (b sqrt(pi t))

    Closed-form counterpart of theta(p) = energy / (b sqrt(p)). Used as the
    reference case for validating the inversion.
    """
    t = np.asarray(t, dtype=float)
    return energy / (b * np.sqrt(np.pi * t))
