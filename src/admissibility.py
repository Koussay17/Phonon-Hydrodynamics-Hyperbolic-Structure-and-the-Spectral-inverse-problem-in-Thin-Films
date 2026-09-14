"""Admissibility of the non-Fourier conduction laws.

Two criteria are usually conflated and are in fact independent: whether a law
admits a convex entropy with non-negative production, and whether it
propagates at finite speed. The results below show that all three laws pass
the first and that only one passes the second.

Dispersion
----------
Combining the constitutive law with the energy balance gives, for a plane wave
in a homogeneous medium,

    sigma^2 = p (1 + tau_R p) / [ a (1 + tau_l p) ],    p = i omega,

with sigma = i k. The high-frequency behaviour of k decides the propagation
speed:

    Fourier            k ~ sqrt(omega)     phase velocity unbounded
    Cattaneo           k ~ omega           phase velocity -> sqrt(a / tau_R)
    Guyer-Krumhansl    k ~ sqrt(omega)     phase velocity unbounded

The nonlocal term is diffusive in the flux itself. Above the frequency where
it dominates the term in tau_R, the principal part of the system becomes
parabolic again and the finite speed that Cattaneo had restored is lost.

Entropy
-------
Extended irreversible thermodynamics takes the flux as an independent state
variable,

    s(u, q) = s_eq(u) - (tau_R / 2 lambda T0^2) q^2 ,

which is concave provided tau_R > 0. With the classical entropy flux q / T the
production of Cattaneo is

    sigma_s = q^2 / (lambda T0^2) .

Guyer-Krumhansl requires an extended entropy flux,

    J_s = q / T + (3 ell^2 / lambda T0^2) q dq/dz ,

without which the production is not a sum of squares. With it,

    sigma_s = [ q^2 + 3 ell^2 (dq/dz)^2 ] / (lambda T0^2) .

Both are non-negative for lambda > 0 and ell^2 >= 0.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "wavenumber",
    "phase_velocity",
    "cattaneo_speed",
    "characteristic_matrix",
    "characteristic_speeds",
    "entropy_flux_coefficient",
    "entropy_production",
    "admissibility_table",
]


# --------------------------------------------------------------------------
# Dispersion
# --------------------------------------------------------------------------

def wavenumber(omega, a: float, tau_r: float = 0.0, tau_l: float = 0.0):
    """Complex wavenumber of a plane thermal wave.

    Returns k such that the field varies as exp(i k z), with
    sigma = i k and sigma^2 = p (1 + tau_R p) / [a (1 + tau_l p)].
    """
    p = 1j * np.asarray(omega, dtype=float)
    sigma = np.sqrt(p * (1.0 + tau_r * p) / (a * (1.0 + tau_l * p)))
    return -1j * sigma


def phase_velocity(omega, a: float, tau_r: float = 0.0, tau_l: float = 0.0):
    """Phase velocity, omega divided by the real part of the wavenumber."""
    k = wavenumber(omega, a, tau_r, tau_l)
    return np.asarray(omega, dtype=float) / np.abs(np.real(k))


def cattaneo_speed(a: float, tau_r: float) -> float:
    """Limiting propagation speed of the Cattaneo law, sqrt(a / tau_R).

    Also called the second sound velocity. Equal to v / sqrt(3) in the Debye
    model, since a = v^2 tau_R / 3.
    """
    if tau_r <= 0.0:
        return np.inf
    return float(np.sqrt(a / tau_r))


# --------------------------------------------------------------------------
# Characteristic analysis of the first-order system
# --------------------------------------------------------------------------

def characteristic_matrix(a: float, tau_r: float) -> np.ndarray:
    """Principal part of the Cattaneo system in the variables (u, q).

        du/dt + dq/dz = 0
        tau_R dq/dt + (a / tau_R is carried here) du/dz = 0

    The matrix is not symmetric but has two real distinct eigenvalues, so the
    system is strictly hyperbolic. Godunov and Mock's theorem then relates
    that property to the existence of a convex entropy.
    """
    if tau_r <= 0.0:
        raise ValueError("tau_R must be strictly positive for the system to be hyperbolic")
    return np.array([[0.0, 1.0], [a / tau_r, 0.0]])


def characteristic_speeds(a: float, tau_r: float) -> np.ndarray:
    """Eigenvalues of the principal part, plus and minus sqrt(a / tau_R)."""
    return np.linalg.eigvals(characteristic_matrix(a, tau_r))


# --------------------------------------------------------------------------
# Entropy
# --------------------------------------------------------------------------

def entropy_flux_coefficient(lam: float, t0: float, ell_sq: float) -> float:
    """Coefficient of the extra entropy flux term q dq/dz.

    Equal to 3 ell^2 / (lambda T0^2). It is what makes the production of
    Guyer-Krumhansl a sum of squares; with the classical flux q / T alone a
    cross term of indefinite sign survives.
    """
    return 3.0 * ell_sq / (lam * t0 ** 2)


def entropy_production(q, dq_dz, lam: float, t0: float, ell_sq: float = 0.0):
    """Entropy production of Cattaneo (ell_sq = 0) or Guyer-Krumhansl.

        sigma_s = [ q^2 + 3 ell^2 (dq/dz)^2 ] / (lambda T0^2)

    Non-negative whenever lambda > 0 and ell^2 >= 0, whatever the fields.
    """
    q = np.asarray(q, dtype=float)
    dq_dz = np.asarray(dq_dz, dtype=float)
    return (q ** 2 + 3.0 * ell_sq * dq_dz ** 2) / (lam * t0 ** 2)


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------

def admissibility_table() -> str:
    """The comparison the two criteria produce, as readable text."""
    return (
        "law               convex entropy      second law            propagation\n"
        "Fourier           yes                 yes, lambda > 0       infinite\n"
        "Cattaneo          yes, tau_R > 0      yes, lambda > 0       finite, "
        "sqrt(a/tau_R)\n"
        "Guyer-Krumhansl   yes, tau_R > 0      yes, lambda > 0,      infinite\n"
        "                                      ell^2 >= 0\n"
        "\n"
        "Thermodynamic admissibility and finite propagation speed are "
        "independent criteria.\n"
        "All three laws pass the first; only Cattaneo passes the second. The "
        "hierarchy\n"
        "Fourier, Cattaneo, Guyer-Krumhansl is therefore not a monotone "
        "refinement."
    )
