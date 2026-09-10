"""Thermal quadrupoles for one-dimensional heat conduction.

Conventions
-----------
All quantities are per unit area. The flux is a heat flux density in W/m^2,
so the cross-section factor S appearing in Maillet et al. (2000) is set to 1.
See notes/conventions.md.

State vector ordering, front face to back face:

    [theta_0]       [A  B] [theta_1]
    [ phi_0 ]   =   [C  D] [ phi_1 ]

For a stack, the total matrix is the ordered product of the individual
matrices, front to back. The determinant of every quadrupole is 1, a
consequence of Abel's theorem applied to the Liouville normal form.

Reference for the homogeneous wall: Maillet, Andre, Batsale, Degiovanni,
Moyne, *Thermal Quadrupoles*, Wiley (2000), section 1.3.4 page 10.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "diffusivity",
    "effusivity",
    "from_diffusivity_effusivity",
    "homogeneous_wall",
    "semi_infinite_impedance",
    "compose",
    "front_face_temperature",
    "argument_magnitude",
]

# Above this value cosh and sinh overflow in double precision.
_OVERFLOW_GUARD = 700.0


# --------------------------------------------------------------------------
# Material properties
# --------------------------------------------------------------------------

def diffusivity(lam: float, rho_c: float) -> float:
    """Thermal diffusivity a = lambda / (rho c), in m^2/s."""
    return lam / rho_c


def effusivity(lam: float, rho_c: float) -> float:
    """Thermal effusivity b = sqrt(lambda * rho c), in W s^(1/2) m^-2 K^-1."""
    return np.sqrt(lam * rho_c)


def from_diffusivity_effusivity(a: float, b: float) -> tuple[float, float]:
    """Recover (lambda, rho c) from (a, b).

    Inverse of the two definitions above:  lambda = b sqrt(a),  rho c = b / sqrt(a).
    """
    return b * np.sqrt(a), b / np.sqrt(a)


# --------------------------------------------------------------------------
# Quadrupoles
# --------------------------------------------------------------------------

def homogeneous_wall(p, lam: float, rho_c: float, e: float) -> np.ndarray:
    """Quadrupole of a homogeneous slab of thickness `e`.

    Parameters
    ----------
    p : complex or array-like
        Spectral parameter. Use p = 1j * omega for the harmonic regime, or a
        real positive value for the Laplace regime.
    lam, rho_c : float
        Thermal conductivity and volumetric heat capacity.
    e : float
        Slab thickness in metres.

    Returns
    -------
    ndarray of shape (2, 2) or (..., 2, 2)
        The quadrupole matrix, per unit area.

    Notes
    -----
    With k = sqrt(p / a) and b the effusivity, the identity lambda * k = b sqrt(p)
    holds, which is what makes the two possible parameterisations agree.
    """
    p = np.asarray(p, dtype=complex)
    a = diffusivity(lam, rho_c)
    b = effusivity(lam, rho_c)

    k = np.sqrt(p / a)
    ke = k * e

    if np.any(np.abs(ke) > _OVERFLOW_GUARD):
        raise OverflowError(
            f"|k*e| exceeds {_OVERFLOW_GUARD:.0f}; cosh and sinh overflow in double "
            "precision. Reduce the frequency range or split the layer."
        )

    ch = np.cosh(ke)
    sh = np.sinh(ke)
    bsp = b * np.sqrt(p)

    m = np.empty(p.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = ch
    m[..., 0, 1] = sh / bsp
    m[..., 1, 0] = bsp * sh
    m[..., 1, 1] = ch
    return m if p.shape else m.reshape(2, 2)


def semi_infinite_impedance(p, b: float):
    """Thermal impedance of a semi-infinite homogeneous medium, per unit area.

    Z = 1 / (b sqrt(p)).

    Cross-checked against Maillet et al. section 1.4.2 remark 2 page 15, and
    Krapez (2018) equation (29).
    """
    p = np.asarray(p, dtype=complex)
    return 1.0 / (b * np.sqrt(p))


def compose(*matrices: np.ndarray) -> np.ndarray:
    """Ordered product of quadrupoles, front face to back face."""
    if not matrices:
        raise ValueError("compose() requires at least one matrix")
    out = matrices[0]
    for m in matrices[1:]:
        out = out @ m
    return out


def front_face_temperature(m: np.ndarray, z, power=1.0, h: float = 0.0):
    """Front-face temperature of a stack closed by an impedance.

        theta_0 = power * (A Z + B) / (C Z + D + h (A Z + B))

    `h` is the heat exchange coefficient on the free front surface; h = 0
    corresponds to the adiabatic case. Krapez (2018) equation (29).
    """
    a_, b_ = m[..., 0, 0], m[..., 0, 1]
    c_, d_ = m[..., 1, 0], m[..., 1, 1]
    num = a_ * z + b_
    return power * num / (c_ * z + d_ + h * num)


def argument_magnitude(p, lam: float, rho_c: float, e: float):
    """Return |k e|, the argument that drives numerical overflow.

    Useful to check a frequency range before building the matrices.
    """
    p = np.asarray(p, dtype=complex)
    return np.abs(np.sqrt(p / diffusivity(lam, rho_c)) * e)
