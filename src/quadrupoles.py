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
    "xi_from_thickness",
    "linear_effusivity_profile",
    "graded_linear_layer",
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


# --------------------------------------------------------------------------
# Graded layers: linear metaproperty profiles
# --------------------------------------------------------------------------

def xi_from_thickness(e: float, a: float) -> float:
    """Liouville coordinate of a slab of thickness `e` and constant diffusivity.

    xi = e / sqrt(a), in s^(1/2).
    """
    return e / np.sqrt(a)


def linear_effusivity_profile(xi, b0: float, b1: float, xi1: float,
                              form: str = "T"):
    """Effusivity profile of a 'linear' graded layer, Krapez (2018) eq. (18).

    The metaproperty s = b^(+1/2) (form 'T') or b^(-1/2) (form 'phi') varies
    linearly with the Liouville coordinate between the two faces.
    """
    sign = 1.0 if form == "T" else -1.0
    xi = np.asarray(xi, dtype=float)
    t = xi / xi1
    s = b0 ** (sign / 2) * (1.0 - t) + b1 ** (sign / 2) * t
    return s ** (2.0 / sign)


def graded_linear_layer(p, b0: float, b1: float, xi1: float,
                        form: str = "T") -> np.ndarray:
    """Quadrupole of a graded layer with a linear metaproperty profile.

    Krapez (2018) equation (27), per unit area.

    Parameters
    ----------
    p : complex or array-like
        Spectral parameter.
    b0, b1 : float
        Effusivity at the front and back faces.
    xi1 : float
        Liouville thickness of the layer, in s^(1/2).
    form : {'T', 'phi'}
        Temperature form uses s = b^(+1/2); flux form uses s = b^(-1/2)
        followed by the pseudo-permutation of equation (A-6).

    Notes
    -----
    Reduces exactly to `homogeneous_wall` when b1 == b0, and the determinant
    is 1 for any parameters. Both properties are verified analytically.
    """
    if form not in ("T", "phi"):
        raise ValueError("form must be 'T' or 'phi'")

    p = np.asarray(p, dtype=complex)
    sign = 1.0 if form == "T" else -1.0
    s0 = b0 ** (sign / 2)
    s1 = b1 ** (sign / 2)
    x = s1 / s0

    sp_ = np.sqrt(p)
    arg = sp_ * xi1

    if np.any(np.abs(arg) > _OVERFLOW_GUARD):
        raise OverflowError(
            f"|sqrt(p)*xi1| exceeds {_OVERFLOW_GUARD:.0f}; cosh and sinh overflow."
        )

    ch = np.cosh(arg)
    # sinh(z)/z, continuous at z = 0
    sn = np.where(np.abs(arg) < 1e-8, 1.0 + arg ** 2 / 6.0, np.sinh(arg) / arg)

    u = (x - 1.0) * (1.0 - 1.0 / x)          # recurring coefficient of row C

    a_ = x * ch + (1.0 - x) * sn
    b_ = (xi1 / (s0 * s1)) * sn
    c_ = (s0 * s1 / xi1) * (u * ch + (p * xi1 ** 2 - u) * sn)
    d_ = (1.0 / x) * ch + (1.0 - 1.0 / x) * sn

    m = np.empty(p.shape + (2, 2), dtype=complex)
    if form == "T":
        m[..., 0, 0], m[..., 0, 1] = a_, b_
        m[..., 1, 0], m[..., 1, 1] = c_, d_
    else:
        # Pseudo-permutation, Krapez equation (A-6)
        m[..., 0, 0], m[..., 0, 1] = d_, c_ / p
        m[..., 1, 0], m[..., 1, 1] = p * b_, a_

    return m if p.shape else m.reshape(2, 2)
