"""Axisymmetric FDTR digital twin with Fourier, Cattaneo and Guyer-Krumhansl layers.

Synthetic-experiment generator. It never substitutes for laboratory data.
Conventions follow notes/conventions.md: exp(+i omega t), p = i omega, z
positive into the sample, heat flux positive towards +z, quantities per unit
area, beam radii at 1/e^2, response in kelvin per absorbed watt. Absorption is
at the surface; there is no optical penetration, electron-phonon coupling or
volume source.

Constitutive laws
-----------------
Energy (every layer):   C dT/dt + div q = 0
Fourier:                q = -Lambda grad T,   Lambda = diag(lam_r, lam_r, lam_z)
Cattaneo:               tau_R dq/dt + q = -Lambda grad T
Guyer-Krumhansl:        tau_R dq/dt + q = -lam grad T + ell^2 [lap q + alpha grad div q]
                        (isotropic lam; alpha = 1/3 conserving, alpha = 2 historical)
The one-dimensional reduction is the lambda_eff quadrupole of forward_model.py
with tau_l = (1 + alpha) ell^2 C / lam (note 14, section 6).

Hankel-space layer system
-------------------------
T = int T(k,z) J0(kr) k dk,  q_z = int Q_z J0(kr) k dk,  q_r = int Q_r J1(kr) k dk.
With g = 1 + tau_R p and D = k Q_r + Q_z' (Hankel divergence, D = -C p T):
    C p T + k Q_r + Q_z' = 0
    g Q_z = -lam T' + ell^2 [Q_z'' - k^2 Q_z + alpha D']
    g Q_r =  lam k T + ell^2 [Q_r'' - k^2 Q_r - alpha k D]
The determinant of the symbol factorises into a longitudinal and a solenoidal
factor, giving the characteristic roots
    gamma_P^2 = k^2 + C p g / Lambda_P,  Lambda_P = lam + (1 + alpha) ell^2 C p
    gamma_S^2 = k^2 + g / ell^2
A GK layer has four modes (P and S, decaying downwards or upwards); a Fourier
or Cattaneo layer has two, with gamma^2 = (C p g + lam_r k^2) / lam_z.
P modes carry temperature; S modes have T = 0 and zero divergence. At k = 0
the two families decouple and the P family is the 1D lambda_eff solution.

Boundary conditions
-------------------
Always imposed: unit absorbed flux Q_z = 1 at z = 0; normal-flux continuity
and a temperature jump at each interface (resistance R below each layer);
decay in the semi-infinite substrate. Each face of a GK layer needs one more
condition. The declared family is `GKBoundary`:

* tangential="slip": q_r = C_s ell dq_r/dn_in, n_in pointing into the GK
  layer (Navier slip: Beardo et al. 2020, eq. A7; Hennessy and Myers 2021,
  eq. 5). C_s = 0 no slip, C_s = 1 diffuse walls, C_s = inf free slip. With
  this sign the slip term dissipates.
* tangential="normal_gradient": dq_z/dz = 0. Declared to document it, not
  recommended: at k = 0 it constrains the P family twice and leaves the S
  family free (singular system). For k > 0 the kernel tends to a finite limit
  that differs from the 1D model, carried by S modes whose amplitude grows
  like 1/k; the dissipation identity then balances through boundary terms of
  indefinite sign (see `dissipation_balance`).
* temperature: the jump is theta_above - theta_below = R q_z with
  theta = T - c_n dq_z/dz - c_t div_t(q_t) on GK sides, theta = T elsewhere.
    "local"     c_n = c_t = 0. Convention of the 1D quadrupoles of this project.
    "entropic"  c_n = (1+alpha) ell^2/lam, c_t = alpha ell^2/lam. Natural
                condition of the per-wavenumber dissipation identity; at k = 0
                theta = (1 + tau_l p) T.
    "kinetic"   c_n, c_t supplied, e.g. Beardo et al. (2020) eq. (A6):
                c_n = gamma^-1 (beta - chi_nn), c_t = gamma^-1 (beta - chi_tt).

Observables and passivity
-------------------------
The returned kernel is the surface temperature T(0) per unit absorbed flux.
The per-wavenumber dissipation identity (see `dissipation_balance`) bounds
Re theta_e(0), theta_e the entropic temperature, which equals T(0) when the
top layer is Fourier or Cattaneo (a transducer). Hence:
* Fourier/Cattaneo top layer, entropic temperature, C_s >= 0: Re Z >= 0
  (derived; every term of the identity is nonnegative).
* local temperature: the interface term has no sign; no Re Z < 0 was found
  in adversarial searches with a transducer, but nothing guarantees it.
* kinetic temperature: not passive in general (Re Z/|Z| = -0.94 found with a
  50 nm Au transducer). Gaussian averages check passivity on their k nodes
  (`passivity="warn"`, `"raise"` or `"ignore"`).
* GK top layer (no transducer): T(0) is not the passive variable and Re T(0)
  can be negative even for the entropic temperature; the identity bounds
  theta_e(0), available as `observable="theta_e"`. A warning flags T(0) on
  such stacks in the Gaussian averages.

Numerics
--------
Each mode is referenced at the face where it is largest, so only decaying
exponentials exp(-gamma d) are evaluated. The interface and surface conditions
form a small linear system per (k, p) (7 unknowns for Au/GK-film/sapphire),
equilibrated by row and column scaling and solved in batches. The Gaussian
probe average uses the trapezoidal rule in x = ln(sqrt(eta) k), in which the
integrand is analytic in a strip; convergence is exponential in 1/step. The
step is chosen per frequency from the branch points that lie inside the
Gaussian window, and batches are capped in memory.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import warnings

import numpy as np

__all__ = [
    "Material", "GKBoundary", "Layer", "Stack",
    "kernel", "solve_kernel", "KernelSolution", "dissipation_balance",
    "gaussian_response", "gaussian_response_adaptive", "plane_response",
    "instrument", "synthetic_measurement", "au_aln_sapphire",
    "trapezoid_nodes", "quadrature_strip", "passivity_margin",
    "PassivityWarning", "PassivityError", "NonPassiveObservableWarning",
]


class PassivityWarning(UserWarning):
    """Re of the passive variable is negative somewhere on the evaluated (f, k) nodes."""


class PassivityError(ValueError):
    """Raised instead of PassivityWarning when passivity="raise"."""


class NonPassiveObservableWarning(UserWarning):
    """T(0) of a stack with a GK top layer is not the variable bounded by the identity."""


# --------------------------------------------------------------------------
# Description of the sample
# --------------------------------------------------------------------------

def _finite_positive(name, value):
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _finite_nonnegative(name, value):
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")


@dataclass(frozen=True)
class Material:
    """Homogeneous medium.

    The law follows from the parameters: Fourier when tau_r = ell = 0,
    Cattaneo when ell = 0 < tau_r, Guyer-Krumhansl when ell > 0. GK layers must
    be isotropic (lam_r = lam_z). alpha is the coefficient of grad(div q); it
    only matters for GK. tau_l is the nonlocal time of the 1D reduction.
    """

    lam_z: float
    rho_c: float
    lam_r: float | None = None
    tau_r: float = 0.0
    ell: float = 0.0
    alpha: float = 1.0 / 3.0

    def __post_init__(self):
        if self.lam_r is None:
            object.__setattr__(self, "lam_r", self.lam_z)
        for name in ("lam_z", "lam_r", "rho_c"):
            _finite_positive(name, getattr(self, name))
        for name in ("tau_r", "ell"):
            _finite_nonnegative(name, getattr(self, name))
        if not np.isfinite(self.alpha) or self.alpha <= -1.0 / 3.0:
            raise ValueError("alpha must exceed -1/3 for a nonnegative bulk dissipation")
        if self.ell > 0 and not np.isclose(self.lam_r, self.lam_z, rtol=1e-14, atol=0):
            raise ValueError("Guyer-Krumhansl layers are isotropic here: lam_r must equal lam_z")

    @property
    def law(self) -> str:
        if self.ell > 0:
            return "gk"
        return "cattaneo" if self.tau_r > 0 else "fourier"

    @property
    def diffusivity(self) -> float:
        return self.lam_z / self.rho_c

    @property
    def tau_l(self) -> float:
        """Nonlocal time of the 1D reduction, (1 + alpha) ell^2 C / lam."""
        return (1.0 + self.alpha) * self.ell ** 2 / self.diffusivity


@dataclass(frozen=True)
class GKBoundary:
    """Declared extra conditions at the faces of a Guyer-Krumhansl layer.

    slip        Navier coefficient C_s >= 0 (inf allowed) for tangential="slip".
    tangential  "slip" or "normal_gradient" (the latter is singular at k = 0).
    temperature "local", "entropic" or "kinetic".
    c_n, c_t    coefficients of the kinetic temperature, in m^3 K/W.
    Ignored for Fourier and Cattaneo media.
    """

    slip: float = 1.0
    tangential: str = "slip"
    temperature: str = "local"
    c_n: float = 0.0
    c_t: float = 0.0

    def __post_init__(self):
        if self.tangential not in ("slip", "normal_gradient"):
            raise ValueError("tangential must be 'slip' or 'normal_gradient'")
        if self.temperature not in ("local", "entropic", "kinetic"):
            raise ValueError("temperature must be 'local', 'entropic' or 'kinetic'")
        if np.isnan(self.slip) or self.slip < 0:
            raise ValueError("slip coefficient must be nonnegative (inf allowed)")
        for name in ("c_n", "c_t"):
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if self.temperature != "kinetic" and (self.c_n != 0 or self.c_t != 0):
            raise ValueError("c_n and c_t are used only with temperature='kinetic'")

    @classmethod
    def beardo(cls, gamma_inv, beta, chi_nn, chi_tt, slip=1.0):
        """Kinetic temperature of Beardo et al. (2020), eq. (A6).

        T_transducer - T_GK = R q.n - gamma_inv (beta div q - chi : grad q), n into
        the GK medium, chi diagonal. In the present notation c_n =
        gamma_inv (beta - chi_nn), c_t = gamma_inv (beta - chi_tt).
        """
        return cls(slip=slip, temperature="kinetic",
                   c_n=gamma_inv * (beta - chi_nn), c_t=gamma_inv * (beta - chi_tt))

    def coefficients(self, material: Material):
        """(c_n, c_t) of the temperature entering the jump condition."""
        if material.law != "gk" or self.temperature == "local":
            return 0.0, 0.0
        if self.temperature == "entropic":
            w = material.ell ** 2 / material.lam_z
            return (1.0 + material.alpha) * w, material.alpha * w
        return self.c_n, self.c_t

    def slip_weights(self):
        """(a, b) with a Q_r - b ell dQ_r/dn_in = 0, a + b = 1."""
        if np.isinf(self.slip):
            return 0.0, 1.0
        return 1.0 / (1.0 + self.slip), self.slip / (1.0 + self.slip)


@dataclass(frozen=True)
class Layer:
    """Finite layer; `resistance_below` is the contact resistance at its bottom."""

    material: Material
    thickness: float
    resistance_below: float = 0.0
    boundary: GKBoundary = field(default_factory=GKBoundary)

    def __post_init__(self):
        if not isinstance(self.material, Material):
            raise TypeError("layer material must be a Material")
        _finite_positive("thickness", self.thickness)
        _finite_nonnegative("resistance_below", self.resistance_below)


@dataclass(frozen=True)
class Stack:
    """Layers from the heated face down, on a semi-infinite substrate."""

    layers: tuple
    substrate: Material
    substrate_boundary: GKBoundary = field(default_factory=GKBoundary)

    def __post_init__(self):
        object.__setattr__(self, "layers", tuple(self.layers))
        for layer in self.layers:
            if not isinstance(layer, Layer):
                raise TypeError("layers must be Layer instances")
        if not isinstance(self.substrate, Material):
            raise TypeError("substrate must be a Material")

    @property
    def media(self):
        return [l.material for l in self.layers] + [self.substrate]

    @property
    def boundaries(self):
        return [l.boundary for l in self.layers] + [self.substrate_boundary]

    @property
    def top_is_gk(self) -> bool:
        """True without a Fourier/Cattaneo transducer: T(0) is then not the passive variable."""
        return self.media[0].law == "gk"

    @property
    def unknowns(self) -> int:
        """Size of the global linear system per (k, p)."""
        n = sum(4 if l.material.law == "gk" else 2 for l in self.layers)
        return n + (2 if self.substrate.law == "gk" else 1)


def au_aln_sapphire(aln: Material, aln_thickness=500e-9, au_thickness=80e-9,
                    au=None, sapphire=None, r_au_aln=1e-8, r_aln_sapphire=2e-8,
                    boundary=None):
    """Illustrative Au/AlN/sapphire stack (values of note 16, not measurements).

    Defaults: metal 150 W/(m K), 2.49e6 J/(m^3 K), 80 nm; sapphire 35 W/(m K),
    3.03e6 J/(m^3 K); R = 1e-8 and 2e-8 m^2 K/W.
    """
    au = Material(150.0, 2.49e6) if au is None else au
    sapphire = Material(35.0, 3.03e6) if sapphire is None else sapphire
    boundary = GKBoundary() if boundary is None else boundary
    return Stack((Layer(au, au_thickness, r_au_aln),
                  Layer(aln, aln_thickness, r_aln_sapphire, boundary)),
                 sapphire)


# --------------------------------------------------------------------------
# Modes of one medium
# --------------------------------------------------------------------------

def _modes(m: Material, k, p, semi_infinite: bool):
    """Mode list [(s, decays_downward, (T, Q_z, Q_r, D_r))] at the reference face.

    Down-decaying modes (Re s < 0) are referenced at the top face, up-decaying
    ones at the bottom face. Components are per unit amplitude.
    """
    g = 1.0 + m.tau_r * p
    C = m.rho_c
    out = []
    if m.law == "gk":
        lam, ell, al = m.lam_z, m.ell, m.alpha
        lam_p = lam + (1.0 + al) * ell * ell * C * p
        gp = np.sqrt(k * k + C * p * g / lam_p)
        gs = np.sqrt(k * k + g / (ell * ell))
        one = np.ones_like(gp)
        zero = np.zeros_like(gp)
        y = lam_p / g
        out.append((-gp, True, (one, gp * y, k * y, -gp * k * y)))
        out.append((-gs, True, (zero, k / gs, one + zero, -gs)))
        if not semi_infinite:
            out.append((gp, False, (one, -gp * y, k * y, gp * k * y)))
            out.append((gs, False, (zero, -k / gs, one + zero, gs)))
    else:
        gam = np.sqrt((C * p * g + m.lam_r * k * k) / m.lam_z)
        one = np.ones_like(gam)
        qz = m.lam_z * gam / g
        qr = m.lam_r * k / g + 0 * gam
        out.append((-gam, True, (one, qz, qr, -gam * qr)))
        if not semi_infinite:
            out.append((gam, False, (one, -qz, qr, gam * qr)))
    return out


def _face_factor(s, down, thickness, face):
    """exp(s (zeta - zeta_ref)) at a face; never exceeds one in modulus."""
    if face == "top":
        return 1.0 if down else np.exp(-s * thickness)
    if not np.isfinite(thickness):
        raise ValueError("a semi-infinite medium has no bottom face")
    return np.exp(s * thickness) if down else 1.0


# --------------------------------------------------------------------------
# Global system for one wavenumber and frequency
# --------------------------------------------------------------------------

@dataclass
class KernelSolution:
    """Mode amplitudes of every medium for broadcast (k, p).

    impedance: T(0) per unit absorbed flux. theta_e: entropic temperature
    T - (ell^2/lam)(Q_z' + alpha D) at z = 0, the variable bounded by the
    dissipation identity; equal to impedance unless the top layer is GK.
    """

    stack: Stack
    k: np.ndarray
    p: np.ndarray
    modes: list
    amplitudes: list
    impedance: np.ndarray
    condition: np.ndarray | None = None
    theta_e: np.ndarray | None = None

    def medium_fields(self, j, zeta, derivatives=False):
        """Fields of medium j at local depths zeta (0 at its top), scalar (k, p).

        Returns (T, Q_z, Q_r, D_r), followed by their z derivatives when
        `derivatives` is true; derivatives come from the mode exponents.
        """
        zeta = np.atleast_1d(np.asarray(zeta, dtype=float))
        thick = self.stack.layers[j].thickness if j < len(self.stack.layers) else np.inf
        out = np.zeros((8 if derivatives else 4, zeta.size), dtype=complex)
        for (s, down, comp), a in zip(self.modes[j], self.amplitudes[j]):
            s = complex(s)
            e = np.exp(s * (zeta - (0.0 if down else thick)))
            for c in range(4):
                v = complex(a * comp[c]) * e
                out[c] += v
                if derivatives:
                    out[4 + c] += s * v
        return out

    def fields(self, z):
        """(T, Q_z, Q_r, D_r) at depths z for scalar k and p (per unit absorbed flux).

        An interface depth returns the medium below it.
        """
        z = np.atleast_1d(np.asarray(z, dtype=float))
        tops = np.concatenate([[0.0], np.cumsum([l.thickness for l in self.stack.layers])])
        out = np.zeros((4, z.size), dtype=complex)
        idx = np.clip(np.searchsorted(tops, z, side="right") - 1, 0, len(self.modes) - 1)
        for j in np.unique(idx):
            sel = idx == j
            out[:, sel] = self.medium_fields(j, z[sel] - tops[j])
        return out

    def face(self, j, face, derivatives=False):
        """Field components of medium j at 'top' or 'bottom' (see medium_fields)."""
        thick = self.stack.layers[j].thickness if j < len(self.stack.layers) else np.inf
        res = [0j] * (8 if derivatives else 4)
        for (s, down, comp), a in zip(self.modes[j], self.amplitudes[j]):
            f = _face_factor(s, down, thick, face)
            for c in range(4):
                v = a * comp[c] * f
                res[c] = res[c] + v
                if derivatives:
                    res[4 + c] = res[4 + c] + s * v
        return tuple(complex(np.asarray(v)) for v in res)


def _row_values(m, bc, modes, thick, face, quantity, k, p):
    """Values of a linear functional for every mode of one medium at one face."""
    vals = []
    for s, down, (T, Qz, Qr, Dr) in modes:
        f = _face_factor(s, down, thick, face)
        if quantity == "Qz":
            v = Qz
        elif quantity == "T":
            v = T
        elif quantity == "theta":
            c_n, c_t = bc.coefficients(m)
            # dQ_z/dz = -C p T - k Q_r exactly (energy equation)
            v = T + c_n * (m.rho_c * p * T + k * Qr) - c_t * k * Qr
        elif quantity == "extra":
            if bc.tangential == "slip":
                a, b = bc.slip_weights()
                sign = 1.0 if face == "top" else -1.0   # d/dn_in
                v = a * Qr - sign * b * m.ell * Dr
            else:
                v = (m.rho_c * p * T + k * Qr) / (m.rho_c * np.abs(p))
        else:
            raise ValueError(quantity)
        vals.append(v * f)
    return vals


def solve_kernel(freq, k, stack: Stack, diagnostics=False):
    """Solve the layered problem for broadcast arrays of frequency and wavenumber.

    Returns a KernelSolution whose `impedance` is the surface temperature per
    unit absorbed flux, Z(k, p) = T(0) / Q_z(0), in m^2 K/W.
    """
    f = np.asarray(freq, dtype=float)
    kk = np.asarray(k, dtype=float)
    if np.any(~np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("frequencies must be finite and positive")
    if np.any(~np.isfinite(kk)) or np.any(kk < 0):
        raise ValueError("radial wavenumbers must be finite and nonnegative")
    f, kk = np.broadcast_arrays(f, kk)
    p = 2j * np.pi * f
    shape = f.shape
    media = stack.media
    bcs = stack.boundaries
    nlay = len(stack.layers)
    thick = [l.thickness for l in stack.layers] + [np.inf]
    modes = [_modes(m, kk, p, j == nlay) for j, m in enumerate(media)]
    sizes = [len(md) for md in modes]
    offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    n = int(offs[-1])
    M = np.zeros(shape + (n, n), dtype=complex)
    rhs = np.zeros(shape + (n,), dtype=complex)

    def put(row, j, face, quantity, scale=1.0):
        vals = _row_values(media[j], bcs[j], modes[j], thick[j], face, quantity, kk, p)
        for i, v in enumerate(vals):
            M[..., row, offs[j] + i] += scale * v

    row = 0
    put(row, 0, "top", "Qz")
    rhs[..., row] = 1.0
    row += 1
    if media[0].law == "gk":
        put(row, 0, "top", "extra")
        row += 1
    for j in range(nlay):
        r = stack.layers[j].resistance_below
        put(row, j, "bottom", "Qz")
        put(row, j + 1, "top", "Qz", -1.0)
        row += 1
        put(row, j, "bottom", "theta")
        if r:
            put(row, j, "bottom", "Qz", -r)
        put(row, j + 1, "top", "theta", -1.0)
        row += 1
        if media[j].law == "gk":
            put(row, j, "bottom", "extra")
            row += 1
        if media[j + 1].law == "gk":
            put(row, j + 1, "top", "extra")
            row += 1
    if row != n:
        raise RuntimeError("boundary-condition count does not match the mode count")

    # Row and column equilibration, then batched LU solve.
    col = np.max(np.abs(M), axis=-2, keepdims=True)
    col = np.where(col > 0, col, 1.0)
    Ms = M / col
    rw = np.max(np.abs(Ms), axis=-1, keepdims=True)
    rw = np.where(rw > 0, rw, 1.0)
    Ms = Ms / rw
    y = np.linalg.solve(Ms, (rhs / rw[..., 0])[..., None])[..., 0]
    x = y / col[..., 0, :]
    amps = [x[..., offs[j]:offs[j + 1]] for j in range(len(media))]
    T0 = 0
    for i, (s, down, comp) in enumerate(modes[0]):
        T0 = T0 + amps[0][..., i] * comp[0] * _face_factor(s, down, thick[0], "top")
    T0 = np.asarray(T0)
    theta0 = T0
    m0 = media[0]
    if m0.law == "gk":
        Qr0 = dQz0 = 0
        for i, (s, down, comp) in enumerate(modes[0]):
            v = amps[0][..., i] * _face_factor(s, down, thick[0], "top")
            Qr0 = Qr0 + v * comp[2]
            dQz0 = dQz0 + v * s * comp[1]
        D0 = kk * Qr0 + dQz0
        theta0 = T0 - (m0.ell ** 2 / m0.lam_z) * (dQz0 + m0.alpha * D0)
    cond = np.linalg.cond(Ms) if diagnostics else None
    return KernelSolution(stack, kk, p, modes,
                          [[a[..., i] for i in range(a.shape[-1])] for a in amps],
                          T0, cond, np.asarray(theta0))


def kernel(freq, k, stack: Stack, observable="T"):
    """Surface kernel in m^2 K/W per unit absorbed flux, broadcast over (freq, k).

    observable="T" (default): T(0)/Q_z(0), the probe temperature. With a GK top
    layer (no transducer) this is not the passive variable: Re T(0) can be
    negative. observable="theta_e": the entropic temperature at z = 0, which the
    dissipation identity bounds (equal to T(0) for Fourier/Cattaneo top layers).
    """
    sol = solve_kernel(freq, k, stack)
    if observable == "T":
        return sol.impedance
    if observable == "theta_e":
        return sol.theta_e
    raise ValueError("observable must be 'T' or 'theta_e'")


def passivity_margin(freq, k, stack: Stack):
    """min over (freq, k) of Re(P)/|P|, P the passive variable (theta_e at z = 0).

    Nonnegative when the stack satisfies the dissipation identity with
    nonnegative terms (entropic temperature, C_s >= 0); a negative value is a
    counterexample to passivity for the declared boundary conditions.
    """
    th = solve_kernel(freq, k, stack).theta_e
    return float(np.min(th.real / np.abs(th)))


def _passivity_policy(margin, policy, where):
    if policy == "ignore" or margin >= -1e-8:
        return
    msg = (f"{where}: Re of the passive variable is negative (min Re/|.| = {margin:.3g}); "
           "the declared boundary conditions are not dissipative for this stack. Only the "
           "entropic temperature with C_s >= 0 is guaranteed passive.")
    if policy == "raise":
        raise PassivityError(msg)
    warnings.warn(msg, PassivityWarning, stacklevel=3)


# --------------------------------------------------------------------------
# Energy and dissipation diagnostics
# --------------------------------------------------------------------------

_TINY_FLUX = 1e-250   # fluxes below this fraction of the unit input are treated as zero


def _gauss_panels(breaks, n):
    xg, wg = np.polynomial.legendre.leggauss(n)
    h = np.diff(breaks)
    nodes = 0.5 * h[:, None] * xg + 0.5 * (breaks[1:] + breaks[:-1])[:, None]
    return nodes.ravel(), (0.5 * h[:, None] * wg).ravel()


def _cap_panels(breaks, modes, thickness):
    """Split panels wider than half a wavelength of a mode still alive on them.

    modes: (gamma, decays_downward). A mode referenced at the top (bottom) face is
    alive on [a, b] while Re(gamma) a (Re(gamma) (d - b)) < 40; its half
    wavelength is pi/|Im gamma|. Short-wavelength modes confined to thin
    boundary layers therefore do not refine the interior.
    """
    out = [breaks[0]]
    for a, b in zip(breaks[:-1], breaks[1:]):
        cap = np.inf
        for s, down in modes:
            if s.imag == 0:
                continue
            dist = a if down else thickness - b
            if abs(s.real) * dist < 40.0:
                cap = min(cap, np.pi / abs(s.imag))
        m = 1 if not np.isfinite(cap) else max(1, int(math.ceil((b - a) / cap)))
        out.extend(a + (b - a) * np.arange(1, m + 1) / m)
    return np.asarray(out)


def _depth_rule(thickness, modes, n=16):
    """Gauss-Legendre panels for a finite layer.

    Graded geometrically towards both faces down to 1/(8 max|gamma|), which
    resolves thin boundary layers and short wavelengths alike, then split by
    _cap_panels so that oscillatory fields are resolved across the layer.
    """
    rate = max(abs(s) for s, _ in modes)
    half = thickness / 2.0
    s0 = min(1.0 / (8.0 * rate), half / 8.0)
    m = int(math.ceil(math.log2(half / s0))) + 1
    left = np.concatenate([[0.0], np.geomspace(s0, half, m)])
    breaks = np.unique(np.concatenate([left, thickness - left[::-1]]))
    return _gauss_panels(_cap_panels(breaks, modes, thickness), n)


def _halfspace_rule(modes, n=16):
    """Same for the substrate, up to the depth where single fields fall below exp(-46)."""
    decay_min = min(abs(s.real) for s, _ in modes)
    rate = max(abs(s) for s, _ in modes)
    length = 46.0 / decay_min
    s0 = 1.0 / (8.0 * rate)
    m = int(math.ceil(math.log2(max(length / s0, 2.0)))) + 1
    breaks = np.concatenate([[0.0], np.geomspace(s0, length, m)])
    return _gauss_panels(_cap_panels(breaks, modes, np.inf), n)


def _entropic_theta(m: Material, T, Qz, Qr, dQz, k):
    """theta_e = T - (ell^2/lam)(Q_z' + alpha D), D = k Q_r + Q_z' (T for non-GK)."""
    if m.law != "gk":
        return T
    D = k * Qr + dQz
    return T - (m.ell ** 2 / m.lam_z) * (dQz + m.alpha * D)


def dissipation_balance(freq, k, stack: Stack, n=16):
    """Per-wavenumber energy and dissipation balance of a solved kernel.

    For a solution at real (k, omega), the Hankel-space equations give
        Re theta_e(0) - surface_slip = sum(bulk) + sum(interface terms),
    with bulk = int |Q|^2/lam + (ell^2/lam)(|Q_z'|^2 + |Q_r'|^2 + k^2|Q|^2 + alpha|D|^2) dz
    (|Q_z|^2/lam_z + |Q_r|^2/lam_r for Fourier and Cattaneo media), and per interface
        Re[conj(q)(theta_e,above - theta_e,below)] - (ell^2/lam)Re(conj(Q_r) Q_r')|above
                                                    + (ell^2/lam)Re(conj(Q_r) Q_r')|below.
    theta_e is the entropic temperature, equal to T for Fourier top layers, so
    the left side is then Re Z. With the entropic temperature and Navier slip
    every term is nonnegative: Re Z >= 0. Integrals use graded Gauss-Legendre
    quadrature of the reconstructed fields, independent of the linear solve,
    with panels graded by |gamma| and limited to half a wavelength of every mode
    alive on them. Also returns the integrated energy balance of every medium.
    Relative residuals use max(scale, 1e-250 x input flux) as denominator, so
    fields that underflow to zero deep in a stack give zero residuals.
    """
    sol = solve_kernel(float(freq), float(k), stack)
    p = complex(sol.p)
    kk = float(k)
    media = stack.media
    nlay = len(stack.layers)
    bulk, energy = [], []
    for j, m in enumerate(media):
        md = [(complex(np.asarray(s)), down) for s, down, _ in sol.modes[j]]
        if j < nlay:
            z, w = _depth_rule(stack.layers[j].thickness, md, n)
        else:
            z, w = _halfspace_rule(md, n)
        T, Qz, Qr, Dr, dT, dQz, dQr, dDr = sol.medium_fields(j, z, derivatives=True)
        Q2 = np.abs(Qz) ** 2 + np.abs(Qr) ** 2
        if m.law == "gk":
            D = kk * Qr + dQz
            dens = Q2 / m.lam_z + (m.ell ** 2 / m.lam_z) * (
                np.abs(dQz) ** 2 + np.abs(dQr) ** 2 + kk * kk * Q2 + m.alpha * np.abs(D) ** 2)
        else:
            dens = np.abs(Qz) ** 2 / m.lam_z + np.abs(Qr) ** 2 / m.lam_r
        bulk.append(float(np.sum(w * dens)))
        top = sol.face(j, "top")
        qbot = sol.face(j, "bottom")[1] if j < nlay else 0.0
        storage = m.rho_c * p * np.sum(w * T)
        lateral = kk * np.sum(w * Qr)
        scale = abs(storage) + abs(lateral) + abs(top[1]) + abs(qbot)
        energy.append(abs(storage + lateral + qbot - top[1]) / max(scale, _TINY_FLUX))
    interfaces = []
    for j in range(nlay):
        mu, ml = media[j], media[j + 1]
        U = sol.face(j, "bottom", derivatives=True)
        L = sol.face(j + 1, "top", derivatives=True)
        q = U[1]
        thU = _entropic_theta(mu, U[0], U[1], U[2], U[5], kk)
        thL = _entropic_theta(ml, L[0], L[1], L[2], L[5], kk)
        R = stack.layers[j].resistance_below
        slipU = -(mu.ell ** 2 / mu.lam_z) * (np.conj(U[2]) * U[3]).real if mu.law == "gk" else 0.0
        slipL = (ml.ell ** 2 / ml.lam_z) * (np.conj(L[2]) * L[3]).real if ml.law == "gk" else 0.0
        interfaces.append({
            "resistive": R * abs(q) ** 2,
            "nonequilibrium": (np.conj(q) * (thU - thL - R * q)).real,
            "tangential_above": float(slipU), "tangential_below": float(slipL),
            "flux_jump": abs(U[1] - L[1]) / max(abs(U[1]), abs(L[1]), _TINY_FLUX),
        })
    S = sol.face(0, "top", derivatives=True)
    m0 = media[0]
    surface_input = _entropic_theta(m0, S[0], S[1], S[2], S[5], kk).real
    surface_slip = (m0.ell ** 2 / m0.lam_z) * (np.conj(S[2]) * S[3]).real if m0.law == "gk" else 0.0
    total = sum(bulk) + sum(i["resistive"] + i["nonequilibrium"] + i["tangential_above"] + i["tangential_below"]
                            for i in interfaces)
    lhs = surface_input - surface_slip
    return {
        "re_impedance": float(np.real(sol.impedance)),
        "surface_input": float(surface_input),
        "surface_slip": float(surface_slip),
        "bulk": bulk,
        "interfaces": interfaces,
        "energy_residuals": energy,
        "residual": abs(lhs - total) / (abs(lhs) + sum(bulk)),
    }


# --------------------------------------------------------------------------
# Gaussian pump and probe
# --------------------------------------------------------------------------

def quadrature_strip(freq, stack: Stack, eta=None, window=8.0):
    """Half-width of the analyticity strip of the integrand in x = ln k, per frequency.

    Branch points of gamma(k) lie at k*^2 = -c with c = C p g / lam_r (Fourier,
    Cattaneo), C p g / Lambda_P (GK, P modes) or g / ell^2 (GK, S modes), i.e. at
    |Im x| = (pi - arg c)/2: pi/4 for Fourier, narrower for a weakly damped wave.
    With eta given, branch points with sqrt(eta)|k*| > window lie where the
    Gaussian factor is below exp(-window^2) and are ignored. The Gaussian factor
    itself limits the strip to pi/4.
    """
    w = 2 * np.pi * np.atleast_1d(np.asarray(freq, dtype=float))
    p = 1j * w
    width = np.full(w.shape, np.pi / 4)

    def narrow(c, width):
        strip = (np.pi - np.angle(c)) / 2
        if eta is not None:
            strip = np.where(np.sqrt(eta * np.abs(c)) > window, np.pi / 4, strip)
        return np.minimum(width, strip)

    for m in stack.media:
        g = 1 + m.tau_r * p
        if m.law == "gk":
            lam_p = m.lam_z + (1 + m.alpha) * m.ell ** 2 * m.rho_c * p
            width = narrow(m.rho_c * p * g / lam_p, width)
            width = narrow(g / m.ell ** 2, width)
        else:
            width = narrow(m.rho_c * p * g / m.lam_r, width)
    return width


def trapezoid_nodes(freq, stack: Stack, pump_radius, probe_radius, step=None,
                    x_low=None, x_high=2.0):
    """Nodes and weights of the log-variable trapezoidal rule.

    H = (1/(2 pi eta)) int exp(2x - exp(2x)) Z(exp(x)/sqrt(eta)) dx.
    The default step is 0.06 * (strip / (pi/4)) with the strip of the given
    frequencies inside the Gaussian window; the lower cut-off sits 16 e-folds
    below the smallest thermal crossover of the stack.
    """
    eta = (pump_radius ** 2 + probe_radius ** 2) / 8.0
    f = np.atleast_1d(np.asarray(freq, dtype=float))
    if step is None:
        strip = float(np.min(quadrature_strip(f, stack, eta=eta)))
        step = 0.06 * max(strip, 1e-3) / (np.pi / 4)
    if x_low is None:
        uc = min(np.sqrt(eta * m.rho_c * 2 * np.pi * f.min() / max(m.lam_r, m.lam_z))
                 for m in stack.media)
        x_low = min(-18.0, math.log(uc) - 16.0)
    n = int(math.ceil((x_high - x_low) / step)) + 1
    x = np.linspace(x_low, x_high, n)
    h = x[1] - x[0]
    w = h * np.exp(2 * x - np.exp(2 * x)) / (2 * np.pi * eta)
    w[0] *= 0.5
    w[-1] *= 0.5
    return np.exp(x) / np.sqrt(eta), w


def _check_beams(pump_radius, probe_radius):
    for name, v in (("pump_radius", pump_radius), ("probe_radius", probe_radius)):
        _finite_positive(name, v)


def _check_options(observable, passivity):
    if observable not in ("T", "theta_e"):
        raise ValueError("observable must be 'T' or 'theta_e'")
    if passivity not in ("warn", "raise", "ignore"):
        raise ValueError("passivity must be 'warn', 'raise' or 'ignore'")


def _warn_gk_top(stack, observable, where):
    if stack.top_is_gk and observable == "T":
        warnings.warn(f"{where}: the top layer is Guyer-Krumhansl (no transducer); T(0) is not "
                      "the variable bounded by the dissipation identity and Re T(0) may be "
                      "negative. Use observable='theta_e' for the entropic surface temperature.",
                      NonPassiveObservableWarning, stacklevel=3)


def gaussian_response(freq, stack: Stack, pump_radius, probe_radius, step=None,
                      x_low=None, chunk=64, observable="T", passivity="warn",
                      max_bytes=2 ** 28):
    """Probe-averaged complex temperature per absorbed watt, K/W.

    Pump q(r) = 2P/(pi wp^2) exp(-2 r^2/wp^2); probe weight 2/(pi wr^2)
    exp(-2 r^2/wr^2). H = (1/2pi) int k exp(-eta k^2) Z(k,p) dk,
    eta = (wp^2 + wr^2)/8 (note 16). The phase is arg H, not a mean phase.

    Each frequency gets the step of its own strip (rounded down on a half-octave
    ladder so that frequencies share rules); batches are limited to about
    max_bytes of system matrices. observable: "T" (probe temperature) or
    "theta_e" (see `kernel`). passivity: the passive variable theta_e(0) is
    checked on every evaluated node and a negative real part triggers a
    PassivityWarning ("warn"), a PassivityError ("raise") or nothing ("ignore").
    """
    _check_beams(pump_radius, probe_radius)
    _check_options(observable, passivity)
    f = np.atleast_1d(np.asarray(freq, dtype=float))
    if f.ndim != 1 or not f.size or np.any(~np.isfinite(f)) or np.any(f <= 0):
        raise ValueError("frequencies must be a nonempty positive vector")
    eta = (pump_radius ** 2 + probe_radius ** 2) / 8.0
    if step is None:
        hf = 0.06 * np.maximum(quadrature_strip(f, stack, eta=eta), 1e-3) / (np.pi / 4)
        hq = 0.06 * 2.0 ** (-np.ceil(2 * np.log2(0.06 / hf) - 1e-9) / 2)
    else:
        hq = np.full(f.size, float(step))
    n2 = stack.unknowns ** 2
    out = np.empty(f.size, dtype=complex)
    margin = np.inf
    for h in np.unique(hq):
        idx = np.nonzero(hq == h)[0]
        kn, w = trapezoid_nodes(f[idx], stack, pump_radius, probe_radius, h, x_low)
        batch = max(1, min(chunk, int(max_bytes // (4 * 16 * n2 * kn.size))))
        for i in range(0, idx.size, batch):
            sel = idx[i:i + batch]
            sol = solve_kernel(f[sel][:, None], kn[None, :], stack)
            z = sol.impedance if observable == "T" else sol.theta_e
            out[sel] = z @ w
            margin = min(margin, float(np.min(sol.theta_e.real / np.abs(sol.theta_e))))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite probe-averaged response")
    _passivity_policy(margin, passivity, "gaussian_response")
    _warn_gk_top(stack, observable, "gaussian_response")
    return out


def gaussian_response_adaptive(freq, stack: Stack, pump_radius, probe_radius,
                               epsrel=1e-11, upper=10.0, observable="T", passivity="warn"):
    """Same observable by adaptive quadrature in u = sqrt(eta) k (cross-check)."""
    from scipy.integrate import quad_vec
    _check_beams(pump_radius, probe_radius)
    _check_options(observable, passivity)
    f = np.atleast_1d(np.asarray(freq, dtype=float))
    eta = (pump_radius ** 2 + probe_radius ** 2) / 8.0
    margin = [np.inf]

    def integrand(u):
        sol = solve_kernel(f, u / np.sqrt(eta), stack)
        z = sol.impedance if observable == "T" else sol.theta_e
        margin[0] = min(margin[0], float(np.min(sol.theta_e.real / np.abs(sol.theta_e))))
        return u * np.exp(-u * u) * z / (2 * np.pi * eta)

    pts = []
    for m in stack.media:
        pts.extend(np.sqrt(eta * m.rho_c * 2 * np.pi * f / m.lam_r))
    pts = np.unique([x for x in pts if 0 < x < upper])
    val, err, info = quad_vec(integrand, 0.0, upper, epsabs=0.0, epsrel=epsrel,
                              points=pts if len(pts) else None, full_output=True,
                              limit=20000)
    if not info.success or not np.all(np.isfinite(val)):
        raise RuntimeError(f"adaptive radial integration failed: {info.message}")
    _passivity_policy(margin[0], passivity, "gaussian_response_adaptive")
    _warn_gk_top(stack, observable, "gaussian_response_adaptive")
    return val


def plane_response(freq, stack: Stack, pump_radius, probe_radius, observable="T"):
    """1D approximation of the same observable: Z(0, p) * 2/(pi (wp^2 + wr^2))."""
    _check_beams(pump_radius, probe_radius)
    f = np.atleast_1d(np.asarray(freq, dtype=float))
    return kernel(f, 0.0, stack, observable) * 2.0 / (np.pi * (pump_radius ** 2 + probe_radius ** 2))


# --------------------------------------------------------------------------
# Instrument and synthetic data
# --------------------------------------------------------------------------

def instrument(freq, response, log_gain=0.0, phase_offset_deg=0.0, delay=0.0):
    """H_obs = exp(g + i phi0 - i 2 pi f t_d) H (note 16); delay in seconds."""
    f = np.asarray(freq, dtype=float)
    return np.asarray(response) * np.exp(
        log_gain + 1j * (np.deg2rad(phase_offset_deg) - 2 * np.pi * f * delay))


def synthetic_measurement(freq, stack: Stack, pump_radius, probe_radius,
                          log_gain=0.0, phase_offset_deg=0.0, delay=0.0,
                          sigma_log_amplitude=0.01, sigma_phase_deg=0.1, seed=None,
                          passivity="warn"):
    """Synthetic lock-in data: SYNTHETIC ONLY, never a measurement.

    Noise is Gaussian, independent across frequencies: sigma on ln|H| and on
    the phase in degrees (note 16). Returns a dict with the noiseless and
    noisy complex signals, amplitude, phase in degrees and the settings.
    """
    f = np.atleast_1d(np.asarray(freq, dtype=float))
    truth = instrument(f, gaussian_response(f, stack, pump_radius, probe_radius,
                                            passivity=passivity),
                       log_gain, phase_offset_deg, delay)
    rng = np.random.default_rng(seed)
    noise = (sigma_log_amplitude * rng.standard_normal(f.size)
             + 1j * np.deg2rad(sigma_phase_deg) * rng.standard_normal(f.size))
    observed = truth * np.exp(noise)
    return {
        "status": "SYNTHETIC ONLY",
        "frequency_Hz": f,
        "truth": truth,
        "observed": observed,
        "amplitude": np.abs(observed),
        "phase_deg": np.degrees(np.angle(observed)),
        "sigma_log_amplitude": sigma_log_amplitude,
        "sigma_phase_deg": sigma_phase_deg,
        "seed": seed,
        "instrument": {"log_gain": log_gain, "phase_offset_deg": phase_offset_deg,
                       "delay_s": delay},
    }
