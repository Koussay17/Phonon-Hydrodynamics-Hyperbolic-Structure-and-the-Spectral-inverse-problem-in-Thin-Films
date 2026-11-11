"""Independent forward model for workstream A (no code shared with src/twin_fdtr.py).

Formulation: thermal quadrupoles in Hankel space (impedance recursion from the substrate upwards), not the
twin's global mode-matching linear system. Each Fourier/Cattaneo layer at radial wavenumber k:
    T' = -g_z Q_z / lam_z,   Q_z' = -(C p + lam_r k^2 / g_r) T,    g_z = 1 + tau_z p,  g_r = 1 + tau_r p,
    gamma^2 = g_z (C p + lam_r k^2 / g_r) / lam_z,   Y = lam_z gamma / g_z,
    Z_top = (Z_b + th/Y) / (1 + Y Z_b th),  th = tanh(gamma d) evaluated as (1 - e^{-2 gamma d})/(1 + e^{-2 gamma d}),
    Z_b = R_below + Z_top(next),  half-space Z = 1/Y.
The in-plane relaxation time tau_r is a separate parameter here (the twin uses one tau_R for both directions).
1D Guyer-Krumhansl (k = 0, 'local' temperature convention): lam_eff = lam (1 + tau_l p)/(1 + tau_R p),
gamma = sqrt(C p / lam_eff), Y = sqrt(lam_eff C p).
Gaussian average: H = 1/(2 pi eta) int_0^U u exp(-u^2) Z(u/sqrt(eta)) du by composite Gauss-Legendre in u
(different variable and rule from the twin's trapezoid in ln k). Derivatives: exact forward-mode autodiff (JAX).
"""
from __future__ import annotations

import os

os.environ.setdefault("XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=3")
os.environ.setdefault("OMP_NUM_THREADS", "3")
import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402

# parameter vector order (log of the physical value unless stated)
PNAMES = ("ln_tauR", "ln_taur_inplane", "ln_lamz", "ln_lamr", "ln_C", "ln_R1", "ln_R2", "ln_lamAu", "ln_CAu",
          "ln_dAu", "ln_dAlN", "ln_lamS", "ln_CS", "ln_wscale", "g", "phi0", "td_ns", "taul_lin")


def _tanh_pos(z):
    e = jnp.exp(-2.0 * z)
    return (1.0 - e) / (1.0 + e)


def _layer_up(Zb, p, k2, lamz, lamr, C, tz, tr, d):
    gz = 1.0 + tz * p
    gr = 1.0 + tr * p
    gam = jnp.sqrt(gz * (C * p + lamr * k2 / gr) / lamz)
    Y = lamz * gam / gz
    th = _tanh_pos(gam * d)
    return (Zb + th / Y) / (1.0 + Y * Zb * th)


def surface_impedance(p, k2, P):
    """Z(k, p) for Au / R1 / AlN (anisotropic Cattaneo) / R2 / sapphire. p, k2 broadcastable."""
    gs = jnp.sqrt((P["CS"] * p + P["lamS"] * k2) / P["lamS"])
    Zs = 1.0 / (P["lamS"] * gs)
    Z = _layer_up(Zs + P["R2"], p, k2, P["lamz"], P["lamr"], P["C"], P["tauR"], P["taur_ip"], P["dAlN"])
    Z = _layer_up(Z + P["R1"], p, k2, P["lamAu"], P["lamAu"], P["CAu"], 0.0, 0.0, P["dAu"])
    return Z


def impedance_1d_gk(p, P):
    """k = 0 impedance with a 1D GK (lam_eff) AlN layer, local temperature convention."""
    lam_eff = P["lamz"] * (1.0 + P["taul"] * p) / (1.0 + P["tauR"] * p)
    gam = jnp.sqrt(P["C"] * p / lam_eff)
    Y = lam_eff * gam
    th = _tanh_pos(gam * P["dAlN"])
    Zs = 1.0 / jnp.sqrt(P["lamS"] * P["CS"] * p)
    Zb = Zs + P["R2"]
    Z = (Zb + th / Y) / (1.0 + Y * Zb * th)
    Z = Z + P["R1"]
    gAu = jnp.sqrt(P["CAu"] * p / P["lamAu"])
    YAu = P["lamAu"] * gAu
    thAu = _tanh_pos(gAu * P["dAu"])
    return (Z + thAu / YAu) / (1.0 + YAu * Z * thAu)


def gl_rule(U=7.0, panels=56, n=16):
    x, w = np.polynomial.legendre.leggauss(n)
    edges = np.linspace(0.0, U, panels + 1)
    a, b = edges[:-1, None], edges[1:, None]
    u = (0.5 * (b - a) * x + 0.5 * (a + b)).ravel()
    wu = (0.5 * (b - a) * w).ravel()
    return u, wu


class Indep:
    """Forward map theta (vector in PNAMES order, increments) -> [ln|H_obs|, phase] for one base point."""

    def __init__(self, base: dict, f, geometry="beam", rule=None, law="cattaneo"):
        self.base = dict(base)
        self.f = np.asarray(f, float)
        self.geometry = geometry
        self.law = law
        self.u, self.wu = gl_rule() if rule is None else rule
        self._fn = jax.jit(self._lnH)
        self._jac = jax.jit(jax.jacfwd(self._lnH_real))

    def phys(self, th):
        b = self.base
        e = jnp.exp
        P = {"tauR": b["tauR"] * e(th[0]), "taur_ip": b["taur_ip"] * e(th[1]) * e(th[0]),
             "lamz": b["lamz"] * e(th[2]), "lamr": b["lamr"] * e(th[3]), "C": b["C"] * e(th[4]),
             "R1": b["R1"] * e(th[5]), "R2": b["R2"] * e(th[6]), "lamAu": b["lamAu"] * e(th[7]),
             "CAu": b["CAu"] * e(th[8]), "dAu": b["dAu"] * e(th[9]), "dAlN": b["dAlN"] * e(th[10]),
             "lamS": b["lamS"] * e(th[11]), "CS": b["CS"] * e(th[12]),
             "wp": b["wp"] * e(th[13]), "wr": b["wr"] * e(th[13]),
             "g": th[14], "phi0": th[15], "td_ns": th[16],
             "taul": b.get("taul", 0.0) + th[17] * b["tauR"]}
        return P

    def _lnH(self, th):
        P = self.phys(th)
        f = jnp.asarray(self.f)
        p = 2j * jnp.pi * f
        eta = (P["wp"] ** 2 + P["wr"] ** 2) / 8.0
        if self.geometry == "1d":
            if self.law == "gk":
                Z0 = impedance_1d_gk(p, P)
            else:
                Z0 = surface_impedance(p, 0.0, P)
            H = Z0 * 2.0 / (jnp.pi * (P["wp"] ** 2 + P["wr"] ** 2))
        else:
            u = jnp.asarray(self.u)
            k2 = (u / jnp.sqrt(eta)) ** 2
            Z = surface_impedance(p[:, None], k2[None, :], P)
            H = (Z * (u * jnp.exp(-u * u) * jnp.asarray(self.wu))[None, :]).sum(axis=1) / (2.0 * jnp.pi * eta)
        lnH = jnp.log(H)
        return lnH + P["g"] + 1j * (P["phi0"] - 2 * jnp.pi * f * P["td_ns"] * 1e-9)

    def _lnH_real(self, th):
        z = self._lnH(th)
        return jnp.concatenate([z.real, z.imag])

    def lnH(self, th=None):
        th = jnp.zeros(len(PNAMES)) if th is None else jnp.asarray(th)
        return np.asarray(self._fn(th))

    def jacobian(self, names, th=None):
        """Complex d ln H_obs / d theta for the named parameters (exact autodiff)."""
        th = jnp.zeros(len(PNAMES)) if th is None else jnp.asarray(th)
        Jr = np.asarray(self._jac(th))            # (2 nf, nparam)
        nf = self.f.size
        Jc = Jr[:nf] + 1j * Jr[nf:]
        idx = [PNAMES.index(n) for n in names]
        return Jc[:, idx]


def base_from_model(m, **over):
    """Base dict for Indep from a common.Model (same physical values)."""
    from common import AU_LAM, AU_C, AU_D, SAP_LAM, SAP_C, R1, R2, C_ALN
    b = {"tauR": m.tauR, "taur_ip": m.tauR, "lamz": m.lamz, "lamr": m.lamr if m.law != "gk" else m.lamz,
         "C": C_ALN, "R1": R1, "R2": R2, "lamAu": AU_LAM, "CAu": AU_C, "dAu": AU_D, "dAlN": m.L_nm * 1e-9,
         "lamS": SAP_LAM, "CS": SAP_C, "wp": m.wp, "wr": m.wr, "taul": m.taul}
    b.update(over)
    return b
