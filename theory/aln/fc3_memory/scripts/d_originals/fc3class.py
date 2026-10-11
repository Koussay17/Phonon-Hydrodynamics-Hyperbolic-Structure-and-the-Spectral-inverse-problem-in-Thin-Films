"""fc3-parametrised rate class on the fixed AlN event set: g_alpha = W_alpha |Phi_alpha(fc3)|^2.

W_alpha (fc3-independent) = conv * g0w / (4 s_p s_a s_b)  (distinct daughters)
                          = conv * g0w / (8 s_p s_a^2)     (repeated daughter),  s = sinh(x/2), x = eps/kT,
with g0w the Gaussian weight delta_sigma(w0 - w1 - w2) exported by p3_events.py and conv the CollisionMatrix
unit conversion; pp = |Phi|^2 * Interaction.unit_conversion is phono3py's interaction strength.
Rates are tied on the 24-operation orbits (6mm x time reversal) by orbit averaging, as for the memory-gap
reference. Units: g in THz; tau[ps] = (N/K)/(4 pi); kappa via aln_geometry.kappa_factor.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True

from aln_geometry import load as load_geom, kappa_factor  # noqa: E402
from fastops import FastOps  # noqa: E402
from symmetry import orbits  # noqa: E402

FOURPI = 4.0 * np.pi


class FC3Class:
    def __init__(self, events_npz, pp_ref_npz):
        self.geom = geom = load_geom(events_npz)
        self.tie = orbits(geom, geom.labels["maps"])
        self.orb = self.tie["ev_orb"]
        self.osize = np.bincount(self.orb)
        self.n_orb = len(self.osize)
        self.ops = FastOps(geom)
        P = np.load(pp_ref_npz)
        self.conv = float(P["conv"])
        self.g0w = P["g0w"]
        self.tp = P["tp"]
        self.bands = P["bands"]
        ev = geom.events
        x = geom.eps / geom.kT
        s = np.sinh(x / 2)
        rep = ev[:, 2] < 0
        sp, sa = s[ev[:, 0]], s[ev[:, 1]]
        sb = np.where(rep, sa, s[np.where(rep, 0, ev[:, 2])])
        self.rep = rep
        self.W = np.where(rep, self.conv * self.g0w / (8 * sp * sa * sa), self.conv * self.g0w / (4 * sp * sa * sb))
        self.kf = kappa_factor(geom)
        self.bnorm2 = (geom.b ** 2).sum(axis=0)
        self.pp_ref_raw = P["pp"]

    # ------------------------------------------------------------------ rates
    def orbit_avg(self, g):
        return (np.bincount(self.orb, weights=g, minlength=self.n_orb) / self.osize)[self.orb]

    def rates(self, pp, tie=True):
        g = self.W * pp
        return self.orbit_avg(g) if tie else g

    # ------------------------------------------------------------------ observables
    def metrics(self, g, cols=(0, 2), want_lmin=True, want_x=False):
        ops = self.ops
        C = ops.C(g)
        e = self.geom.ehat
        scale = np.trace(C) / self.geom.n
        Creg = C + scale * np.outer(e, e)
        cf = sla.cho_factor(Creg, lower=True, check_finite=False) if len(cols) else None
        out = {"r": np.diag(C).copy()}
        for c in cols:
            b = self.geom.b[:, c]
            xx = sla.cho_solve(cf, b, check_finite=False)
            K = float(b @ xx)
            N = float(xx @ xx)
            name = "xyz"[c]
            out[f"K_{name}"] = K
            out[f"kappa_{name}"] = self.kf * K
            out[f"tau_{name}"] = N / K / FOURPI
            out[f"tauCS_{name}"] = K / self.bnorm2[c] / FOURPI
            if want_x:
                out[f"x_{name}"] = xx
        if want_lmin:
            w = np.linalg.eigvalsh(Creg)
            # the regularised direction e has eigenvalue 'scale' (> lambda_min of C|H for these operators)
            wH = np.sort(w)
            # remove the eigenvalue belonging to e: identify via Rayleigh quotient of e
            out["lam_min_H"] = float(wH[0]) if wH[0] < scale * (1 - 1e-9) else float(wH[1])
            out["ceiling_F1_ps"] = (1.0 / out["lam_min_H"] / FOURPI) if out["lam_min_H"] > 0 else float("inf")
        return out


def summary(m):
    return {k: v for k, v in m.items() if not isinstance(v, np.ndarray)}
