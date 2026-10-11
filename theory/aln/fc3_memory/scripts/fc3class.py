"""fc3-parametrised rate class on the fixed AlN event set: g_alpha = W_alpha |Phi_alpha(fc3)|^2 (C: mirror copy).

W_alpha = conv * g0w / (4 s_p s_a s_b) (distinct daughters) or conv * g0w / (8 s_p s_a^2) (repeated daughter),
s = sinh(x/2); pp = |Phi|^2 * unit_conversion (phono3py interaction strength). Orbit tying by averaging over the
24-operation orbits (memory-gap convention). tau[ps] = (N/K)/(4 pi).
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
        self.g0w, self.tp, self.bands = P["g0w"], P["tp"], P["bands"]
        ev = geom.events
        s = np.sinh(geom.eps / geom.kT / 2)
        rep = ev[:, 2] < 0
        sp_, sa = s[ev[:, 0]], s[ev[:, 1]]
        sb = np.where(rep, sa, s[np.where(rep, 0, ev[:, 2])])
        self.W = np.where(rep, self.conv * self.g0w / (8 * sp_ * sa * sa), self.conv * self.g0w / (4 * sp_ * sa * sb))
        self.kf = kappa_factor(geom)
        self.bnorm2 = (geom.b ** 2).sum(axis=0)

    def orbit_avg(self, g):
        return (np.bincount(self.orb, weights=g, minlength=self.n_orb) / self.osize)[self.orb]

    def rates(self, pp, tie=True):
        g = self.W * pp
        return self.orbit_avg(g) if tie else g

    def metrics(self, g, cols=(0, 2), want_lmin=True):
        C = self.ops.C(g)
        e = self.geom.ehat
        scale = np.trace(C) / self.geom.n
        Creg = C + scale * np.outer(e, e)
        out = {"r": np.diag(C).copy()}
        if len(cols):
            cf = sla.cho_factor(Creg, lower=True, check_finite=False)
            for c in cols:
                b = self.geom.b[:, c]
                x = sla.cho_solve(cf, b, check_finite=False)
                K, N = float(b @ x), float(x @ x)
                nm = "xyz"[c]
                out[f"K_{nm}"], out[f"kappa_{nm}"] = K, self.kf * K
                out[f"tau_{nm}"], out[f"tauCS_{nm}"] = N / K / FOURPI, K / self.bnorm2[c] / FOURPI
        if want_lmin:
            w = np.sort(np.linalg.eigvalsh(Creg))
            out["lam_min_H"] = float(w[0]) if w[0] < scale * (1 - 1e-9) else float(w[1])
            out["ceiling_F1_ps"] = (1.0 / out["lam_min_H"] / FOURPI) if out["lam_min_H"] > 0 else float("inf")
        return out


def summary(m):
    return {k: v for k, v in m.items() if not isinstance(v, np.ndarray)}
