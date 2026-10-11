"""Numerical audit helpers: an event-operator implementation written from the stated definitions only
(02-assumptions.md A/B; aln_events.py docstring), sharing no code with fastops / memgap / aln_geometry / fc3class.

Definitions reproduced here:
  modes      full-grid (q, j) with f > cutoff (the three acoustic Gamma modes are excluded)
  d_mu       sqrt(n_mu (1 + n_mu)), n = Bose occupation at T
  event      p -> a + b; conserving stoichiometry s = -l_p e_p + l_d (e_a + e_b), l_p = sqrt((f_a + f_b)/f_p),
             l_d = 1/l_p (repeated daughter: s = -l_p e_p + 2 l_d e_a); event vector a = D^-1 s
  rate       g = conv pp g0w / (4 s_p s_a s_b) (distinct) or / (8 s_p s_a^2) (repeated), s = sinh(f/(2 kT))
  operator   C = sum_alpha g_alpha a_alpha a_alpha^T;  e = D f (energy null vector);  b = D (f v)
  response   K = b^T C^+ b, N = |C^+ b|^2, tau[ps] = N / K / (4 pi);  kappa = kf K
  tying      orbit average of g over the 24 operations (rotations x time reversal) of rot_maps, computed here as the
             GROUP average (1/24) sum_R g(R alpha), which equals the orbit average for a group action.
C^+ b is evaluated by a full symmetric eigendecomposition with the null pair removed explicitly (no regularisation).
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "4")

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

FOURPI = 4.0 * np.pi
RES = r"C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results"
MESHES = {
    "m553": (RES + "/m553/events_m553_s0.1.npz", RES + "/m553"),
    "m775": (RES + "/m775x/events_m775_s0.1.npz", RES + "/m775x"),
    "m996": (RES + "/m996/events_m996_s0.1.npz", RES + "/m996"),
}


class Operator:
    def __init__(self, events_path):
        Z = np.load(events_path)
        f = Z["freqs"]
        self.Nq, self.nb = Nq, nb = f.shape
        self.mesh = Z["mesh"].tolist()
        self.sigma = float(Z["sigma"])
        valid = (f > float(Z["cutoff"])).ravel()
        newidx = -np.ones(Nq * nb, dtype=np.int64)
        newidx[valid] = np.arange(int(valid.sum()))
        self.valid, self.newidx = valid, newidx
        self.full_index = np.nonzero(valid)[0]
        self.f = f.ravel()[valid]
        self.n = n = int(valid.sum())
        self.T = T = float(Z["T"])
        self.kT = float(Z["KB"]) * T / float(Z["THzToEv"])
        x = self.f / self.kT
        nocc = 1.0 / np.expm1(x)
        self.d = np.sqrt(nocc * (1.0 + nocc))
        self.sh = np.sinh(x / 2.0)
        self.v = Z["gv"].reshape(Nq * nb, 3)[valid]
        self.q_of = self.full_index // nb
        self.band_of = self.full_index % nb
        self.addr = Z["addr"]
        self.lattice = Z["lattice"]
        # events (full indices kept for the symmetry lookup)
        self.ev_full = np.stack([Z["ev_p"], Z["ev_a"], Z["ev_b"]], axis=1).astype(np.int64)
        P, A, B = (newidx[self.ev_full[:, k]] for k in range(3))
        assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
        self.P, self.A, self.B = P, A, B
        self.rep = rep = A == B
        self.m = len(P)
        fp, fa, fb = self.f[P], self.f[A], self.f[B]
        lp = np.sqrt((fa + fb) / fp)
        ld = 1.0 / lp
        self.slots = np.stack([P, A, B], axis=1)
        s = np.stack([-lp, np.where(rep, 2.0 * ld, ld), np.where(rep, 0.0, ld)], axis=1)
        self.acoef = s / self.d[self.slots]
        self.dl = Z["ev_dl"]
        self.ev_g = Z["ev_g"]
        self.e = self.d * self.f
        self.ehat = self.e / np.linalg.norm(self.e)
        self.b = (self.d * self.f)[:, None] * self.v
        self.kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (
            2.0 * Nq * float(Z["KB"]) * T ** 2)
        self.rot = Z["rot_maps"]
        self._img = None

    # ------------------------------------------------------------------ rates
    def rates_from_pp(self, pp, g0w, conv):
        sp_, sa, sb = self.sh[self.P], self.sh[self.A], self.sh[self.B]
        return np.where(self.rep, conv * pp * g0w / (8.0 * sp_ * sa * sa), conv * pp * g0w / (4.0 * sp_ * sa * sb))

    def images(self):
        """(24, m) event index of R(alpha) for every symmetry operation R (full-index lookup, canonical a <= b)."""
        if self._img is not None:
            return self._img
        nb, M = self.nb, self.Nq * self.nb
        ev = self.ev_full
        key = (ev[:, 0] * M + ev[:, 1]) * M + ev[:, 2]
        order = np.argsort(key)
        skey = key[order]
        img = np.empty((len(self.rot), self.m), dtype=np.int64)
        for r, R in enumerate(self.rot):
            mp = R[ev // nb] * nb + ev % nb
            a2, b2 = np.minimum(mp[:, 1], mp[:, 2]), np.maximum(mp[:, 1], mp[:, 2])
            k2 = (mp[:, 0] * M + a2) * M + b2
            pos = np.searchsorted(skey, k2)
            pos = np.minimum(pos, len(skey) - 1)
            if not np.array_equal(skey[pos], k2):
                raise ValueError("event set not closed under the symmetry operations")
            img[r] = order[pos]
        self._img = img
        return img

    def group_average(self, g):
        img = self.images()
        return g[img].mean(axis=0)

    # ------------------------------------------------------------------ operator
    def C(self, g):
        n = self.n
        sl, ac = self.slots, self.acoef
        rows = np.repeat(sl, 3, axis=1).ravel()            # (m, 9): k slowest
        cols = np.tile(sl, (1, 3)).ravel()
        vals = (g[:, None] * (np.repeat(ac, 3, axis=1) * np.tile(ac, (1, 3)))).ravel()
        Cm = np.bincount(rows * n + cols, weights=vals, minlength=n * n).reshape(n, n)
        return 0.5 * (Cm + Cm.T)

    def diagC(self, g):
        return np.bincount(self.slots.ravel(), weights=(g[:, None] * self.acoef ** 2).ravel(), minlength=self.n)

    def ATx(self, x):
        return np.sum(self.acoef * x[self.slots], axis=1)


def eig_response(op, C, cols=(0, 2), keep_vectors=False):
    """K, N, tau from the full eigendecomposition; null pair identified by its overlap with e_hat."""
    lam, Q = sla.eigh(C, overwrite_a=False, check_finite=False, driver="evd")
    ov = Q.T @ op.ehat
    k0 = int(np.argmax(np.abs(ov)))
    mask = np.ones(len(lam), dtype=bool)
    mask[k0] = False
    beta = Q.T @ op.b[:, list(cols)]                       # (n, ncols)
    out = {"lam_null": float(lam[k0]), "null_overlap": float(abs(ov[k0])), "lam_max": float(lam[-1]),
           "n_negative_H": int(np.sum(lam[mask] <= 0)), "lam_sorted_H": None, "cols": {}}
    lH = lam[mask]
    out["lam_min_H"] = float(lH.min())
    out["lam_H_lowest"] = [float(v) for v in np.sort(lH)[:10]]
    out["cond_H"] = float(lam[-1] / lH.min())
    out["beta_on_null"] = [float(v) for v in beta[k0]]
    for i, c in enumerate(cols):
        bt = beta[mask, i]
        o = np.argsort(lH)
        cK = (bt ** 2 / lH)[o]
        cN = (bt ** 2 / lH ** 2)[o]
        K, N = float(cK.sum()), float(cN.sum())
        cumN, cumK = np.cumsum(cN) / N, np.cumsum(cK) / K
        rec = {"K": K, "N": N, "tau_ps": N / K / FOURPI, "kappa": op.kf * K,
               "N_frac_lowest": {str(k): float(cumN[k - 1]) for k in (1, 5, 20, 100, 500)},
               "K_frac_lowest": {str(k): float(cumK[k - 1]) for k in (1, 5, 20, 100, 500)},
               "n_pairs_for_N_50_90_99": [int(np.searchsorted(cumN, t) + 1) for t in (0.5, 0.9, 0.99)],
               "n_pairs_for_K_50_90_99": [int(np.searchsorted(cumK, t) + 1) for t in (0.5, 0.9, 0.99)],
               "top_N_contributors": [], }
        top = np.argsort(-cN)[:5]
        for t in top:
            rec["top_N_contributors"].append({"rank_by_lambda": int(t), "lambda": float(lH[o][t]),
                                              "frac_N": float(cN[t] / N), "frac_K": float(cK[t] / K)})
        out["cols"]["xyz"[c]] = rec
    if keep_vectors:
        out["_lam"], out["_Q"], out["_k0"] = lam, Q, k0
    return out


def chol_response(op, C, gamma, cols=(0, 2)):
    Creg = C + gamma * np.outer(op.ehat, op.ehat)
    cf = sla.cho_factor(Creg, lower=True, check_finite=False)
    res = {}
    for c in cols:
        b = op.b[:, c]
        x = sla.cho_solve(cf, b, check_finite=False)
        r = Creg @ x - b
        K, N = float(b @ x), float(x @ x)
        res["xyz"[c]] = {"K": K, "N": N, "tau_ps": N / K / FOURPI, "kappa": op.kf * K,
                         "rel_residual": float(np.linalg.norm(r) / np.linalg.norm(b)),
                         "x_dot_ehat_rel": float(abs(x @ op.ehat) / np.linalg.norm(x))}
    return res
