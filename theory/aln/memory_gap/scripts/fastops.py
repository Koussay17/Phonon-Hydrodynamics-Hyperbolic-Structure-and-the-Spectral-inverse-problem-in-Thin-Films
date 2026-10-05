"""Fast assembly for three-mode event geometries (no sparse-matrix temporaries).

Each event alpha touches at most 3 modes (slot value 0 and index n for an absent slot), so
  C(g)            = sum_alpha g_alpha a_alpha a_alpha^T        (9-pair bincount),
  W g             = diag C(g)                                   (3-slot bincount),
  W diag(g) W^T   = sum_alpha g_alpha w_alpha w_alpha^T         (9-pair bincount, w = a**2),
  A^T x           = per-event gather.
Tied variables gamma (one per orbit): g = gamma[orb]; P v = bincount(orb, v).
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla


class FastOps:
    def __init__(self, geom, tie=None):
        n, m = geom.n, geom.m
        self.n, self.m = n, m
        idx = geom.events.copy()
        vals = np.zeros((m, 3))
        for s in range(3):
            ok = idx[:, s] >= 0
            vals[ok, s] = geom.coef[ok, s] / geom.dscale[idx[ok, s]]
            idx[~ok, s] = n                       # dummy slot
        self.idx = idx
        self.vals = vals
        self.vals2 = vals * vals
        n1 = n + 1
        self.n1 = n1
        # 9 (k,l) pairs flattened indices into (n+1)^2
        self.pair_flat = np.stack([idx[:, k] * n1 + idx[:, l] for k in range(3) for l in range(3)], axis=1)
        self.pair_val = np.stack([vals[:, k] * vals[:, l] for k in range(3) for l in range(3)], axis=1)
        self.ehat = geom.ehat
        self.b = geom.b
        if tie is None:
            self.orb = None
            self.m_var = m
            self.rows = np.arange(n)
        else:
            self.orb = np.asarray(tie["ev_orb"], dtype=np.int64)
            self.m_var = int(self.orb.max()) + 1
            self.rows = np.asarray(tie["reps"], dtype=np.int64)
        self.n_rows = len(self.rows)
        # row map: full mode index -> position among constraint rows (or -1)
        rowpos = -np.ones(n1, dtype=np.int64)
        rowpos[self.rows] = np.arange(self.n_rows)
        self.rowpos = rowpos
        # constraint "column" structure in the variable space
        if self.orb is not None:
            # dense W_gamma (n_rows x m_var) if affordable
            Wg = np.zeros((self.n_rows, self.m_var))
            for s in range(3):
                rp = rowpos[idx[:, s]]
                ok = rp >= 0
                np.add.at(Wg, (rp[ok], self.orb[ok]), self.vals2[ok, s])
            self.Wg_dense = Wg
        else:
            self.Wg_dense = None
            nr1 = self.n_rows + 1
            rp = np.where(rowpos[idx] >= 0, rowpos[idx], self.n_rows)
            self.wpair_flat = np.stack([rp[:, k] * nr1 + rp[:, l] for k in range(3) for l in range(3)], axis=1)
            self.wpair_val = np.stack([self.vals2[:, k] * self.vals2[:, l] for k in range(3) for l in range(3)], axis=1)
            self.rp = rp

    # ------------------------------------------------------------ variable maps
    def g_of(self, gam):
        return gam if self.orb is None else gam[self.orb]

    def P_dot(self, v):
        return v if self.orb is None else np.bincount(self.orb, weights=v, minlength=self.m_var)

    def gamma_of(self, g):
        if self.orb is None:
            return g.copy()
        sizes = np.bincount(self.orb, minlength=self.m_var)
        return self.P_dot(g) / sizes

    # ------------------------------------------------------------ operator
    def C(self, g):
        n1 = self.n1
        w = (g[:, None] * self.pair_val).ravel()
        Cf = np.bincount(self.pair_flat.ravel(), weights=w, minlength=n1 * n1)
        return Cf.reshape(n1, n1)[: self.n, : self.n]

    def diag_full(self, g):
        d = np.bincount(self.idx.ravel(), weights=(g[:, None] * self.vals2).ravel(), minlength=self.n1)
        return d[: self.n]

    def ATx(self, x):
        xe = np.append(x, 0.0)
        return np.sum(self.vals * xe[self.idx], axis=1)

    # ------------------------------------------------------------ constraints (variable space)
    def Wv(self, gam_like):
        """W_gamma @ v for v in variable space (constraint rows)."""
        if self.Wg_dense is not None:
            return self.Wg_dense @ gam_like
        d = np.bincount(self.rp.ravel(), weights=(gam_like[:, None] * self.vals2).ravel(),
                        minlength=self.n_rows + 1)
        return d[: self.n_rows]

    def WTv(self, lam):
        """W_gamma^T @ lam (variable space)."""
        if self.Wg_dense is not None:
            return self.Wg_dense.T @ lam
        le = np.append(lam, 0.0)
        return np.sum(self.vals2 * le[self.rp], axis=1)

    def WGW(self, gam):
        """W_gamma diag(gam) W_gamma^T (n_rows x n_rows)."""
        if self.Wg_dense is not None:
            return (self.Wg_dense * gam[None, :]) @ self.Wg_dense.T
        nr1 = self.n_rows + 1
        w = (gam[:, None] * self.wpair_val).ravel()
        M = np.bincount(self.wpair_flat.ravel(), weights=w, minlength=nr1 * nr1).reshape(nr1, nr1)
        return M[: self.n_rows, : self.n_rows]

    # ------------------------------------------------------------ response
    def response(self, gam, col=0, want_grad=False):
        g = self.g_of(gam)
        C = self.C(g)
        gamma_reg = max(np.trace(C) / self.n, 1e-300)
        Creg = C + gamma_reg * np.outer(self.ehat, self.ehat)
        cf = sla.cho_factor(Creg, lower=True, check_finite=False)
        b = self.b[:, col]
        x = sla.cho_solve(cf, b, check_finite=False)
        K = float(b @ x)
        N = float(x @ x)
        out = {"K": K, "N": N, "tau": N / K, "x": x}
        if want_grad:
            u = sla.cho_solve(cf, x, check_finite=False)
            p = self.ATx(x)
            q = self.ATx(u)
            out["dK"] = self.P_dot(-(p * p))
            out["dN"] = self.P_dot(-2.0 * p * q)
            out["p2"] = self.P_dot(p * p)
            out["p"] = p
        return out
