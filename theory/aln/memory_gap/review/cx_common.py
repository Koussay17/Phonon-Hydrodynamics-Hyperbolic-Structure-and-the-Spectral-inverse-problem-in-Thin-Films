"""Independent re-implementation (counterexample review) of the AlN event-cone objects.

Written from the definitions in 02-assumptions.md only; it does not import memgap.py, fastops.py,
symmetry.py or optim*.py. The single exception is `their_ev_orb`, used ONLY to decode the orbit
variables stored in the campaign's witness files (their orbit numbering); every decoded witness is
re-checked here for orbit-constancy against the orbits computed independently in this module.

Objects (entropy coordinates):
  modes mu: full-grid (q, j) with f > cutoff;  eps = f [THz];  x_mu = h f / k T
  D_mu = sqrt(N(1+N)) = 1/(2 sinh(x/2));  e = D eps;  b = D (eps v)  (3 columns)
  event alpha = (p, a, b): s = -l_p e_p + l_d (e_a + e_b), l_p = sqrt((f_a+f_b)/f_p), l_d = 1/l_p
  a_alpha = D^{-1} s_alpha ;  C(g) = sum g a a^T ;  r = diag C(g)
  K = b^T C^+ b ; N = |C^+ b|^2 ; tau = N/K ; tau[ps] = tau / (4 pi)
Temperature dependence (T-independent |Phi|^2 delta):
  g_alpha(T) = w_alpha * Bf_alpha(T),  Bf = 1/(4 sinh(x_p/2) sinh(x_a/2) sinh(x_b/2))   (a != b)
                                        Bf = 1/(8 sinh(x_p/2) sinh(x_a/2)^2)             (a == b)
  (aln_events.py: g = conv*pp*g0 * Bf, conv*pp*g0 is temperature independent)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

CAMP = Path(r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap")
EVENTS = CAMP / "results" / "aln" / "events_m553_s0.1.npz"
SCAN = CAMP / "results" / "aln" / "scan_m553_s0.1"
PHYSOP = CAMP / "results" / "aln" / "physop_m553_s0.1.npz"
FOURPI = 4.0 * np.pi


class AlN:
    def __init__(self, path=EVENTS, T=None):
        Z = np.load(path)
        self.Z = Z
        f = Z["freqs"]
        Nq, nb = f.shape
        self.Nq, self.nb = Nq, nb
        cutoff = float(Z["cutoff"])
        valid = (f > cutoff).ravel()
        newidx = -np.ones(Nq * nb, dtype=np.int64)
        newidx[valid] = np.arange(valid.sum())
        self.valid, self.newidx = valid, newidx
        self.eps = f.ravel()[valid].astype(float)
        self.n = len(self.eps)
        self.T0 = float(Z["T"])
        self.THzToEv, self.KB = float(Z["THzToEv"]), float(Z["KB"])
        self.vel = Z["gv"].reshape(Nq * nb, 3)[valid]
        P, A, B = newidx[Z["ev_p"]], newidx[Z["ev_a"]], newidx[Z["ev_b"]]
        assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
        self.P, self.A, self.B = P, A, B
        self.rep = A == B
        self.m = len(P)
        self.g_raw = Z["ev_g"].astype(float)
        self.delta = Z["ev_dl"].astype(float)
        self.umk = Z["ev_umk"].astype(bool)
        e = self.eps
        lp = np.sqrt((e[A] + e[B]) / e[P])
        self.lp, self.ld = lp, 1.0 / lp
        # stoichiometric coefficients (T independent)
        self.sP = -lp
        self.sA = np.where(self.rep, 2.0 * self.ld, self.ld)
        self.sB = np.where(self.rep, 0.0, self.ld)
        self.unit_to_WmK = float(Z["unit_to_WmK"])
        self.volume = float(Z["volume"])
        self.set_T(self.T0 if T is None else T)
        self._orbits = None

    # ------------------------------------------------------------------ temperature
    def xk(self, T):
        return self.eps * self.THzToEv / (self.KB * T)

    def set_T(self, T):
        self.T = float(T)
        x = self.xk(T)
        self.sh = np.sinh(x / 2.0)
        self.D = 1.0 / (2.0 * self.sh)                      # sqrt(N(1+N))
        self.e = self.D * self.eps
        self.ehat = self.e / np.linalg.norm(self.e)
        self.b = (self.D * self.eps)[:, None] * self.vel
        # event matrix A_T (n x m), columns a_alpha = D^{-1} s_alpha
        Dinv = 1.0 / self.D
        rows = np.concatenate([self.P, self.A, self.B[~self.rep]])
        cols = np.concatenate([np.arange(self.m), np.arange(self.m), np.nonzero(~self.rep)[0]])
        vals = np.concatenate([self.sP * Dinv[self.P], self.sA * Dinv[self.A], (self.sB * Dinv[self.B])[~self.rep]])
        self.Amat = sp.csc_matrix((vals, (rows, cols)), shape=(self.n, self.m))
        self.AT = self.Amat.T.tocsr()
        self.W = self.Amat.multiply(self.Amat).tocsr()     # diag C(g) = W g
        sh = self.sh
        self.Bf = np.where(self.rep, 1.0 / (8 * sh[self.P] * sh[self.A] ** 2),
                           1.0 / (4 * sh[self.P] * sh[self.A] * sh[self.B]))
        return self

    def kappa_factor(self):
        return self.unit_to_WmK / self.volume * self.THzToEv ** 2 / (2 * self.Nq * self.KB * self.T ** 2)

    # ------------------------------------------------------------------ symmetry (independent)
    def orbits(self):
        if self._orbits is not None:
            return self._orbits
        rot = self.Z["rot_maps"]
        nb = self.nb
        mode_maps = []
        for R in rot:
            full = (R[:, None] * nb + np.arange(nb)[None, :]).ravel()
            mp = self.newidx[full][self.valid]
            assert (mp >= 0).all()
            mode_maps.append(mp)
        mode_maps = np.array(mode_maps)
        # energies invariant
        for mp in mode_maps:
            assert np.allclose(self.eps[mp], self.eps, rtol=1e-6, atol=1e-6)
        n = self.n
        key = lambda p, a, b: (p * (n + 1) + np.minimum(a, b)) * (n + 1) + np.maximum(a, b)
        k0 = key(self.P, self.A, self.B)
        order = np.argsort(k0)
        ks = k0[order]
        assert np.all(np.diff(ks) > 0), "duplicate events"
        imgs = np.empty((len(mode_maps), self.m), dtype=np.int64)
        for i, mp in enumerate(mode_maps):
            ki = key(mp[self.P], mp[self.A], mp[self.B])
            pos = np.searchsorted(ks, ki)
            assert np.all(pos < len(ks)) and np.all(ks[pos] == ki), "event set not invariant"
            imgs[i] = order[pos]
        ev_rep = imgs.min(axis=0)
        # orbit = closure; for a group the image set is the orbit, check consistency
        assert np.all(ev_rep[imgs] == ev_rep[None, :]), "orbit representative not consistent (not a group?)"
        uniq, ev_orb = np.unique(ev_rep, return_inverse=True)
        mode_rep = mode_maps.min(axis=0)
        assert np.all(mode_rep[mode_maps] == mode_rep[None, :])
        muniq, mode_orb = np.unique(mode_rep, return_inverse=True)
        self._orbits = dict(ev_orb=ev_orb, n_ev_orb=len(uniq), mode_orb=mode_orb, n_mode_orb=len(muniq),
                            mode_reps=muniq, mode_maps=mode_maps, ev_imgs=imgs)
        return self._orbits

    def orbit_average(self, g):
        o = self.orbits()
        s = np.bincount(o["ev_orb"], weights=g, minlength=o["n_ev_orb"])
        c = np.bincount(o["ev_orb"], minlength=o["n_ev_orb"])
        return (s / c)[o["ev_orb"]]

    def g_ref(self):
        return self.orbit_average(self.g_raw)

    # ------------------------------------------------------------------ operator
    def C(self, g):
        Ag = self.Amat @ sp.diags(g)
        return (Ag @ self.AT).toarray()

    def diag(self, g):
        return self.W @ g

    def solve(self, g, cols=(0, 2), C=None):
        """x = C^+ b via Cholesky of C + gam e e^T (valid for admissible g, b orthogonal to e)."""
        if C is None:
            C = self.C(g)
        gam = np.trace(C) / self.n
        cf = sla.cho_factor(C + gam * np.outer(self.ehat, self.ehat), lower=True, check_finite=False)
        Bm = self.b[:, list(cols)]
        X = sla.cho_solve(cf, Bm, check_finite=False)
        return C, cf, X

    def moments(self, g, cols=(0, 2)):
        C, cf, X = self.solve(g, cols)
        out = {}
        for k, c in enumerate(cols):
            b = self.b[:, c]
            x = X[:, k]
            K = float(b @ x)
            N = float(x @ x)
            out[c] = dict(K=K, N=N, tau=N / K, tau_ps=N / K / FOURPI, x=x)
        return out

    def moments_eig(self, g, cols=(0, 2)):
        """Second, independent route: orthonormal basis of H = e^perp, eigendecomposition of Q^T C Q."""
        C = self.C(g)
        Q = basis_H(self.ehat)
        CH = Q.T @ C @ Q
        CH = 0.5 * (CH + CH.T)
        w, U = np.linalg.eigh(CH)
        out = {"lam_min_H": float(w[0]), "lam_max": float(w[-1]), "cond_H": float(w[-1] / w[0])}
        for c in cols:
            bh = U.T @ (Q.T @ self.b[:, c])
            K = float(np.sum(bh ** 2 / w))
            N = float(np.sum(bh ** 2 / w ** 2))
            out[c] = dict(K=K, N=N, tau=N / K, tau_ps=N / K / FOURPI)
        out["w"], out["U"], out["Q"] = w, U, Q
        return out


def basis_H(ehat):
    """Orthonormal basis of e^perp via a Householder reflection (exact up to rounding)."""
    n = len(ehat)
    v = ehat.copy()
    v[0] += np.sign(v[0]) if v[0] != 0 else 1.0
    v /= np.linalg.norm(v)
    Hh = np.eye(n) - 2.0 * np.outer(v, v)        # Hh e_1 = -sign * ehat
    return Hh[:, 1:]


def their_ev_orb():
    """Orbit numbering used by the campaign's witness files (decoding only)."""
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(CAMP / "scripts"))
    from aln_geometry import load       # noqa
    from symmetry import orbits          # noqa
    geom = load(str(EVENTS))
    tie = orbits(geom, geom.labels["maps"])
    return np.asarray(tie["ev_orb"]), geom


def load_witness(name, ev_orb_theirs):
    gam = np.load(SCAN / f"{name}.npz")["gamma"]
    return gam[ev_orb_theirs] if len(gam) != len(ev_orb_theirs) else gam
