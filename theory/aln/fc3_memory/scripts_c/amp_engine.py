"""Independent NumPy implementation of the three-phonon amplitude Phi(q0 j0, q1 j1, q2 j2) from a compact fc3
(C: mirror copy of the validated D: version).

Re-implements phono3py 4.5.0's real-to-reciprocal transform (single leg, or three-leg r0 average with the
all_shortest shortcut) and the reciprocal-to-normal contraction (no conjugation of eigenvectors, division by
sqrt(m_i m_j m_k) and sqrt(f0 f1 f2), zero below the cutoff frequency). Pre-phase: exp(2 pi i G.r_i) with the FULL
G = q0 + q1 + q2 (BZ addresses) in every component; this reproduces the Rust kernel for every event (7e-12);
phono3py's pure-Python prototype zeroes components unless all three addresses are non-zero and disagrees for
38 of 1,613 Umklapp triplets at 5x5x3.
Phi is linear in fc3: event_functionals() returns L with Phi_m = L[m] . fc3_compact.ravel().
"""
from __future__ import annotations

import numpy as np


class AmpEngine:
    def __init__(self, static_npz, r0avg=True, unit_conversion=1.0):
        S = np.load(static_npz)
        self.f = np.array(S["freqs_bz"])
        self.E = np.array(S["eigv_bz"])
        self.addr = np.array(S["bz_addresses"])
        self.Dd = np.array(S["D_diag"]).astype(float)
        svecs, multi = S["svecs"], S["multi"]
        p2s, s2p = S["p2s"], S["s2p"]
        self.masses = np.array(S["masses"])
        self.pos = np.array(S["prim_scaled"])
        self.allsh = np.array(S["all_shortest"]).astype(bool)
        self.cutoff = float(S["cutoff_frequency"])
        self.r0avg = bool(r0avg)
        self.uc = float(unit_conversion)
        self.npa = npa = len(p2s)
        self.nsat = nsat = len(s2p)
        self.ptype = np.array([int(np.nonzero(p2s == s)[0][0]) for s in s2p])
        self.groups = [np.nonzero(self.ptype == t)[0] for t in range(npa)]
        q = self.addr / self.Dd
        PH = np.zeros((npa, nsat, len(q)), dtype=complex)
        for a in range(npa):
            for J in range(nsat):
                mlt, adr = int(multi[J, a, 0]), int(multi[J, a, 1])
                PH[a, J] = np.exp(2j * np.pi * svecs[adr:adr + mlt] @ q.T).sum(axis=0) / mlt
        self.PH = PH
        self.msq = np.repeat(np.sqrt(self.masses), 3)

    def _prephase(self, tps):
        sumq = self.addr[tps].sum(axis=1) / self.Dd
        return np.exp(2j * np.pi * sumq @ self.pos.T)

    def _bilinear(self, F, qA, qB):
        T = len(qA)
        npa = self.npa
        S = np.zeros((T, npa, npa, npa, 27), dtype=complex)
        for a in range(npa):
            PA = self.PH[a][:, qA].T
            for t2 in range(npa):
                K = self.groups[t2]
                X = np.tensordot(self.PH[a][K][:, qB].T, F[a][:, K, :].transpose(1, 0, 2), axes=([1], [0]))
                X *= PA[:, :, None]
                for t1 in range(npa):
                    S[:, a, t1, t2, :] = X[:, self.groups[t1], :].sum(axis=1)
        return S

    def fc3_reciprocal(self, fc3c, tps):
        npa, nsat = self.npa, self.nsat
        F = np.asarray(fc3c, dtype=float).reshape(npa, nsat, nsat, 27)
        pre = self._prephase(tps)
        T = len(tps)
        if not self.r0avg:
            return (pre[:, :, None, None, None] * self._bilinear(F, tps[:, 1], tps[:, 2])).reshape(T, npa, npa, npa, 3, 3, 3)
        w1 = np.where(self.allsh, 3.0, 1.0)[..., None]
        w23 = np.where(self.allsh, 0.0, 1.0)[..., None]
        S1 = self._bilinear(F * w1, tps[:, 1], tps[:, 2]).reshape(T, npa, npa, npa, 3, 3, 3)
        S2 = self._bilinear(F * w23, tps[:, 0], tps[:, 2]).reshape(T, npa, npa, npa, 3, 3, 3)
        S3 = self._bilinear(F * w23, tps[:, 1], tps[:, 0]).reshape(T, npa, npa, npa, 3, 3, 3)
        rec = pre[:, :, None, None, None, None, None] * S1
        rec += pre[:, None, :, None, None, None, None] * S2.transpose(0, 2, 1, 3, 5, 4, 6)
        rec += pre[:, None, None, :, None, None, None] * S3.transpose(0, 3, 2, 1, 6, 5, 4)
        return rec / 3.0

    def normal(self, rec, tps):
        T = len(tps)
        nb = 3 * self.npa
        R = rec.transpose(0, 1, 4, 2, 5, 3, 6).reshape(T, nb, nb, nb)
        R = R / (self.msq[None, :, None, None] * self.msq[None, None, :, None] * self.msq[None, None, None, :])
        X = np.einsum("tabc,tai->tibc", R, self.E[tps[:, 0]], optimize=True)
        X = np.einsum("tibc,tbj->tijc", X, self.E[tps[:, 1]], optimize=True)
        X = np.einsum("tijc,tck->tijk", X, self.E[tps[:, 2]], optimize=True)
        f0, f1, f2 = self.f[tps[:, 0]], self.f[tps[:, 1]], self.f[tps[:, 2]]
        ok = ((f0 > self.cutoff)[:, :, None, None] & (f1 > self.cutoff)[:, None, :, None]
              & (f2 > self.cutoff)[:, None, None, :])
        fff = np.sqrt(np.abs(f0[:, :, None, None] * f1[:, None, :, None] * f2[:, None, None, :]))
        return np.where(ok, X / np.where(ok, fff, 1.0), 0.0)

    def phi_events(self, fc3c, tp, bands, chunk=400):
        utp, inv = np.unique(tp, axis=0, return_inverse=True)
        inv = inv.ravel()
        out = np.zeros(len(tp), dtype=complex)
        for s in range(0, len(utp), chunk):
            sl = slice(s, min(s + chunk, len(utp)))
            Phi = self.normal(self.fc3_reciprocal(fc3c, utp[sl]), utp[sl])
            sel = np.nonzero((inv >= sl.start) & (inv < sl.stop))[0]
            out[sel] = Phi[inv[sel] - sl.start, bands[sel, 0], bands[sel, 1], bands[sel, 2]]
        return out

    def event_functionals(self, tps, bands):
        npa, nsat = self.npa, self.nsat
        m = len(tps)
        pre = self._prephase(tps)
        ar = np.arange(m)
        e0 = self.E[tps[:, 0]][ar, :, bands[:, 0]]
        e1 = self.E[tps[:, 1]][ar, :, bands[:, 1]]
        e2 = self.E[tps[:, 2]][ar, :, bands[:, 2]]
        f0 = self.f[tps[:, 0], bands[:, 0]]
        f1 = self.f[tps[:, 1], bands[:, 1]]
        f2 = self.f[tps[:, 2], bands[:, 2]]
        ok = (f0 > self.cutoff) & (f1 > self.cutoff) & (f2 > self.cutoff)
        scale = np.where(ok, 1.0 / np.sqrt(np.abs(f0 * f1 * f2) + (~ok)), 0.0)
        u0 = (e0 / self.msq[None, :]).reshape(m, npa, 3)
        u1 = (e1 / self.msq[None, :]).reshape(m, npa, 3)
        u2 = (e2 / self.msq[None, :]).reshape(m, npa, 3)
        ty = self.ptype
        L = np.zeros((m, npa, nsat, nsat, 3, 3, 3), dtype=complex)
        ein = "m,mJ,mK,ma,mJb,mKc->mJKabc"
        for a in range(npa):
            P1 = self.PH[a][:, tps[:, 1]].T
            P2 = self.PH[a][:, tps[:, 2]].T
            if not self.r0avg:
                L[:, a] = np.einsum(ein, pre[:, a], P1, P2, u0[:, a], u1[:, ty], u2[:, ty], optimize=True)
                continue
            w1 = np.where(self.allsh[a], 3.0, 1.0)[None, :, :, None, None, None]
            w23 = np.where(self.allsh[a], 0.0, 1.0)[None, :, :, None, None, None]
            Q0 = self.PH[a][:, tps[:, 0]].T
            L1 = np.einsum(ein, pre[:, a], P1, P2, u0[:, a], u1[:, ty], u2[:, ty], optimize=True) * w1
            L2 = np.einsum(ein, pre[:, a], Q0, P2, u1[:, a], u0[:, ty], u2[:, ty], optimize=True) * w23
            L3 = np.einsum(ein, pre[:, a], P1, Q0, u2[:, a], u1[:, ty], u0[:, ty], optimize=True) * w23
            L[:, a] = (L1 + L2 + L3) / 3.0
        L *= scale[:, None, None, None, None, None, None]
        return L.reshape(m, -1)
