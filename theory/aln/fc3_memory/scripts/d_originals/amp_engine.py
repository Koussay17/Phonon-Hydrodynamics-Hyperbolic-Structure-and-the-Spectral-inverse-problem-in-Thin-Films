"""Independent NumPy implementation of the three-phonon amplitude Phi(q0 j0, q1 j1, q2 j2) from a compact fc3.

Re-implements phono3py 4.5.0's real-to-reciprocal transform (single-leg, or three-leg r0 average with the
all_shortest shortcut) and the reciprocal-to-normal contraction (no conjugation of eigenvectors,
division by sqrt(m_i m_j m_k) and sqrt(f0 f1 f2), zero below the cutoff frequency), vectorised over triplets.
pp = |Phi|^2 * unit_conversion reproduces phono3py's interaction strength (checked in t5_engine_check.py).

Linearity: Phi is linear in fc3, Phi(fc3) = sum_e L_e fc3_e. event_functionals() returns L for a set of
events as dense vectors over the compact fc3 elements (used for the linearisation in linearize.py).

Inputs: static_*.npz written by p3_events.py --static (phonons on phono3py's BZ grid, smallest vectors,
all_shortest, ...). Compact fc3 layout: (n_patom, n_satom, n_satom, 3, 3, 3), first index = primitive atom.
"""
from __future__ import annotations

import numpy as np


class AmpEngine:
    def __init__(self, static_npz, r0avg=True, unit_conversion=1.0, gamma_context="parent1"):
        S = np.load(static_npz)
        self.f = np.array(S["freqs_bz"])
        self.E = np.array(S["eigv_bz"])
        gp = int(S["gp_Gamma"])
        if gamma_context == "parentGamma":
            self.f[gp] = S["gamma_f_parentGamma"]
            self.E[gp] = S["gamma_e_parentGamma"]
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
                vs = svecs[adr:adr + mlt]
                PH[a, J] = np.exp(2j * np.pi * vs @ q.T).sum(axis=0) / mlt
        self.PH = PH
        self.msq = np.repeat(np.sqrt(self.masses), 3)

    # ---------------------------------------------------------------- real -> reciprocal
    def _prephase(self, tps):
        # G = q0 + q1 + q2 (BZ addresses) for every component. phono3py's pure-Python prototype
        # (RealToReciprocal) zeroes components unless all three addresses are non-zero; its Rust/C
        # kernels do not: with the Python rule 38 of 1613 Umklapp triplets (5x5x3) disagree with the
        # Rust interaction strength by up to a factor 400; with the full G all events agree to 7e-12.
        a3 = self.addr[tps]                                     # (T, 3 q's, 3 comps)
        sumq = a3.sum(axis=1) / self.Dd
        return np.exp(2j * np.pi * sumq @ self.pos.T)           # (T, npa)

    def _bilinear(self, F, qA, qB):
        """S[t, a, t1, t2, 27] = sum_{J in t1, K in t2} F[a, J, K] PH[a, J, qA_t] PH[a, K, qB_t]."""
        T = len(qA)
        npa, nsat = self.npa, self.nsat
        S = np.zeros((T, npa, npa, npa, 27), dtype=complex)
        for a in range(npa):
            PA = self.PH[a][:, qA].T                            # (T, nsat)
            for t2 in range(npa):
                K = self.groups[t2]
                Bm = self.PH[a][K][:, qB]                       # (|K|, T)
                Fa = F[a][:, K, :]                              # (nsat, |K|, 27)
                X = np.tensordot(Bm.T, Fa.transpose(1, 0, 2), axes=([1], [0]))   # (T, nsat, 27)
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
            S1 = self._bilinear(F, tps[:, 1], tps[:, 2])
            rec = pre[:, :, None, None, None] * S1
            return rec.reshape(T, npa, npa, npa, 3, 3, 3)
        w1 = np.where(self.allsh, 3.0, 1.0)[..., None]
        w23 = np.where(self.allsh, 0.0, 1.0)[..., None]
        S1 = self._bilinear(F * w1, tps[:, 1], tps[:, 2]).reshape(T, npa, npa, npa, 3, 3, 3)
        S2 = self._bilinear(F * w23, tps[:, 0], tps[:, 2]).reshape(T, npa, npa, npa, 3, 3, 3)
        S3 = self._bilinear(F * w23, tps[:, 1], tps[:, 0]).reshape(T, npa, npa, npa, 3, 3, 3)
        # leg 1: anchor i -> S1[t, i, j, k, a, b, c]
        rec = pre[:, :, None, None, None, None, None] * S1
        # leg 2: anchor j, elements [alpha_j, beta_i, gamma_k] at S2[t, j, i, k] -> transpose cart (1,0,2)
        rec += (pre[:, None, :, None, None, None, None]
                * S2.transpose(0, 2, 1, 3, 5, 4, 6))
        # leg 3: anchor k, elements [alpha_k, beta_j, gamma_i] at S3[t, k, j, i] -> transpose cart (2,1,0)
        rec += (pre[:, None, None, :, None, None, None]
                * S3.transpose(0, 3, 2, 1, 6, 5, 4))
        return rec / 3.0

    # ---------------------------------------------------------------- reciprocal -> normal
    def normal(self, rec, tps):
        """Phi[t, b0, b1, b2] (all band triples) from fc3_reciprocal at triplets tps."""
        T = len(tps)
        npa = self.npa
        nb = 3 * npa
        R = rec.transpose(0, 1, 4, 2, 5, 3, 6).reshape(T, nb, nb, nb)
        R = R / (self.msq[None, :, None, None] * self.msq[None, None, :, None] * self.msq[None, None, None, :])
        e0, e1, e2 = self.E[tps[:, 0]], self.E[tps[:, 1]], self.E[tps[:, 2]]
        X = np.einsum("tabc,tai->tibc", R, e0, optimize=True)
        X = np.einsum("tibc,tbj->tijc", X, e1, optimize=True)
        X = np.einsum("tijc,tck->tijk", X, e2, optimize=True)
        f0, f1, f2 = self.f[tps[:, 0]], self.f[tps[:, 1]], self.f[tps[:, 2]]
        ok = ((f0 > self.cutoff)[:, :, None, None] & (f1 > self.cutoff)[:, None, :, None]
              & (f2 > self.cutoff)[:, None, None, :])
        fff = np.sqrt(np.abs(f0[:, :, None, None] * f1[:, None, :, None] * f2[:, None, None, :]))
        return np.where(ok, X / np.where(ok, fff, 1.0), 0.0)

    def phi_events(self, fc3c, tp, bands, chunk=400):
        """Complex Phi for events given by triplet BZ indices tp (m,3) and band triples bands (m,3)."""
        utp, inv = np.unique(tp, axis=0, return_inverse=True)
        inv = inv.ravel()
        out = np.zeros(len(tp), dtype=complex)
        for s in range(0, len(utp), chunk):
            sl = slice(s, min(s + chunk, len(utp)))
            rec = self.fc3_reciprocal(fc3c, utp[sl])
            Phi = self.normal(rec, utp[sl])
            sel = np.nonzero((inv >= sl.start) & (inv < sl.stop))[0]
            out[sel] = Phi[inv[sel] - sl.start, bands[sel, 0], bands[sel, 1], bands[sel, 2]]
        return out

    def phi_blocks(self, fc3c, utp, chunk=400):
        """Full Phi[t, 12, 12, 12] at the given unique triplets."""
        res = []
        for s in range(0, len(utp), chunk):
            rec = self.fc3_reciprocal(fc3c, utp[s:s + chunk])
            res.append(self.normal(rec, utp[s:s + chunk]))
        return np.concatenate(res, axis=0)

    # ---------------------------------------------------------------- linear functionals
    def event_functionals(self, tps, bands):
        """L[m, e]: Phi_m = sum_e L[m, e] fc3c_flat[e] for events m (dense, (m, npa*nsat*nsat*27))."""
        npa, nsat = self.npa, self.nsat
        nb = 3 * npa
        m = len(tps)
        pre = self._prephase(tps)                                    # (m, npa)
        e0 = self.E[tps[:, 0]][np.arange(m), :, bands[:, 0]]         # (m, 12 comps)
        e1 = self.E[tps[:, 1]][np.arange(m), :, bands[:, 1]]
        e2 = self.E[tps[:, 2]][np.arange(m), :, bands[:, 2]]
        f0 = self.f[tps[:, 0], bands[:, 0]]
        f1 = self.f[tps[:, 1], bands[:, 1]]
        f2 = self.f[tps[:, 2], bands[:, 2]]
        ok = (f0 > self.cutoff) & (f1 > self.cutoff) & (f2 > self.cutoff)
        scale = np.where(ok, 1.0 / np.sqrt(np.abs(f0 * f1 * f2) + (~ok)), 0.0)
        u0 = (e0 / self.msq[None, :]).reshape(m, npa, 3)             # mass-scaled eigenvector components
        u1 = (e1 / self.msq[None, :]).reshape(m, npa, 3)
        u2 = (e2 / self.msq[None, :]).reshape(m, npa, 3)
        ty = self.ptype
        L = np.zeros((m, npa, nsat, nsat, 3, 3, 3), dtype=complex)
        if not self.r0avg:
            for a in range(npa):
                P1 = self.PH[a][:, tps[:, 1]].T                      # (m, nsat) phase of J with q1
                P2 = self.PH[a][:, tps[:, 2]].T
                # Phi += pre_a sum_{J,K} f[a,J,K,al,be,ga] P1_J P2_K u0[a,al] u1[type J, be] u2[type K, ga]
                L[:, a] = np.einsum("m,mJ,mK,ma,mJb,mKc->mJKabc", pre[:, a], P1, P2, u0[:, a],
                                    u1[:, ty], u2[:, ty], optimize=True)
        else:
            w1 = np.where(self.allsh, 3.0, 1.0)
            w23 = np.where(self.allsh, 0.0, 1.0)
            for a in range(npa):
                # leg 1: anchor i = a; J ~ j with q1, K ~ k with q2; element [al_i, be_j, ga_k]
                P1 = self.PH[a][:, tps[:, 1]].T
                P2 = self.PH[a][:, tps[:, 2]].T
                L1 = np.einsum("m,mJ,mK,ma,mJb,mKc->mJKabc", pre[:, a], P1, P2, u0[:, a],
                               u1[:, ty], u2[:, ty], optimize=True) * w1[a][None, :, :, None, None, None]
                # leg 2: anchor j = a; J ~ i with q0, K ~ k with q2; element [al_j, be_i, ga_k]
                Q0 = self.PH[a][:, tps[:, 0]].T
                L2 = np.einsum("m,mJ,mK,ma,mJb,mKc->mJKabc", pre[:, a], Q0, P2, u1[:, a],
                               u0[:, ty], u2[:, ty], optimize=True) * w23[a][None, :, :, None, None, None]
                # leg 3: anchor k = a; J ~ j with q1, K ~ i with q0; element [al_k, be_j, ga_i]
                L3 = np.einsum("m,mJ,mK,ma,mJb,mKc->mJKabc", pre[:, a], P1, Q0, u2[:, a],
                               u1[:, ty], u0[:, ty], optimize=True) * w23[a][None, :, :, None, None, None]
                L[:, a] = (L1 + L2 + L3) / 3.0
        L *= scale[:, None, None, None, None, None, None]
        return L.reshape(m, -1)
