"""Exact-rational checks of the cone statements T1, T2, P3, P4, P6 (13-final-report.md, section 2)
and of the remark "the cone bound is sharp only without the lifetime constraint".

Setting reproduced in exact arithmetic (fractions.Fraction):
  integer energies eps_mu > 0, integer stoichiometry s_alpha (parent -1, daughters +1 or +2),
  s_alpha . eps = 0 exactly; positive rational entropy scales d_mu (the statements are algebraic and
  do not use the Bose form of d_mu); a_alpha = D^{-1} s_alpha, e = D eps, H = e^perp;
  b in H; C(g) = sum g a a^T; x = C^+ b computed as the solution of (C + e e^T) x = b
  (valid for admissible g and b in H); K = b.x, N = x.x, tau = N/K.
Basic solutions y_I (|I| = n - 2): [a_I^T; e^T; b^T] y = (0,...,0,0,1), when nonsingular.
M^2 = max_I |y_I|^2 is an exact rational.
Cauchy-Binet vectors: V_I with V_I . y = det[a_I, y, e] for all y (= |e| v_I in note 18's notation).

Run:  python -B check_cone_bounds.py  (writes check_cone_bounds.json next to this file)
"""
from __future__ import annotations

import itertools
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
OUT = Path(__file__).with_suffix(".json")


# ----------------------------------------------------------------------------- exact linear algebra
def det(Mx):
    A = [list(r) for r in Mx]
    n = len(A)
    d = F(1)
    for c in range(n):
        piv = next((i for i in range(c, n) if A[i][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[c][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return d


def solve(A, b):
    n = len(A)
    Mx = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        piv = next((i for i in range(c, n) if Mx[i][c] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular")
        Mx[c], Mx[piv] = Mx[piv], Mx[c]
        for i in range(n):
            if i != c and Mx[i][c] != 0:
                f = Mx[i][c] / Mx[c][c]
                Mx[i] = [a - f * bb for a, bb in zip(Mx[i], Mx[c])]
    return [Mx[i][n] / Mx[i][i] for i in range(n)]


def rank(rows):
    A = [[F(v) for v in r] for r in rows]          # exact: never fall back to float division
    nr, nc = len(A), len(A[0])
    rk = 0
    for c in range(nc):
        piv = next((i for i in range(rk, nr) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(nr):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / A[rk][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[rk])]
        rk += 1
    return rk


def nullvec(rows, n):
    """One nonzero vector in the null space of rows (k x n, k = n - 1, full rank)."""
    A = [list(r) for r in rows]
    k = len(A)
    piv_cols = []
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, k) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [a / pv for a in A[r]]
        for i in range(k):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
        if r == k:
            break
    free = [c for c in range(n) if c not in piv_cols]
    assert len(free) == 1
    fc = free[0]
    v = [F(0)] * n
    v[fc] = F(1)
    for i, c in enumerate(piv_cols):
        v[c] = -A[i][fc]
    return v


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


# ----------------------------------------------------------------------------- model
def all_events(eps):
    n = len(eps)
    ev = []
    for p in range(n):
        for a in range(n):
            for b in range(a, n):
                if p in (a, b):
                    continue
                if eps[p] == eps[a] + eps[b]:
                    ev.append((p, a, b))
    return ev


def stoich(ev, n):
    p, a, b = ev
    s = [0] * n
    s[p] -= 1
    s[a] += 1
    s[b] += 1          # a == b gives the repeated-daughter coefficient +2
    return s


class Instance:
    def __init__(self, eps, events, dsc, b=None, w_seed=0):
        self.n = n = len(eps)
        self.eps = [F(x) for x in eps]
        self.events = list(events)
        self.m = len(events)
        self.D = [F(x) for x in dsc]
        self.S = [stoich(e, n) for e in events]                        # m rows
        self.A = [[F(s[i]) / self.D[i] for i in range(n)] for s in self.S]
        self.e = [self.D[i] * self.eps[i] for i in range(n)]
        assert all(dot(a, self.e) == 0 for a in self.A), "energy conservation violated"
        if b is None:
            rng = random.Random(w_seed)
            w = [F(rng.randint(-9, 9)) for _ in range(n)]
            c = dot(w, self.e) / dot(self.e, self.e)
            b = [wi - c * ei for wi, ei in zip(w, self.e)]
        self.b = list(b)
        assert dot(self.b, self.e) == 0

    def C(self, g):
        n = self.n
        C = [[F(0)] * n for _ in range(n)]
        for gk, a in zip(g, self.A):
            if gk == 0:
                continue
            nz = [i for i in range(n) if a[i] != 0]
            for i in nz:
                for j in nz:
                    C[i][j] += gk * a[i] * a[j]
        return C

    def admissible(self, g):
        return rank(self.C(g)) == self.n - 1

    def solve(self, g):
        C = self.C(g)
        n = self.n
        Creg = [[C[i][j] + self.e[i] * self.e[j] for j in range(n)] for i in range(n)]
        x = solve(Creg, self.b)
        assert [dot(C[i], x) for i in range(n)] == self.b, "b not in range C(g)"
        return x, dot(self.b, x), dot(x, x)

    def basic_solutions(self):
        n = self.n
        out = {}
        rhs = [F(0)] * (n - 1) + [F(1)]
        for I in itertools.combinations(range(self.m), n - 2):
            B = [self.A[k] for k in I] + [self.e, self.b]
            try:
                out[I] = solve(B, rhs)
            except ZeroDivisionError:
                continue
        return out

    def V(self, I):
        """V_I with V_I . y = det[a_I, y, e] for all y (columns a_I, y, e)."""
        rows = [self.A[k] for k in I] + [self.e]
        if rank(rows) < self.n - 1:
            return [F(0)] * self.n
        v = nullvec(rows, self.n)
        cols = [self.A[k] for k in I] + [v] + [self.e]
        Mx = [[cols[j][i] for j in range(self.n)] for i in range(self.n)]
        lam = det(Mx) / dot(v, v)
        return [lam * vi for vi in v]


def make_instance(seed, n_target, m_target, eps_pool):
    rng = random.Random(seed)
    eps = list(eps_pool[:n_target])
    ev_all = all_events(eps)
    for _ in range(5000):
        ev = rng.sample(ev_all, min(m_target, len(ev_all)))
        S = [stoich(e, len(eps)) for e in ev]
        if rank(S) == len(eps) - 1 and all(any(s[i] != 0 for s in S) for i in range(len(eps))):
            break
    else:
        raise RuntimeError("no spanning event set found")
    dsc = [F(rng.randint(2, 12), rng.randint(2, 7)) for _ in eps]
    return Instance(eps, ev, dsc, w_seed=seed + 1)


def rand_rates(rng, m, lo=1, hi=9):
    return [F(rng.randint(lo, hi), rng.randint(1, 5)) for _ in range(m)]


# ----------------------------------------------------------------------------- checks
def check_T1_identity(inst, g, Vs):
    x, K, N = inst.solve(g)
    num = [F(0)] * inst.n
    den = F(0)
    n_terms = 0
    for I, VI in Vs.items():
        cI = F(1)
        for k in I:
            cI *= g[k]
        if cI == 0:
            continue
        vb = dot(VI, inst.b)
        if vb == 0:
            continue
        num = [a + cI * vb * v for a, v in zip(num, VI)]
        den += cI * vb * vb
        n_terms += 1
    rhs = [a / den for a in num]
    lhs = [xi / K for xi in x]
    return lhs == rhs, n_terms


def main():
    report = {}
    rng = random.Random(12345)

    # ---------------- instance A: random spanning event set with distinct and repeated daughters
    inst = make_instance(seed=3, n_target=7, m_target=10, eps_pool=(1, 1, 2, 2, 3, 3, 4))
    report["A_instance"] = {"n": inst.n, "m": inst.m, "eps": [str(v) for v in inst.eps],
                            "events": inst.events, "d": [str(v) for v in inst.D],
                            "b": [str(v) for v in inst.b]}
    ys = inst.basic_solutions()
    norms2 = {I: dot(y, y) for I, y in ys.items()}
    Istar = max(norms2, key=lambda I: norms2[I])
    M2 = norms2[Istar]
    report["A_M2"] = {"M2": str(M2), "M2_float": float(M2), "Istar": list(Istar), "n_bases": len(ys)}
    Vs = {I: inst.V(I) for I in itertools.combinations(range(inst.m), inst.n - 2)}
    # consistency: y_I = V_I / (V_I . b) for every basis
    cons = all([a == v / dot(Vs[I], inst.b) for a, v in zip(y, Vs[I])] == [True] * inst.n
               for I, y in ys.items())
    report["A_yI_equals_VI_over_VIb"] = cons
    # V_I . b != 0  <=>  I is a basis
    report["A_basis_iff_VIb_nonzero"] = all((dot(Vs[I], inst.b) != 0) == (I in ys) for I in Vs)

    t1 = []
    for trial in range(5):
        g = rand_rates(rng, inst.m)
        assert inst.admissible(g)
        ok, nt = check_T1_identity(inst, g, Vs)
        x, K, N = inst.solve(g)
        tau = N / K
        t1.append({"identity_exact": ok, "n_terms": nt, "tau_over_M2K": float(tau / (M2 * K)),
                   "bound_holds": tau <= M2 * K})
    # a sparse admissible g (some rates zero)
    for trial in range(3):
        g = rand_rates(rng, inst.m)
        for k in rng.sample(range(inst.m), 2):
            g[k] = F(0)
        if not inst.admissible(g):
            continue
        ok, nt = check_T1_identity(inst, g, Vs)
        x, K, N = inst.solve(g)
        t1.append({"identity_exact": ok, "n_terms": nt, "tau_over_M2K": float((N / K) / (M2 * K)),
                   "bound_holds": N / K <= M2 * K, "zeros": True})
    report["A_T1"] = t1

    # P3: rational-kernel form of |y_I|^2 for every basis (phi_I from ker S_I^T, orthogonal to eps)
    p3_ok = True
    for I, y in ys.items():
        rows = [[F(inst.S[k][i]) for i in range(inst.n)] for k in I] + [list(inst.eps)]
        phi = nullvec(rows, inst.n)                      # s_alpha.phi = 0 (alpha in I) and eps.phi = 0
        Dphi = [inst.D[i] * phi[i] for i in range(inst.n)]
        c = dot(Dphi, inst.e) / dot(inst.e, inst.e)
        PH = [a - c * b for a, b in zip(Dphi, inst.e)]
        val = dot(PH, PH) / dot(inst.b, Dphi) ** 2
        if val != dot(y, y):
            p3_ok = False
    report["A_P3_all_bases_exact"] = p3_ok

    # T2: family g_eps = 1 on I*, eps elsewhere; (|x|/K)^2 = tau/K -> M^2, deficit O(eps), K ~ 1/eps
    fam = []
    for k in range(1, 9):
        epsl = F(1, 10 ** k)
        g = [F(1) if j in Istar else epsl for j in range(inst.m)]
        x, K, N = inst.solve(g)
        ratio = N / K ** 2
        fam.append({"eps": f"1e-{k}", "M2_minus_ratio": float(M2 - ratio),
                    "deficit_over_eps": float((M2 - ratio) / epsl), "K_times_eps": float(K * epsl),
                    "rescale_lambda_for_K0=1": float(K)})
    report["A_T2_family"] = fam

    # P4(a), P4(b), P6 on random rates
    p4 = []
    bb = dot(inst.b, inst.b)
    for trial in range(6):
        g = rand_rates(rng, inst.m)
        x, K, N = inst.solve(g)
        tau = N / K
        Cg = inst.C(g)
        r = [Cg[i][i] for i in range(inst.n)]
        lb_b = F(0)
        for mu in range(inst.n):
            rho2 = max(dot(a, a) / a[mu] ** 2 for a in inst.A if a[mu] != 0)
            lb_b = max(lb_b, inst.b[mu] ** 2 / (rho2 * r[mu] ** 2))
        p4.append({"P4a_tau_ge_K_over_b2": tau >= K / bb, "P4a_ratio": float(tau / (K / bb)),
                   "P4b_N_ge_single_mode": N >= lb_b, "P6": tau / (K / bb) <= M2 * bb})
    report["A_P4_P6"] = p4

    # ---------------- instance B: lifetime-constrained sharpness holds when b is parallel to an event vector
    beta = 0
    instB = Instance(inst.eps, inst.events, inst.D, b=inst.A[beta])
    ysB = instB.basic_solutions()
    n2B = {I: dot(y, y) for I, y in ysB.items()}
    JB = max(n2B, key=lambda I: n2B[I])
    M2B = n2B[JB]
    basis = set(JB) | {beta}
    gB = [rand_rates(rng, 1)[0] if j in basis else F(0) for j in range(instB.m)]
    admB = instB.admissible(gB)
    xB, KB, NB = instB.solve(gB)
    CgB = instB.C(gB)
    report["B_counterexample_sharp_with_lifetimes"] = {
        "b_equals_event_vector": beta, "beta_in_Jstar": beta in JB, "M2": str(M2B), "J_star": list(JB),
        "support_of_g": sorted(basis), "admissible": admB,
        "r_data": [str(CgB[i][i]) for i in range(instB.n)], "K0": str(KB), "tau": str(NB / KB),
        "M2_K0": str(M2B * KB), "tau_equals_M2K0_exactly": NB / KB == M2B * KB}
    eq_all = []
    for trial in range(4):
        g2 = [rand_rates(rng, 1)[0] if j in basis else F(0) for j in range(instB.m)]
        x2, K2, N2 = instB.solve(g2)
        eq_all.append(N2 / K2 == M2B * K2)
    report["B_equality_on_whole_face"] = eq_all
    # with all events active (generic positive rates) the inequality is strict
    g3 = rand_rates(rng, instB.m)
    x3, K3, N3 = instB.solve(g3)
    report["B_full_support_ratio_tau_over_M2K"] = float((N3 / K3) / (M2B * K3))

    # ---------------- instance C: generic b, lifetimes fixed: dense scan of a 2-parameter polygon P(r)
    found = None
    for seed in range(400):
        try:
            I2 = make_instance(seed=1000 + seed, n_target=6, m_target=8, eps_pool=(1, 1, 2, 2, 3, 4))
        except RuntimeError:
            continue
        W = [[I2.A[k][i] ** 2 for k in range(I2.m)] for i in range(I2.n)]
        if rank(W) == I2.n and I2.m - I2.n == 2:
            found = (I2, W)
            break
    if found is not None:
        from scipy.optimize import linprog
        I2, W = found
        ys2 = I2.basic_solutions()
        M22 = max(dot(y, y) for y in ys2.values())
        g0 = rand_rates(rng, I2.m)
        x0, K0, N0 = I2.solve(g0)
        Wf = np.array([[float(v) for v in row] for row in W])
        gf0 = np.array([float(v) for v in g0])
        _, s, Vt = np.linalg.svd(Wf)
        Z = Vt[-2:].T
        Af = np.array([[float(v) for v in a] for a in I2.A])
        ef = np.array([float(v) for v in I2.e])
        bf = np.array([float(v) for v in I2.b])
        outer = np.einsum("mi,mj->mij", Af, Af)
        EE = np.outer(ef, ef)
        box = []
        for i in range(2):
            lohi = []
            for sgn in (1, -1):
                c = np.zeros(2)
                c[i] = sgn
                res = linprog(c, A_ub=-Z, b_ub=gf0, bounds=[(None, None)] * 2, method="highs")
                lohi.append(sgn * res.fun)
            box.append(sorted(lohi))
        Kf0 = float(K0)
        nz = 1501
        z1s = np.linspace(box[0][0], box[0][1], nz)
        z2s = np.linspace(box[1][0], box[1][1], nz)
        Kg = np.full((nz, nz), np.nan)
        Ng = np.full((nz, nz), np.nan)
        for i, z1 in enumerate(z1s):
            G = gf0[None, :] + (Z @ np.vstack([np.full(nz, z1), z2s])).T
            ok = np.min(G, axis=1) > 1e-12 * np.max(gf0)
            if not ok.any():
                continue
            Cb = np.einsum("pm,mij->pij", G[ok], outer) + EE[None]
            X = np.linalg.solve(Cb, np.broadcast_to(bf, (int(ok.sum()), len(bf)))[..., None])[..., 0]
            Kg[i, ok] = X @ bf
            Ng[i, ok] = np.sum(X * X, axis=1)
        taus = []
        for i in range(nz):
            for j in range(nz - 1):
                for (k1, k2, n1, n2) in ((Kg[i, j], Kg[i, j + 1], Ng[i, j], Ng[i, j + 1]),
                                         (Kg[j, i], Kg[j + 1, i], Ng[j, i], Ng[j + 1, i])):
                    if np.isfinite(k1) and np.isfinite(k2) and (k1 - Kf0) * (k2 - Kf0) <= 0 and k1 != k2:
                        t = (Kf0 - k1) / (k2 - k1)
                        taus.append(((1 - t) * n1 + t * n2) / Kf0)
        taus = np.array(taus)
        report["C_generic_lifetime_slice"] = {
            "n": I2.n, "m": I2.m, "dim_P(r)": 2, "M2": float(M22), "K0": Kf0,
            "M2_K0": float(M22) * Kf0, "tau_ref": float(N0 / K0),
            "n_level_set_points": int(len(taus)),
            "tau_max_on_scanned_F": float(taus.max()) if len(taus) else None,
            "tau_min_on_scanned_F": float(taus.min()) if len(taus) else None,
            "ratio_taumax_to_M2K0": float(taus.max() / (float(M22) * Kf0)) if len(taus) else None,
            "K_range_on_P(r)": [float(np.nanmin(Kg)), float(np.nanmax(Kg))],
            "note": "numerical scan of the polygon P(r) (1501^2 grid); illustrative, not a proof"}
    OUT.write_text(json.dumps(report, indent=1, default=str))
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
