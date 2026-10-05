"""Exact checks of L3 (K >= K_RTA / s for events touching at most s modes) - 13-final-report.md, L3.

Questions audited:
  (1) does the bound hold with repeated-daughter events (stoichiometry -e_p + 2 e_a, entries not +-1)?
  (2) what is the sharp constant?  Claim checked here: the sharp constant is s_max, the largest number of
      DISTINCT modes touched by an active event (3 for three-phonon events, 2 if every active event has a
      repeated daughter), independently of the stoichiometric coefficients.
  (3) can K / K_RTA = 1/3 be attained exactly by an admissible, energy-conserving three-phonon network?

Exact arithmetic (fractions.Fraction). Entropy scales d_mu are arbitrary positive rationals: L3 is algebraic.
Equality construction (derived in the audit): if psi is a vector with s_{alpha,mu} psi_mu = c_alpha for every
mode mu in the support of every active event alpha, then with y = D psi one has a_alpha.y = |supp alpha| c_alpha,
C y = s R y when all supports have size s, and b := s R y gives K = b.y = s y^T R y, K_RTA = s^2 y^T R y,
so K = K_RTA / s exactly.  For s = 3 this needs: in each event the two daughters carry psi = c and the parent
psi = -c (a "bipartite" sign structure).  For repeated daughters (-1, 2): psi_a = c/2, psi_p = -c.

Run: python -B check_rta_bound.py   (writes check_rta_bound.json)
"""
from __future__ import annotations

import itertools
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
OUT = Path(__file__).with_suffix(".json")


# ---------------------------------------------------------------- exact linear algebra
def rank(rows):
    A = [[F(v) for v in r] for r in rows]          # exact: never fall back to float division
    nr, nc = len(A), len(A[0]) if A else 0
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


def solve(A, b):
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [a - f * bb for a, bb in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


# ---------------------------------------------------------------- model
def stoich(ev, n):
    p, a, b = ev
    s = [0] * n
    s[p] -= 1
    s[a] += 1
    s[b] += 1
    return s


class Net:
    def __init__(self, eps, events, dsc):
        self.n = len(eps)
        self.eps = [F(x) for x in eps]
        self.events = events
        self.D = [F(x) for x in dsc]
        self.A = [[F(s) / self.D[i] for i, s in enumerate(stoich(ev, self.n))] for ev in events]
        self.e = [self.D[i] * self.eps[i] for i in range(self.n)]
        assert all(dot(a, self.e) == 0 for a in self.A), "energy conservation"

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

    def response(self, g, b):
        C = self.C(g)
        n = self.n
        Creg = [[C[i][j] + self.e[i] * self.e[j] for j in range(n)] for i in range(n)]
        x = solve(Creg, b)
        Cx = [dot(C[i], x) for i in range(n)]
        assert Cx == list(b), "b not in range of C (non-admissible?)"
        K = dot(b, x)
        r = [C[i][i] for i in range(n)]
        KR = sum(bi * bi / ri for bi, ri in zip(b, r))
        return K, KR, r

    def project_H(self, w):
        ee = dot(self.e, self.e)
        c = dot(w, self.e) / ee
        return [wi - c * ei for wi, ei in zip(w, self.e)]


def all_events(eps, allow_repeated=True):
    n = len(eps)
    out = []
    for p in range(n):
        for a in range(n):
            for b in range(a, n):
                if p in (a, b) or (a == b and not allow_repeated):
                    continue
                if eps[p] == eps[a] + eps[b]:
                    out.append((p, a, b))
    return out


def main():
    rng = random.Random(2026)
    rep = {}

    # (1) random tests: bound with constant 3, events with and without repeated daughters
    viol3, viol2_wrong, ratios = 0, 0, []
    eps = [1, 1, 2, 2, 2, 3, 3, 4, 4, 5]
    ev_all = all_events(eps)
    n_tests = 0
    while n_tests < 300:
        ev = rng.sample(ev_all, 14)
        if rank([stoich(e, len(eps)) for e in ev]) != len(eps) - 1:
            continue
        dsc = [F(rng.randint(1, 20), rng.randint(1, 9)) for _ in eps]
        net = Net(eps, ev, dsc)
        g = [F(rng.randint(1, 30), rng.randint(1, 7)) for _ in ev]
        if not net.admissible(g):
            continue
        w = [F(rng.randint(-20, 20)) for _ in eps]
        b = net.project_H(w)
        K, KR, r = net.response(g, b)
        smax = max(sum(1 for v in a if v != 0) for a in net.A)
        ratios.append(float(K / KR))
        if K < KR / smax:
            viol3 += 1
        n_tests += 1
    rep["random_tests"] = {"n": n_tests, "violations_of_K_ge_KRTA_over_3": viol3,
                           "min_K_over_KRTA": min(ratios), "max_K_over_KRTA": max(ratios),
                           "events_include_repeated_daughters": True}

    # (2) exact equality K = K_RTA/3 with support-3 events (bipartite sign structure).
    # Reduction used: with phi = sigma * phi', the invariants of a bipartite network are the solutions of the
    # zero-sum system phi'_p + phi'_a + phi'_b = 0 on its triples.  Conversely, any 3-uniform hypergraph whose
    # zero-sum space is spanned by a nowhere-zero integer w gives an admissible bipartite network with
    # eps = |w|, sigma = sign(w), parent = the vertex of the odd sign (|w_p| = |w_a| + |w_b|).
    # (A naive random search over (eps, sigma) found none in 2e4 trials; the zero-sum route finds one at once.)
    found = None
    rngz = random.Random(7)
    for trial in range(20000):
        n = rngz.randint(6, 12)
        vals = [v for v in range(-6, 7) if v != 0]
        w = [rngz.choice(vals) for _ in range(n)]
        trip = [t for t in itertools.combinations(range(n), 3) if w[t[0]] + w[t[1]] + w[t[2]] == 0]
        if len(trip) < n - 1:
            continue
        ev = []
        for t in trip:
            sg = [1 if w[i] > 0 else -1 for i in t]
            odd = [i for i, s_ in zip(t, sg) if sg.count(s_) == 1][0]
            da = sorted(i for i in t if i != odd)
            ev.append((odd, da[0], da[1]))
        eps = [abs(v) for v in w]
        S = [stoich(e, n) for e in ev]
        if any(all(s[i] == 0 for s in S) for i in range(n)):
            continue
        if rank(S) != n - 1:
            continue
        sig = [1 if v > 0 else -1 for v in w]
        found = (eps, sig, ev)
        break
    assert found is not None
    eps, sig, ev = found
    assert all(sig[e[1]] == sig[e[2]] == -sig[e[0]] for e in ev)
    dsc = [F(rng.randint(1, 20), rng.randint(1, 9)) for _ in eps]
    net = Net(eps, ev, dsc)
    g = [F(rng.randint(1, 30), rng.randint(1, 7)) for _ in ev]
    assert net.admissible(g)
    C = net.C(g)
    R = [C[i][i] for i in range(net.n)]
    y = [net.D[i] * sig[i] for i in range(net.n)]
    b = [3 * R[i] * y[i] for i in range(net.n)]
    K, KR, r = net.response(g, b)
    rep["sharp_support3"] = {"eps": eps, "sign_class": sig, "events": ev, "d": [str(v) for v in dsc],
                             "g": [str(v) for v in g], "b_dot_e": str(dot(b, net.e)),
                             "admissible": True, "K": str(K), "K_RTA": str(KR), "K_over_KRTA": str(K / KR)}

    # same network, random rates: still exactly 1/3 for b = 3 R(g) D sigma (b depends on g)
    eq = []
    for t in range(5):
        g = [F(rng.randint(1, 30), rng.randint(1, 7)) for _ in ev]
        C = net.C(g)
        R = [C[i][i] for i in range(net.n)]
        b = [3 * R[i] * y[i] for i in range(net.n)]
        K, KR, r = net.response(g, b)
        eq.append(str(K / KR))
    rep["sharp_support3_other_rates"] = eq

    # (3) repeated-daughter-only chain: constant 2 attained (incidence (-1, 2))
    epsc = [1, 2, 4, 8, 16]
    evc = [(1, 0, 0), (2, 1, 1), (3, 2, 2), (4, 3, 3)]
    dscc = [F(rng.randint(1, 20), rng.randint(1, 9)) for _ in epsc]
    netc = Net(epsc, evc, dscc)
    gc = [F(rng.randint(1, 30), rng.randint(1, 7)) for _ in evc]
    assert netc.admissible(gc)
    psi = [F(1), F(-2), F(4), F(-8), F(16)]
    Cc = netc.C(gc)
    Rc = [Cc[i][i] for i in range(netc.n)]
    yc = [netc.D[i] * psi[i] for i in range(netc.n)]
    bc = [2 * Rc[i] * yc[i] for i in range(netc.n)]
    Kc, KRc, _ = netc.response(gc, bc)
    rep["repeated_only_chain"] = {"K_over_KRTA": str(Kc / KRc),
                                  "note": "support-2 events: sharp constant 2; the coefficient 2 is irrelevant"}

    # (4) the would-be 'entry-bounded' version fails: z z^T <= 3 diag(z^2) is a support statement.
    #     For z = (-1, 2): 2 diag(z^2) - z z^T = [[1, 2], [2, 4]] (PSD, singular); 3 diag(z^2) - z z^T PSD.
    z = (F(-1), F(2))
    M2 = [[2 * z[0] ** 2 - z[0] * z[0], -z[0] * z[1]], [-z[1] * z[0], 2 * z[1] ** 2 - z[1] * z[1]]]
    rep["repeated_daughter_matrix_inequality"] = {
        "2diag(z^2)-zz^T": [[str(v) for v in row] for row in M2],
        "det": str(M2[0][0] * M2[1][1] - M2[0][1] * M2[1][0]), "trace": str(M2[0][0] + M2[1][1])}

    OUT.write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
