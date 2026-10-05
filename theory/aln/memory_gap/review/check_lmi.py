"""P5 (13-final-report.md): "{g in P(r): K(g) <= K0} = {g in P(r): [[C(g), b], [b^T, K0]] >= 0}".
With K(g) defined only for ADMISSIBLE g (02-assumptions A: ker C(g) = span(e)), the right-hand side is larger:
it contains every g in P(r) with b in range C(g) and b^T C(g)^+ b <= K0, admissible or not.

Exact example (fractions): two energy-conserving sub-networks A (modes 0-3) and B (modes 4-7), each with
energies (1, 2, 3, 4), joined by "bridge" events.  g0 = intra-network rates only (bridges off): ker C(g0) =
span(e_A, e_B) (two invariants, NOT admissible), but b chosen orthogonal to e_A and e_B lies in range C(g0), so
K(g0) := b^T C(g0)^+ b is finite and [[C, b], [b^T, K(g0)]] = [I; x^T] C [I, x] is PSD exactly.
r := diag C(g0); an admissible g1 in P(r) with bridges on is constructed by compensating the bridge diagonal
with intra rates.  Along g_theta = (1 - theta) g0 + theta g1 (admissible for theta > 0), K(g_theta) -> K(g0).
Run: python -B check_lmi.py  (writes check_lmi.json)
"""
from __future__ import annotations

import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from check_rta_bound import Net, rank, solve, dot, stoich  # noqa: E402


def main():
    rng = random.Random(5)
    eps = [1, 2, 3, 4, 1, 2, 3, 4]
    intra = [(1, 0, 0), (2, 0, 1), (3, 0, 2), (3, 1, 1),
             (5, 4, 4), (6, 4, 5), (7, 4, 6), (7, 5, 5)]
    bridges = [(2, 1, 4), (6, 0, 5), (3, 2, 4), (7, 0, 6)]
    events = intra + bridges
    n = len(eps)
    dsc = [F(rng.randint(2, 15), rng.randint(2, 7)) for _ in eps]
    net = Net(eps, events, dsc)
    SA = [stoich(e, n) for e in intra]
    out = {"rank_intra_S": rank(SA), "rank_all_S": rank([stoich(e, n) for e in events]), "n": n}
    eA = [net.e[i] if i < 4 else F(0) for i in range(n)]
    eB = [net.e[i] if i >= 4 else F(0) for i in range(n)]
    # b orthogonal to e_A and e_B
    w = [F(rng.randint(-9, 9)) for _ in range(n)]
    for ev in (eA, eB):
        c = dot(w, ev) / dot(ev, ev)
        w = [wi - c * vi for wi, vi in zip(w, ev)]
    b = w
    m = len(events)
    g0 = [F(rng.randint(1, 9), rng.randint(1, 4)) for _ in intra] + [F(0)] * len(bridges)
    C0 = net.C(g0)
    out["rank_C(g0)"] = rank(C0)
    out["g0_admissible"] = rank(C0) == n - 1
    # pseudo-inverse solve on the range: (C0 + eA eA^T + eB eB^T) x = b, then C0 x = b
    Creg = [[C0[i][j] + eA[i] * eA[j] + eB[i] * eB[j] for j in range(n)] for i in range(n)]
    x0 = solve(Creg, b)
    assert [dot(C0[i], x0) for i in range(n)] == b
    K0 = dot(b, x0)
    out["K(g0)_finite"] = str(K0)
    # LMI at g0 with t = K0 is a Gram matrix: [[C, b],[b^T, K0]] = [I; x^T] C [I, x]  (exact identity check)
    Mx = [[C0[i][j] for j in range(n)] + [b[i]] for i in range(n)] + [list(b) + [K0]]
    # direct verification: L C L^T with L = [I; x0^T]
    L = [[F(1) if i == k else F(0) for k in range(n)] for i in range(n)] + [list(x0)]
    LCLT = [[sum(L[i][k] * sum(C0[k][l] * L[j][l] for l in range(n)) for k in range(n)) for j in range(n + 1)]
            for i in range(n + 1)]
    out["LMI_at_g0_is_gram_matrix_exactly"] = LCLT == Mx
    # r and an admissible point g1 of P(r)
    W = [[a[i] ** 2 for a in net.A] for i in range(n)]               # n x m
    r = [dot(W[i], g0) for i in range(n)]
    Wi = [[W[i][k] for k in range(len(intra))] for i in range(n)]    # n x 8
    out["rank_W_intra"] = rank(Wi)
    tb = F(1, 50)
    db = [tb] * len(bridges)
    rhs = [-sum(W[i][len(intra) + k] * db[k] for k in range(len(bridges))) for i in range(n)]
    dI = solve(Wi, rhs)
    g1 = [g0[k] + dI[k] for k in range(len(intra))] + db
    assert min(g1) > 0
    assert [dot(W[i], g1) for i in range(n)] == r
    out["g1_in_P(r)"] = True
    out["g1_admissible"] = net.admissible(g1)
    Ks = []
    for th in (F(1, 10), F(1, 100), F(1, 1000), F(1, 10 ** 5)):
        gt = [(1 - th) * a + th * c for a, c in zip(g0, g1)]
        assert net.admissible(gt)
        K, KR, _ = net.response(gt, b)
        Ks.append({"theta": str(th), "K": float(K), "K_minus_K(g0)": float(K - K0)})
    out["K_along_segment_to_g0"] = Ks
    (HERE / "check_lmi.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
