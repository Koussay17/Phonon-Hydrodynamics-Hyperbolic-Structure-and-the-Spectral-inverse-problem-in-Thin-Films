"""Attack 5: search for counterexamples to L3,  K(g) >= K_RTA(r)/s_max,  r = diag C(g).

Random energy-conserving event networks (exact rational energies, integer or conserving stoichiometry):
  three-phonon p -> a + b (s = 3), repeated daughter p -> 2a (s = 2), four-phonon p + q -> a + b (s = 4),
  repeated four-phonon 2p -> a + b and p + q -> 2a (s = 3), 'split-merge' with mixed coefficients,
  mixtures; random log-normal rates spanning many decades (near-singular C), random currents
  b orthogonal to e (also adversarial: b chosen to minimise K/K_RTA via the generalised eigenproblem).
Also: the variant with the SELF-ENERGY linewidth r_SE (probe-partner counting: a repeated slot
contributes coefficient |s| instead of s^2) instead of diag C, which is what an RTA code reports.
Each candidate is checked in float and the minimum ratio is re-checked in exact rational arithmetic.
Output: review/cx_attack5_L3.json
"""
import itertools
import json
import sys
from fractions import Fraction

import numpy as np
import scipy.linalg as sla
from hypothesis import given, settings, strategies as st, HealthCheck

sys.dont_write_bytecode = True
OUT = r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review\cx_attack5_L3.json"


def random_network(rng, n, kinds, m_target):
    """Integer energies 1..E; events from the allowed kinds with exact energy conservation."""
    eps = np.sort(rng.integers(1, 3 * n, size=n)).astype(float)
    eps += rng.random(n) * 1e-3 * 0  # keep exact integers
    events = []
    tries = 0
    while len(events) < m_target and tries < 200 * m_target:
        tries += 1
        k = kinds[rng.integers(len(kinds))]
        if k == "3":
            a, b = rng.choice(n, 2, replace=False)
            s = eps[a] + eps[b]
            ps = np.nonzero(eps == s)[0]
            if len(ps):
                p = rng.choice(ps)
                v = np.zeros(n); v[p] -= 1; v[a] += 1; v[b] += 1
                events.append(v)
        elif k == "2a":
            a = rng.integers(n)
            ps = np.nonzero(eps == 2 * eps[a])[0]
            if len(ps):
                p = rng.choice(ps)
                v = np.zeros(n); v[p] -= 1; v[a] += 2
                events.append(v)
        elif k == "4":
            p, q, a = rng.choice(n, 3, replace=False)
            s = eps[p] + eps[q] - eps[a]
            bs = np.nonzero(eps == s)[0]
            bs = bs[(bs != p) & (bs != q) & (bs != a)]
            if len(bs):
                b = rng.choice(bs)
                v = np.zeros(n); v[p] -= 1; v[q] -= 1; v[a] += 1; v[b] += 1
                events.append(v)
        elif k == "4r":  # 2p -> a + b
            p, a = rng.choice(n, 2, replace=False)
            s = 2 * eps[p] - eps[a]
            bs = np.nonzero(eps == s)[0]
            bs = bs[(bs != p) & (bs != a)]
            if len(bs):
                b = rng.choice(bs)
                v = np.zeros(n); v[p] -= 2; v[a] += 1; v[b] += 1
                events.append(v)
        elif k == "mix":  # p -> a + b with conserving non-integer coefficients (tolerance-type)
            a, b = rng.choice(n, 2, replace=False)
            p = rng.integers(n)
            if p in (a, b) or eps[p] <= 0:
                continue
            lp = np.sqrt((eps[a] + eps[b]) / eps[p])
            v = np.zeros(n); v[p] -= lp; v[a] += 1 / lp; v[b] += 1 / lp
            events.append(v)
    if not events:
        return None
    S = np.unique(np.array(events), axis=0)
    return eps, S


def l3_ratio(eps, S, g, b, kT=None, rse=False):
    n = len(eps)
    # entropy scaling with random positive D (any positive diagonal scaling is allowed)
    A = S.T                       # n x m (stoichiometric columns, D = I here; D is absorbed in eps/b)
    C = (A * g) @ A.T
    e = eps / np.linalg.norm(eps)
    w, U = np.linalg.eigh(C)
    keep = np.abs(U.T @ e) < 0.999999
    if np.sum(w[keep] < 1e-12 * w.max()) > 0:
        return None               # not admissible (extra kernel)
    x = sla.lstsq(C, b)[0]
    x -= e * (e @ x)
    if np.linalg.norm(C @ x - b) > 1e-8 * np.linalg.norm(b):
        return None
    K = float(b @ x)
    if rse:
        r = np.abs(A) @ g  # self-energy-type counting: |s| instead of s^2 for repeated slots
        r = np.array([np.sum(g * np.where(np.abs(A[i]) > 0, np.abs(A[i]) * np.minimum(np.abs(A[i]), 1.0), 0)) for i in range(n)])
    else:
        r = (A * A) @ g
    ok = r > 0
    KR = float(np.sum(b[ok] ** 2 / r[ok]))
    smax = int(np.max(np.sum(np.abs(S) > 0, axis=1)[g > 0]))
    return K * smax / KR, K / KR, smax


def adversarial_b(eps, S, g):
    """Minimise K/K_RTA over b orthogonal to e: generalised eigenproblem on H."""
    n = len(eps)
    A = S.T
    C = (A * g) @ A.T
    r = (A * A) @ g
    e = eps / np.linalg.norm(eps)
    Q = sla.null_space(e[None, :])
    CH = Q.T @ C @ Q
    try:
        CHinv = np.linalg.inv(CH)
    except np.linalg.LinAlgError:
        return None
    Rinv = Q.T @ np.diag(1 / r) @ Q
    # minimise y^T CHinv y / y^T Rinv y
    try:
        w, V = sla.eigh(0.5 * (CHinv + CHinv.T), 0.5 * (Rinv + Rinv.T))
    except Exception:
        return None
    return Q @ V[:, 0]


def exact_check(eps, S, g, b):
    """Exact rational evaluation of K*smax/K_RTA for integer S, rational g, b (b projected exactly)."""
    import sympy as spy
    n, m = len(eps), len(S)
    Sx = spy.Matrix([[Fraction(int(round(v))) for v in row] for row in S])
    gq = [spy.Rational(Fraction(float(v)).limit_denominator(10 ** 12)) for v in g]
    A = Sx.T
    C = spy.zeros(n, n)
    for k in range(m):
        col = A[:, k]
        C += gq[k] * (col * col.T)
    eq = spy.Matrix([spy.Rational(int(v)) for v in eps])
    bq = spy.Matrix([spy.Rational(Fraction(float(v)).limit_denominator(10 ** 9)) for v in b])
    bq = bq - eq * (eq.dot(bq) / eq.dot(eq))
    Creg = C + eq * eq.T
    x = Creg.LUsolve(bq)
    K = bq.dot(x)
    r = [sum(gq[k] * A[i, k] ** 2 for k in range(m)) for i in range(n)]
    KR = sum(bq[i] ** 2 / r[i] for i in range(n) if r[i] != 0)
    smax = max(sum(1 for v in S[k] if v != 0) for k in range(m) if g[k] > 0)
    return float(K * smax / KR), str(spy.nsimplify(K * smax / KR)) if False else None


def main():
    rng = np.random.default_rng(20261004)
    kinds_sets = {"3": ["3"], "3+2a": ["3", "2a"], "4": ["4"], "3+4": ["3", "4"], "4r": ["4r", "4"],
                  "all": ["3", "2a", "4", "4r"], "mix": ["mix", "3"]}
    out = {"random": {}, "adversarial": {}, "rse_variant": {}}
    worst_global = (np.inf, None)
    for name, kinds in kinds_sets.items():
        best = np.inf
        best_adv = np.inf
        best_rse = np.inf
        count = 0
        for trial in range(1500):
            n = int(rng.integers(4, 12))
            net = random_network(rng, n, kinds, int(rng.integers(n, 4 * n)))
            if net is None:
                continue
            eps, S = net
            if np.linalg.matrix_rank(S) < n - 1:
                continue  # extra invariants: not admissible with full support
            g = np.exp(rng.normal(0, rng.choice([0.5, 3.0, 8.0]), len(S)))
            e = eps / np.linalg.norm(eps)
            b = rng.normal(size=n); b -= e * (e @ b)
            res = l3_ratio(eps, S, g, b)
            if res is None:
                continue
            count += 1
            best = min(best, res[0])
            if res[0] < worst_global[0]:
                worst_global = (res[0], (name, eps.tolist(), S.tolist(), g.tolist(), b.tolist()))
            ba = adversarial_b(eps, S, g)
            if ba is not None:
                ra = l3_ratio(eps, S, g, ba)
                if ra is not None:
                    best_adv = min(best_adv, ra[0])
                    if ra[0] < worst_global[0]:
                        worst_global = (ra[0], (name, eps.tolist(), S.tolist(), g.tolist(), ba.tolist()))
                    rs = l3_ratio(eps, S, g, ba, rse=True)
                    if rs is not None:
                        best_rse = min(best_rse, rs[0])
        out["random"][name] = {"n_admissible": count, "min_K_smax_over_KRTA": best}
        out["adversarial"][name] = {"min_K_smax_over_KRTA": best_adv}
        out["rse_variant"][name] = {"min_K_smax_over_KRTA_with_self_energy_linewidth": best_rse}
        print(name, count, best, best_adv, best_rse, flush=True)
    out["worst_case_ratio_float"] = worst_global[0]
    # exact re-check of the worst case if integer stoichiometry
    if worst_global[1] is not None:
        name, eps, S, g, b = worst_global[1]
        S = np.array(S)
        if np.allclose(S, np.round(S)):
            out["worst_case_exact_ratio"] = exact_check(np.array(eps), S, np.array(g), np.array(b))[0]
        out["worst_case_family"] = name

    # explicit repeated-daughter chain (energies 1,2,4,8,16) with the self-energy linewidth
    eps = np.array([1, 2, 4, 8, 16], float)
    S = np.array([[2, -1, 0, 0, 0], [0, 2, -1, 0, 0], [0, 0, 2, -1, 0], [0, 0, 0, 2, -1]], float)
    g = np.ones(4)
    ba = adversarial_b(eps, S, g)
    out["chain_diagC"] = l3_ratio(eps, S, g, ba)
    out["chain_selfenergy_linewidth"] = l3_ratio(eps, S, g, ba, rse=True)
    print("chain diagC (ratio*smax, ratio, smax):", out["chain_diagC"], " self-energy:", out["chain_selfenergy_linewidth"])

    # Hypothesis property test (float), three- and four-phonon networks with random rates
    fails = []

    @settings(max_examples=400, deadline=None, suppress_health_check=list(HealthCheck))
    @given(st.integers(0, 2 ** 32 - 1), st.sampled_from(list(kinds_sets.keys())), st.floats(0.1, 10.0))
    def prop(seed, name, spread):
        r = np.random.default_rng(seed)
        n = int(r.integers(4, 10))
        net = random_network(r, n, kinds_sets[name], int(r.integers(n, 3 * n)))
        if net is None:
            return
        eps, S = net
        if np.linalg.matrix_rank(S) < n - 1:
            return
        g = np.exp(r.normal(0, spread, len(S)))
        ba = adversarial_b(eps, S, g)
        if ba is None:
            return
        res = l3_ratio(eps, S, g, ba)
        if res is None:
            return
        if res[0] < 1 - 1e-7:
            fails.append((seed, name, spread, res))
    prop()
    out["hypothesis_failures"] = fails
    print("hypothesis failures:", len(fails))
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
