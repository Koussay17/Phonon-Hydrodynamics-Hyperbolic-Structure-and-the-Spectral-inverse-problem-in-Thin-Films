"""Independent checks of Lemma 2 (K0 >= K_RTA / varsigma_max) and of its sharpness for varsigma_max = 2, 3, 4.

1. Matrix inequality C <= varsigma_max diag C, exactly (rational char. polynomial sign test), on random
   energy-conserving networks with distinct daughters, repeated daughters and four-phonon events.
2. The variational bound K0 >= 2 j^T y - y^T C y at y = R^{-1} j / varsigma_max, and K0 >= K_RTA/varsigma_max.
3. Equality condition: K0 = K_RTA/s for some j in H iff there is y with C y = s R y (all active
   supports of size s). Exact equality examples built here (not taken from the sources):
     s = 2: chain of repeated-daughter events;
     s = 3: 3-uniform zero-sum network (p -> a + b);
     s = 4: 4-uniform zero-sum network of 1 -> 3 events (p -> a + b + c).
4. Mixed 3-phonon + 2<->2 networks (varsigma_max = 4): the ratio approaches 1/4 but is never attained
   (lambda_max(R^-1/2 C R^-1/2) < 4 strictly).

Usage: python -B audit_lemma2.py > audit_lemma2.out.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools
import random
import sympy as sp
import numpy as np

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


def is_psd_exact(M):
    """Real symmetric rational M is PSD iff the coefficients of det(lam I - M) alternate weakly in sign."""
    lam = sp.Symbol("lam")
    c = sp.Poly(M.charpoly(lam).as_expr(), lam).all_coeffs()  # c[0] lam^n + c[1] lam^(n-1) + ...
    return all((-1) ** k * c[k] >= 0 for k in range(len(c)))


def build(energies, events, W, g):
    """events: list of integer stoichiometric vectors s (s^T eps = 0); returns C (n x n), e, zetas."""
    nm = len(energies)
    eps = sp.Matrix(energies)
    e = W * eps
    zetas = [W.inv() * sp.Matrix(s) for s in events]
    C = sp.zeros(nm, nm)
    for gg, z in zip(g, zetas):
        C += gg * z * z.T
    return C, e, zetas


def K0_of(C, e, j):
    # K0 = j^T C^+ j for j in H; C^+ via (C + e e^T/|e|^2)^{-1} restricted to H (exact)
    Pe = e * e.T / (e.T * e)[0]
    x = (C + Pe).LUsolve(j)
    x = x - Pe * x
    return (j.T * x)[0], x


def supp(s):
    return sum(1 for c in s if c != 0)


# ---------------------------------------------------------------------------- 1-2: random networks
random.seed(11)
okM, okV, okK, nnet = True, True, True, 0
minratio = None
for trial in range(150):
    nm = random.randint(4, 7)
    energies = sorted(random.sample(range(1, 12), nm))
    cand = []
    for p_, a_, b_ in itertools.product(range(nm), repeat=3):          # p -> a + b, a <= b (incl. a = b)
        if a_ <= b_ and p_ not in (a_, b_) and energies[p_] == energies[a_] + energies[b_]:
            s = [0] * nm; s[p_] -= 1; s[a_] += 1; s[b_] += 1; cand.append(s)
    for p_, q_, a_, b_ in itertools.product(range(nm), repeat=4):      # p + q -> a + b (four-phonon)
        if p_ < q_ and a_ < b_ and len({p_, q_, a_, b_}) == 4 and p_ < a_ and energies[p_] + energies[q_] == energies[a_] + energies[b_]:
            s = [0] * nm; s[p_] -= 1; s[q_] -= 1; s[a_] += 1; s[b_] += 1; cand.append(s)
    if len(cand) < nm - 1:
        continue
    events = random.sample(cand, min(len(cand), random.randint(nm - 1, nm + 4)))
    W = sp.diag(*[sp.Rational(random.randint(2, 9), random.randint(1, 5)) for _ in range(nm)])
    g = [sp.Rational(random.randint(1, 20), random.randint(1, 7)) for _ in events]
    C, e, zetas = build(energies, events, W, g)
    if C.rank() != nm - 1:
        continue
    nnet += 1
    sig = max(supp(s) for s in events)
    R = sp.diag(*[C[k, k] for k in range(nm)])
    okM = okM and is_psd_exact(sig * R - C)
    for _ in range(3):
        jr = sp.Matrix([sp.Rational(random.randint(-9, 9), random.randint(1, 4)) for _ in range(nm)])
        j = jr - e * (e.T * jr)[0] / (e.T * e)[0]
        if j.is_zero_matrix:
            continue
        K0, x = K0_of(C, e, j)
        KR = sum(j[k] ** 2 / R[k, k] for k in range(nm))
        y = R.inv() * j / sig
        okV = okV and K0 >= 2 * (j.T * y)[0] - (y.T * C * y)[0]
        okK = okK and K0 >= KR / sig
        rr = K0 * sig / KR
        minratio = rr if minratio is None else min(minratio, rr)
check(f"C <= varsigma_max diag C exactly on {nnet} admissible random networks (3-phonon, repeated daughters, 2<->2)", okM and nnet >= 10)
check("variational inequality K0 >= 2 j^T y - y^T C y at y = R^-1 j / varsigma_max", okV)
check("K0 >= K_RTA / varsigma_max (exact rationals)", okK, f"min K0*varsigma/K_RTA = {float(minratio):.6f}")


# ---------------------------------------------------------------------------- 3: exact equality examples
def equality_example(energies, events, W, g, sval):
    C, e, zetas = build(energies, events, W, g)
    nm = len(energies)
    admissible = C.rank() == nm - 1
    R = sp.diag(*[C[k, k] for k in range(nm)])
    # psi with s_{gamma,nu} psi_nu = c_gamma on supp: here psi = -sign pattern (built by caller via energies)
    return C, e, R, admissible


def ratio_for_y(C, e, R, y, sval):
    j = sval * R * y
    orth = sp.simplify((e.T * j)[0]) == 0
    K0, _ = K0_of(C, e, j)
    KR = sum(j[k] ** 2 / R[k, k] for k in range(len(j)))
    return sp.nsimplify(K0 / KR), orth, (C * y - sval * R * y).is_zero_matrix


# s = 2: chain eps = 1,2,4,8,16, events 2->1+1, 4->2+2, 8->4+4, 16->8+8; psi = (1,-2,4,-8,16)
en2 = [1, 2, 4, 8, 16]
ev2 = []
for k in range(1, 5):
    s = [0] * 5; s[k] = -1; s[k - 1] = 2; ev2.append(s)
W2 = sp.diag(sp.Rational(3, 2), sp.Rational(5, 7), 2, sp.Rational(1, 3), sp.Rational(9, 4))
g2 = [sp.Rational(2), sp.Rational(1, 3), sp.Rational(7, 5), sp.Rational(4)]
C2, e2, R2, adm2 = equality_example(en2, ev2, W2, g2, 2)
psi2 = sp.Matrix([1, -2, 4, -8, 16])
rat2, orth2, eig2 = ratio_for_y(C2, e2, R2, W2 * psi2, 2)
check("s=2 exact example (repeated-daughter chain, non-Bose weights, unequal rates): admissible, C y = 2 R y, K0/K_RTA = 1/2",
      adm2 and orth2 and eig2 and rat2 == sp.Rational(1, 2), str(rat2))


def zero_sum_network(wv, k):
    """Edges: k-subsets with zero sum of wv (all entries nonzero): the parent is the odd-sign vertex."""
    nm = len(wv)
    edges = []
    for S in itertools.combinations(range(nm), k):
        if sum(wv[i] for i in S) == 0:
            sg = [1 if wv[i] > 0 else -1 for i in S]
            if sg.count(1) == 1 or sg.count(-1) == 1:          # exactly one vertex of the minority sign
                par = [S[t] for t in range(k) if sg.count(sg[t]) == 1][0]
                s = [0] * nm
                s[par] = -1
                for i in S:
                    if i != par:
                        s[i] += 1
                edges.append(s)
    return edges


def select_rank(edges, target):
    chosen, M = [], sp.zeros(0, len(edges[0]))
    for s in edges:
        M2 = M.col_join(sp.Matrix([s]))
        if M2.rank() > M.rank():
            M, chosen = M2, chosen + [s]
        if M.rank() == target:
            break
    return chosen, M.rank()


def find_network(kk, nm, seed, tries=4000):
    """Search for an admissible equality network: parents (w < 0) and daughters (w > 0) with energies in a
    small range, all resonant (1 parent, kk-1 daughters) events enumerated; accept if the stoichiometric
    rank is nm-1 (energy is the only invariant)."""
    rnd = random.Random(seed)
    for _ in range(tries):
        npar = rnd.randint(2, nm // 2)
        nd = nm - npar
        dmax = 3 if kk == 3 else 2
        dE = [rnd.randint(1, dmax) for _ in range(nd)]
        pE = [rnd.randint(kk - 1, dmax * (kk - 1)) for _ in range(npar)]
        wv = [-x for x in pE] + dE
        edges = zero_sum_network(wv, kk)
        if len(edges) < nm - 1:
            continue
        if sp.Matrix(edges).rank() == nm - 1:
            chosen, rk = select_rank(edges, nm - 1)
            return wv, chosen, rk
    return None, None, None


for kk, nm0 in ((3, 8), (4, 9)):
    wv, chosen, rk = find_network(kk, nm0, seed=100 + kk)
    if wv is None and kk == 4:
        # deterministic fallback: two parents of energy 3, six daughters of energy 1, all 1 -> 3 events
        wv = [-3, -3, 1, 1, 1, 1, 1, 1]
        chosen, rk = select_rank(zero_sum_network(wv, 4), len(wv) - 1)
    if wv is None:
        check(f"s={kk}: zero-sum network found", False)
        continue
    nm = len(wv)
    energies = [abs(x) for x in wv]
    # energy conservation of every chosen event
    cons = all(sum(s[i] * energies[i] for i in range(nm)) == 0 for s in chosen)
    random.seed(kk)
    W = sp.diag(*[sp.Rational(random.randint(2, 9), random.randint(1, 5)) for _ in range(nm)])
    g = [sp.Rational(random.randint(1, 9), random.randint(1, 4)) for _ in chosen]
    C, e, R, adm = equality_example(energies, chosen, W, g, kk)
    psi = sp.Matrix([-1 if x > 0 else 1 for x in wv])  # psi = -sign(w)
    # sign convention check: s_{gamma nu} psi_nu equal on the support of every event
    signok = all(len({s[i] * psi[i] for i in range(nm) if s[i] != 0}) == 1 for s in chosen)
    rat, orth, eig = ratio_for_y(C, e, R, W * psi, kk)
    lab = "p -> a + b" if kk == 3 else "p -> a + b + c (1 <-> 3 four-phonon events)"
    check(f"s={kk} exact example ({lab}): {nm} modes, {len(chosen)} events, energy-conserving, admissible, "
          f"C y = {kk} R y, K0/K_RTA = 1/{kk}",
          cons and rk == nm - 1 and adm and signok and orth and eig and rat == sp.Rational(1, kk), str(rat))
    print("      energies:", energies)
    print("      events (stoichiometry):", chosen)

# ---------------------------------------------------------------------------- 4: mixed 3-phonon + 2<->2 networks
# 2<->2 sub-network with psi in {+-1} attaining equality, plus 3-phonon events of rate delta that
# break phonon-number conservation; varsigma_max = 4.
energies = [1, 2, 3, 4, 5, 6]
psi_m = [1, 1, -1, -1, 1, -1]
c22 = []
for p_, q_, a_, b_ in itertools.product(range(6), repeat=4):
    if p_ < q_ and a_ < b_ and len({p_, q_, a_, b_}) == 4 and energies[p_] + energies[q_] == energies[a_] + energies[b_]:
        s = [0] * 6; s[p_] -= 1; s[q_] -= 1; s[a_] += 1; s[b_] += 1
        if len({s[i] * psi_m[i] for i in range(6) if s[i] != 0}) == 1 and (tuple(-x for x in s) not in [tuple(t) for t in c22]):
            c22.append(s)
c3 = []
for p_, a_, b_ in itertools.product(range(6), repeat=3):
    if a_ < b_ and p_ not in (a_, b_) and energies[p_] == energies[a_] + energies[b_]:
        s = [0] * 6; s[p_] -= 1; s[a_] += 1; s[b_] += 1; c3.append(s)
Wm = sp.eye(6)
ok_mixed = True
vals = []
for delta in (sp.Rational(1, 10), sp.Rational(1, 100), sp.Rational(1, 1000), sp.Rational(1, 10 ** 5)):
    ev = c22 + c3
    gg = [1] * len(c22) + [delta] * len(c3)
    C, e, zetas = build(energies, ev, Wm, gg)
    if C.rank() != 5:
        ok_mixed = False
        break
    R = sp.diag(*[C[k, k] for k in range(6)])
    Rm = np.diag([1 / float(R[k, k]) ** 0.5 for k in range(6)])
    Cn = np.array(C.evalf(30).tolist(), dtype=float)
    lmax = float(np.max(np.linalg.eigvalsh(Rm @ Cn @ Rm)))
    vals.append((float(delta), lmax, 4 / lmax))
    ok_mixed = ok_mixed and lmax < 4
check("mixed 3-phonon + 2<->2 network (varsigma_max = 4): lambda_max(R^-1/2 C R^-1/2) < 4 for delta > 0, -> 4 as delta -> 0",
      ok_mixed and vals[-1][1] > 3.999 and len(c22) > 0,
      "; ".join(f"delta={d:g}: min K0*4/K_RTA = {r_:.6f}" for d, lm, r_ in vals))
# pure 2<->2 is never admissible (phonon number is an extra invariant)
C22, e22, _ = build(energies, c22, Wm, [1] * len(c22))
Nvec = Wm.inv() * sp.ones(6, 1)  # phonon number in entropy coordinates: W^{-1}... check s^T 1 = 0 for every 2<->2 event
check("pure 2<->2 network is not admissible (s^T 1 = 0: phonon number conserved, rank < n-1)",
      all(sum(s) == 0 for s in c22) and C22.rank() < 5)

print()
print(f"{sum(res)} of {len(res)} checks passed")
