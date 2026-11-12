"""Independent checks of Supplemental Sec. S2 (sharpness of K0 >= K_RTA / varsigma_max, Lemma 2).

Routes deliberately differ from checks/verify_manuscript_algebra.py and checks/audit/audit_lemma2.py:
  * rank n-1 is checked through the explicit null space of S^T (it must be span{eps});
  * the equality identity is checked event by event, zeta (zeta^T y) = s diag(zeta o zeta) y, with
    symbolic weights (rates factor out of each event term);
  * K0 is computed from the bordered system [[C, e], [e^T, 0]] (exact rationals) and, with Bose weights
    at several temperatures, in 50-digit arithmetic;
  * the generalized eigenvalue varsigma of (C, R) is certified by det(varsigma R - C) = 0 together with
    positive semidefiniteness of varsigma R - C (exact sign test of the characteristic polynomial);
  * for the mixed family: weight independence of y^T C y and y^T R y for y = W psi (symbolic weights),
    strictness by det(4R - C) > 0 for delta > 0, certified enclosures of 4/lambda_max(G) from exact real-root
    isolation of det(lambda R - C), and the slope d(4/lambda_max)/d delta at delta = 0 by implicit
    differentiation;
  * min_j K0/K_RTA = 1/lambda_max(G) against the explicit minimiser and random currents;
  * the equality criterion (existence of y with zeta_{g,nu} y_nu constant on every support, all supports of
    size s) against det(sR - C) = 0 on random uniform networks;
  * the self-energy statement that the main text attributes to Sec. S2 (pending there): chain network,
    self-coupling of the repeated daughter counted with half weight.

Usage (from this directory): python -B audit_s2_networks.py > audit_s2_networks.out.txt
Single thread; < 200 MB.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import itertools
import random
import sympy as sp
import mpmath as mp

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""), flush=True)


def stoich(n, parent, daughters):
    s = [0] * n
    s[parent - 1] -= 1
    for d in daughters:
        s[d - 1] += 1
    return s


def psd_exact(M):
    """Rational symmetric M is PSD iff the coefficients of det(x I - M) alternate weakly in sign."""
    x = sp.Symbol("x")
    c = sp.Poly(M.charpoly(x).as_expr(), x).all_coeffs()
    return all((-1) ** k * c[k] >= 0 for k in range(len(c)))


def K0_bordered(C, e, j):
    """K0 = j^T x with C x = j, e^T x = 0, from the nonsingular bordered system (exact)."""
    n = C.shape[0]
    Bd = sp.zeros(n + 1, n + 1)
    Bd[:n, :n] = C
    Bd[:n, n] = e
    Bd[n, :n] = e.T
    rhs = sp.Matrix(list(j) + [0])
    sol = Bd.LUsolve(rhs)
    x = sol[:n, 0]
    return (j.T * x)[0], x, sol[n]


# ================================================================ Table S1 networks
nets = {
    2: ([1, 2, 4, 8, 16], [(k + 1, [k, k]) for k in range(1, 5)], [1, -2, 4, -8, 16]),
    3: ([6, 6, 6, 3, 3, 3, 3, 3],
        [(1, [4, 5]), (1, [4, 6]), (1, [4, 7]), (1, [4, 8]), (1, [5, 6]), (2, [4, 5]), (3, [4, 5])],
        [1, 1, 1, -1, -1, -1, -1, -1]),
    4: ([6, 6, 6, 2, 2, 2, 2, 2, 2],
        [(1, [4, 5, 6]), (1, [4, 5, 7]), (1, [4, 5, 8]), (1, [4, 5, 9]), (1, [4, 6, 7]), (1, [5, 6, 7]),
         (2, [4, 5, 6]), (3, [4, 5, 6])],
        [1, 1, 1, -1, -1, -1, -1, -1, -1]),
}
print("=== Table S1 networks")
for vs, (eps, evdef, psi) in nets.items():
    n = len(eps)
    S = [stoich(n, p_, ds) for p_, ds in evdef]
    St = sp.Matrix(S)  # rows = events
    ns = St.nullspace()  # kernel of the map x -> S x, i.e. vectors with s_g^T x = 0 for all g
    ker_is_eps = len(ns) == 1 and sp.Matrix.hstack(ns[0], sp.Matrix(eps)).rank() == 1
    sizes = {sum(1 for c in s if c != 0) for s in S}
    const = all(len({s[i] * psi[i] for i in range(n) if s[i] != 0}) == 1 for s in S)
    # event-wise identity with symbolic weights: zeta (zeta^T y) == vs * diag(zeta o zeta) y, y = W psi
    ws = sp.symbols(f"w1:{n + 1}", positive=True)
    y = sp.Matrix([ws[i] * psi[i] for i in range(n)])
    evwise = True
    for s in S:
        z = sp.Matrix([sp.Integer(s[i]) / ws[i] for i in range(n)])
        lhs = z * (z.T * y)[0]
        rhs = vs * sp.Matrix([z[i] ** 2 * y[i] for i in range(n)])
        evwise = evwise and sp.simplify(lhs - rhs).is_zero_matrix
    check(f"S1 net varsigma={vs}: n={n}, {len(S)} events all of size {vs}, energy-conserving, ker S^T = span(eps) "
          f"(rank n-1), s*psi constant on supports, event-wise identity for symbolic weights",
          sizes == {vs} and all(sum(s[i] * eps[i] for i in range(n)) == 0 for s in S) and ker_is_eps and const and evwise)
    # exact rational instances: bordered-system K0, generalized eigenvalue certified
    rnd = random.Random(1000 + vs)
    okr = True
    for trial in range(4):
        w = [sp.Rational(rnd.randint(1, 30), rnd.randint(1, 30)) for _ in range(n)]
        g = [sp.Rational(rnd.randint(1, 50), rnd.randint(1, 13)) for _ in S]
        C = sp.zeros(n, n)
        for s, gg in zip(S, g):
            z = sp.Matrix([sp.Rational(s[i]) / w[i] for i in range(n)])
            C += gg * z * z.T
        e = sp.Matrix([w[i] * eps[i] for i in range(n)])
        R = sp.diag(*[C[i, i] for i in range(n)])
        yv = sp.Matrix([w[i] * psi[i] for i in range(n)])
        j = C * yv
        K0, x, lag = K0_bordered(C, e, j)
        KR = sum(j[i] ** 2 / R[i, i] for i in range(n))
        detz = (vs * R - C).det() == 0
        psd = psd_exact(vs * R - C)
        okr = okr and (e.T * j)[0] == 0 and lag == 0 and K0 / KR == sp.Rational(1, vs) and detz and psd \
            and C.rank() == n - 1
    check(f"S1 net varsigma={vs}: 4 random rational weight/rate sets: K0/K_RTA = 1/{vs} (bordered system), "
          f"det({vs}R - C) = 0 and {vs}R - C PSD (so lambda_max(G) = {vs} exactly)", okr)
    # Bose weights at three temperatures, 50 digits
    mp.mp.dps = 50
    okb, worst = True, mp.mpf(0)
    for T in (mp.mpf("0.7"), mp.mpf(3), mp.mpf(20)):
        f = [1 / (mp.e ** (mp.mpf(ei) / T) - 1) for ei in eps]
        w = [mp.sqrt(fi * (1 + fi)) for fi in f]
        g = [mp.mpf(k + 1) ** 2 / 7 for k in range(len(S))]
        C = mp.matrix(n, n)
        for s, gg in zip(S, g):
            for a in range(n):
                for b in range(n):
                    C[a, b] += gg * s[a] * s[b] / (w[a] * w[b])
        e = mp.matrix([w[i] * eps[i] for i in range(n)])
        yv = mp.matrix([w[i] * psi[i] for i in range(n)])
        j = C * yv
        Bd = mp.matrix(n + 1, n + 1)
        for a in range(n):
            for b in range(n):
                Bd[a, b] = C[a, b]
            Bd[a, n] = e[a]
            Bd[n, a] = e[a]
        sol = mp.lu_solve(Bd, mp.matrix(list(j) + [0]))
        K0 = sum(j[a] * sol[a] for a in range(n))
        KR = sum(j[a] ** 2 / C[a, a] for a in range(n))
        dev = abs(K0 / KR - mp.mpf(1) / vs)
        worst = max(worst, dev)
    check(f"S1 net varsigma={vs}: Bose weights at T = 0.7, 3, 20 (energy units), unequal rates: |K0/K_RTA - 1/{vs}| < 1e-40",
          worst < mp.mpf("1e-40"), f"max deviation {mp.nstr(worst, 3)}")

# ================================================================ equality criterion on random uniform networks
print("=== equality criterion vs det(sR - C) on random uniform networks")


def has_equal_product_vector(S, w):
    """Exists y != 0 with zeta_{g,nu} y_nu equal on supp of every event (zeta = W^-1 s)? Linear system."""
    n = len(w)
    rows = []
    for s in S:
        sup = [i for i in range(n) if s[i] != 0]
        for a, b in zip(sup[:-1], sup[1:]):
            r = [0] * n
            r[a] = sp.Rational(s[a]) / w[a]
            r[b] = -sp.Rational(s[b]) / w[b]
            rows.append(r)
    A = sp.Matrix(rows)
    return A.rank() < n


rnd = random.Random(7)
agree, tested, n_pos = True, 0, 0
for trial in range(400):
    n = rnd.randint(4, 7)
    kk = rnd.choice([2, 3])
    eps = [rnd.randint(1, 4) for _ in range(n)]
    cand = []
    if kk == 3:
        for p_, a_, b_ in itertools.product(range(n), repeat=3):
            if a_ < b_ and p_ not in (a_, b_) and eps[p_] == eps[a_] + eps[b_]:
                s = [0] * n; s[p_] = -1; s[a_] += 1; s[b_] += 1; cand.append(s)
    else:
        for p_, a_ in itertools.product(range(n), repeat=2):
            if p_ != a_ and eps[p_] == 2 * eps[a_]:
                s = [0] * n; s[p_] = -1; s[a_] = 2; cand.append(s)
    if len(cand) < n - 1:
        continue
    S = rnd.sample(cand, min(len(cand), rnd.randint(n - 1, n + 3)))
    if sp.Matrix(S).rank() != n - 1:
        continue
    w = [sp.Rational(rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(n)]
    g = [sp.Rational(rnd.randint(1, 9), rnd.randint(1, 5)) for _ in S]
    C = sp.zeros(n, n)
    for s, gg in zip(S, g):
        z = sp.Matrix([sp.Rational(s[i]) / w[i] for i in range(n)])
        C += gg * z * z.T
    R = sp.diag(*[C[i, i] for i in range(n)])
    crit = has_equal_product_vector(S, w)
    detz = (kk * R - C).det() == 0
    tested += 1
    n_pos += int(crit)
    agree = agree and (crit == detz)
check(f"equality criterion: 'exists y with constant zeta*y on all supports' <=> det(sR - C) = 0, on {tested} random "
      f"admissible uniform networks ({n_pos} with equality)", agree and tested >= 40 and n_pos >= 3)

# ================================================================ mixed family
print("=== mixed family (energies 1..6; 2+5->3+4 rate 1; six three-phonon events rate delta)")
eps = [1, 2, 3, 4, 5, 6]
S22 = [[0, -1, 1, 1, -1, 0]]
S3 = [stoich(6, 3, [1, 2]), stoich(6, 4, [1, 3]), stoich(6, 5, [1, 4]), stoich(6, 5, [2, 3]),
      stoich(6, 6, [1, 5]), stoich(6, 6, [2, 4])]
Smx = S22 + S3
dl = sp.Symbol("delta", positive=True)
ns = sp.Matrix(Smx).nullspace()
ns3 = sp.Matrix(S3).nullspace()
check("mixed: all events conserve energy; ker S^T = span(eps) (rank 5); the three-phonon events alone already have rank 5",
      all(sum(s[i] * eps[i] for i in range(6)) == 0 for s in Smx)
      and len(ns) == 1 and sp.Matrix.hstack(ns[0], sp.Matrix(eps)).rank() == 1
      and len(ns3) == 1 and sp.Matrix.hstack(ns3[0], sp.Matrix(eps)).rank() == 1)
# symbolic weights: y = W psi gives weight-independent quadratic forms
ws = sp.symbols("w1:7", positive=True)
Cw = sp.zeros(6, 6)
for s, gg in zip(Smx, [1] + [dl] * 6):
    z = sp.Matrix([sp.Integer(s[i]) / ws[i] for i in range(6)])
    Cw += gg * z * z.T
psi = sp.Matrix([0, 1, -1, -1, 1, 0])
yw = sp.Matrix([ws[i] * psi[i] for i in range(6)])
Rw = sp.diag(*[Cw[i, i] for i in range(6)])
qC = sp.simplify((yw.T * Cw * yw)[0])
qR = sp.simplify((yw.T * Rw * yw)[0])
check("mixed: for y = W psi, psi = (0,1,-1,-1,1,0), y^T C y = 16 + 10 delta and y^T R y = 4 + 12 delta for ALL weights "
      "(the bound (S2) is weight independent)", sp.expand(qC - 16 - 10 * dl) == 0 and sp.expand(qR - 4 - 12 * dl) == 0)
# unit weights: strictness and certified minima
C = sp.zeros(6, 6)
for s, gg in zip(Smx, [1] + [dl] * 6):
    C += gg * sp.Matrix(s) * sp.Matrix(s).T
R = sp.diag(*[C[i, i] for i in range(6)])
print("    diag(C) =", [sp.factor(R[i, i]) for i in range(6)])
D4 = sp.factor((4 * R - C).det())
print("    det(4R - C) =", D4)
D4p = sp.Poly(sp.expand(D4), dl)
minors_ok = True
for kmin in range(1, 7):
    mnr = sp.Poly(sp.expand((4 * R - C)[:kmin, :kmin].det()), dl)
    # positive for every delta > 0: all coefficients >= 0 and not identically zero is sufficient
    minors_ok = minors_ok and all(c >= 0 for c in mnr.all_coeffs()) and not mnr.is_zero
check("mixed: 4R - C positive definite for every delta > 0 (all leading principal minors are polynomials in delta "
      "with nonnegative coefficients), hence lambda_max(G) < 4 and min 4K0/K_RTA > 1 strictly", minors_ok)
lam = sp.Symbol("lam")
P = sp.expand((lam * R - C).det())
Pf = sp.factor(P)
print("    det(lam R - C) =", Pf)
stated = {sp.Rational(1, 10): "1.2198", sp.Rational(1, 100): "1.0234", sp.Rational(1, 1000): "1.00236",
          sp.Rational(1, 10 ** 5): "1.000024"}
ok_enc, prev = True, None
for dv, txt in stated.items():
    Pd = sp.Poly(P.subs(dl, dv), lam)
    ivs = Pd.intervals(eps=sp.Rational(1, 10 ** 30))
    lo, hi = max(iv[0] for iv in ivs)  # largest real root enclosure
    m_lo, m_hi = 4 / hi, 4 / lo
    bound = (16 + 48 * dv) / (16 + 10 * dv)
    digits = len(txt.split(".")[1])
    rounded_ok = abs(sp.N(m_lo, 40) - sp.Float(txt, 40)) <= sp.Rational(1, 2) * sp.Rational(1, 10 ** digits)
    ok_enc = ok_enc and 1 < m_lo and m_hi <= bound and rounded_ok and (prev is None or m_hi < prev)
    prev = m_lo
    print(f"    delta = {float(dv):g}: 4/lambda_max in [{sp.N(m_lo, 22)}, {sp.N(m_hi, 22)}]; bound {sp.N(bound, 12)}; stated {txt}")
check("mixed: certified enclosures of min 4K0/K_RTA = 4/lambda_max(G) (exact root isolation of det(lam R - C)) "
      "agree with the stated 1.2198, 1.0234, 1.00236, 1.000024, lie in (1, (16+48d)/(16+10d)] and decrease", ok_enc)
# slope at delta = 0: remove the delta^k factor and differentiate implicitly at (lam, delta) = (4, 0)
Q = sp.Poly(P, dl)
kpow = min(m[0] for m in Q.monoms())
Qr = sp.expand(sp.cancel(P / dl ** kpow))
Q40 = Qr.subs({lam: 4, dl: 0})
Ql = sp.diff(Qr, lam).subs({lam: 4, dl: 0})
Qd = sp.diff(Qr, dl).subs({lam: 4, dl: 0})
slope_lam = -Qd / Ql if Ql != 0 else None
c_slope = -slope_lam / 4 if slope_lam is not None else None
print(f"    det(lam R - C) = delta^{kpow} Q(lam, delta); Q(4,0) = {Q40}, Q_lam(4,0) = {Ql}; "
      f"d lam_max/d delta (0) = {slope_lam}; 4/lam_max = 1 + {c_slope} delta + O(delta^2)")
check("mixed: lambda_max(0+) = 4 is a simple root of the reduced determinant; 4/lambda_max = 1 + c delta + O(delta^2) with "
      "c below the slope 38/16 of the bound (S2)",
      Q40 == 0 and Ql != 0 and c_slope is not None and 0 < c_slope < sp.Rational(38, 16), f"c = {c_slope} = {sp.N(c_slope, 12)}")

# ================================================================ min_j K0/K_RTA = 1/lambda_max(G)
print("=== min over j in H of K0/K_RTA equals 1/lambda_max(G)")
mp.mp.dps = 40
okmin = True
rnd = random.Random(3)
for dv in (sp.Rational(1, 10), sp.Rational(1, 1000)):
    Cd = C.subs(dl, dv)
    Rd = sp.diag(*[Cd[i, i] for i in range(6)])
    e = sp.Matrix(eps)
    Pd = sp.Poly(P.subs(dl, dv), lam)
    lo, hi = max(iv[0] for iv in Pd.intervals(eps=sp.Rational(1, 10 ** 30)))
    lmax = (lo + hi) / 2
    # explicit minimiser: generalized eigenvector y* (C - lam R) y = 0, j* = C y*
    Mn = mp.matrix([[mp.mpf(sp.N(Cd[a, b] - lmax * Rd[a, b], 45)) for b in range(6)] for a in range(6)])
    U, Sv, V = mp.svd_r(Mn)
    ystar = mp.matrix([V[5, i] for i in range(6)])
    Cm = mp.matrix([[mp.mpf(sp.N(Cd[a, b], 45)) for b in range(6)] for a in range(6)])
    jstar = Cm * ystar
    # K0(j*) via bordered system
    Bd = mp.matrix(7, 7)
    for a in range(6):
        for b in range(6):
            Bd[a, b] = Cm[a, b]
        Bd[a, 6] = eps[a]
        Bd[6, a] = eps[a]
    sol = mp.lu_solve(Bd, mp.matrix(list(jstar) + [0]))
    K0s = sum(jstar[a] * sol[a] for a in range(6))
    KRs = sum(jstar[a] ** 2 / Cm[a, a] for a in range(6))
    okmin = okmin and abs(K0s / KRs - 1 / mp.mpf(sp.N(lmax, 45))) < mp.mpf("1e-25")
    # random currents in H
    worst = mp.mpf(10)
    for _ in range(300):
        jr = mp.matrix([mp.mpf(rnd.gauss(0, 1)) for _ in range(6)])
        jr = jr - mp.matrix(eps) * (sum(jr[a] * eps[a] for a in range(6)) / sum(x * x for x in eps))
        sol = mp.lu_solve(Bd, mp.matrix(list(jr) + [0]))
        K0r = sum(jr[a] * sol[a] for a in range(6))
        KRr = sum(jr[a] ** 2 / Cm[a, a] for a in range(6))
        worst = min(worst, K0r / KRr)
    okmin = okmin and worst >= 1 / mp.mpf(sp.N(lmax, 45)) - mp.mpf("1e-30")
    print(f"    delta = {float(dv):g}: 1/lambda_max = {mp.nstr(1 / mp.mpf(sp.N(lmax, 45)), 15)}, "
          f"K0/K_RTA at j* = {mp.nstr(K0s / KRs, 15)}, min over 300 random j = {mp.nstr(worst, 15)}")
check("min_j K0/K_RTA = 1/lambda_max(G): attained at j* = C y* (y* top generalized eigenvector), never undercut by random j in H",
      okmin)

# ================================================================ self-energy linewidths (main text cites Sec. S2; pending there)
print("=== self-energy linewidth with half weight on the repeated daughter (chain 1,2,4,8,16)")
eps, evdef, psi = nets[2]
n = 5
S = [stoich(n, p_, ds) for p_, ds in evdef]
rnd = random.Random(99)
okse = True
for trial in range(3):
    w = [sp.Rational(rnd.randint(1, 20), rnd.randint(1, 20)) for _ in range(n)]
    g = [sp.Rational(rnd.randint(1, 20), rnd.randint(1, 7)) for _ in S]
    C = sp.zeros(n, n)
    rt = [0] * n
    for s, gg in zip(S, g):
        z = sp.Matrix([sp.Rational(s[i]) / w[i] for i in range(n)])
        C += gg * z * z.T
        for i in range(n):
            if s[i] != 0:
                rt[i] += gg * z[i] ** 2 * (sp.Rational(1, 2) if s[i] == 2 else 1)
    e = sp.Matrix([w[i] * eps[i] for i in range(n)])
    j = C * sp.Matrix([w[i] * psi[i] for i in range(n)])
    K0, _, _ = K0_bordered(C, e, j)
    KR = sum(j[i] ** 2 / C[i, i] for i in range(n))
    KRt = sum(j[i] ** 2 / rt[i] for i in range(n))
    ratio = K0 / (KRt / 2)
    okse = okse and K0 == KR / 2 and ratio < 1
    print(f"    trial {trial}: K0/(K_RTA/2) = {K0 / (KR / 2)}, K0/(K~_RTA/2) = {sp.N(ratio, 12)}")
check("self-energy: with the self-coupling of the repeated daughter at half weight, K0 < K~_RTA/varsigma_max at the "
      "equality current of the chain (the inequality fails, as the main text states)", okse)
# a general fact: any linewidth r~ <= r with r~ < r on a mode that carries current breaks the inequality at an
# equality current. Check on the varsigma = 3 net with a 1% reduction of one lifetime-rate.
eps3, evdef3, psi3 = nets[3]
S3n = [stoich(8, p_, ds) for p_, ds in evdef3]
w = [sp.Rational(k + 2, k + 1) for k in range(8)]
g = [sp.Rational(k + 1, 3) for k in range(len(S3n))]
C = sp.zeros(8, 8)
for s, gg in zip(S3n, g):
    z = sp.Matrix([sp.Rational(s[i]) / w[i] for i in range(8)])
    C += gg * z * z.T
j = C * sp.Matrix([w[i] * psi3[i] for i in range(8)])
K0, _, _ = K0_bordered(C, sp.Matrix([w[i] * eps3[i] for i in range(8)]), j)
rt = [C[i, i] * (sp.Rational(99, 100) if i == 3 else 1) for i in range(8)]
KRt = sum(j[i] ** 2 / rt[i] for i in range(8))
check("any linewidth r~ <= diag C with r~ < diag C on a current-carrying mode breaks K0 >= K~_RTA/varsigma at an equality "
      "current (varsigma = 3 net, one rate reduced by 1%)", K0 < KRt / 3, f"K0/(K~_RTA/3) = {sp.N(K0 / (KRt / 3), 12)}")

print()
print(f"{sum(res)} of {len(res)} checks passed")
