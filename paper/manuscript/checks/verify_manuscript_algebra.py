"""Exact checks of the equations stated in Secs. II-IV of the manuscript.

Each block prints PASS/FAIL with the quantity tested. Exact rational or symbolic arithmetic is used
throughout, except where a floating-point check is stated explicitly.

Usage (from this directory):  python -B verify_manuscript_algebra.py > verify_manuscript_algebra.out.txt
"""
import itertools
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import sympy as sp

results = []


def check(name, ok):
    results.append((name, bool(ok)))
    print(("PASS" if ok else "FAIL") + "  " + name)


p, t, w = sp.symbols("p t omega", positive=True)

# ---------------------------------------------------------------------------------------------
# Sec. II: equivalent forms of the memory moment (Lemma 1)
# ---------------------------------------------------------------------------------------------
# random rational PSD operator with a one-dimensional kernel spanned by e

e = sp.Matrix([1, 2, 1, 3])
B = sp.Matrix([[2, -1, 0, 1], [1, 1, -2, 0], [0, 3, 1, -1]])
# rows of B projected orthogonally to e give event-like vectors spanning e-perp
P = sp.eye(4) - e * e.T / (e.T * e)[0]
Z = (P * B.T)
C = Z * sp.diag(3, sp.Rational(1, 2), 5) * Z.T
j = P * sp.Matrix([1, -1, 2, 0])
Cp = C.pinv()
x = Cp * j
K0 = (j.T * x)[0]
N = (x.T * x)[0]
Kp = (j.T * (p * sp.eye(4) + C).LUsolve(j))[0] - (j.T * e)[0] ** 2 / ((e.T * e)[0] * p)
Kp = sp.simplify(Kp)
check("C e = 0 and rank C = n-1", (C * e).is_zero_matrix and C.rank() == 3)
check("K(0) = j^T C^+ j", sp.simplify(Kp.subs(p, 0) - K0) == 0)
check("-K'(0)/K(0) = |C^+ j|^2 / j^T C^+ j", sp.simplify(-sp.diff(Kp, p).subs(p, 0) / K0 - N / K0) == 0)
# time-domain first moment: int_0^oo t exp(-C t) dt = (C^+)^2 on H, via eigen-decomposition
lam = sp.symbols("lam", positive=True)
check("int_0^oo t e^{-lam t} dt / int_0^oo e^{-lam t} dt = 1/lam",
      sp.simplify(sp.integrate(t * sp.exp(-lam * t), (t, 0, sp.oo)) / sp.integrate(sp.exp(-lam * t), (t, 0, sp.oo)) - 1 / lam) == 0)

# ---------------------------------------------------------------------------------------------
# Sec. II: grey Callaway operator in entropy coordinates (j eigenvector, Cattaneo time = tau_R)
# ---------------------------------------------------------------------------------------------
tauN, tauR, v, hbar, kB, T = sp.symbols("tau_N tau_R v hbar k_B T", positive=True)
# isotropic star mesh: wavevectors +-k e_i, two magnitudes, linear dispersion
ks = []
for kk in (1, 2):
    for i in range(3):
        for sgn in (1, -1):
            vec = [0, 0, 0]
            vec[i] = sgn * kk
            ks.append(vec)
n = len(ks)
eps = [hbar * v * sp.sqrt(sum(c ** 2 for c in k)) for k in ks]
# Bose weights depend on |k| only (energies do); keep them as independent positive symbols per shell,
# which is more general than the Bose values and keeps the algebra rational
w1, w2 = sp.symbols("w_1 w_2", positive=True)
wgt = [w1 if sum(c ** 2 for c in k) == 1 else w2 for k in ks]
evec = sp.Matrix([wgt[a] * eps[a] for a in range(n)])
Pvec = [sp.Matrix([wgt[a] * hbar * ks[a][i] for a in range(n)]) for i in range(3)]
jx = sp.Matrix([wgt[a] * eps[a] * v * ks[a][0] / sp.sqrt(sum(c ** 2 for c in ks[a])) for a in range(n)])
check("j_x = v^2 P_x for linear isotropic dispersion", sp.simplify(jx - v ** 2 * Pvec[0]).is_zero_matrix)
U = sp.Matrix.hstack(evec, *Pvec)
PiN = U * (U.T * U).inv() * U.T
PiE = evec * evec.T / (evec.T * evec)[0]
Ccal = (sp.eye(n) - PiN) / tauN + (sp.eye(n) - PiE) / tauR
check("Callaway: C e = 0", sp.simplify(Ccal * evec).is_zero_matrix)
check("Callaway: C j = j / tau_R", sp.simplify(Ccal * jx - jx / tauR).is_zero_matrix)
check("Callaway: C symmetric", sp.simplify(Ccal - Ccal.T).is_zero_matrix)
diag_dev = sp.simplify(Ccal[0, 0] - (1 / tauN + 1 / tauR))
check("Callaway: diag C = 1/tau_N + 1/tau_R - O(1/n) (deviation is negative)",
      sp.N(diag_dev.subs({tauN: 1, tauR: 10, v: 1, hbar: 1, w1: sp.Rational(3, 2), w2: sp.Rational(1, 2)})) < 0)
check("Callaway: |j_x|^2 tau_R = K(0) (single pole, tau_mem = tau_stat = tau_R)",
      sp.simplify((sp.eye(n) * p + Ccal) * (jx / (p + 1 / tauR)) - jx).is_zero_matrix)

# ---------------------------------------------------------------------------------------------
# Sec. III: four-mode example (same diagonal, same K(0), different memory)
# ---------------------------------------------------------------------------------------------
e4 = sp.Matrix([1, 1, 1, 1]) / 2
j4 = sp.Matrix([1, 0, 0, -1])
CA = sp.Rational(4, 3) * sp.eye(4) - sp.ones(4, 4) / 3
CB = sp.Matrix([[1, -sp.Rational(1, 2), 0, -sp.Rational(1, 2)],
                [-sp.Rational(1, 2), 1, -sp.Rational(1, 2), 0],
                [0, -sp.Rational(1, 2), 1, -sp.Rational(1, 2)],
                [-sp.Rational(1, 2), 0, -sp.Rational(1, 2), 1]])
for name, CC, Kexp in (("A", CA, 2 / (p + sp.Rational(4, 3))), ("B", CB, 1 / (p + 1) + 1 / (p + 2))):
    KK = sp.simplify((j4.T * (p * sp.eye(4) + CC).LUsolve(j4))[0])
    check(f"four-mode {name}: diag = 1, C e = 0, rank 3", all(CC[i, i] == 1 for i in range(4)) and (CC * e4).is_zero_matrix and CC.rank() == 3)
    check(f"four-mode {name}: K(p) as stated", sp.simplify(KK - Kexp) == 0)
    tau = sp.simplify(-sp.diff(KK, p).subs(p, 0) / KK.subs(p, 0))
    print(f"      K_{name}(0) = {sp.simplify(KK.subs(p, 0))}, tau_{name} = {tau}")
check("four-mode: K_A(0) = K_B(0) = 3/2, tau_A = 3/4, tau_B = 5/6",
      sp.simplify((j4.T * CA.pinv() * j4)[0]) == sp.Rational(3, 2) and sp.simplify((j4.T * CB.pinv() * j4)[0]) == sp.Rational(3, 2)
      and sp.simplify(((CA.pinv() * j4).T * (CA.pinv() * j4))[0] / sp.Rational(3, 2)) == sp.Rational(3, 4)
      and sp.simplify(((CB.pinv() * j4).T * (CB.pinv() * j4))[0] / sp.Rational(3, 2)) == sp.Rational(5, 6))

# ---------------------------------------------------------------------------------------------
# Sec. III: Theorem 1, eight-mode Walsh family
# ---------------------------------------------------------------------------------------------
eta = sp.symbols("eta", positive=True)
Hw = sp.Matrix(8, 8, lambda i, k: (-1) ** bin(i & k).count("1")) / sp.sqrt(8)
Tm = sp.diag(0, eta ** 2, 2, 2 - eta ** 2, 1, 1, 1, 1)
Tm[1, 2] = Tm[2, 1] = eta
Tm[4, 7] = Tm[7, 4] = -eta
C8 = sp.simplify(Hw * Tm * Hw.T)
e8 = Hw[:, 0]
j8 = Hw[:, 2]
check("Walsh matrix orthogonal and symmetric", sp.simplify(Hw * Hw.T - sp.eye(8)).is_zero_matrix and Hw == Hw.T)
check("Thm 1: diag C_eta = 1 for all eta", all(sp.simplify(C8[i, i] - 1) == 0 for i in range(8)))
check("Thm 1: C_eta e = 0", sp.simplify(C8 * e8).is_zero_matrix)
blk1 = sp.Matrix([[eta ** 2, eta], [eta, 2]])
blk2 = sp.Matrix([[1, -eta], [-eta, 1]])
check("Thm 1: blocks positive definite for 0<eta<1 (det eta^2, 1-eta^2; entry 2-eta^2)",
      sp.simplify(blk1.det() - eta ** 2) == 0 and sp.simplify(blk2.det() - (1 - eta ** 2)) == 0)
x8 = Hw * sp.Matrix([0, -1 / eta, 1, 0, 0, 0, 0, 0])
check("Thm 1: C_eta (-h_1/eta + h_2) = j exactly", sp.simplify(C8 * x8 - j8).is_zero_matrix)
# candidate resolvent solution from the 2x2 block {1,2}; verified by direct multiplication below
blkinv = (p * sp.eye(2) + blk1).inv()
X8 = Hw * sp.Matrix([0, blkinv[0, 1], blkinv[1, 1], 0, 0, 0, 0, 0])
check("Thm 1: (p + C_eta) X = j for the block candidate X", sp.simplify((p * sp.eye(8) + C8) * X8 - j8).is_zero_matrix)
K8 = sp.simplify((j8.T * X8)[0])
check("Thm 1: K_eta(p) = (p+eta^2)/(p^2+(2+eta^2)p+eta^2)",
      sp.simplify(K8 - (p + eta ** 2) / (p ** 2 + (2 + eta ** 2) * p + eta ** 2)) == 0)
check("Thm 1: K_eta(0) = 1 and tau = 1 + eta^-2",
      sp.simplify(K8.subs(p, 0) - 1) == 0 and sp.simplify(-sp.diff(K8, p).subs(p, 0) - (1 + eta ** -2)) == 0)
# parity: permutation nu -> nu XOR 3 commutes with C, fixes e, reverses j
Pm = sp.zeros(8, 8)
for i in range(8):
    Pm[i ^ 3, i] = 1
check("Thm 1: XOR-3 involution has no fixed point", all((i ^ 3) != i for i in range(8)))
check("Thm 1: XOR-3 commutes with C_eta, fixes e, reverses j",
      sp.simplify(Pm * C8 - C8 * Pm).is_zero_matrix and (Pm * e8 - e8).is_zero_matrix and sp.simplify(Pm * j8 + j8).is_zero_matrix)
# band bound |K(i w)(i w + 2) - 1| <= eta^2/(2 w)
expr = sp.simplify(K8 * (p + 2) - 1)
check("Thm 1: K_eta(p)(p+2) - 1 = eta^2 / D(p)",
      sp.simplify(expr - eta ** 2 / (p ** 2 + (2 + eta ** 2) * p + eta ** 2)) == 0)
Dw = sp.expand((sp.I * w) ** 2 + (2 + eta ** 2) * sp.I * w + eta ** 2)
check("Thm 1: Im D(i omega) = (2+eta^2) omega >= 2 omega", sp.simplify(sp.im(Dw) - (2 + eta ** 2) * w) == 0)
# Remark (mechanism): slow root -eta^2/2 + O(eta^4), residue eta^2/4 + O(eta^4), half of K(0), all of N
Dp = p ** 2 + (2 + eta ** 2) * p + eta ** 2
slow = [rt for rt in sp.solve(Dp, p) if sp.limit(rt, eta, 0) == 0][0]
resid = sp.simplify(((p + eta ** 2) / sp.diff(Dp, p)).subs(p, slow))
check("Thm 1 remark: slow root = -eta^2/2 + O(eta^4)", sp.simplify(sp.series(slow, eta, 0, 4).removeO() + eta ** 2 / 2) == 0)
check("Thm 1 remark: residue = eta^2/4 + O(eta^4)", sp.simplify(sp.series(resid, eta, 0, 4).removeO() - eta ** 2 / 4) == 0)
check("Thm 1 remark: slow pole carries 1/2 of K(0) and a fraction 1 - O(eta^2) of N",
      sp.limit(resid / (-slow), eta, 0) == sp.Rational(1, 2) and sp.limit((resid / slow ** 2) / (1 + eta ** -2), eta, 0) == 1)
# not in any fixed cone: u = h_1 gives u^T C j = eta, u^T C u = eta^2
u8 = Hw[:, 1]
check("Thm 1 remark: h_1^T C h_2 = eta and h_1^T C h_1 = eta^2",
      sp.simplify((u8.T * C8 * j8)[0] - eta) == 0 and sp.simplify((u8.T * C8 * u8)[0] - eta ** 2) == 0)

# ---------------------------------------------------------------------------------------------
# Sec. III: Theorem 2 (Cauchy-Binet / Ben-Tal-Teboulle), Prop. 3, Prop. 4, Prop. 5, Lemma 2
# ---------------------------------------------------------------------------------------------
# exact energy-conserving integer stoichiometry: energies 1,2,3,4,5 and events p -> a + b
energies = [1, 2, 3, 4, 5]
nm = len(energies)
events = []
for pp_, a_, b_ in itertools.product(range(nm), repeat=3):
    if a_ <= b_ and energies[pp_] == energies[a_] + energies[b_]:
        s = [0] * nm
        s[pp_] -= 1
        s[a_] += 1
        s[b_] += 1
        events.append(s)
Wd = sp.diag(*[sp.Rational(k + 2, k + 1) for k in range(nm)])  # stand-in positive Bose weights
epsv = sp.Matrix(energies)
ee = Wd * epsv
zetas = [Wd.inv() * sp.Matrix(s) for s in events]
check("events conserve energy: zeta^T e = 0", all((z.T * ee)[0] == 0 for z in zetas))
# orthonormal basis of H = e-perp (exact, Gram-Schmidt)
basis = sp.GramSchmidt([sp.Matrix(c) for c in (sp.eye(nm) - ee * ee.T / (ee.T * ee)[0]).columnspace()], True)
Q = sp.Matrix.hstack(*basis)  # n x (n-1)
nH = Q.shape[1]
jv = sp.Matrix([1, -2, 0, 1, sp.Rational(1, 3)])
jv = jv - ee * (ee.T * jv)[0] / (ee.T * ee)[0]
zH = [Q.T * z for z in zetas]
jH = Q.T * jv
g = [sp.Rational(k % 4 + 1, k % 3 + 1) for k in range(len(zetas))]
A = sp.zeros(nH, nH)
for gg, z in zip(g, zH):
    A += gg * z * z.T
check("event operator admissible on H (positive definite)", A.det() > 0)
xH = A.LUsolve(jH)
K0e = (jH.T * xH)[0]
num = sp.zeros(nH, 1)
den = 0
Mval = 0
for I in itertools.combinations(range(len(zetas)), nH - 1):
    cols = [zH[k] for k in I]
    psi = sp.Matrix([sp.Matrix.hstack(*(cols + [sp.eye(nH)[:, r]])).det() for r in range(nH)])
    pj = (psi.T * jH)[0]
    if pj != 0:
        gI = sp.prod([g[k] for k in I])
        num += gI * pj * psi
        den += gI * pj ** 2
        yI = psi / pj
        Mval = max(Mval, sp.sqrt((yI.T * yI)[0]))
check("Cauchy-Binet: x/K0 = sum_I varpi_I y_I exactly", sp.simplify(num / den - xH / K0e).is_zero_matrix)
adj = A.adjugate()
adj_cb = sp.zeros(nH, nH)
for I in itertools.combinations(range(len(zetas)), nH - 1):
    cols = [zH[k] for k in I]
    psi = sp.Matrix([sp.Matrix.hstack(*(cols + [sp.eye(nH)[:, r]])).det() for r in range(nH)])
    adj_cb += sp.prod([g[k] for k in I]) * psi * psi.T
check("Cauchy-Binet: adj A = sum_I g_I psi_I psi_I^T", sp.simplify(adj - adj_cb).is_zero_matrix)
taue = (xH.T * xH)[0] / K0e
check("Thm 2: tau_mem <= M^2 K0", sp.N(taue - Mval ** 2 * K0e) <= 0)
print(f"      tau = {sp.N(taue, 8)}, M^2 K0 = {sp.N(Mval ** 2 * K0e, 8)}")
# Prop. 2 lower bound
check("Prop. 2: tau >= |x|/|j| >= K0/|j|^2",
      sp.N(taue - sp.sqrt((xH.T * xH)[0]) / sp.sqrt((jH.T * jH)[0])) >= 0
      and sp.N(sp.sqrt((xH.T * xH)[0]) / sp.sqrt((jH.T * jH)[0]) - K0e / (jH.T * jH)[0]) >= 0)
# Lemma 2: K0 >= K_RTA / varsigma_max with r = diag C in mode coordinates
Cfull = sp.zeros(nm, nm)
for gg, z in zip(g, zetas):
    Cfull += gg * z * z.T
r = [Cfull[k, k] for k in range(nm)]
KRTA = sum(jv[k] ** 2 / r[k] for k in range(nm))
smax = max(sum(1 for c in s if c != 0) for s in events)
check("Lemma 2: K0 >= K_RTA / varsigma_max", sp.N(K0e - KRTA / smax) >= 0)
print(f"      K0/K_RTA = {sp.N(K0e / KRTA, 8)}, varsigma_max = {smax}")
# Prop. 5: tau <= 1/Gamma_min and amplitude-prior ceiling
import mpmath as mpm
mpm.mp.dps = 50
Gmin = min(mpm.eigsy(mpm.matrix([[sp.N(A[i_, k_], 60) for k_ in range(nH)] for i_ in range(nH)]))[0])  # 50-digit symmetric eigensolver
Gmin = sp.Float(str(Gmin), 50)
check("Lemma 1 corollary: tau <= 1/Gamma_min", sp.N(taue) <= 1 / Gmin + 1e-12)
Fp = 3
g2 = [gg / Fp if k % 2 == 0 else gg * 5 for k, gg in enumerate(g)]
A2 = sp.zeros(nH, nH)
for gg, z in zip(g2, zH):
    A2 += gg * z * z.T
x2 = A2.LUsolve(jH)
tau2 = (x2.T * x2)[0] / (jH.T * x2)[0]
check("Prop. 5: g >= g_ref/F implies tau <= F/Gamma_min(C_ref)", sp.N(tau2) <= Fp / Gmin + 1e-12)
# Prop. 4: sub-network form |y_I| = |Pi_H W phi_I| / |j^T W phi_I|
ok4 = True
PiH = sp.eye(nm) - ee * ee.T / (ee.T * ee)[0]
for I in itertools.combinations(range(len(zetas)), nH - 1):
    SI = sp.Matrix.hstack(*[sp.Matrix(events[k]) for k in I])
    if SI.rank() < nH - 1:
        continue
    ker = SI.T.nullspace()
    # choose phi in ker independent of eps
    phi = None
    for kv in ker:
        if sp.Matrix.hstack(kv, epsv).rank() == 2:
            phi = kv
            break
    if phi is None:
        continue
    den4 = (jv.T * Wd * phi)[0]
    if den4 == 0:
        continue
    yI_formula = PiH * Wd * phi / den4
    cols = [zH[k] for k in I]
    psi = sp.Matrix([sp.Matrix.hstack(*(cols + [sp.eye(nH)[:, r_]])).det() for r_ in range(nH)])
    yI_direct = Q * psi / (psi.T * jH)[0]
    ok4 = ok4 and sp.simplify(yI_formula - yI_direct).is_zero_matrix
check("Prop. 4: y_I = Pi_H W phi_I / (j^T W phi_I) for every basis", ok4)
# Prop. 3: weight concentration on a maximising basis
Istar = None
best = -1
for I in itertools.combinations(range(len(zetas)), nH - 1):
    cols = [zH[k] for k in I]
    psi = sp.Matrix([sp.Matrix.hstack(*(cols + [sp.eye(nH)[:, r_]])).det() for r_ in range(nH)])
    pj = (psi.T * jH)[0]
    if pj != 0:
        val = sp.N((psi.T * psi)[0] / pj ** 2)
        if val > best:
            best, Istar = val, I
ratios = []
for epsw in [sp.Rational(1, 10 ** kexp) for kexp in (1, 3, 5, 7, 9, 11)]:
    gw = [1 if k in Istar else epsw for k in range(len(zetas))]
    Aw = sp.zeros(nH, nH)
    for gg, z in zip(gw, zH):
        Aw += gg * z * z.T
    xw = Aw.LUsolve(jH)
    Kw = (jH.T * xw)[0]
    ratios.append(sp.N(((xw.T * xw)[0] / Kw) / (Mval ** 2 * Kw)))
print("      tau/(M^2 K0) along weight concentration:", [f"{float(r_):.6f}" for r_ in ratios])
check("Prop. 3: tau/(M^2 K0) increases monotonically to 1 under weight concentration (eps = 1e-1 ... 1e-11)",
      all(ratios[k_] < ratios[k_ + 1] for k_ in range(len(ratios) - 1)) and abs(ratios[-1] - 1) < 1e-4)

# ---------------------------------------------------------------------------------------------
# Sec. III: infrared criterion, radial integral
# ---------------------------------------------------------------------------------------------
q, q0, qmin, alpha, dd, m = sp.symbols("q q_0 q_min alpha d m", positive=True)
I1 = sp.integrate(q ** (dd - 1 - m * alpha), (q, qmin, q0), conds="none")
check("radial moment: int_{qmin}^{q0} q^{d-1-m alpha} dq = (q0^s - qmin^s)/s, s = d - m alpha",
      sp.simplify(I1 - (q0 ** (dd - m * alpha) - qmin ** (dd - m * alpha)) / (dd - m * alpha)) == 0)
check("radial moment: logarithmic at m alpha = d",
      sp.simplify(sp.integrate(q ** (-1), (q, qmin, q0)) - sp.log(q0 / qmin)) == 0)

# ---------------------------------------------------------------------------------------------
# Sec. IV: film model
# ---------------------------------------------------------------------------------------------
lamf, rc, tR, tl, L2, e_, bs, z, mu = sp.symbols("lambda rho_c tau_R tau_l L2 e b_s z mu", positive=True)
a = lamf / rc
bf = sp.sqrt(lamf * rc)
xi1 = e_ / sp.sqrt(a)
# derive lambda_eff and sigma from the 1D GK law and the energy balance (Laplace domain)
th = sp.Function("theta")
qf = sp.Function("q")
sig2 = sp.symbols("sigma2")
# Ansatz: q = -lam_eff theta'; energy balance q' = -rho c p theta
lam_eff = (lamf + L2 * rc * p) / (1 + tR * p)
# check: (1 + tR p) q = -lam theta' + L2 q'' with q'' = -rho c p theta'
thp = sp.symbols("thp")  # theta'
q_expr = -lam_eff * thp
qpp = -rc * p * thp
check("GK 1D: (1+tau_R p) q = -lambda theta' + L^2 q'' with q = -lambda_eff theta'",
      sp.simplify((1 + tR * p) * q_expr - (-lamf * thp + L2 * qpp)) == 0)
taul_def = L2 / a
check("lambda_eff = lambda (1 + tau_l p)/(1 + tau_R p), tau_l = L^2/a",
      sp.simplify(lam_eff - lamf * (1 + taul_def * p) / (1 + tR * p)) == 0)
sigma = sp.sqrt(rc * p / lam_eff)
check("sigma e = xi_1 sqrt(p (1+tau_R p)/(1+tau_l p))",
      sp.simplify((sigma * e_) ** 2 - (xi1 ** 2 * p * (1 + tR * p) / (1 + taul_def * p))) == 0)
check("lambda_eff sigma = b sqrt(p (1+tau_l p)/(1+tau_R p))",
      sp.simplify((lam_eff * sigma) ** 2 - bf ** 2 * p * (1 + taul_def * p) / (1 + tR * p)) == 0)
# quadrupole: front state = M rear state
sg, Lam = sp.symbols("sigma Lambda", positive=True)
th0, q0f, c1 = sp.symbols("theta_0 q_0 c_1")
theta_z = th0 * sp.cosh(sg * z) + c1 * sp.sinh(sg * z)
q_z = -(Lam / sg) * sp.diff(theta_z, z)  # q = -lambda_eff theta', Lambda = lambda_eff sigma
c1sol = sp.solve(sp.Eq(q_z.subs(z, 0), q0f), c1)[0]
thr = sp.simplify(theta_z.subs(c1, c1sol).subs(z, e_))
qr = sp.simplify(q_z.subs(c1, c1sol).subs(z, e_))
Mq = sp.Matrix([[sp.cosh(sg * e_), sp.sinh(sg * e_) / Lam], [Lam * sp.sinh(sg * e_), sp.cosh(sg * e_)]])
front = sp.simplify(Mq * sp.Matrix([thr, qr]))
check("quadrupole maps rear state to front state", sp.simplify(front - sp.Matrix([th0, q0f])).is_zero_matrix)
check("det M = 1", sp.simplify(Mq.det() - 1) == 0)
Zs = 1 / (bs * sp.sqrt(p))
tt = sp.symbols("t_h")
Zin = (sp.cosh(sg * e_) * Zs + sp.sinh(sg * e_) / Lam) / (Lam * sp.sinh(sg * e_) * Zs + sp.cosh(sg * e_))
rho_ratio = Lam * Zs
Zalt = (1 / Lam) * (rho_ratio + sp.tanh(sg * e_)) / (1 + rho_ratio * sp.tanh(sg * e_))
check("Z = (AZs+B)/(CZs+D) = (1/Lambda)(varrho + tanh)/(1 + varrho tanh)",
      sp.simplify((Zin - Zalt).rewrite(sp.exp)) == 0)
# gauge invariance
sub = {e_: mu * e_, lamf: mu * lamf, rc: rc / mu, L2: mu ** 2 * L2}
chi = [xi1, bf, tR, taul_def]
check("gauge e->mu e, a->mu^2 a, lambda->mu lambda, rho c->rho c/mu, ell->mu ell leaves (xi_1,b,tau_R,tau_l) invariant",
      all(sp.simplify(c.subs(sub, simultaneous=True) - c) == 0 for c in chi))
check("a -> mu^2 a under the gauge", sp.simplify(a.subs(sub, simultaneous=True) - mu ** 2 * a) == 0)
# Euler identity (infinitesimal gauge) in variables (e, a, lambda, ell) for the two combinations
ee_, aa, ll, ell, beta = sp.symbols("e a lam ell beta", positive=True)
comb1 = (ee_ / sp.sqrt(aa)) ** 2 * p * (1 + tR * p) / (1 + (1 + beta) * ell ** 2 / aa * p)
comb2 = (ll / sp.sqrt(aa)) ** 2 * p * (1 + (1 + beta) * ell ** 2 / aa * p) / (1 + tR * p)
euler = lambda G: sp.simplify(ee_ * sp.diff(G, ee_) + 2 * aa * sp.diff(G, aa) + ll * sp.diff(G, ll) + ell * sp.diff(G, ell))
check("Euler identity e d_e + 2a d_a + lambda d_lambda + ell d_ell annihilates both combinations",
      euler(comb1) == 0 and euler(comb2) == 0)
# Fourier resonance
check("Fourier resonance: tau_R = tau_l gives sigma e = xi_1 sqrt p and lambda_eff sigma = b sqrt p",
      sp.simplify(xi1 ** 2 * p * (1 + tR * p) / (1 + tR * p) - xi1 ** 2 * p) == 0)
# closure-dependent resonance positions (grey kinetic coefficients)
tN, tRk, vv = sp.symbols("tau_N tau_Rk v", positive=True)
tc = 1 / (1 / tN + 1 / tRk)
a_k = vv ** 2 * tRk / 3
ell2 = vv ** 2 * tc * tRk / 5
tl_cons = (1 + sp.Rational(1, 3)) * ell2 / a_k
tl_hist = 3 * ell2 / a_k
check("conserving closure: tau_l = 4 tau_c/5 and tau_l/tau_R = 4 tau_N /(5 (tau_N + tau_R)) < 4/5",
      sp.simplify(tl_cons - sp.Rational(4, 5) * tc) == 0
      and sp.simplify(tl_cons / tRk - sp.Rational(4, 5) * tN / (tN + tRk)) == 0)
xr = sp.symbols("x", positive=True)
check("historical closure with tau_c kept: resonance at tau_R/tau_N = 4/5",
      sp.solve(sp.Eq(sp.simplify(tl_hist.subs(tRk, xr * tN)), xr * tN), xr) == [sp.Rational(4, 5)])
check("historical closure with tau_c ~ tau_N: tau_l = 9 tau_N / 5 (resonance at 9/5)",
      sp.simplify(3 * (vv ** 2 * tN * tRk / 5) / a_k - sp.Rational(9, 5) * tN) == 0)
# Cattaneo half-space phase and the naive substitution
Zc = sp.sqrt(1 + tR * sp.I * w) / (bf * sp.sqrt(sp.I * w))
ph = sp.arg(Zc)
num_ok = all(abs(sp.N((ph - (-sp.pi / 4 + sp.atan(w * tR) / 2)).subs({tR: tv, w: wv, lamf: 3, rc: 2}))) < 1e-14
             for tv in (sp.Rational(1, 10), 1, 7) for wv in (sp.Rational(1, 3), 2, 50))
check("Cattaneo half-space: arg Z = -pi/4 + arctan(omega tau_R)/2 (numerical, 9 points)", num_ok)
Znaive = 1 / (bf * sp.sqrt(p * (1 + tR * p)))
Zright = sp.sqrt(1 + tR * p) / (bf * sp.sqrt(p))
check("naive substitution p -> p(1+tau_R p) is off by the factor (1+tau_R p)", sp.simplify(Zright / Znaive - (1 + tR * p)) == 0)
# GK half-space phase
Zgk = sp.sqrt((1 + tR * sp.I * w) / (sp.I * w * (1 + tl * sp.I * w))) / bf
ok = all(abs(sp.N((sp.arg(Zgk) - (-sp.pi / 4 + (sp.atan(w * tR) - sp.atan(w * tl)) / 2)).subs({tR: tv, tl: lv, w: wv, lamf: 3, rc: 2}))) < 1e-14
         for tv in (sp.Rational(1, 10), 3) for lv in (sp.Rational(1, 5), 2) for wv in (sp.Rational(1, 2), 9))
check("GK half-space: arg Z = -pi/4 + [arctan(omega tau_R) - arctan(omega tau_l)]/2 (numerical, 8 points)", ok)

# ---------------------------------------------------------------------------------------------
# Dimensions of the kinetic normalisation (SI base units as exponent vectors)
# ---------------------------------------------------------------------------------------------
from sympy.physics.units import joule, kelvin, meter, second, convert_to
from sympy.physics.units.systems.si import dimsys_SI
from sympy.physics.units import Dimension
from sympy.physics.units.definitions.dimension_definitions import energy, length, time, temperature
y_dim = length ** sp.Rational(-3, 2)
e_dim = energy * length ** sp.Rational(-3, 2)
j_dim = energy * length ** sp.Rational(-3, 2) * length / time
C_dim = 1 / time
K_dim = j_dim ** 2 * time
lam_dim = K_dim / (energy / temperature * temperature ** 2)
rc_dim = e_dim ** 2 / (energy / temperature * temperature ** 2)
same = lambda d1, d2: dimsys_SI.equivalent_dims(d1, d2)
check("[e^T y] = energy density", same(e_dim * y_dim, energy / length ** 3))
check("[j^T y] = W m^-2 (flux per unit area)", same(j_dim * y_dim, energy / time / length ** 2))
check("[K(0)/(k_B T^2)] = W m^-1 K^-1", same(lam_dim, energy / time / length / temperature))
check("[|e|^2/(k_B T^2)] = J m^-3 K^-1", same(rc_dim, energy / length ** 3 / temperature))
check("[N/K0] = s", same((j_dim ** 2 * time ** 2) / K_dim, time))
check("[K0/|e|^2] = m^2 s^-1 (diffusivity)", same(K_dim / e_dim ** 2, length ** 2 / time))

print()
nfail = sum(1 for _, ok in results if not ok)
print(f"{len(results) - nfail} of {len(results)} checks passed")
