"""Independent checks of Example 1, Theorem 1 (i)-(v), Remark 1 (mechanism) and Remark 3 (no fixed cone).

Methods differ from checks/verify_manuscript_algebra.py on purpose:
  * K_eta(p) from the full 8x8 resolvent (no reduction to the {1,2} block);
  * spectral data (slow root, residues) from the exact characteristic polynomial of C_eta and,
    independently, from 60-digit eigen-decompositions at several eta;
  * the time-domain first moment from Phi(t) = j^T exp(-C t) j by numerical quadrature;
  * the parity involution from its permutation matrix, not from Walsh characters.

Usage: python -B audit_thm1_walsh.py > audit_thm1_walsh.out.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


p = sp.symbols("p", positive=True)
eta = sp.symbols("eta", positive=True)

# ------------------------------------------------------------------ Example 1
e4 = sp.Matrix([1, 1, 1, 1]) / 2
j4 = sp.Matrix([1, 0, 0, -1])
A4 = sp.Matrix([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]])  # cycle 1-2-3-4-1
CA = sp.Rational(4, 3) * sp.eye(4) - sp.ones(4, 4) / 3
CB = sp.eye(4) - A4 / 2
for nm, CC, Kst, tst in (("A", CA, 2 / (p + sp.Rational(4, 3)), sp.Rational(3, 4)),
                         ("B", CB, 1 / (p + 1) + 1 / (p + 2), sp.Rational(5, 6))):
    # spectral route: eigen-decomposition of the symmetric rational matrix
    Kspec = 0
    for val, mult, vecs in CC.eigenvects():
        if val == 0:
            continue
        Q, _ = sp.Matrix.hstack(*vecs).QRdecomposition()
        for c in range(Q.shape[1]):
            u = Q[:, c]
            Kspec += (u.T * j4)[0] ** 2 / (p + val)
    Kspec = sp.simplify(Kspec)
    lap = all(sum(CC[i, :]) == 0 for i in range(4)) and all(CC[i, k] <= 0 for i in range(4) for k in range(4) if i != k)
    check(f"Ex.1 C_{nm}: unit diagonal, Ce=0, rank 3, graph Laplacian (zero row sums, nonpositive off-diagonal)",
          all(CC[i, i] == 1 for i in range(4)) and (CC * e4).is_zero_matrix and CC.rank() == 3 and lap)
    check(f"Ex.1 K_{nm}(p) from eigen-decomposition equals Eq. (17)", sp.simplify(Kspec - Kst) == 0, str(Kspec))
    K0 = sp.simplify(Kspec.subs(p, 0))
    tau = sp.simplify(-sp.diff(Kspec, p).subs(p, 0) / K0)
    check(f"Ex.1 K_{nm}(0) = 3/2 and tau_mem = {tst}", K0 == sp.Rational(3, 2) and tau == tst)

# ------------------------------------------------------------------ Theorem 1 family
# Exact rational form: H = H0/sqrt(8) with H0 the integer +-1 Walsh matrix, so C = H0 T H0^T / 8.
n = 8
H0 = sp.Matrix(n, n, lambda i, k: (-1) ** bin(i & k).count("1"))
T = sp.zeros(n, n)
for i, v in enumerate([0, eta ** 2, 2, 2 - eta ** 2, 1, 1, 1, 1]):
    T[i, i] = v
T[1, 2] = T[2, 1] = eta
T[4, 7] = T[7, 4] = -eta
C = (H0 * T * H0.T / 8).applyfunc(sp.expand)
e0 = H0[:, 0]          # sqrt(8) e
j0 = H0[:, 2]          # sqrt(8) j ; ratios below are normalization-free after dividing by 8
check("H0 H0^T = 8 I and H0 symmetric", (H0 * H0.T - 8 * sp.eye(n)).is_zero_matrix and H0 == H0.T)
lam = sp.Symbol("lam")
# (i) characteristic polynomial computed from C itself (Berkowitz, polynomial entries)
charC = sp.expand(C.charpoly(lam).as_expr())
expect = lam * (lam ** 2 - (2 + eta ** 2) * lam + eta ** 2) * (lam - (2 - eta ** 2)) * (lam - 1) ** 2 * ((lam - 1) ** 2 - eta ** 2)
check("(i) char. polynomial of C_eta = lam (lam^2-(2+eta^2)lam+eta^2)(lam-2+eta^2)(lam-1)^2((lam-1)^2-eta^2)",
      sp.expand(charC - expect) == 0)
check("(i) C_eta e = 0; the zero eigenvalue is simple (nonzero roots have product eta^2(2-eta^2)(1-eta^2) > 0 on (0,1))",
      (C * e0).is_zero_matrix and sp.factor(sp.expand(expect / lam).subs(lam, 0)) != 0)
# symmetric matrix with nonnegative roots: PSD; positivity of each root factor on (0,1)
check("(i) every nonzero root positive on (0,1): sum 2+eta^2>0, product eta^2>0, 2-eta^2>0, 1-eta>0",
      all(float(v) > 0 for et in (0.01, 0.5, 0.99) for v in (2 + et ** 2, et ** 2, 2 - et ** 2, 1 - et)))
# (ii) diagonal
check("(ii) diag C_eta = 1 identically in eta", all(sp.expand(C[i, i] - 1) == 0 for i in range(n)))
# (iii) K from the full 8x8 resolvent (p>0, p+C invertible on R^8); K = j^T (p+C)^{-1} j = j0^T(...)j0/8
from sympy.polys.matrices import DomainMatrix
FF = sp.QQ.frac_field(p, eta)


def dsolve(Msym, bsym):
    """Exact solve over the fraction field QQ(p, eta)."""
    Md = DomainMatrix.from_Matrix(Msym).convert_to(FF)
    bd = DomainMatrix.from_Matrix(bsym).convert_to(FF)
    return Md.lu_solve(bd).to_Matrix()


Kfull = sp.factor(sp.cancel((j0.T * dsolve(p * sp.eye(n) + C, j0))[0] / 8))
D = p ** 2 + (2 + eta ** 2) * p + eta ** 2
check("(iii) K_eta(p) from the full 8x8 resolvent = (p+eta^2)/D(p)", sp.cancel(Kfull - (p + eta ** 2) / D) == 0, str(Kfull))
# x = C^+ j by the regularized solve (C + e e^T/|e|^2) x = j, then projection on H (exact)
Pe = e0 * e0.T / 8
x0 = dsolve(C + Pe, j0)
x0 = (x0 - Pe * x0).applyfunc(sp.cancel)
K0 = sp.cancel((j0.T * x0)[0] / 8)
N = sp.cancel((x0.T * x0)[0] / 8)
check("(iii) via the regularized solve: K0 = 1, N = 1 + eta^-2, tau = N/K0 = 1 + eta^-2",
      sp.cancel(K0 - 1) == 0 and sp.cancel(N - 1 - eta ** -2) == 0)
# H, e, j in normalized form for later blocks
H = H0 / sp.sqrt(8)
e = H[:, 0]
j = H[:, 2]
# (iv) parity involution from its permutation matrix
P = sp.zeros(n, n)
for i in range(n):
    P[i ^ 3, i] = 1
check("(iv) i -> i XOR 3 is an involution without fixed points", (P * P - sp.eye(n)).is_zero_matrix and all((i ^ 3) != i for i in range(n)))
check("(iv) P C_eta = C_eta P, P e = e, P j = -j", (P * C - C * P).applyfunc(sp.expand).is_zero_matrix
      and (P * e0 - e0).is_zero_matrix and (P * j0 + j0).is_zero_matrix)
# the claim in the proof: P h_k = (-1)^<3,k> h_k for every k
check("(iv) proof step: P h_k = (-1)^popcount(3&k) h_k for k = 0..7",
      all((P * H0[:, k] - (-1) ** bin(3 & k).count("1") * H0[:, k]).is_zero_matrix for k in range(n)))
# energies of partners: e_i = e_{i xor 3} trivially (e constant); j odd: j_i = -j_{i xor 3}
check("(iv) j_{i xor 3} = -j_i for all i (current odd, energy even)", all(j0[i ^ 3] + j0[i] == 0 for i in range(n)))

# (v) band bound, symbolic identity then dense numerical test of the inequality
w = sp.symbols("omega", positive=True)
check("(v) K(p)(p+2) - 1 = eta^2/D(p)", sp.cancel(Kfull * (p + 2) - 1 - eta ** 2 / D) == 0)
worst = 0.0
for et in (mp.mpf("0.9"), mp.mpf("0.5"), mp.mpf("0.1"), mp.mpf("0.01")):
    for om in [mp.mpf(10) ** (k / mp.mpf(8)) for k in range(-40, 41)]:
        Kv = (1j * om + et ** 2) / ((1j * om) ** 2 + (2 + et ** 2) * 1j * om + et ** 2)
        lhs = abs(Kv * (1j * om + 2) - 1)
        worst = max(worst, float(lhs / (et ** 2 / (2 * om))))
check("(v) |K(i w)(i w+2)-1| <= eta^2/(2 w) on 81 x 4 points (max ratio LHS/RHS < 1)", worst < 1, f"max ratio {worst:.6f}")

# ------------------------------------------------------------------ Remark 1 (mechanism), independent routes
sfree = sp.Symbol("s")  # unrestricted symbol: both roots of D are negative
Ds = D.subs(p, sfree)
roots = sp.solve(Ds, sfree)
check("Rem.1 D has two real negative roots for 0<eta<1", len(roots) == 2 and all(float(r.subs(eta, sp.Rational(1, 3))) < 0 for r in roots))
slow = [r for r in roots if sp.limit(r, eta, 0) == 0][0]
fast = [r for r in roots if r != slow][0]
res_s = sp.simplify(((sfree + eta ** 2) / sp.diff(Ds, sfree)).subs(sfree, slow))
res_f = sp.simplify(((sfree + eta ** 2) / sp.diff(Ds, sfree)).subs(sfree, fast))
ser_slow = sp.series(slow, eta, 0, 6).removeO()
ser_res = sp.series(res_s, eta, 0, 6).removeO()
check("Rem.1 slow root = -eta^2/2 + eta^4/8 + O(eta^6)", sp.simplify(ser_slow - (-eta ** 2 / 2 + eta ** 4 / 8)) == 0, str(ser_slow))
check("Rem.1 residue of slow pole = eta^2/4 + O(eta^4)", sp.simplify(sp.series(res_s, eta, 0, 4).removeO() - eta ** 2 / 4) == 0, str(ser_res))
check("Rem.1 residues sum to |j|^2 = 1", sp.simplify(res_s + res_f - 1) == 0)
fracK = sp.limit(res_s / (-slow), eta, 0)
fracN_div = sp.limit((res_s / slow ** 2) - eta ** -2, eta, 0)
check("Rem.1 slow pole carries K0 fraction -> 1/2 and N = eta^-2 + O(1) (all of the divergence)",
      fracK == sp.Rational(1, 2) and fracN_div.is_finite, f"N_slow - eta^-2 -> {fracN_div}")
# numerical eigen-decomposition at small eta (60 digits) reproduces the residues as current weights (u^T j)^2
ok_num = True
for et in (mp.mpf("0.2"), mp.mpf("0.05"), mp.mpf("0.01")):
    Cn = mp.matrix([[mp.mpf(sp.Rational(C[i, k].subs(eta, sp.Rational(str(et)))).p) / mp.mpf(sp.Rational(C[i, k].subs(eta, sp.Rational(str(et)))).q) for k in range(n)] for i in range(n)])
    ev, U = mp.eigsy(Cn)
    jn = mp.matrix([mp.mpf(int(j0[i])) / mp.sqrt(8) for i in range(n)])
    pos = [k for k in range(n) if ev[k] > mp.mpf(10) ** -40]
    kslow = min(pos, key=lambda k: ev[k])
    wslow = sum(U[i, kslow] * jn[i] for i in range(n)) ** 2
    s_ex = -mp.mpf(sp.N(slow.subs(eta, sp.Rational(str(et))), 70))
    r_ex = mp.mpf(sp.N(res_s.subs(eta, sp.Rational(str(et))), 70))
    ok_num = ok_num and abs(ev[kslow] - s_ex) < mp.mpf(10) ** -45 and abs(wslow - r_ex) < mp.mpf(10) ** -45
    # weight of the other (fast, j-coupled) mode
check("Rem.1 60-digit eigen-decomposition: slowest rate and its current weight equal -slow root and its residue (eta=0.2,0.05,0.01)", ok_num)
# time domain: int t Phi / int Phi = 1 + eta^-2 (Lemma 1) by quadrature
et = mp.mpf("0.3")
gs = -mp.mpf(sp.N(slow.subs(eta, sp.Rational(3, 10)), 70)); gf = -mp.mpf(sp.N(fast.subs(eta, sp.Rational(3, 10)), 70))
ws = mp.mpf(sp.N(res_s.subs(eta, sp.Rational(3, 10)), 70)); wf = mp.mpf(sp.N(res_f.subs(eta, sp.Rational(3, 10)), 70))
Phi = lambda t: ws * mp.e ** (-gs * t) + wf * mp.e ** (-gf * t)
I0 = mp.quad(Phi, [0, 1, 10, 100, mp.inf]); I1 = mp.quad(lambda t: t * Phi(t), [0, 1, 10, 100, mp.inf])
check("Lemma 1 time-domain form at eta=0.3: int t Phi / int Phi = 1 + eta^-2 (quadrature, 60 digits)",
      abs(I1 / I0 - (1 + et ** -2)) < mp.mpf(10) ** -40 and abs(I0 - 1) < mp.mpf(10) ** -40, mp.nstr(I1 / I0, 20))
# Consequence: eta = (tau-1)^(-1/2) lies in (0,1) for every tau in (2, oo) and reproduces tau
tt = sp.symbols("tau", positive=True)
eta_of_tau = (tt - 1) ** sp.Rational(-1, 2)
check("consequence: tau_mem(eta(tau)) = tau and eta(tau) in (0,1) for tau > 2",
      sp.simplify(1 + eta_of_tau ** -2 - tt) == 0 and all(0 < float(eta_of_tau.subs(tt, v)) < 1 for v in (2.0001, 3, 10, 1e6)))

# ------------------------------------------------------------------ Remark 3
u0 = H0[:, 1]
check("Rem.3 h_1^T C h_2 = eta, h_1^T C h_1 = eta^2", sp.expand((u0.T * C * j0)[0] / 8 - eta) == 0 and sp.expand((u0.T * C * u0)[0] / 8 - eta ** 2) == 0)
# the obstruction inequality |u^T C j| <= R u^T C u on random cones (exact rationals)
import random
random.seed(7)
okR = True
for trial in range(200):
    m = random.randint(3, 9)
    Z = [sp.Matrix([sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for _ in range(5)]) for _ in range(m)]
    g = [sp.Rational(random.randint(0, 9), random.randint(1, 5)) for _ in range(m)]
    uu = sp.Matrix([sp.Rational(random.randint(-3, 3), 1) for _ in range(5)])
    jj = sp.Matrix([sp.Rational(random.randint(-3, 3), 1) for _ in range(5)])
    Cg = sp.zeros(5, 5)
    for gg, z in zip(g, Z):
        Cg += gg * z * z.T
    rat = [abs((jj.T * z)[0]) / abs((uu.T * z)[0]) for z in Z if (uu.T * z)[0] != 0]
    R = max(rat) if rat else 0
    lhs = abs((uu.T * Cg * jj)[0]); rhs = R * (uu.T * Cg * uu)[0]
    okR = okR and lhs <= rhs
check("Rem.3 |u^T C j| <= R u^T C u on 200 random rational cones", okR)

print()
print(f"{sum(res)} of {len(res)} checks passed")
