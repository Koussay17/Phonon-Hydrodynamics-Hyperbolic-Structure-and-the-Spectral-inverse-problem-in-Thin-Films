"""Independent checks of Theorem 2 (adj(A) via the matrix determinant lemma), the extension to
non-admissible rates, Proposition 3 (det A = Theta(eps), K0 ~ c/eps), Remark 5 (current parallel to an
event vector), Proposition 4 (sub-network form) and the LMI identity (27).

Everything is computed in mode space R^n with rational arithmetic (no orthonormal basis, no square roots):
for y in H, det_H[zeta_I, y] in an orthonormal basis of H equals det_{R^n}[zeta_I, y, e]/|e|, so
psi~_I := (det[zeta_I, delta_r, e])_r = |e| psi_I lies in H, and Eq. (24) becomes
    det_H(A) C^+ = |e|^{-2} sum_I g_I psi~_I psi~_I^T,   det_H(A) = det(C + e e^T/|e|^2).
Network: n = 6 modes, energies 1..6, all resonant p -> a + b events (distinct and repeated daughters),
non-uniform rational weights W, rational rates with some zeros.

Usage: python -B audit_thm2_cone.py > audit_thm2_cone.out.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools
import random
import sympy as sp

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


# ------------------------------------------------------------------ matrix determinant lemma (symbolic, 3x3)
a = sp.Matrix(3, 3, sp.symbols("a0:9"))
yv = sp.Matrix(sp.symbols("y0:3"))
zv = sp.Matrix(sp.symbols("z0:3"))
check("matrix determinant lemma det(A + y z^T) = det A + z^T adj(A) y for a general (possibly singular) 3x3 A",
      sp.expand((a + yv * zv.T).det() - a.det() - (zv.T * a.adjugate() * yv)[0]) == 0)

# ------------------------------------------------------------------ network
energies = [1, 2, 3, 4, 5, 6]
nm = len(energies)
events = []
for p_, a_, b_ in itertools.product(range(nm), repeat=3):
    if a_ <= b_ and p_ not in (a_, b_) and energies[p_] == energies[a_] + energies[b_]:
        s = [0] * nm; s[p_] -= 1; s[a_] += 1; s[b_] += 1; events.append(s)
m = len(events)
W = sp.diag(sp.Rational(7, 5), sp.Rational(2, 3), 3, sp.Rational(5, 4), sp.Rational(1, 2), 2)
eps = sp.Matrix(energies)
e = W * eps
e2 = (e.T * e)[0]
Pe = sp.eye(nm) - e * e.T / e2
zeta = [W.inv() * sp.Matrix(s) for s in events]
nH = nm - 1
check("events conserve energy (zeta in e-perp)", all((zz.T * e)[0] == 0 for zz in zeta))
check(f"Z spans H ({m} events, rank {sp.Matrix.hstack(*zeta).rank()} = n_H = {nH})", sp.Matrix.hstack(*zeta).rank() == nH)


def psi_t(I):
    cols = [zeta[k] for k in I]
    return sp.Matrix([sp.Matrix.hstack(*(cols + [sp.eye(nm)[:, r], e])).det() for r in range(nm)])


subsets = list(itertools.combinations(range(m), nH - 1))
PSI = {I: psi_t(I) for I in subsets}
check("psi~_I lies in H for every subset (psi~^T e = det[zeta_I, e, e] = 0)", all((PSI[I].T * e)[0] == 0 for I in subsets))


def Cop(gv):
    C = sp.zeros(nm, nm)
    for gg, zz in zip(gv, zeta):
        C += gg * zz * zz.T
    return C


def solve_H(C, j):
    """x = C^+ j for j in H (admissible C): (C + e e^T/|e|^2) x = j, then x is already in H."""
    x = (C + e * e.T / e2).LUsolve(j)
    return Pe * x


jr = sp.Matrix([1, -2, 0, 1, sp.Rational(1, 3), -1])
j_gen = Pe * jr
j_spec = zeta[0] + 2 * zeta[1]          # lies in the span of two events
random.seed(3)
g = [sp.Rational(random.randint(1, 9), random.randint(1, 4)) if k % 5 else 0 for k in range(m)]  # zeros at k = 0, 5
C = Cop(g)
check("g admissible despite zero rates: rank C(g) = n-1", C.rank() == nm - 1)
detH = (C + e * e.T / e2).det()
Cplus = (C + e * e.T / e2).inv() - e * e.T / e2
rhs = sp.zeros(nm, nm)
for I in subsets:
    rhs += sp.prod([g[k] for k in I]) * PSI[I] * PSI[I].T
rhs = rhs / e2
nb_dep = sum(1 for I in subsets if PSI[I].is_zero_matrix)
nb_nb = sum(1 for I in subsets if not PSI[I].is_zero_matrix and (PSI[I].T * j_spec)[0] == 0)
check("Eq. (24): det_H(A) C^+ = sum over ALL (n_H-1)-subsets g_I psi_I psi_I^T (exact)",
      (detH * Cplus - rhs).is_zero_matrix,
      f"{len(subsets)} subsets, {nb_dep} dependent (psi=0), {nb_nb} independent with psi^T j_spec = 0")
check("statement gap: for j_spec some independent subsets are not bases but have psi_I != 0 and enter adj(A)", nb_nb > 0)


def basic_solution(I, j):
    """Solve Eq. (21) directly: zeta_g^T y = 0 (g in I), e^T y = 0, j^T y = 1 (n equations)."""
    Mx = sp.Matrix.vstack(*[zeta[k].T for k in I], e.T, j.T)
    if Mx.det() == 0:
        return None
    return Mx.LUsolve(sp.Matrix([0] * (nH - 1) + [0, 1]))


for lab, jj in (("generic j", j_gen), ("j in span of two events", j_spec)):
    x = solve_H(C, jj)
    K0 = (jj.T * x)[0]
    num = sp.zeros(nm, 1); den = 0; okb = True; M2 = 0
    for I in subsets:
        pj = (PSI[I].T * jj)[0]
        yI = basic_solution(I, jj)
        okb = okb and ((pj != 0) == (yI is not None))          # basis iff psi_I^T j != 0
        if pj != 0:
            okb = okb and (PSI[I] / pj - yI).is_zero_matrix     # psi_I/(psi_I^T j) = y_I
            gI = sp.prod([g[k] for k in I])
            num += gI * pj * PSI[I]; den += gI * pj ** 2
            M2 = max(M2, (yI.T * yI)[0])
    check(f"[{lab}] basis iff psi_I^T j != 0, and then psi_I/(psi_I^T j) = y_I (solved directly from Eq. 21)", okb)
    check(f"[{lab}] Eq. (22): x/K0 = sum_I varpi_I y_I exactly", (num / den - x / K0).is_zero_matrix)
    N = (x.T * x)[0]
    tau = N / K0
    check(f"[{lab}] Eq. (23): |x|^2 <= M^2 K0^2 and tau <= M^2 K0", N <= M2 * K0 ** 2 and tau <= M2 * K0,
          f"tau = {float(tau):.6g}, M^2 K0 = {float(M2 * K0):.6g}")

# ------------------------------------------------------------------ extension: non-admissible g with j in range C(g)
sub = [0, 1, 3]
gna = [sp.Rational(k + 2, 3) if k in sub else 0 for k in range(m)]
Cna = Cop(gna)
jna = zeta[0] - zeta[1] + 3 * zeta[3]
inrange = Cna.rank() == Cna.row_join(jna).rank()
check(f"non-admissible g: rank C(g) = {Cna.rank()} < n-1 = {nm - 1}, j in range C(g)", Cna.rank() < nm - 1 and inrange)
Rm = sp.Matrix.hstack(*Cna.columnspace())
PiR = Rm * (Rm.T * Rm).inv() * Rm.T                 # projector onto range C(g)
PiQ = Pe - PiR                                       # projector onto ker C(g) cap H
check("H = range C(g) (+) (ker C(g) cap H): PiR + PiQ = Pi_H, PiQ idempotent, Cna PiQ = 0",
      (PiQ * PiQ - PiQ).is_zero_matrix and (Cna * PiQ).is_zero_matrix)
w = Cna.pinv() * jna
w = PiR * w
K0na = (jna.T * w)[0]
E = Cop([1] * m)
Qm = sp.Matrix.hstack(*[c for c in PiQ.columnspace()])
# u_inf = -(Pi_Q E Pi_Q)^{-1} Pi_Q E w, inverse taken on Q: solve (Qm^T E Qm) c = -Qm^T E w, u = Qm c
cvec = (Qm.T * E * Qm).LUsolve(-Qm.T * E * w)
uinf = Qm * cvec
rows = []
for th in (sp.Rational(1, 10 ** 2), sp.Rational(1, 10 ** 4), sp.Rational(1, 10 ** 6)):
    xth = solve_H(Cop([gg + th for gg in gna]), jna)
    wth, uth = PiR * xth, PiQ * xth
    rows.append((float(th), float(sp.sqrt(((wth - w).T * (wth - w))[0])), float(sp.sqrt(((uth - uinf).T * (uth - uinf))[0])),
                 float((jna.T * xth)[0] - K0na)))
check("extension: w_theta -> C(g)^+ j, u_theta -> -(Pi_Q E Pi_Q)^{-1} Pi_Q E w, K0(g_theta) -> K0(g), errors O(theta)",
      rows[-1][1] < 1e-4 and rows[-1][2] < 1e-4 and abs(rows[-1][3]) < 1e-4 and rows[-1][1] < rows[0][1],
      "; ".join(f"theta={t:g}: |dw|={a_:.1e}, |du|={b_:.1e}, dK={c_:.1e}" for t, a_, b_, c_ in rows))
M2na = 0
for I in subsets:
    yI = basic_solution(I, jna)
    if yI is not None:
        M2na = max(M2na, (yI.T * yI)[0])
check("extension: w orthogonal to u_inf, and |w|^2 <= |w + u_inf|^2 <= M^2 K0(g)^2",
      (w.T * uinf)[0] == 0 and (w.T * w)[0] <= ((w + uinf).T * (w + uinf))[0] <= M2na * K0na ** 2)

# ------------------------------------------------------------------ Proposition 3: det A = Theta(eps), K0 ~ c/eps
best, Istar = -1, None
for I in subsets:
    yI = basic_solution(I, j_gen)
    if yI is not None and (yI.T * yI)[0] > best:
        best, Istar = (yI.T * yI)[0], I
rows = []
for ex in (2, 6, 10, 14, 18):
    ep = sp.Rational(1, 10 ** ex)
    ge = [1 if k in Istar else ep for k in range(m)]
    Ce = Cop(ge)
    dA = (Ce + e * e.T / e2).det()
    xe = solve_H(Ce, j_gen)
    Ke = (j_gen.T * xe)[0]
    rows.append((float(ep), float(dA / ep), float(Ke * ep), float(((xe.T * xe)[0] / Ke) / (best * Ke))))
check("Prop. 3: det A/eps and K0*eps converge to positive constants; tau/(M^2 K0) -> 1",
      abs(rows[-1][1] / rows[-2][1] - 1) < 1e-3 and abs(rows[-1][2] / rows[-2][2] - 1) < 1e-3 and rows[-1][1] > 0
      and abs(rows[-1][3] - 1) < 1e-4,
      "; ".join(f"eps={a_:g}: detA/eps={b_:.6g}, K0*eps={c_:.6g}, tau/(M^2K0)={d_:.6f}" for a_, b_, c_, d_ in rows))
okb = any(sp.Matrix.hstack(*([zeta[k] for k in Istar + (gm,)] + [e])).det() != 0 for gm in range(m) if gm not in Istar)
check("Prop. 3 proof: some I* U {gamma} (gamma not in I*) is a basis of H, so det A >= eps det(Z_S)^2 > 0", okb)

# ------------------------------------------------------------------ Remark 5: current parallel to an event vector
beta = 2
jb = zeta[beta]
contain_beta_nonbasis = all(basic_solution(I, jb) is None for I in subsets if beta in I)
cand = []
for I in subsets:
    if beta not in I:
        yI = basic_solution(I, jb)
        if yI is not None:
            cand.append(((yI.T * yI)[0], I, yI))
cand.sort(key=lambda t: t[0], reverse=True)
Mb2, Jst, yJ = cand[0]
okR5 = contain_beta_nonbasis
for trial in range(4):
    gb = [sp.Rational(random.randint(1, 9), random.randint(1, 5)) if k in Jst + (beta,) else 0 for k in range(m)]
    Cb = Cop(gb)
    if Cb.rank() != nm - 1:
        okR5 = False
        break
    xb = solve_H(Cb, jb)
    Kb_ = (jb.T * xb)[0]
    okR5 = okR5 and (xb / Kb_ - yJ).is_zero_matrix and ((xb.T * xb)[0] / Kb_) == Mb2 * Kb_
check("Rem. 5: j = zeta_beta; no subset containing beta is a basis; g supported on J* U {beta} gives x/K0 = y_J*, tau = M^2 K0",
      okR5, f"M^2 = {float(Mb2):.6g}")

# ------------------------------------------------------------------ Proposition 4 on every basis (non-uniform W)
ok4, nb = True, 0
for I in subsets:
    SI = sp.Matrix.hstack(*[sp.Matrix(events[k]) for k in I])
    yI = basic_solution(I, j_gen)
    if yI is None:
        continue
    nb += 1
    ker = SI.T.nullspace()
    phi = [kv for kv in ker if sp.Matrix.hstack(kv, eps).rank() == 2]
    ok4 = ok4 and len(ker) == 2 and len(phi) > 0
    if phi:
        den = (j_gen.T * W * phi[0])[0]
        ok4 = ok4 and den != 0 and (Pe * W * phi[0] / den - yI).is_zero_matrix
check(f"Prop. 4: ker S_I^T = span(eps, phi_I) and y_I = Pi_H W phi_I / (j^T W phi_I) on all {nb} bases", ok4 and nb > 0)

# ------------------------------------------------------------------ LMI (27)
lam = sp.Symbol("lam")


def psd(Mx):
    c = sp.Poly(Mx.charpoly(lam).as_expr(), lam).all_coeffs()
    return all((-1) ** k * c[k] >= 0 for k in range(len(c)))


x_ = solve_H(C, j_gen)
K0f = (j_gen.T * x_)[0]
okL = True
for Ks in (K0f * sp.Rational(99, 100), K0f, K0f * sp.Rational(101, 100)):
    Mx = C.row_join(j_gen).col_join(j_gen.T.row_join(sp.Matrix([[Ks]])))
    okL = okL and (psd(Mx) == (K0f <= Ks))
# non-admissible point with j in range: LMI holds at the extended value K0(g)
Mx = Cna.row_join(jna).col_join(jna.T.row_join(sp.Matrix([[K0na]])))
okL = okL and psd(Mx)
check("LMI (27): [[C, j],[j^T, K*]] PSD iff K0(g) <= K* (admissible g), and PSD at K* = K0(g) for the non-admissible g", okL)

print()
print(f"{sum(res)} of {len(res)} checks passed")
