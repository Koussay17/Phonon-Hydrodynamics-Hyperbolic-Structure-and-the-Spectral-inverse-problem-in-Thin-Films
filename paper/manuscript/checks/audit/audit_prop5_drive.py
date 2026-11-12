"""Checks of Proposition 5 (per-event prior, F > 0, converse) and of Eq. (6) (zero-wave-vector drive).

Prop. 5: exact rational networks; Gamma_min by exact real-root isolation of the characteristic polynomial
of C on H (sympy RootOf / interval refinement), so the inequalities are certified, not floating point.
Eq. (6): a random time-reversal-symmetric system; the H-projected plane-wave problem at wave vector k
with an external energy source maintaining theta exp(ikx+pt); the current response converges to
j^T (p + C)^{-1} j as k -> 0 with an O(k^2) error (the O(k) term vanishes by time reversal).

Usage: python -B audit_prop5_drive.py > audit_prop5_drive.out.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import itertools
import random
import numpy as np
import sympy as sp

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


lam = sp.Symbol("lam")
energies = [1, 2, 3, 4, 5, 6, 7]
nm = len(energies)
events = []
for p_, a_, b_ in itertools.product(range(nm), repeat=3):
    if a_ <= b_ and p_ not in (a_, b_) and energies[p_] == energies[a_] + energies[b_]:
        s = [0] * nm; s[p_] -= 1; s[a_] += 1; s[b_] += 1; events.append(s)
m = len(events)
W = sp.diag(*[sp.Rational(k + 3, 2 * k + 1) for k in range(nm)])
eps = sp.Matrix(energies)
e = W * eps
e2 = (e.T * e)[0]
Pe = sp.eye(nm) - e * e.T / e2
zeta = [W.inv() * sp.Matrix(s) for s in events]


def Cop(gv):
    C = sp.zeros(nm, nm)
    for gg, zz in zip(gv, zeta):
        C += gg * zz * zz.T
    return C


def gamma_min_interval(C):
    """Certified enclosure [lo, hi] of the smallest eigenvalue of C on H (C e = 0, rank n-1)."""
    poly = sp.Poly(sp.cancel(C.charpoly(lam).as_expr() / lam), lam)
    rts = poly.intervals(eps=sp.Rational(1, 10 ** 30))
    lo, hi = min(rts, key=lambda iv: iv[0][0])[0]
    return lo, hi


def tau_of(C, j):
    x = Pe * (C + e * e.T / e2).LUsolve(j)
    return (x.T * x)[0] / (j.T * x)[0]


random.seed(5)
jr = sp.Matrix([random.randint(-5, 5) for _ in range(nm)])
j = Pe * jr
gref = [sp.Rational(random.randint(1, 9), random.randint(1, 4)) if k % 4 else sp.Integer(0) for k in range(m)]  # some zero reference rates
Cref = Cop(gref)
check("g_ref admissible (with zero reference rates on some events)", Cref.rank() == nm - 1)
lo, hi = gamma_min_interval(Cref)
check("Gamma_min(C_ref) enclosed to 1e-30 by exact root isolation", hi - lo < sp.Rational(1, 10 ** 29) and lo > 0, f"{float(lo):.15g}")
ok_bound, ok_adm, nF = True, True, 0
for F in (sp.Rational(1, 3), sp.Rational(1, 2), sp.Integer(1), sp.Integer(3), sp.Integer(10)):  # includes F < 1
    for trial in range(6):
        g = [gr / F + (sp.Rational(random.randint(0, 6), random.randint(1, 3)) if random.random() < 0.5 else sp.Integer(0)) for gr in gref]
        C = Cop(g)
        ok_adm = ok_adm and C.rank() == nm - 1
        t = tau_of(C, j)
        ok_bound = ok_bound and t <= F / hi                    # certified: tau <= F/hi <= F/Gamma_min
        nF += 1
check(f"Prop. 5: every g >= g_ref/F is admissible ({nF} cases, F in {{1/3,1/2,1,3,10}})", ok_adm)
check("Prop. 5: tau(g) <= F / Gamma_min(C_ref), certified (F < 1 included)", ok_bound)
# converse: rho = min over supp(g_ref) of g/g_ref; tau(g) <= 1/(rho Gamma_min); rho = 0 makes the converse trivial
ok_conv, n0 = True, 0
for trial in range(40):
    g = [sp.Rational(random.randint(0, 30), random.randint(1, 5)) for _ in range(m)]
    C = Cop(g)
    if C.rank() != nm - 1:
        continue
    rho = min(g[k] / gref[k] for k in range(m) if gref[k] != 0)
    t = tau_of(C, j)
    if rho == 0:
        n0 += 1
        continue
    # tau(g) > tau  =>  rho < 1/(tau Gamma_min): equivalent to tau(g) <= 1/(rho Gamma_min) for all g with rho > 0
    ok_conv = ok_conv and t <= 1 / (rho * hi)
check("Prop. 5 converse: tau(g) <= 1/(rho Gamma_min(C_ref)) with rho = min over supp(g_ref) of g/g_ref", ok_conv,
      f"{n0} admissible cases with rho = 0 (converse trivially true)")
# the ratio g/g_ref is undefined where g_ref = 0: the minimum must be restricted to supp(g_ref)
check("definitional gap: g_ref has zero entries, so min over all events of g/g_ref is undefined", any(gr == 0 for gr in gref))

# ------------------------------------------------------------------ Eq. (6): k -> 0 limit of the plane-wave problem
rng = np.random.default_rng(1)
half = 6
n = 2 * half
# time-reversal pairs (nu, nu+half): equal energies/weights, opposite velocities
ev = rng.uniform(0.5, 2.0, half); ev = np.concatenate([ev, ev])
vx = rng.uniform(-1, 1, half); vx = np.concatenate([vx, -vx])
Pm = np.zeros((n, n))
for i in range(n):
    Pm[(i + half) % n, i] = 1
X = rng.standard_normal((n, n)); X = X @ X.T
S = X + Pm @ X @ Pm.T
PiH = np.eye(n) - np.outer(ev, ev) / (ev @ ev)
C = PiH @ S @ PiH
jx = vx * ev   # V_x e
Vx = np.diag(vx)
check("drive setup: C symmetric PSD, C e = 0, P C P = C, j = V_x e in H",
      np.allclose(C, C.T) and np.min(np.linalg.eigvalsh(C)) > -1e-12 and np.allclose(C @ ev, 0) and np.allclose(Pm @ C @ Pm.T, C)
      and abs(jx @ ev) < 1e-12)
p = 0.37 + 0.81j
Hb = np.linalg.svd(PiH)[0][:, :n - 1]  # orthonormal basis of H
Kp = jx @ Hb @ np.linalg.solve(p * np.eye(n - 1) + Hb.T @ C @ Hb, Hb.T @ jx)
errs = []
for kk in (1e-1, 1e-2, 1e-3):
    # unknowns: y in H (coords c), external source amplitude s; theta = 1 prescribed
    # full equation: (p + i k V_x + C)(e theta + Hb c) = s e
    Mfull = p * np.eye(n) + 1j * kk * Vx + C
    A = np.column_stack([Mfull @ Hb, -ev])           # n x n system for (c, s)
    rhs = -Mfull @ ev
    sol = np.linalg.solve(A, rhs)
    c, sext = sol[:n - 1], sol[n - 1]
    y = Hb @ c
    q = jx @ y                                       # heat flux (j orthogonal to e)
    lam_k = -q / (1j * kk)                           # flux / (-gradient), gradient = i k theta
    errs.append(abs(lam_k - Kp))
    # e-projection is the energy balance with source: |e|^2 s = p |e|^2 theta + i k j^T y
    eb = abs(sext * (ev @ ev) - (p * (ev @ ev) + 1j * kk * q))
order = np.log10(errs[0] / errs[1]), np.log10(errs[1] / errs[2])
check("Eq. (6): -q/(grad T) -> j^T (p + C)^{-1} j as k -> 0, error O(k^2) (time reversal kills the O(k) term)",
      errs[-1] < 1e-5 and all(1.9 < o < 2.1 for o in order), f"errors {errs[0]:.2e}, {errs[1]:.2e}, {errs[2]:.2e}; orders {order[0]:.3f}, {order[1]:.3f}")
check("Eq. (6): the e-projection is the energy balance p rho c theta + i k q = external source", eb < 1e-10)
# without the external source the gradient cannot be prescribed at fixed p: the free problem has the diffusive pole
check("Eq. (6) remark: external source needed (s_ext != 0 at k -> 0 and p != 0)", abs(sext) > 1e-3, f"|s_ext| = {abs(sext):.3e}")
# time domain: q(t) = -G int_0^t Phi for a step gradient G (Eq. 7), from integrating Eq. (6)
from scipy.integrate import solve_ivp
G = 0.7
Ch = Hb.T @ C @ Hb; jh = Hb.T @ jx
sol = solve_ivp(lambda t, yy: -Ch @ yy - jh * G, (0, 3.0), np.zeros(n - 1), rtol=1e-11, atol=1e-13, dense_output=True)
ev_, U_ = np.linalg.eigh(Ch)
Phi_int = lambda t: np.sum((U_.T @ jh) ** 2 * (1 - np.exp(-ev_ * t)) / ev_)
okt = all(abs(jh @ sol.sol(t) + G * Phi_int(t)) < 1e-8 for t in (0.1, 1.0, 3.0))
check("Eq. (7): q(t) = -G int_0^t Phi(s) ds for a step gradient (sign and kernel)", okt)

print()
print(f"{sum(res)} of {len(res)} checks passed")
