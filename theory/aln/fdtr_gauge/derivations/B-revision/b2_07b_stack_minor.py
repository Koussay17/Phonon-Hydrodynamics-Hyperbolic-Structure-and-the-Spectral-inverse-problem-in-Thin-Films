"""Second-system check (SymPy) of the leading three-frequency coefficient for the transducer stack (Theorem T4(iv)).

W = R1 + N Mobius(M; L2/N) (b_s = 1, R = R2), local (N = 1) and entropic (N = 1 + tl s^2).
Computes G'_12345 = det[grad c_1; ...; grad c_5] over (xi, b, tR, tl, R2) and compares with the Wolfram Language results
recorded in b2_03_reduce_record.md (local: 2 (tl - tR) xi^7 P/(135 b^2); entropic on Delta, R2 = 0: -(4/45)(3 + b^2) xi^10 t^2),
and the frequency factor t_012345(lam, mu) of the Cauchy-Binet expansion for omega_j = (lam_j sigma)^2.
Output: b2_07b_stack_minor.out
"""
import os
import sympy as sp

xi, b, tR, tl, R, s, lam, mu, t = sp.symbols("xi b tau_R tau_l R s lambda mu t")
n = 7


def trunc(e):
    P = sp.Poly(sp.expand(e), s)
    return sum(c * s ** m[0] for m, c in zip(P.monoms(), P.coeffs()) if m[0] < n)


def coeffs(entropic):
    ginv = sum((-tl * s ** 2) ** m for m in range(n // 2 + 1))
    g = trunc((1 + tR * s ** 2) * ginv)
    z = trunc(xi ** 2 * s ** 2 * g)
    C = 0
    Sz = 0
    zp = 1
    for m in range(n // 2 + 2):
        C += zp / sp.factorial(2 * m)
        Sz += zp / sp.factorial(2 * m + 1)
        zp = trunc(zp * z)
    C, Sz = trunc(C), trunc(Sz)
    X = trunc(xi / b * g * Sz)
    N = 1 + (tl * s ** 2 if entropic else 0)
    Nh = trunc(N * (C * (1 + R * s) + N * X * s))
    D1 = trunc(N * C + xi * b * s * Sz * (1 + R * s))
    nk = [sp.expand(Nh).coeff(s, k) for k in range(n)]
    dk = [sp.expand(D1).coeff(s, k) for k in range(n)]
    q = []
    for k in range(n):
        acc = nk[k] - sum(dk[j] * q[k - j] for j in range(1, k + 1))
        q.append(sp.expand(sp.cancel(acc / dk[0])))
    return q          # q[k] = c_{k-1}


lines = []
P_ref = (90 * b * R * tl ** 2 - 90 * b * R * tl * tR + 15 * b ** 4 * R ** 4 * xi - 30 * b ** 2 * R ** 2 * tl * xi
         - 45 * b ** 2 * tl ** 2 * xi + 30 * b ** 2 * R ** 2 * tR * xi + 45 * tl * tR * xi + 45 * b ** 2 * tl * tR * xi
         - 45 * tR ** 2 * xi + 24 * b ** 3 * R ** 3 * xi ** 2 - 30 * b ** 5 * R ** 3 * xi ** 2 - 6 * b * R * tl * xi ** 2
         + 60 * b ** 3 * R * tl * xi ** 2 + 6 * b * R * tR * xi ** 2 - 60 * b ** 3 * R * tR * xi ** 2 + 11 * b ** 2 * R ** 2 * xi ** 3
         - 31 * b ** 4 * R ** 2 * xi ** 3 - 3 * tl * xi ** 3 + 3 * b ** 2 * tl * xi ** 3 + 3 * tR * xi ** 3 - 3 * b ** 2 * tR * xi ** 3
         + 2 * b * R * xi ** 4 - 4 * b ** 3 * R * xi ** 4 + 2 * b ** 5 * R * xi ** 4)
for entropic in (False, True):
    q = coeffs(entropic)
    cs = q[2:7]                     # c_1 .. c_5
    J = sp.Matrix([[sp.diff(c, v) for v in (xi, b, tR, tl, R)] for c in cs])
    G = sp.factor(J.det(method="berkowitz"))
    if not entropic:
        ref = 2 * (tl - tR) * xi ** 7 * P_ref / (135 * b ** 2)
        lines.append(f"local: G'_12345 - Wolfram = {sp.simplify(G - ref)}")
    else:
        GD = sp.factor(G.subs({tR: t, tl: t, R: 0}))
        lines.append(f"entropic, Delta, R2=0: G'_12345 = {GD}")
        lines.append(f"  minus Wolfram = {sp.simplify(GD + sp.Rational(4, 45) * (3 + b ** 2) * xi ** 10 * t ** 2)}")
rows = []
for L in (1, lam, mu):
    rows.append([sp.cos(k * sp.pi / 4) * L ** k for k in range(6)])
    rows.append([sp.sin(k * sp.pi / 4) * L ** k for k in range(6)])
tK = sp.factor(sp.simplify(sp.Matrix(rows).det()))
tref = (lam - 1) ** 2 * lam * (1 + lam ** 2) * (lam - mu) ** 2 * (mu - 1) ** 2 * mu * (1 + mu ** 2) * (lam ** 2 + mu ** 2) / (2 * sp.sqrt(2))
lines.append(f"t_012345(lam, mu) = {tK}")
lines.append(f"  minus Wolfram = {sp.simplify(tK - tref)}")
out = "\n".join(lines)
print(out)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_07b_stack_minor.out"), "w") as fh:
    fh.write(out + "\n")
