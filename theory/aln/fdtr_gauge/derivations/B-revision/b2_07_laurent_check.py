"""Second-system check (SymPy) of the low-frequency coefficients of Proposition LF and of the moment minor G_0123.

Film block Phi = N0 * Mobius(M; L/N), N0 = 1 + to s^2, N = 1 + tn s^2, L = R + 1/(bs s); s Phi = N0 Nh/D1 with
Nh = C (1 + R bs s) + N X bs s,  D1 = bs N C + xi b s Sz (1 + R bs s)  (truncated power series, as in B-analytic laurent.py).
Compares c_{-1}, c_0, c_1, c_2 (R = 0) with the closed forms obtained independently by Wolfram Language series
(derivations/B-revision/b2_03_reduce_record.md), and recomputes G_0123 for the read-out T(0) with tn = eps tl.
Output: b2_07_laurent_check.out
"""
import os
import sympy as sp

xi, b, bs, tR, tl, R, tn, to, eps, s = sp.symbols("xi b b_s tau_R tau_l R tau_n tau_o epsilon s")
n = 7


def trunc(e):
    P = sp.Poly(sp.expand(e), s)
    return sum(c * s ** m[0] for m, c in zip(P.monoms(), P.coeffs()) if m[0] < n)


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
N = 1 + tn * s ** 2
N0 = 1 + to * s ** 2
Nh = trunc(N0 * (C * (1 + R * bs * s) + N * X * bs * s))
D1 = trunc(bs * N * C + xi * b * s * Sz * (1 + R * bs * s))
nk = [sp.expand(Nh).coeff(s, k) for k in range(n)]
dk = [sp.expand(D1).coeff(s, k) for k in range(n)]
q = []
for k in range(n):
    acc = nk[k] - sum(dk[j] * q[k - j] for j in range(1, k + 1))
    q.append(sp.expand(sp.cancel(acc / dk[0])))
cm1, c0, c1, c2 = q[0], q[1], q[2], q[3]

ref = {
    "c_-1": 1 / bs,
    "c_0": R + xi * (bs ** 2 - b ** 2) / (b * bs ** 2),
    "c_1": -(bs ** 2 * tn - bs ** 2 * to + 2 * b * bs ** 2 * R * xi - b ** 2 * xi ** 2 + bs ** 2 * xi ** 2) / bs ** 3,
    "c_2(R=0)": (-(tl * xi) / b + (2 * b * tn * xi) / bs ** 2 - (b - bs) * (b + bs) * to * xi / (b * bs ** 2) + (tR * xi) / b
                 - (b - bs) * (b + bs) * (3 * b ** 2 - bs ** 2) * xi ** 3 / (3 * b * bs ** 4)),
}
got = {"c_-1": cm1, "c_0": c0, "c_1": c1, "c_2(R=0)": c2.subs(R, 0)}
lines = []
for k in ref:
    d = sp.simplify(got[k] - ref[k])
    lines.append(f"{k}: SymPy - Wolfram = {d}")

# G_0123 for the read-out T(0), tn = eps tl, R = 0, bs = 1
cs = [sp.expand(c.subs({to: 0, tn: eps * tl, R: 0, bs: 1})) for c in (q[1], q[2], q[3], q[4])]
J = sp.Matrix([[sp.diff(c, v) for v in (xi, b, tR, tl)] for c in cs])
G = sp.factor(J.det())
Gref = -2 * eps * xi ** 3 * (3 * (tR - tl) + 3 * eps * (1 + b ** 2) * tl + b ** 2 * (b ** 2 - 1) * xi ** 2) / (3 * b ** 3)
lines.append(f"G_0123 (T read-out, tied eps, R=0, bs=1) = {G}")
lines.append(f"G_0123 - Wolfram = {sp.simplify(G - Gref)}")
out = "\n".join(lines)
print(out)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_07_laurent_check.out"), "w") as fh:
    fh.write(out + "\n")
