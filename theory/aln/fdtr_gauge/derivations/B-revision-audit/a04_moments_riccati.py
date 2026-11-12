"""Audit of Proposition LF and Theorem T4(iv) by an independent route: Laurent coefficients of the film block from
the Riccati equation of the temperature impedance (no closed-form quadrupole, no series of cosh/sinh).

In Liouville depth x = z/sqrt(a) (0 <= x <= xi), with w_T = T/q:   w_T' = b p w_T^2 - g/b,   w_T(xi) = L/N,
L = R + 1/(bs s), N = 1 + tau_n p, p = s^2; the film block is Phi = N0 w_T(0), N0 = 1 + tau_o p.
(From dT/dx = -(g/b) q and dq/dx = -b p T, which follow from q = -lam_eff dT/dz, C p T + dq/dz = 0.)
Write w_T = 1/(bs s) + sum_k v_k(x) s^k; then
   v_k' = [b/bs^2 - g/b]_k + (2 b/bs) v_{k-1} + b sum_{i+j=k-2} v_i v_j,
   v_k(xi) = [R/N + (1/N - 1)/(bs s)]_k,
and c_k(Phi) = v_k(0) + tau_o v_{k-2}(0) + delta_{k1} tau_o/bs.
Checks:
 (1) c_{-1}, c_0, c_1, c_2|_{R=0} of Proposition LF for general tau_n, tau_o;
 (2) the stack (tau_o = tau_n = eps tau_l, W = R1 + Phi): G'_12345 = det[d c_k/d(xi, b, tR, tl, R2)]_{k=1..5}, b_s = 1,
     local (eps = 0) and entropic (eps = 1), against the polynomials P, Q of b2_03_reduce_record.md;
 (3) t(lam, mu) = det of the 6x6 frequency matrix for K = {0..5}, against the record;
 (4) the restrictions G'|_{local, R2=0} and G'|_{entropic, Delta, R2=0}.
Output: a04_moments_riccati.out
"""
import os
import time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


x, s = sp.symbols("x s")
xi, b, bs = sp.symbols("xi b b_s", positive=True)
tR, tl, R = sp.symbols("tau_R tau_l R", nonnegative=True)
tn, to, eps = sp.symbols("tau_n tau_o epsilon", real=True)
KMAX = 6


def gcoef(k):
    """s^k coefficient of g = (1 + tR s^2)/(1 + tl s^2)."""
    if k % 2:
        return sp.Integer(0)
    m = k // 2
    return sp.Integer(1) if m == 0 else (tR - tl) * (-tl) ** (m - 1)


def bcoef(k, tau_n, Rv, bsv):
    """s^k coefficient of R/N + (1/N - 1)/(bs s), N = 1 + tau_n s^2."""
    if k % 2 == 0:
        return Rv * (-tau_n) ** (k // 2)
    j = (k + 1) // 2
    return (-tau_n) ** j / bsv


def moments(tau_n, tau_o, Rv=R, bsv=bs, kmax=KMAX):
    v = []
    for k in range(kmax + 1):
        rhs = (b / bsv ** 2 if k == 0 else 0) - gcoef(k) / b
        if k >= 1:
            rhs += 2 * b / bsv * v[k - 1]
        rhs += b * sum(v[i] * v[k - 2 - i] for i in range(0, k - 1)) if k >= 2 else 0
        rhs = sp.expand(rhs)
        prim = sp.integrate(rhs, x)
        vk = sp.expand(bcoef(k, tau_n, Rv, bsv) - (prim.subs(x, xi) - prim))
        v.append(vk)
    v0 = [sp.expand(vk.subs(x, 0)) for vk in v]
    c = {}
    for k in range(kmax + 1):
        ck = v0[k] + (tau_o * v0[k - 2] if k >= 2 else 0) + (tau_o / bsv if k == 1 else 0)
        c[k] = sp.expand(ck)
    return c


t0 = time.time()
# ------------------------------------------------------------------ (1) Proposition LF
c = moments(tn, to, kmax=2)
LF = {0: R + xi * (bs ** 2 - b ** 2) / (b * bs ** 2),
      1: xi ** 2 * (b ** 2 - bs ** 2) / bs ** 3 - 2 * b * R * xi / bs - (tn - to) / bs,
      2: xi / b * ((tR - tl) + 2 * b ** 2 / bs ** 2 * tn - (b ** 2 - bs ** 2) / bs ** 2 * to)
      - xi ** 3 * (b ** 2 - bs ** 2) * (3 * b ** 2 - bs ** 2) / (3 * b * bs ** 4)}
log("(1) c_0 - LF:", sp.simplify(c[0] - LF[0]))
log("(1) c_1 - LF:", sp.simplify(c[1] - LF[1]))
log("(1) c_2|_{R=0} - LF:", sp.simplify((c[2] - LF[2]).subs(R, 0)))
log("(1) c_2 with R (for the record):", sp.factor(c[2] - c[2].subs(R, 0)))
# delta_eps of LF(b): tau_o = tau_n = eps tau_l
d2 = sp.expand(c[2].subs({R: 0, tn: eps * tl, to: eps * tl}))
time_part = sp.expand(d2 - d2.subs({tR: 0, tl: 0}))
log("(1b) time part of c_2 (stack, R=0) * b/xi:", sp.factor(sp.expand(time_part * b / xi)),
    "  vs tR - tl + eps (1 + b^2/bs^2) tl:", sp.simplify(time_part * b / xi - (tR - tl + eps * (1 + b ** 2 / bs ** 2) * tl)))
log("    time", round(time.time() - t0, 1), "s")

# ------------------------------------------------------------------ (2) G'_12345 for the stack, b_s = 1
R2 = sp.Symbol("R2", nonnegative=True)
P_rec = sp.sympify("90*b*R*tl**2 - 90*b*R*tl*tR + 15*b**4*R**4*xi - 30*b**2*R**2*tl*xi - 45*b**2*tl**2*xi + 30*b**2*R**2*tR*xi"
                   " + 45*tl*tR*xi + 45*b**2*tl*tR*xi - 45*tR**2*xi + 24*b**3*R**3*xi**2 - 30*b**5*R**3*xi**2 - 6*b*R*tl*xi**2"
                   " + 60*b**3*R*tl*xi**2 + 6*b*R*tR*xi**2 - 60*b**3*R*tR*xi**2 + 11*b**2*R**2*xi**3 - 31*b**4*R**2*xi**3"
                   " - 3*tl*xi**3 + 3*b**2*tl*xi**3 + 3*tR*xi**3 - 3*b**2*tR*xi**3 + 2*b*R*xi**4 - 4*b**3*R*xi**4 + 2*b**5*R*xi**4",
                   locals={"b": b, "R": R2, "tl": tl, "tR": tR, "xi": xi})
Q_rec = sp.sympify("-1350*b*R*tl**2*tR + 1350*b*R*tl*tR**2 - 45*b**4*R**4*tl*xi + 630*b**4*R**2*tl**2*xi + 675*b**4*tl**3*xi"
                   " + 225*b**4*R**4*tR*xi - 450*b**2*R**2*tl*tR*xi - 450*b**4*R**2*tl*tR*xi - 675*tl**2*tR*xi"
                   " + 2025*b**2*tl**2*tR*xi - 675*b**4*tl**2*tR*xi + 450*b**2*R**2*tR**2*xi + 1350*tl*tR**2*xi"
                   " - 2025*b**2*tl*tR**2*xi - 675*tR**3*xi + 90*b**5*R**5*xi**2 - 90*b**3*R**3*tl*xi**2"
                   " + 540*b**5*R**3*tl*xi**2 + 630*b**3*R*tl**2*xi**2 + 450*b**5*R*tl**2*xi**2 + 360*b**3*R**3*tR*xi**2"
                   " - 450*b**5*R**3*tR*xi**2 - 90*b*R*tl*tR*xi**2 + 540*b**3*R*tl*tR*xi**2 - 450*b**5*R*tl*tR*xi**2"
                   " - 900*b**3*R*tR**2*xi**2 + 150*b**4*R**4*xi**3 - 15*b**2*R**2*tl*xi**3 + 675*b**4*R**2*tl*xi**3"
                   " + 45*b**2*tl**2*xi**3 + 45*b**4*tl**2*xi**3 + 135*b**2*R**2*tR*xi**3 - 435*b**4*R**2*tR*xi**3"
                   " - 45*tl*tR*xi**3 + 180*b**2*tl*tR*xi**3 + 45*b**4*tl*tR*xi**3 + 45*tR**2*xi**3 + 45*b**2*tR**2*xi**3"
                   " + 90*b**3*R**3*xi**4 + 30*b**5*R**3*xi**4 + 150*b**3*R*tl*xi**4 - 30*b**5*R*tl*xi**4 + 30*b*R*tR*xi**4"
                   " + 30*b**3*R*tR*xi**4 + 60*b**5*R*tR*xi**4 + 24*b**2*R**2*xi**5 + 24*b**4*R**2*xi**5 + 2*b*R*xi**6"
                   " - 2*b**5*R*xi**6", locals={"b": b, "R": R2, "tl": tl, "tR": tR, "xi": xi})
VARS = [xi, b, tR, tl, R2]
for label, epsv, claim in (("local", 0, 2 * (tl - tR) * xi ** 7 * P_rec / (135 * b ** 2)),
                           ("entropic", 1, -2 * xi ** 7 * Q_rec / (2025 * b ** 2))):
    t1 = time.time()
    cc = moments(epsv * tl, epsv * tl, Rv=R2, bsv=1, kmax=5)
    Gm = sp.Matrix([[sp.diff(cc[k], v) for v in VARS] for k in range(1, 6)])
    G = sp.expand(Gm.det(method="berkowitz"))
    diff = sp.expand(G - sp.expand(claim))
    log(f"(2) {label}: G'_12345 - record = {diff}   ({round(time.time() - t1, 1)} s, {len(G.args)} terms)")
    if label == "local":
        log("(4) local, R2 = 0: G' =", sp.factor(G.subs(R2, 0)))
    else:
        log("(4) entropic, Delta, R2 = 0: G' =", sp.factor(G.subs({R2: 0, tR: tl})))
        log("    entropic, R2 = 0 (off Delta):", sp.factor(G.subs(R2, 0)))
        log("    entropic, Delta (R2 free):", sp.factor(G.subs(tR, tl)))
        log("    entropic, corner tR = tl = 0:", sp.factor(G.subs({tR: 0, tl: 0})))

# ------------------------------------------------------------------ (3) t(lam, mu)
lam, mu = sp.symbols("lam mu", positive=True)
rows = []
for L_ in (1, lam, mu):
    rows.append([sp.cos(k * sp.pi / 4) * L_ ** k for k in range(6)])
    rows.append([sp.sin(k * sp.pi / 4) * L_ ** k for k in range(6)])
T6 = sp.Matrix(rows)
tdet = sp.factor(sp.simplify(T6.det(method="berkowitz")))
t_rec = (lam - 1) ** 2 * lam * (1 + lam ** 2) * (lam - mu) ** 2 * (mu - 1) ** 2 * mu * (1 + mu ** 2) * (lam ** 2 + mu ** 2) / (2 * sp.sqrt(2))
log("(3) t(lam,mu) =", tdet)
log("(3) t / t_record =", sp.simplify(tdet / t_rec))
log("total time", round(time.time() - t0, 1), "s")
with open(os.path.join(HERE, "a04_moments_riccati.out"), "w") as fh:
    fh.write("\n".join(out) + "\n")
