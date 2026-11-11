"""B, regime (i): low-frequency Puiseux expansion of Z - Z_s in q = sqrt(p).

1. Symbolic: L_k(theta), k = 0..5, by series in q (SymPy), with the structural claims
     L0 = R_f - C_f/b_s^2,  L1 = -(xi1^2/b_s)(1 - b^2/b_s^2),  L2 = R_f delta + (xi1^3/b_s) K(beta),
     tau_l first enters at q^4 with coefficient -R_f delta, i.e. dL4/dtau_l = -R_f delta,
   R_f = xi1/b, C_f = xi1 b, delta = tau_R - tau_l, beta = b_s/b, K(beta) = 4/(3 beta) - 1/beta^3 - beta/3.
2. Coefficient Jacobians: det d(L0,L1,L2,L4)/d(xi1,b,tau_R,tau_l) and, at b = b_s,
   det d(L0,L2,L3,L4)/d(b,xi1,tau_R,tau_l)  (rank-four certificates).
3. Numerical: mpmath (50 digits) remainder |Z - Z_s - sum_{k<K} L_k q^k| vs omega -> slope K/2.
4. Explicit inversion of (L0, L1, L2, L4) -> theta, tested on random theta.
"""
import json
import os
import sys
import time

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "b2_lowfreq_series.json")

t0 = time.time()
q, x, b, s, tr, tl = sp.symbols("q xi1 b b_s tau_R tau_l", positive=True)
G = (1 + tr * q**2) / (1 + tl * q**2)
g = sp.sqrt(G)
r = s * g / b
T = sp.tanh(x * q * g)
# Z - Z_s = (1/(s q)) (r - 1/r) T / (1 + T/r)
expr = (r - 1 / r) * T / (1 + T / r) / (s * q)
NMAX = 6
ser = sp.series(expr, q, 0, NMAX).removeO()
ser = sp.expand(ser)
L = [sp.simplify(sp.factor(ser.coeff(q, k))) for k in range(NMAX)]
print("series time %.1f s" % (time.time() - t0))

Rf, Cf, delta, beta = x / b, x * b, tr - tl, s / b
K = sp.Rational(4, 3) / beta - 1 / beta**3 - beta / 3
claims = {
    "L0": sp.simplify(L[0] - (Rf - Cf / s**2)),
    "L1": sp.simplify(L[1] - (-(x**2 / s) * (1 - b**2 / s**2))),
    "L2": sp.simplify(L[2] - (Rf * delta + (x**3 / s) * K)),
    # L4 = -R_f tau_l delta + F(xi1, b, delta): (d/dtR + d/dtl)(L4 + R_f tau_l delta) must vanish
    "L4_plus_Rf_tl_delta_depends_on_delta_only": sp.simplify(sp.diff(L[4] + Rf * tl * delta, tr) + sp.diff(L[4] + Rf * tl * delta, tl)),
    # L3 = -(xi1^2/b_s) delta + (xi1^4/b_s) c4(beta)  (tau-dependence through delta, linear)
    "L3_tau_part": sp.simplify(sp.diff(L[3], tr) - (-(x**2 / s))),
    "L3_minus_claim_at_b_eq_bs": sp.simplify(L[3].subs(b, s) - (-(x**2 / s) * (tr - tl))),
    "L2_at_b_eq_bs": sp.simplify(L[2].subs(b, s) - x * (tr - tl) / s),
}
for k, v in claims.items():
    print(k, "residual:", v)

# tau_l-dependence of the lower coefficients (must be through delta only for L2, none for L0, L1)
dep = {
    "dL0_dtr": sp.simplify(sp.diff(L[0], tr)), "dL1_dtr": sp.simplify(sp.diff(L[1], tr)),
    "dL2_dtr+dL2_dtl": sp.simplify(sp.diff(L[2], tr) + sp.diff(L[2], tl)),
    "dL3_dtr+dL3_dtl": sp.simplify(sp.diff(L[3], tr) + sp.diff(L[3], tl)),
}
for k, v in dep.items():
    print(k, "=", v)

# Coefficient Jacobians
vars4 = [x, b, tr, tl]
J1 = sp.Matrix([[sp.diff(Lk, v) for v in vars4] for Lk in (L[0], L[1], L[2], L[4])])
det1 = sp.factor(sp.simplify(J1.det()))
print("det d(L0,L1,L2,L4)/d(xi1,b,tR,tl) =", det1)
claim_det1 = (-(2 * x**2 / (s * b**2)) * (1 - b**2 / s**2)) * (x / b) * (-(x / b) * (tr - tl))
print("  ratio to hand formula:", sp.simplify(det1 / claim_det1))

Lb = [sp.simplify(Lk.subs(b, s)) for Lk in L]
# Jacobian at b = b_s: derivatives must be taken before substitution
J2 = sp.Matrix([[sp.diff(Lk, v).subs(b, s) for v in [b, x, tr, tl]] for Lk in (L[0], L[2], L[3], L[4])])
det2 = sp.factor(sp.simplify(J2.det()))
print("det d(L0,L2,L3,L4)/d(b,xi1,tR,tl) at b=b_s =", det2)

# At delta = 0 (tl = tr): rank of d(L0..L5)/d(xi1,b,tR,tl)
J3 = sp.Matrix([[sp.diff(Lk, v).subs(tl, tr) for v in vars4] for Lk in L])
print("rank at tl = tr (generic b):", J3.rank(simplify=True))
J4 = sp.Matrix([[sp.diff(Lk, v).subs(tl, tr).subs(b, s) for v in vars4] for Lk in L])
print("rank at tl = tr, b = b_s:", J4.rank(simplify=True))

# ---------------------------------------------------------------------------------------------
# Numerical check of the series with mpmath (50 digits)
# ---------------------------------------------------------------------------------------------
mp.mp.dps = 50
fL = [sp.lambdify((x, b, s, tr, tl), Lk, "mpmath") for Lk in L]


def Zm(w, X, Bb, S, TR, TL):
    pp = mp.mpc(0, w)
    Gm = (1 + TR * pp) / (1 + TL * pp)
    gm = mp.sqrt(Gm)
    sq = mp.sqrt(pp)
    A = X * sq * gm
    B = Bb * sq / gm
    Sm = S * sq
    t = mp.tanh(A)
    return (B + Sm * t) / (B * (Sm + B * t)), 1 / (S * sq), sq


cases = {
    # nondimensional: xi1 = 1, times in units of xi1^2, b in units of b_s
    "generic_GK": (mp.mpf(1), mp.mpf("2.13"), mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.004")),
    "cattaneo": (mp.mpf(1), mp.mpf("2.13"), mp.mpf(1), mp.mpf("0.01"), mp.mpf(0)),
    "matched_b": (mp.mpf(1), mp.mpf(1), mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.004")),
    "large_tau": (mp.mpf(1), mp.mpf("0.6"), mp.mpf(1), mp.mpf("3.0"), mp.mpf("0.5")),
}
num = {}
for name, (X, Bb, S, TR, TL) in cases.items():
    Lv = [f(X, Bb, S, TR, TL) for f in fL]
    om = [mp.mpf(10) ** (-k) for k in range(3, 11)]
    rows = []
    for w in om:
        Z, Zs, sq = Zm(w, X, Bb, S, TR, TL)
        D = Z - Zs
        rem = []
        acc = mp.mpc(0)
        for k in range(NMAX):
            rem.append(abs(D - acc))
            acc += Lv[k] * sq**k
        rem.append(abs(D - acc))
        rows.append([float(w)] + [float(v) for v in rem])
    # slopes d ln rem_K / d ln omega between the two smallest omegas
    slopes = []
    for K in range(NMAX + 1):
        r1, r2 = rows[-2][1 + K], rows[-1][1 + K]
        w1, w2 = rows[-2][0], rows[-1][0]
        slopes.append(float(mp.log(r2 / r1) / mp.log(w2 / w1)) if r1 > 0 and r2 > 0 else None)
    num[name] = {"L": [complex(v) if isinstance(v, mp.mpc) else float(v) for v in Lv],
                 "rows_omega_then_remainders_K0..K6": rows,
                 "observed_slopes_K0..K6": slopes,
                 "expected_slopes_K0..K6": [-0.5] + [k / 2 for k in range(NMAX)]}
    print(name, "slopes:", ["%.3f" % v if v is not None else None for v in slopes])

# ---------------------------------------------------------------------------------------------
# Explicit inversion (L0, L1, L2, L4) -> (xi1, b, tau_R, tau_l), b != b_s, delta != 0
# ---------------------------------------------------------------------------------------------
L4_free = sp.simplify(L[4] + Rf * delta * tl)       # tau_l-independent part? check
print("d(L4 + R_f delta tl)/dtl =", sp.simplify(sp.diff(L4_free, tl)))
L4_free_f = sp.lambdify((x, b, s, tr, tl), L4_free, "mpmath")


def invert(Lv, S):
    L0, L1, L2, L4 = Lv[0], Lv[1], Lv[2], Lv[4]
    Cf_ = -S * L1 / L0
    Rf_ = L0 + Cf_ / S**2
    X = mp.sqrt(Rf_ * Cf_)
    Bb = mp.sqrt(Cf_ / Rf_)
    bet = S / Bb
    Kb = mp.mpf(4) / 3 / bet - 1 / bet**3 - bet / 3
    d = (L2 - X**3 / S * Kb) / Rf_
    # L4 = L4_free(xi1, b, s, tR, tl) - R_f delta tl, with L4_free depending on (tR, tl) through delta only
    # solve for tl with tR = tl + d (L4_free evaluated at tl=0, tR=d since it depends only on delta)
    l4f = L4_free_f(X, Bb, S, d, mp.mpf(0))
    TL = -(L4 - l4f) / (Rf_ * d)
    return X, Bb, TL + d, TL


import random
random.seed(20261111)
inv = []
for i in range(20):
    X = mp.mpf(random.uniform(0.5, 2.0))
    Bb = mp.mpf(random.choice([random.uniform(0.2, 0.9), random.uniform(1.1, 5.0)]))
    S = mp.mpf(1)
    TL = mp.mpf(random.uniform(0.0, 0.5))
    TR = TL + mp.mpf(random.choice([-1, 1])) * mp.mpf(random.uniform(0.01, 0.5))
    if TR < 0:
        TR = -TR
    Lv = [f(X, Bb, S, TR, TL) for f in fL]
    Xr, Br, TRr, TLr = invert(Lv, S)
    err = max(abs(Xr - X) / X, abs(Br - Bb) / Bb, abs(TRr - TR) / max(TR, 1e-30), abs(TLr - TL) / max(TL, 1e-3))
    inv.append(float(err))
print("explicit inversion, max relative error over 20 random theta:", max(inv))

out = {
    "coefficients_latex": [sp.latex(Lk) for Lk in L],
    "coefficients_str": [str(Lk) for Lk in L],
    "claims_residuals": {k: str(v) for k, v in claims.items()},
    "dependence_checks": {k: str(v) for k, v in dep.items()},
    "det_L0L1L2L4": str(det1),
    "det_L0L2L3L4_at_b_eq_bs": str(det2),
    "rank_at_delta0": int(J3.rank(simplify=True)),
    "rank_at_delta0_b_eq_bs": int(J4.rank(simplify=True)),
    "numerical": num,
    "inversion_max_rel_error": max(inv),
    "runtime_s": time.time() - t0,
}
with open(OUT, "w") as fh:
    json.dump(out, fh, indent=1, default=str)
print("written", OUT, "in %.1f s" % (time.time() - t0))
