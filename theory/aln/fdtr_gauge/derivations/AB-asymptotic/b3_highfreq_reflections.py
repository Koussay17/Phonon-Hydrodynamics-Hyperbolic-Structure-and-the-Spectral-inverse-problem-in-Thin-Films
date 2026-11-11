"""B, regimes (ii)-(iv): reflection series, high-frequency expansions (GK and Cattaneo).

Checks (all at 60-digit precision, nondimensional units xi1 = 1 unless stated, b_s = 1):
 R1. Reflection series Z = (1/B)(1 + 2 sum_n Gamma^n e^{-2nA}): |Gamma| < 1 and Re A > 0 on random samples;
     partial-sum error ~ (|Gamma| e^{-2 Re A})^{N+1}.
 G1. GK (tau_l > 0): Z b_inf sqrt(p) - 1 - (1/tau_R - 1/tau_l)/(2p) = O(p^-2) + exp. small (slope -2).
 G2. GK reflection: (B (Z - 1/B)/2) / (Gamma e^{-2A}) -> 1, and
     ln|B (Z - 1/B)/2| = ln|Gamma_inf| - xi_inf sqrt(2 omega) + O(omega^-1/2),  xi_inf = xi1 sqrt(tau_R/tau_l).
 G3. Symbolic: the large-p coefficients of ln(g) fix (tau_R, tau_l) iff tau_R != tau_l (Jacobian determinant).
 C1. Cattaneo (tau_l = 0): Z = (sqrt(tau_R)/b) tanh(t_e/(2 tau_R) + i omega t_e) + O(omega^-1/2),  t_e = xi1 sqrt(tau_R).
 C2. Cattaneo reconstruction of (xi1, b, tau_R) from the almost-periodic high-frequency signal
     (period pi/t_e in omega, mean sqrt(tau_R)/b, first harmonic 2 e^{-t_e/tau_R} sqrt(tau_R)/b).
"""
import json
import os
import random

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "b3_highfreq_reflections.json")
mp.mp.dps = 60
res = {}


def parts(w, X, Bb, S, TR, TL):
    pp = mp.mpc(0, w)
    gm = mp.sqrt((1 + TR * pp) / (1 + TL * pp))
    sq = mp.sqrt(pp)
    A = X * sq * gm
    B = Bb * sq / gm
    Sm = S * sq
    t = mp.tanh(A)
    Z = (B + Sm * t) / (B * (Sm + B * t))
    Gam = (B - Sm) / (B + Sm)
    return Z, A, B, Gam, pp, sq


# ---------------------------------------------------------------- R1
random.seed(1)
worst_gamma, worst_reA = 0.0, 1e9
conv = []
for i in range(2000):
    X = mp.mpf(10) ** random.uniform(-1, 1)
    Bb = mp.mpf(10) ** random.uniform(-1.5, 1.5)
    TR = mp.mpf(10) ** random.uniform(-4, 1) * random.choice([0, 1, 1, 1])
    TL = mp.mpf(10) ** random.uniform(-4, 1) * random.choice([0, 1, 1, 1])
    w = mp.mpf(10) ** random.uniform(-6, 6)
    Z, A, B, Gam, pp, sq = parts(w, X, Bb, mp.mpf(1), TR, TL)
    worst_gamma = max(worst_gamma, float(abs(Gam)))
    worst_reA = min(worst_reA, float(mp.re(A) / abs(A)))
    if i < 40:
        E = mp.exp(-2 * A)
        acc = mp.mpc(1)
        errs = []
        for n in range(1, 8):
            acc += 2 * (Gam * E) ** n
            errs.append(float(abs(acc / B - Z) / abs(Z)))
        rho = float(abs(Gam * E))
        conv.append({"rho": rho, "err_N1..N7": errs})
res["R1"] = {"max_abs_Gamma": worst_gamma, "min_ReA_over_absA": worst_reA,
             "partial_sum_examples": conv[:6]}
print("R1 max|Gamma| =", worst_gamma, " min Re A/|A| =", worst_reA)

# ---------------------------------------------------------------- G1, G2
X, Bb, S = mp.mpf(1), mp.mpf("2.13"), mp.mpf(1)
TR, TL = mp.mpf("0.01"), mp.mpf("0.004")
binf = Bb * mp.sqrt(TL / TR)
xinf = X * mp.sqrt(TR / TL)
c1 = (1 / TR - 1 / TL) / 2
Ginf = (binf - S) / (binf + S)
rows = []
for k in np.arange(2.0, 8.01, 0.5):
    w = mp.mpf(10) ** k
    Z, A, B, Gam, pp, sq = parts(w, X, Bb, S, TR, TL)
    alg = Z * binf * sq - 1 - c1 / pp
    # exponentially small part, evaluated without cancellation from the exact representation:
    #   B (Z - 1/B)/2 = Gamma E/(1 - Gamma E),  E = exp(-2A)
    E = mp.exp(-2 * A)
    ln_refl = mp.log(abs(Gam)) - 2 * mp.re(A) - mp.log(abs(1 - Gam * E))
    corrA = 2 * mp.re(A) - xinf * mp.sqrt(2 * w)
    predA = xinf * (TL - TR) / (TR * TL * mp.sqrt(2 * w))
    rows.append(dict(log10_omega=float(k), alg_remainder=float(abs(alg)),
                     two_ReA=float(2 * mp.re(A)), corr_2ReA_over_pred=float(corrA / predA),
                     Gamma_minus_Ginf=float(abs(Gam - Ginf)),
                     ln_refl_over_first=float(ln_refl), ln_refl_pred=float(mp.log(abs(Ginf)) - xinf * mp.sqrt(2 * w) - predA)))
slopes = [float(np.log(rows[i + 1]["alg_remainder"] / rows[i]["alg_remainder"]) /
                ((rows[i + 1]["log10_omega"] - rows[i]["log10_omega"]) * np.log(10))) for i in range(len(rows) - 1)]
gsl = [float(np.log(rows[i + 1]["Gamma_minus_Ginf"] / rows[i]["Gamma_minus_Ginf"]) /
             (0.5 * np.log(10))) for i in range(len(rows) - 1)]
res["G1_G2"] = {"theta": [str(v) for v in (X, Bb, S, TR, TL)], "rows": rows,
                "alg_remainder_slopes": slopes, "Gamma_minus_Ginf_slopes": gsl,
                "Gamma_inf": float(Ginf), "xi_inf": float(xinf), "b_inf": float(binf)}
print("G1 alg remainder slopes (expect -> -2):", ["%.3f" % v for v in slopes[-6:]])
print("G2 (2ReA - xi_inf sqrt(2w)) / prediction (expect -> 1):", ["%.5f" % r["corr_2ReA_over_pred"] for r in rows[-6:]])
print("G2 |Gamma - Gamma_inf| slopes (expect -> -1):", ["%.3f" % v for v in gsl[-6:]])
print("G2 ln|refl| - prediction:", ["%.2e" % (r["ln_refl_over_first"] - r["ln_refl_pred"]) for r in rows[-6:]])

# ---------------------------------------------------------------- G3 (symbolic)
u, al, be, bb = sp.symbols("u alpha beta b", positive=True)   # u = 1/p, alpha = 1/tau_R, beta = 1/tau_l
lng = sp.Rational(1, 2) * (sp.log(1 + al * u) - sp.log(1 + be * u))
ser = sp.series(lng, u, 0, 4).removeO()
k1 = sp.simplify(ser.coeff(u, 1))
k2 = sp.simplify(ser.coeff(u, 2))
lnn0 = sp.Rational(1, 2) * (sp.log(be) - sp.log(al)) - sp.log(bb)   # ln(sqrt(tau_R/tau_l)/b)
J = sp.Matrix([[sp.diff(f, v) for v in (bb, al, be)] for f in (lnn0, k1, k2)])
detJ = sp.factor(sp.simplify(J.det()))
res["G3"] = {"k1": str(k1), "k2": str(k2), "det": str(detJ)}
print("G3 ln g large-p coefficients:", k1, ",", k2, "; det =", detJ)

# ---------------------------------------------------------------- C1, C2 (Cattaneo)
X, Bb, S = mp.mpf(1), mp.mpf("2.13"), mp.mpf(1)
TR, TL = mp.mpf(1), mp.mpf(0)          # epsilon_R = tau_R/xi1^2 = 1: echoes visible (e^{-1})
te = X * mp.sqrt(TR)
c = te / (2 * TR)
crow = []
for k in np.arange(1.0, 7.01, 0.5):
    w0 = mp.mpf(10) ** k
    # max error over one period pi/t_e
    errs = []
    for j in range(64):
        w = w0 + mp.pi / te * j / 64
        Z, A, B, Gam, pp, sq = parts(w, X, Bb, S, TR, TL)
        Zc = mp.sqrt(TR) / Bb * mp.tanh(c + 1j * w * te)
        errs.append(abs(Z - Zc))
    crow.append(dict(log10_omega=float(k), max_err=float(max(errs))))
cs = [float(np.log(crow[i + 1]["max_err"] / crow[i]["max_err"]) / (0.5 * np.log(10))) for i in range(len(crow) - 1)]
print("C1 Cattaneo HF error slopes (expect -> -1/2):", ["%.3f" % v for v in cs[-6:]])

# C2: reconstruction from a window of uniformly spaced frequencies covering M periods
def reconstruct(w0, nper=8, nsamp=1024):
    T = mp.pi / te                      # period in omega
    ws = [w0 + T * nper * j / nsamp for j in range(nsamp)]
    zs = np.array([complex(parts(w, X, Bb, S, TR, TL)[0]) for w in ws])
    F = np.fft.fft(zs) / nsamp
    # e^{-2 i omega t_e} completes t_e W/pi cycles over the window W = nper T; it sits at index N - m
    mag = np.abs(F)
    mag[0] = 0.0
    kk = int(np.argmax(mag))
    m = kk if kk <= nsamp // 2 else nsamp - kk
    W = float(T) * nper
    te_hat = np.pi * m / W
    M = F[0]
    H1 = F[kk]
    ratio = abs(H1) / (2 * abs(M))      # = exp(-t_e/tau_R)
    tr_hat = te_hat / (-np.log(ratio))
    b_hat = np.sqrt(tr_hat) / abs(M)
    x_hat = te_hat / np.sqrt(tr_hat)
    return x_hat, b_hat, tr_hat

recon = []
for k in (2.0, 3.0, 4.0, 5.0):
    xh, bh, th = reconstruct(mp.mpf(10) ** k)
    recon.append(dict(log10_omega0=k, xi1=xh, b=bh, tau_R=th,
                      rel_err=max(abs(xh - 1), abs(bh / 2.13 - 1), abs(th - 1))))
    print("C2 omega0=1e%g: xi1=%.6f b=%.6f tau_R=%.6f" % (k, xh, bh, th))
res["C1_C2"] = {"rows": crow, "slopes": cs, "reconstruction": recon,
                "theta": dict(xi1=1, b=2.13, b_s=1, tau_R=1, tau_l=0)}

# Echo amplitude e^{-t_e/tau_R} = exp(-1/sqrt(eps_R)) for the campaign films
films = {}
for name, e, lam, tau in (("30nm", 30e-9, 78.0, 5.3e-12), ("100nm", 100e-9, 130.0, 11.0e-12),
                          ("300nm", 300e-9, 180.0, 19.5e-12), ("500nm", 500e-9, 200.0, 27.5e-12)):
    a = lam / 2.41e6
    xi2 = e * e / a
    eps = tau / xi2
    films[name] = dict(xi1_sq_ns=xi2 * 1e9, eps_R=eps, Vernotte=np.sqrt(eps), echo_round_trip=float(np.exp(-1 / np.sqrt(eps))),
                       omega_max_xi2=2 * np.pi * 2e8 * xi2, omega_max_tau=2 * np.pi * 2e8 * tau,
                       beta=np.sqrt(35 * 3.03e6) / np.sqrt(lam * 2.41e6))
res["films"] = films
for k, v in films.items():
    print(k, {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items()})

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1, default=str)
print("written", OUT)
