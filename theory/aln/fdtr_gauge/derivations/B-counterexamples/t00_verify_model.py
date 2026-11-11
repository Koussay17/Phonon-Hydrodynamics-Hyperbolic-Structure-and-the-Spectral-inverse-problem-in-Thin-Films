"""Verification of zmodel.py.

1. Z (variant code) against the project's forward model src/forward_model.response (independent code path).
2. Analytic Jacobian against mpmath numerical differentiation at 50 digits (both variants, with R).
3. Double-precision Jacobian against the mp Jacobian.
4. Limits: e -> 0 (Z -> 1/(bs sqrt p) + R), e -> inf (Z -> 1/B), Fourier resonance (tau_R = tau_l).
   The variant "brief" (formula as typed in the branch brief) fails limit e -> 0.
"""
import sys, json, os
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "C:/Users/Koussay/these/src")
import zmodel as zm
import forward_model as fm

rng = np.random.default_rng(20261111)
out = {}

# ---- 1. against the project's forward model -------------------------------------------------
maxrel = 0.0
for trial in range(200):
    lam = 10 ** rng.uniform(0, 2.7)
    C = 10 ** rng.uniform(5.5, 6.7)
    e = 10 ** rng.uniform(-8, -5.5)
    slam = 10 ** rng.uniform(0, 2.5)
    sC = 10 ** rng.uniform(5.5, 6.7)
    tR = 0.0 if trial % 5 == 0 else 10 ** rng.uniform(-13, -7)
    tl = 0.0 if trial % 7 == 0 else 10 ** rng.uniform(-13, -7)
    Rc = 0.0 if trial % 3 == 0 else 10 ** rng.uniform(-10, -7)
    s = fm.Sample(film_lam=lam, film_rho_c=C, thickness=e, sub_lam=slam, sub_rho_c=sC,
                  contact_resistance=Rc, relaxation_time=tR, nonlocal_time=tl)
    xi1 = e / np.sqrt(lam / C)
    b = np.sqrt(lam * C)
    bs = np.sqrt(slam * sC)
    w = np.logspace(3, 11, 60)
    p = 1j * w
    # skip cases where the project code refuses (overflow guard on Re(sigma e))
    try:
        Zref = fm.response(p, s)
    except OverflowError:
        continue
    Z = zm.Z_np(w, xi1, b, tR, tl, bs=bs, R=Rc)
    maxrel = max(maxrel, float(np.max(np.abs(Z / Zref - 1))))
out["max_rel_diff_vs_forward_model"] = maxrel
print("1. max |Z/Z_forward_model - 1| over 200 random samples:", maxrel)

# ---- 2. analytic Jacobian vs numerical derivative at 50 digits ------------------------------
mp.mp.dps = 50
worst = {"code": 0.0, "brief": 0.0}
for variant in ("code", "brief"):
    for trial in range(30):
        xi1 = mp.mpf(10) ** mp.mpf(rng.uniform(-1, 1))
        b = mp.mpf(10) ** mp.mpf(rng.uniform(-1.5, 1.5))
        tR = mp.mpf(10) ** mp.mpf(rng.uniform(-3, 2))
        tl = mp.mpf(10) ** mp.mpf(rng.uniform(-3, 2))
        R = (mp.mpf(10) ** mp.mpf(rng.uniform(-2, 1))) if variant == "code" else mp.mpf(0)
        bs = mp.mpf(10) ** mp.mpf(rng.uniform(-0.5, 0.5))
        for w in (mp.mpf("1e-3"), mp.mpf("0.37"), mp.mpf(5), mp.mpf(300)):
            Z, J = zm.jac_mp(w, xi1, b, tR, tl, bs, R, variant)
            f = lambda lx, lb, a3, a4, a5: zm.Z_mp(w, mp.exp(lx), mp.exp(lb), a3, a4, bs, a5, variant)
            args = [mp.log(xi1), mp.log(b), tR, tl, R]
            for k in range(5 if variant == "code" else 4):
                def fk(x, k=k):
                    a = list(args)
                    a[k] = x
                    return f(*a)
                num = mp.diff(fk, args[k])
                err = abs(num - J[k]) / (abs(J[k]) + abs(Z) * mp.mpf("1e-30"))
                worst[variant] = max(worst[variant], float(err))
out["jacobian_vs_mpdiff_maxrel"] = worst
print("2. analytic vs numerical derivative (50 digits), max relative error:", worst)

# ---- 3. double vs mp Jacobian ---------------------------------------------------------------
mp.mp.dps = 40
worst_np = 0.0
for trial in range(40):
    xi1 = 10 ** rng.uniform(-1, 1); b = 10 ** rng.uniform(-1.5, 1.5)
    tR = 10 ** rng.uniform(-3, 2); tl = 10 ** rng.uniform(-3, 2); R = 10 ** rng.uniform(-2, 1)
    ws = np.logspace(-6, 6, 25)
    Z, J = zm.jac_np(ws, xi1, b, tR, tl, 1.0, R)
    for i, w in enumerate(ws):
        Zm, Jm = zm.jac_mp(mp.mpf(w), mp.mpf(xi1), mp.mpf(b), mp.mpf(tR), mp.mpf(tl), 1, mp.mpf(R))
        for k in range(5):
            rel = abs(complex(Jm[k]) - J[i, k]) / max(abs(complex(Jm[k])), 1e-300)
            # entries that are tiny relative to Z are allowed absolute error ~ eps |Z|
            if abs(complex(Jm[k])) > 1e-8 * abs(complex(Zm)):
                worst_np = max(worst_np, rel)
out["jac_double_vs_mp_maxrel"] = worst_np
print("3. double vs mp Jacobian, max relative error (entries > 1e-8 |Z|):", worst_np)

# ---- 4. limits ------------------------------------------------------------------------------
w = np.logspace(-3, 3, 13)
for variant in ("code", "brief"):
    Zthin = zm.Z_np(w, 1e-9, 2.0, 0.3, 0.1, bs=1.0, variant=variant)
    ref = 1 / np.sqrt(1j * w)
    print(f"4. variant {variant}: thin-film limit max|Z sqrt(p) bs - 1| =",
          float(np.max(np.abs(Zthin / ref - 1))))
    out[f"thin_limit_{variant}"] = float(np.max(np.abs(Zthin / ref - 1)))
Zthick = zm.Z_np(w, 1e4, 2.0, 0.3, 0.1)
d = zm.pieces_np(w, 1e4, 2.0, 0.3, 0.1)
print("4. thick-film limit max|Z B - 1| =", float(np.max(np.abs(Zthick * d['B'] - 1))))
Zres1 = zm.Z_np(w, 1.0, 2.0, 0.7, 0.7)
Zres2 = zm.Z_np(w, 1.0, 2.0, 0.0, 0.0)
print("4. resonance tau_R = tau_l = 0.7 vs Fourier: max rel diff =",
      float(np.max(np.abs(Zres1 / Zres2 - 1))))
out["resonance_vs_fourier"] = float(np.max(np.abs(Zres1 / Zres2 - 1)))
with open(os.path.join(HERE, "data", "t00_verify_model.json"), "w") as fh:
    json.dump(out, fh, indent=1)
