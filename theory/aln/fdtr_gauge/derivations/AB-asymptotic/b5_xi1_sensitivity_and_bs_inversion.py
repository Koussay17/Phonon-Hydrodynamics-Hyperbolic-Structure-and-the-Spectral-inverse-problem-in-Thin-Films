"""B: (i) exact xi1-sensitivity through the reflections, d ln Z/d ln xi1 = -4 A Gamma E/(1 - Gamma^2 E^2), E = e^{-2A},
and its regime asymptotes; (ii) explicit low-frequency inversion on the matched stratum b = b_s, delta != 0:
    L0 = 0  =>  b = b_s;  xi1 = -L3/L2;  delta = b_s L2/xi1;  tau_l from L4 = -R_f tau_l delta + F(xi1, b_s, delta).
"""
import json
import os
import random

import mpmath as mp
import numpy as np
import sympy as sp

from ab_common import film_jacobian, film_reflection

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "b5_xi1_sensitivity_and_bs_inversion.json")
res = {}

# (i) exact formula vs analytic Jacobian column
Om = np.logspace(-4, 6, 400)
p = 1j * Om
worst = 0.0
for (b, tr, tl) in ((2.13, 0.01, 0.004), (1.0, 0.01, 0.004), (0.6, 0.01, 0.0), (2.13, 1.0, 0.0), (0.3, 0.2, 0.05)):
    Z, J = film_jacobian(p, 1.0, b, tr, tl, 1.0)
    Gam, E, A, B = film_reflection(p, 1.0, b, tr, tl, 1.0)
    form = -4 * A * Gam * E / (1 - Gam ** 2 * E ** 2)
    m = np.abs(J[:, 0]) > 1e-250
    worst = max(worst, float(np.max(np.abs(form[m] - J[m, 0]) / np.abs(J[m, 0]))))
res["exact_xi1_formula_max_rel_err"] = worst
print("d ln Z/d ln xi1 = -4 A Gamma E/(1 - Gamma^2 E^2): max rel. err vs analytic Jacobian = %.2e" % worst)

# asymptotes: regime (ii) 4|Gamma0| xi1 sqrt(omega) e^{-xi1 sqrt(2 omega)}; matched b = b_s: Gamma ~ -delta p/4
b, tr, tl = 2.13, 0.001, 0.0004
Om2 = np.logspace(1, 2.6, 9)            # 1/xi1^2 << omega << 1/tau
Z, J = film_jacobian(1j * Om2, 1.0, b, tr, tl, 1.0)
G0 = (b - 1) / (b + 1)
asym = 4 * abs(G0) * np.sqrt(Om2) * np.exp(-np.sqrt(2 * Om2))
res["regime_ii_ratio_exact_over_asymptote"] = (np.abs(J[:, 0]) / asym).tolist()
print("regime (ii) |d ln Z/d ln xi1| / [4|Gamma0| xi1 sqrt(w) e^{-xi1 sqrt(2w)}]:", ["%.3f" % v for v in res["regime_ii_ratio_exact_over_asymptote"]])
b = 1.0
Om3 = np.logspace(-4, -1, 7)
Z, J = film_jacobian(1j * Om3, 1.0, b, tr, tl, 1.0)
asym3 = (tr - tl) * Om3 * 1.0           # |A Gamma E| with Gamma ~ -delta p/4, 4|A| ~ 4 sqrt(omega), E ~ 1 -> delta omega^{3/2}
res["matched_low_freq_xi1_column"] = {"omega": Om3.tolist(), "abs_col": np.abs(J[:, 0]).tolist()}
sl = np.polyfit(np.log(Om3), np.log(np.abs(J[:, 0])), 1)[0]
res["matched_low_freq_slope"] = float(sl)
print("matched b=b_s, low frequency: |d ln Z/d ln xi1| slope vs omega = %.3f (predicted 3/2: delta p/4 x 4A)" % sl)

# (ii) inversion on b = b_s
q, x, bb, s, TR, TL = sp.symbols("q xi1 b b_s tau_R tau_l", positive=True)
G = (1 + TR * q ** 2) / (1 + TL * q ** 2)
g = sp.sqrt(G)
r = s * g / bb
T = sp.tanh(x * q * g)
expr = (r - 1 / r) * T / (1 + T / r) / (s * q)
ser = sp.expand(sp.series(expr, q, 0, 5).removeO())
L = [sp.simplify(ser.coeff(q, k)) for k in range(5)]
Rf, delta = x / bb, TR - TL
F4 = sp.simplify((L[4] + Rf * TL * delta).subs(bb, s))     # depends on (TR, TL) through delta only
fL = [sp.lambdify((x, bb, s, TR, TL), Lk, "mpmath") for Lk in L]
fF4 = sp.lambdify((x, s, TR, TL), F4, "mpmath")
mp.mp.dps = 50
random.seed(7)
errs = []
for i in range(20):
    X = mp.mpf(random.uniform(0.3, 3.0))
    S = mp.mpf(random.uniform(0.5, 2.0))
    tl_ = mp.mpf(random.uniform(0.0, 0.5))
    tr_ = tl_ + mp.mpf(random.uniform(0.01, 0.5))
    Lv = [f(X, S, S, tr_, tl_) for f in fL]     # b = b_s
    assert abs(Lv[0]) < mp.mpf(10) ** -40
    Xh = -Lv[3] / Lv[2]
    dh = S * Lv[2] / Xh
    tlh = -(Lv[4] - fF4(Xh, S, dh, mp.mpf(0))) / ((Xh / S) * dh)
    errs.append(float(max(abs(Xh - X) / X, abs(dh - (tr_ - tl_)) / (tr_ - tl_), abs(tlh - tl_) / max(tl_, mp.mpf("1e-3")))))
res["bs_inversion_max_rel_err"] = max(errs)
print("explicit inversion on b = b_s (20 random theta): max rel. err = %.2e" % max(errs))
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("written", OUT)
