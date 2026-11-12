"""B-2 / B-4, step 4: function-level Jacobian ranks at 50 digits on every stratum of the classification.

For each (model, point): real Jacobian [Re; Im] of log Z over 41 log-spaced frequencies omega in [1e-3, 1e3]
(nondimensional: b_s = 1), central differences (h = 1e-16 at 50 digits, error ~1e-32), column-normalised SVD.
Predicted rank deficiencies must show singular values at the 1e-30 level and predicted kernel vectors must be
annihilated at that level; full-rank claims must show the smallest singular value well above it.
Transducer stack: known Fourier transducer xi_A = 0.2, b_A = 3 (as in review/r3_six_param.py).
Output: b2_04_numeric_rank.json
"""
import os
import sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp  # noqa: E402
from models_mp import jacobian, normalized_svals, kernel_residual  # noqa: E402

mp.mp.dps = 50
WS = [mp.mpf(10) ** (mp.mpf(k) / 20 * 3) for k in range(-20, 21)]   # 1e-3 .. 1e3, 41 points
F = mp.mpf
SQ2 = mp.sqrt(2)

T4 = ["xi", "b", "tR", "tl"]
T5 = T4 + ["R"]
cases = []


def add(label, spec, par, names, rank, kernels=()):
    cases.append(dict(label=label, spec=spec, par=par, names=names, rank=rank, kernels=[list(k) for k in kernels]))


def P(xi=F("0.8"), b=F(2), tR=F("0.5"), tl=F("0.1"), **kw):
    d = dict(xi=F(xi), b=F(b), tR=F(tR), tl=F(tl))
    d.update({k: F(v) for k, v in kw.items()})
    return d


# ------------------------------------------------------------------ bare film, T read at the heated face
loc = dict(obs="T", iface="local")
ent = dict(obs="T", iface="eps", eps=F(1))
add("T/local generic", loc, P(), T4, 4)
add("T/local Delta", loc, P(tR="0.4", tl="0.4"), T4, 3, [(0, 0, 1, 1)])
add("T/local D2", loc, P(b=1, tR="0.4", tl="0.4"), T4, 2, [(1, 0, 0, 0), (0, 0, 1, 1)])
for eps in (F(1), F("0.084"), F("-0.5"), F("2.5")):
    sp_ = dict(obs="T", iface="eps", eps=eps)
    e = mp.nstr(eps, 3)
    add(f"T/eps={e} generic", sp_, P(), T4, 4)
    add(f"T/eps={e} Delta", sp_, P(tR="0.4", tl="0.4"), T4, 4)
    add(f"T/eps={e} Delta b=bs", sp_, P(b=1, tR="0.4", tl="0.4"), T4, 4)
    add(f"T/eps={e} Fourier corner b=2", sp_, P(tR=0, tl=0), T4, 4)
    add(f"T/eps={e} Fourier corner b=bs", sp_, P(b=1, tR=0, tl=0), T4, 3, [(1, 0, 0, 0)])
    add(f"T/eps={e} tR=0", sp_, P(tR=0, tl="0.3"), T4, 4)
    add(f"T/eps={e} tl=0", sp_, P(tR="0.3", tl=0), T4, 4)
    add(f"T/eps={e} +R: Delta R=0", sp_, P(tR="0.4", tl="0.4", R=0), T5, 5)
    add(f"T/eps={e} +R: Delta R=0.3", sp_, P(tR="0.4", tl="0.4", R="0.3"), T5, 5)
    add(f"T/eps={e} +R: corner R=0.3", sp_, P(tR=0, tl=0, R="0.3"), T5, 5)
    add(f"T/eps={e} +R: corner R=0 b=2.5", sp_, P(b="2.5", tR=0, tl=0, R=0), T5, 4,
        [(F("2.5") / (F("2.5") ** 2 - 1), 0, 0, 0, 1)])
    add(f"T/eps={e} +R: corner R=0 b=bs", sp_, P(b=1, tR=0, tl=0, R=0), T5, 4, [(1, 0, 0, 0, 0)])
fre = dict(obs="T", iface="free")
add("T/free generic", fre, P(tn="0.2"), T4 + ["tn"], 5)
add("T/free Delta tn=0.2", fre, P(tR="0.4", tl="0.4", tn="0.2"), T4 + ["tn"], 4, [(0, 0, 1, 1, 0)])
add("T/free Delta tn=0 b=bs", fre, P(b=1, tR="0.4", tl="0.4", tn=0), T4 + ["tn"], 3, [(0, 0, 1, 1, 0), (1, 0, 0, 0, 0)])
add("T/free tn=tl (entropic value)", fre, P(tn="0.1"), T4 + ["tn"], 5)
add("T/free +R: Delta tn=0 R=0 b=2.5", fre, P(b="2.5", tR="0.4", tl="0.4", tn=0, R=0), T4 + ["tn", "R"], 4,
    [(0, 0, 1, 1, 0, 0), (F("2.5") / (F("2.5") ** 2 - 1), 0, 0, 0, 0, 1)])
add("T/free +R: Delta tn=0.2 R=0", fre, P(b="2.5", tR="0.4", tl="0.4", tn="0.2", R=0), T4 + ["tn", "R"], 5,
    [(0, 0, 1, 1, 0, 0)])
add("T/free +R: generic", fre, P(tn="-0.2", R="0.3"), T4 + ["tn", "R"], 6)

# ------------------------------------------------------------------ bare film, theta_e = (1 + tau_l p) T read at the face
bl = F("0.6")
Rx, tx = 2 * (1 - bl ** 2) / bl, 4 * (1 - bl ** 2)          # exotic curve, local interface (xi = 1, b_s = 1)
add("theta_e/local +R: exotic point", dict(obs="theta_e", iface="local"),
    P(xi=1, b=bl, tR=tx, tl=tx, R=Rx), T5, 4, [(1, 0, -tx, tx, (bl ** 2 - 1) / bl)])
add("theta_e/local +R: Delta R>0 off exotic", dict(obs="theta_e", iface="local"),
    P(xi=1, b=bl, tR=tx, tl=tx, R=Rx * F("1.1")), T5, 5)
be = SQ2
add("theta_e/eps=-1 +R: exotic point", dict(obs="theta_e", iface="eps", eps=F(-1)),
    P(xi=1, b=be, tR=2, tl=2, R=SQ2), T5, 4, [(1, 0, -2, 2, 1 / SQ2)])
add("theta_e/eps=1 +R: Delta R>0", dict(obs="theta_e", iface="eps", eps=F(1)),
    P(xi=1, b=be, tR=2, tl=2, R=SQ2), T5, 5)

# ------------------------------------------------------------------ known transducer, unknown R1, R2 (B-4)
T6 = T4 + ["R1", "R2"]


def S(conv, eps=None):
    d = dict(transducer=True, xiA=F("0.2"), bA=F(3), obs="theta_q")
    if conv == "local":
        d["iface"] = "local"
    elif conv == "free":
        d["iface"] = "free"
    else:
        d["iface"] = "eps"
        d["eps"] = F(eps)
    return d


bb = F(2)
kv = bb / (bb ** 2 - 1)
for conv in ("local", "entropic"):
    sp_ = S("local") if conv == "local" else S("eps", 1)
    loc_ = conv == "local"
    add(f"B4/{conv} generic", sp_, P(R1="0.3", R2="0.2"), T6, 6)
    add(f"B4/{conv} b=bs", sp_, P(b=1, R1="0.3", R2="0.2"), T6, 6)
    add(f"B4/{conv} Cattaneo tl=0", sp_, P(tl=0, R1="0.3", R2="0.2"), T6, 6)
    add(f"B4/{conv} tR=0", sp_, P(tR=0, tl="0.3", R1="0.3", R2="0.2"), T6, 6)
    add(f"B4/{conv} R1=0", sp_, P(R1=0, R2="0.2"), T6, 6)
    add(f"B4/{conv} Delta R2>0", sp_, P(tR="0.4", tl="0.4", R1="0.3", R2="0.2"), T6,
        5 if loc_ else 6, [(0, 0, 1, 1, 0, 0)] if loc_ else [])
    add(f"B4/{conv} Delta R2=0", sp_, P(tR="0.4", tl="0.4", R1="0.3", R2=0), T6,
        4 if loc_ else 6, [(0, 0, 1, 1, 0, 0), (kv, 0, 0, 0, 0, 1)] if loc_ else [])
    add(f"B4/{conv} Delta R2=0 b=bs", sp_, P(b=1, tR="0.4", tl="0.4", R1="0.3", R2=0), T6,
        4 if loc_ else 6, [(0, 0, 1, 1, 0, 0), (1, 0, 0, 0, 0, 0)] if loc_ else [])
    add(f"B4/{conv} Fourier corner R2>0", sp_, P(tR=0, tl=0, R1="0.3", R2="0.2"), T6,
        5 if loc_ else 6, [(0, 0, 1, 1, 0, 0)] if loc_ else [])
    add(f"B4/{conv} Fourier corner R2=0 b=2", sp_, P(tR=0, tl=0, R1="0.3", R2=0), T6,
        4 if loc_ else 5, [(0, 0, 1, 1, 0, 0), (kv, 0, 0, 0, 0, 1)] if loc_ else [(kv, 0, 0, 0, 0, 1)])
    add(f"B4/{conv} Fourier corner R2=0 b=bs", sp_, P(b=1, tR=0, tl=0, R1="0.3", R2=0), T6,
        4, [(0, 0, 1, 1, 0, 0), (1, 0, 0, 0, 0, 0)] if loc_ else [(1, 0, 0, 0, 0, 0), (0, 0, -1, 1, 0, 0)])
for eps in (F("0.5"), F("-0.5"), F("0.084")):
    sp_ = S("eps", eps)
    e = mp.nstr(eps, 3)
    add(f"B4/eps={e} Delta R2=0", sp_, P(tR="0.4", tl="0.4", R1="0.3", R2=0), T6, 6)
    add(f"B4/eps={e} Cattaneo tl=0", sp_, P(tl=0, R1="0.3", R2="0.2"), T6, 6)
    add(f"B4/eps={e} corner R2=0 b=bs", sp_, P(b=1, tR=0, tl=0, R1="0.3", R2=0), T6, 4,
        [(1, 0, 0, 0, 0, 0), (0, 0, 1 - 2 * eps, 1, 0, 0)])
spf = S("free")
add("B4/free generic", spf, P(tn="0.15", R1="0.3", R2="0.2"), T4 + ["tn", "R1", "R2"], 7)
add("B4/free Delta", spf, P(tR="0.4", tl="0.4", tn="0.15", R1="0.3", R2="0.2"), T4 + ["tn", "R1", "R2"], 6,
    [(0, 0, 1, 1, 0, 0, 0)])
add("B4/free tn=tl (entropic value)", spf, P(tn="0.1", R1="0.3", R2="0.2"), T4 + ["tn", "R1", "R2"], 7)
add("B4/free Delta tn=0 R2=0", spf, P(tR="0.4", tl="0.4", tn=0, R1="0.3", R2=0), T4 + ["tn", "R1", "R2"], 5,
    [(0, 0, 1, 1, 0, 0, 0), (kv, 0, 0, 0, 0, 0, 1)])

out = []
worst_def, worst_full = mp.mpf(0), mp.mpf(1)
for c in cases:
    J = jacobian(c["par"], c["names"], WS, c["spec"])
    s = normalized_svals(J)
    n = len(c["names"])
    rk = c["rank"]
    res = [kernel_residual(J, k) for k in c["kernels"]]
    # predicted: s[rk-1] >> 1e-25 and s[rk:] < 1e-25
    ok = (s[rk - 1] > mp.mpf("1e-12")) and all(x < mp.mpf("1e-25") for x in s[rk:]) and all(r < mp.mpf("1e-25") for r in res)
    if rk < n:
        worst_def = max(worst_def, max(s[rk:]))
    worst_full = min(worst_full, s[rk - 1])
    row = dict(label=c["label"], n=n, predicted_rank=rk, svals=[mp.nstr(x, 3) for x in s],
               kernel_residuals=[mp.nstr(r, 3) for r in res], consistent=bool(ok))
    out.append(row)
    print(f"{c['label']:42s} n={n} rank={rk} ok={ok}  s={[mp.nstr(x, 2) for x in s]}  ker={[mp.nstr(r, 2) for r in res]}")
summary = dict(cases=len(out), all_consistent=all(r["consistent"] for r in out),
               largest_deficient_singular_value=mp.nstr(worst_def, 3),
               smallest_retained_singular_value=mp.nstr(worst_full, 3))
print(json.dumps(summary, indent=1))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_04_numeric_rank.json"), "w") as fh:
    json.dump(dict(summary=summary, cases=out, frequencies="41 log-spaced in [1e-3, 1e3]", dps=50), fh, indent=1)
