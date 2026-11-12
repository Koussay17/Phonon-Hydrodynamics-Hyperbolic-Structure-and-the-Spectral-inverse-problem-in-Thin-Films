"""T3, free c_n row (7 unknowns xi, b, tR, tl, tn, R1, R2; same tau_n at both film faces): 60-digit function-level
ranks on the coincidence strata tau_n = tau_l, tau_n = tau_R (off Delta), negative tau_n, Delta with tau_n = t,
and Delta with tau_n = 0, R2 = 0 (b != b_s).  Reuses the independent single-valued model of a03_rank_numeric.py.
Output: a10_free_cn_stack.json
"""
import os
import json

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "a03_rank_numeric.py")).read().split("one = mp.mpf(1)")[0]
ns = {"__file__": os.path.join(HERE, "a03_rank_numeric.py")}
exec(compile(src, "a03_head", "exec"), ns)
analyse = ns["analyse"]
mp.mp.dps = 60
NAMES = ["xi", "b", "tR", "tl", "tn", "R1", "R2"]
conv = dict(obs="theta_q", free=True, stack=True)
f = mp.mpf
CASES = [
    ("free stack generic off Delta", dict(xi=f("0.7"), b=f(2), tR=f("0.5"), tl=f("0.1"), tn=f("0.3"), R1=f("0.3"), R2=f("0.2")), [], 7),
    ("free stack tn = tl off Delta, R2 = 0", dict(xi=f("0.7"), b=f(2), tR=f("0.5"), tl=f("0.1"), tn=f("0.1"), R1=f("0.3"), R2=0), [], 7),
    ("free stack tn = tR off Delta, R2 = 0", dict(xi=f("0.7"), b=f(2), tR=f("0.5"), tl=f("0.1"), tn=f("0.5"), R1=f("0.3"), R2=0), [], 7),
    ("free stack tn = tR, tl = 0, R2 = 0", dict(xi=f("0.7"), b=f(2), tR=f("0.5"), tl=0, tn=f("0.5"), R1=f("0.3"), R2=0), [], 7),
    ("free stack tn < 0 off Delta, R2 = 0", dict(xi=f("0.7"), b=f(1), tR=f("0.5"), tl=f("0.1"), tn=f("-0.2"), R1=f("0.3"), R2=0), [], 7),
    ("free stack Delta tn = t, R2 = 0", dict(xi=f("0.7"), b=f(2), tR=f("0.4"), tl=f("0.4"), tn=f("0.4"), R1=f("0.3"), R2=0),
     [[0, 0, 1, 1, 0, 0, 0]], 6),
    ("free stack Delta tn = 0, R2 = 0, b = 2", dict(xi=f("0.7"), b=f(2), tR=f("0.4"), tl=f("0.4"), tn=0, R1=f("0.3"), R2=0),
     [[0, 0, 1, 1, 0, 0, 0], [f(2) / 3, 0, 0, 0, 0, 0, 1]], 5),
]
OUT = []
for name, par, kers, pred in CASES:
    rk, sv, kres, zc = analyse(par, NAMES, conv, kers)
    OUT.append({"case": name, "rank": rk, "predicted": pred, "sv_tail": sv[-3:], "kernel_residuals": kres})
    print(f"{name:42s} rank {rk}/7 (predicted {pred})  sv {sv[-3:]}  ker_res {['%.1e' % x for x in kres]}", flush=True)
with open(os.path.join(HERE, "a10_free_cn_stack.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
