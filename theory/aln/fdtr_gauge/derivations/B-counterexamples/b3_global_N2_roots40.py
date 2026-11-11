"""Store the six exact roots of the N = 2 counterexample (theta = (1, 1, 0.5, 0.1), omega xi1^2 = 0.2, 2.0)
as 40-digit strings, re-polished by Newton in 60-digit arithmetic, with their residuals."""
import sys, os, json
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm
d = json.load(open(os.path.join(HERE, "data", "b3_global_N2.json")))
rec = [r for r in d["N2"] if r["n_solutions"] > 1][0]
theta = [mp.mpf(str(v)) for v in rec["theta"]]
ws = [mp.mpf(str(v)) for v in rec["ws"]]
mp.mp.dps = 60
Zt = [zm.Z_mp(w, *theta) for w in ws]
def F(a, b, c, e):
    th = [mp.exp(a), mp.exp(b), mp.exp(c), mp.exp(e)]
    out = []
    for w, zt in zip(ws, Zt):
        z = zm.Z_mp(w, *th) / zt - 1
        out += [z.real, z.imag]
    return out
roots = []
for s in rec["solutions"]:
    x = mp.findroot(F, [mp.log(mp.mpf(v)) for v in s["theta_p"]], tol=mp.mpf(10) ** -55)
    th = [mp.exp(v) for v in x]
    res = max(abs(v) for v in F(*x))
    roots.append({"theta_p": [mp.nstr(v, 40) for v in th], "residual_60digits": mp.nstr(res, 3)})
    print([mp.nstr(v, 12) for v in th], "residual", mp.nstr(res, 3))
with open(os.path.join(HERE, "data", "b3_global_N2_roots40.json"), "w") as fh:
    json.dump({"theta": [str(v) for v in rec["theta"]], "omega_xi1sq": [str(v) for v in rec["ws"]],
               "bs": 1, "units": "xi1 = bs = 1", "roots": roots}, fh, indent=1)
