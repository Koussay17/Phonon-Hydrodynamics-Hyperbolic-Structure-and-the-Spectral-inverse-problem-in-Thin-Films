"""Addendum to b2_04: transducer stack with different coefficients at the two film faces (50 digits).

Top face (transducer/film) entropic, bottom face (film/substrate) local or kinetic eps < 0: the film block is the
read-out-theta_e model of Proposition C4, so the stack (Lemma K) is predicted rank-deficient by one on the exotic sets
E_0 and E_eps (kernel (1, 0, -t/xi, t/xi, 0, c_R2)) and on Delta x {R2 = 0}; full rank elsewhere (generic mixed values).
Output: b2_04b_mixed_faces.json
"""
import os
import sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp  # noqa: E402
from models_mp import jacobian, normalized_svals, kernel_residual  # noqa: E402

mp.mp.dps = 50
WS = [mp.mpf(10) ** (mp.mpf(k) / 20 * 3) for k in range(-20, 21)]
F = mp.mpf
T6 = ["xi", "b", "tR", "tl", "R1", "R2"]


def spec(eps0, iface, eps=None):
    d = dict(transducer=True, xiA=F("0.2"), bA=F(3), obs="eps0", eps0=F(eps0), iface=iface)
    if eps is not None:
        d["eps"] = F(eps)
    return d


bl = F("0.6")
R0, t0 = 2 * (1 - bl ** 2) / bl, 4 * (1 - bl ** 2)
be, te, Re = mp.sqrt(2), F(2), mp.sqrt(2)
cases = [
    ("top entropic / bottom local: exotic E_0", spec(1, "local"),
     dict(xi=F(1), b=bl, tR=t0, tl=t0, R1=F("0.3"), R2=R0), 5, [(1, 0, -t0, t0, 0, (bl ** 2 - 1) / bl)]),
    ("top entropic / bottom local: Delta, R2 off E_0", spec(1, "local"),
     dict(xi=F(1), b=bl, tR=t0, tl=t0, R1=F("0.3"), R2=R0 * F("1.1")), 6, []),
    ("top entropic / bottom local: Delta, R2 = 0", spec(1, "local"),
     dict(xi=F(1), b=F(2), tR=F("0.4"), tl=F("0.4"), R1=F("0.3"), R2=F(0)), 5, [(F(2) / 3, 0, 0, 0, 0, 1)]),
    ("top entropic / bottom kinetic eps=-1: exotic E_eps", spec(1, "eps", -1),
     dict(xi=F(1), b=be, tR=te, tl=te, R1=F("0.3"), R2=Re), 5, [(1, 0, -te, te, 0, 1 / mp.sqrt(2))]),
    ("top kinetic 0.3 / bottom kinetic 0.7: Delta, R2 = 0", spec("0.3", "eps", "0.7"),
     dict(xi=F("0.8"), b=F(2), tR=F("0.4"), tl=F("0.4"), R1=F("0.3"), R2=F(0)), 6, []),
    ("top kinetic 0.3 / bottom kinetic 0.7: corner R2 = 0, b = 2", spec("0.3", "eps", "0.7"),
     dict(xi=F("0.8"), b=F(2), tR=F(0), tl=F(0), R1=F("0.3"), R2=F(0)), 5, [(F(2) / 3, 0, 0, 0, 0, 1)]),
    ("top kinetic 0.3 / bottom kinetic 0.7: corner R2 = 0, b = bs", spec("0.3", "eps", "0.7"),
     dict(xi=F("0.8"), b=F(1), tR=F(0), tl=F(0), R1=F("0.3"), R2=F(0)), 5, [(1, 0, 0, 0, 0, 0)]),
]
out = []
for label, sp_, par, rk, kers in cases:
    J = jacobian(par, T6, WS, sp_)
    s = normalized_svals(J)
    res = [kernel_residual(J, k) for k in kers]
    ok = (s[rk - 1] > mp.mpf("1e-12")) and all(x < mp.mpf("1e-25") for x in s[rk:]) and all(r < mp.mpf("1e-25") for r in res)
    out.append(dict(label=label, predicted_rank=rk, svals=[mp.nstr(x, 3) for x in s],
                    kernel_residuals=[mp.nstr(r, 3) for r in res], consistent=bool(ok)))
    print(f"{label:58s} rank={rk} ok={ok} s={[mp.nstr(x, 2) for x in s]} ker={[mp.nstr(r, 2) for r in res]}")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_04b_mixed_faces.json"), "w") as fh:
    json.dump(out, fh, indent=1)
