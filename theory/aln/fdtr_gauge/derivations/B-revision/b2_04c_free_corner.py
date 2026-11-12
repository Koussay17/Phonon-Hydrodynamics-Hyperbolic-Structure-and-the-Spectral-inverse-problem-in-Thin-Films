"""Addendum to b2_04: free interface coefficient (tau_n unknown) in the transducer stack at the Fourier corner (50 digits).

Prediction (Lemma K, stack part, r = 0 case): at tau_R = tau_l = tau_n = 0, R2 = 0, b = b_s the kernel is
span{e_xi, (0,0,1,1,0,0,0), (0,0,-2,0,1,0,0)} (coordinates xi, b, tau_R, tau_l, tau_n, R1, R2): rank 4 of 7;
with b != b_s: span{(0,0,1,1,0,0,0), (k_s,0,0,0,0,0,1)}: rank 5 of 7. Output: b2_04c_free_corner.json
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from models_mp import jacobian, normalized_svals, kernel_residual
mp.mp.dps = 50
WS = [mp.mpf(10) ** (mp.mpf(k) / 20 * 3) for k in range(-20, 21)]
F = mp.mpf
spec = dict(transducer=True, xiA=F("0.2"), bA=F(3), obs="theta_q", iface="free")
names = ["xi", "b", "tR", "tl", "tn", "R1", "R2"]
cases = [("corner, b = bs", dict(xi=F("0.8"), b=F(1), tR=F(0), tl=F(0), tn=F(0), R1=F("0.3"), R2=F(0)), 4,
          [(1, 0, 0, 0, 0, 0, 0), (0, 0, 1, 1, 0, 0, 0), (0, 0, -2, 0, 1, 0, 0)]),
         ("corner, b = 2", dict(xi=F("0.8"), b=F(2), tR=F(0), tl=F(0), tn=F(0), R1=F("0.3"), R2=F(0)), 5,
          [(0, 0, 1, 1, 0, 0, 0), (F(2) / 3, 0, 0, 0, 0, 0, 1)]),
         ("Delta t=0.4, b = bs, tn = 0, R2 = 0", dict(xi=F("0.8"), b=F(1), tR=F("0.4"), tl=F("0.4"), tn=F(0), R1=F("0.3"), R2=F(0)), 5,
          [(1, 0, 0, 0, 0, 0, 0), (0, 0, 1, 1, 0, 0, 0)])]
out = []
for label, par, rk, kers in cases:
    J = jacobian(par, names, WS, spec)
    s = normalized_svals(J)
    res = [kernel_residual(J, k) for k in kers]
    ok = (s[rk - 1] > mp.mpf("1e-12")) and all(x < mp.mpf("1e-25") for x in s[rk:]) and all(r < mp.mpf("1e-25") for r in res)
    out.append(dict(label=label, predicted_rank=rk, svals=[mp.nstr(x, 3) for x in s], kernel_residuals=[mp.nstr(r, 3) for r in res], consistent=bool(ok)))
    print(label, "rank", rk, "ok", ok, [mp.nstr(x, 2) for x in s], [mp.nstr(r, 2) for r in res])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_04c_free_corner.json"), "w"), indent=1)
