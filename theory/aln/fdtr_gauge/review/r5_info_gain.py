"""Threshold-free comparison for workstream A: information that FDTR adds to the static-data range.

Prior: relative change x ~ U[-w, w] (admissible half-width w). Linear-Gaussian measurement y = x + n, n ~ N(0, sigma^2)
(the A5 profile shows the response is linear in tau_R at these w). Report, as a function of r = sigma/w:
  * expected posterior variance / prior variance,
  * mutual information I(x; y) in bits (I = 1 bit would mean halving the admissible interval),
  * probability that a 2-sigma test rejects the far end of the range when the truth sits at one end.
sigma values: experiments/A/out (07-experiments.md A4; audit section 7 for the optimal-design and favourable cases).
"""
import json
import numpy as np
from scipy import integrate, stats


def metrics(r):
    s = r
    # density of y for w = 1
    def py(y):
        return (stats.norm.cdf((y + 1) / s) - stats.norm.cdf((y - 1) / s)) / 2.0
    lim = 1 + 12 * s
    def hy_int(y):
        p = py(y)
        return -p * np.log2(p) if p > 0 else 0.0
    hy = integrate.quad(hy_int, -lim, lim, limit=400, points=[-1, 1])[0]
    hn = 0.5 * np.log2(2 * np.pi * np.e * s * s)
    mi = hy - hn
    # expected posterior variance (posterior = N(y, s^2) truncated to [-1, 1])
    def postvar(y):
        a, b = (-1 - y) / s, (1 - y) / s
        tn = stats.truncnorm(a, b, loc=y, scale=s)
        return tn.var()
    epv = integrate.quad(lambda y: py(y) * postvar(y), -lim, lim, limit=400, points=[-1, 1])[0]
    # truth at x = -1 (one end); reject x = +1 at 2 sigma if y < 1 - 2 s
    p_rej = stats.norm.cdf(((1 - 2 * s) - (-1)) / s)
    return {"sigma_over_w": r, "post_var_over_prior": epv / (1 / 3), "MI_bits": mi, "P(2sigma test separates ends)": p_rej}


cases = {
    "500nm N0 0.01deg (w=0.08)": 0.154 / 0.08,
    "500nm N0 0.01deg (w=0.055)": 0.154 / 0.055,
    "500nm N0 0.01deg c-optimal design (w=0.08)": 0.0558 / 0.08,
    "500nm N1 0.01deg (w=0.08)": 0.205 / 0.08,
    "500nm N2 0.01deg (w=0.08)": 0.340 / 0.08,
    "500nm N2 0.01deg c-optimal (w=0.08)": 0.195 / 0.08,
    "500nm N2 0.01deg all favourable (w=0.08)": 0.093 / 0.08,
    "500nm N2 0.1deg (w=0.08)": 2.73 / 0.08,
    "500nm N2 0.1deg all favourable (w=0.08)": 0.159 / 0.08,
    "500nm N3 0.01deg (w=0.08)": 2.25 / 0.08,
    "300nm N0 0.01deg": 0.246 / 0.052,
    "300nm N2 0.01deg": 0.662 / 0.052,
    "300nm N2 0.1deg": 5.36 / 0.052,
    "100nm N2 0.01deg (w=0.037)": 5.68 / 0.037,
    "100nm N2 0.01deg (w=0.08 calibrated)": 5.68 / 0.08,
    "30nm N2 0.01deg": 188 / 0.028,
    "reference r=1": 1.0,
    "reference r=0.5": 0.5,
    "reference r=0.25": 0.25,
}
out = {}
for k, r in cases.items():
    out[k] = metrics(r)
    m = out[k]
    print(f"{k:45s} r={r:9.3f}  postvar/prior={m['post_var_over_prior']:.4f}  MI={m['MI_bits']:.4f} bit  P2s={m['P(2sigma test separates ends)']:.3f}")
with open(__file__.replace(".py", ".json"), "w") as fh:
    json.dump(out, fh, indent=1)
