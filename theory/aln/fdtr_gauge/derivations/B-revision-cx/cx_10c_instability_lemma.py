"""cx_10c: instability for eps < 0, checks of the elementary argument.

Bare film, read-out T(0), tau_n = eps tau_l, eps < 0, tau_R, tau_l > 0.  On real s > 0 all quantities are
real; D(s) = B L tanh A + N is continuous, D(0+) = 1, and D(s) -> -infinity (B L tanh A = O(s), N ~ -|eps| tau_l s^2),
so D has a zero s0 > 0; there the numerator L + N tanh(A)/B = L / cosh^2 A > 0, hence p0 = s0^2 > 0 is a pole of Z.
(1) thick-film limit: p0 -> (1 + b/b_s)/(|eps| tau) (tau_R = tau_l = tau, R = 0) as xi -> infinity;
(2) smallest real pole vs eps for the AlN note-16 film on sapphire (bare, R = 0 and R = R2 of note 16) and the stack:
    growth rate p0/t0 in 1/s;
(3) tau_R = 0 (Cattaneo-free GK limit, g = 1/(1 + tau_l p)): a case where the argument does not apply.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
mp.mp.dps = 40
out = {}


def D_bare(prm, s):
    xi, b, tR, tl, tn, R = prm
    p = s * s
    g = (1 + tR * p) / (1 + tl * p); sg = mp.sqrt(g)
    A = xi * s * sg; B = b * s / sg; N = 1 + tn * p; L = R + 1 / s
    return B * L * mp.tanh(A) + N


def D_stack(prm, s, xiA, bA):
    xi, b, tR, tl, tn, to, R1, R2 = prm
    p = s * s
    g = (1 + tR * p) / (1 + tl * p); sg = mp.sqrt(g)
    A = xi * s * sg; B = b * s / sg; N = 1 + tn * p; N0 = 1 + to * p; L = R2 + 1 / s; T = mp.tanh(A)
    TA = mp.tanh(xiA * s)
    return bA * s * TA * (R1 * (B * L * T + N) + N0 * (L + N * T / B)) + (B * L * T + N)


def first_real_zero(f, smax=1e7, n=4000):
    ss = [mp.mpf(10) ** (mp.mpf(-3) + (mp.log10(smax) + 3) * k / n) for k in range(n + 1)]
    prev = f(ss[0])
    for k in range(1, len(ss)):
        cur = f(ss[k])
        if prev * cur < 0:
            return mp.findroot(f, (ss[k - 1], ss[k]), solver="anderson")
        prev = cur
    return None


asy = []
for eps in [-0.05, -0.2, -0.5, -1.0]:
    for b in [0.5, 2.0]:
        t = mp.mpf("0.3")
        pred = (1 + b) / (abs(eps) * t)
        for xi in [1, 3, 10, 30]:
            s0 = first_real_zero(lambda s: D_bare((xi, b, t, t, eps * t, 0), s))
            asy.append(dict(eps=eps, b=b, tau=0.3, xi=xi, p0=float(s0 ** 2) if s0 is not None else None, predicted=float(pred)))
            print(f"thick film eps={eps} b={b} xi={xi}: p0 = {float(s0**2) if s0 is not None else None}  prediction {float(pred):.6g}", flush=True)
out["thick_film"] = asy

thn, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=10e-12)
xi, b, tR, tl, R1, R2 = [mp.mpf(float(x)) for x in thn]
rows = []
for eps in [-1.0, -0.5, -0.2, -0.1, -0.03, -0.01]:
    r_b0 = first_real_zero(lambda s: D_bare((xi, b, tR, tl, eps * tl, 0), s))
    r_bR = first_real_zero(lambda s: D_bare((xi, b, tR, tl, eps * tl, R2), s))
    r_st = first_real_zero(lambda s: D_stack((xi, b, tR, tl, eps * tl, eps * tl, R1, R2), s, mp.mpf(xiA), mp.mpf(bA)))
    rec = dict(eps=eps, bare_R0_rate=float(r_b0 ** 2 / t0) if r_b0 else None,
               bare_R2_rate=float(r_bR ** 2 / t0) if r_bR else None, stack_rate=float(r_st ** 2 / t0) if r_st else None)
    rows.append(rec)
    print(f"AlN note-16 eps={eps}: growth rate p0/t0 [1/s]: bare R=0 {rec['bare_R0_rate']}, bare R=R2 {rec['bare_R2_rate']}, "
          f"stack {rec['stack_rate']}", flush=True)
out["aln_growth_rates"] = rows

# tau_R = 0
c3 = []
for R in [0.0, 0.5, 5.0]:
    for eps in [-0.1, -1.0]:
        s0 = first_real_zero(lambda s: D_bare((mp.mpf(1), mp.mpf(2), 0, mp.mpf("0.3"), eps * mp.mpf("0.3"), R), s))
        c3.append(dict(R=R, eps=eps, p0=float(s0 ** 2) if s0 else None))
        print(f"tau_R = 0, tau_l = 0.3, R = {R}, eps = {eps}: first real pole p0 = {float(s0**2) if s0 else None}", flush=True)
out["tauR_zero"] = c3
with open(DIR + "cx_10c_instability_lemma.json", "w") as f:
    json.dump(out, f, indent=1)
