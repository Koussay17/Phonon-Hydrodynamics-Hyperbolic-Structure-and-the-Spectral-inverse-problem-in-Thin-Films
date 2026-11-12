"""Referee checks of workstream B (independent implementation, mpmath).

Model: film (xi, b, tauR, taul) on a Fourier half-space of effusivity bs, contact R.
Z = (L cosh A + sinh A / B) / (B L sinh A + cosh A), A = xi sqrt(p g), B = b sqrt(p/g),
g = (1 + tauR p)/(1 + taul p), S = bs sqrt(p), L = R + 1/S.
Interface temperature conventions at the film/substrate face:
  'local'    : T continuous (or jump R q)              -> load L
  'entropic' : theta = (1 + taul p) T on the GK side   -> load L / (1 + taul p)
(the 'entropic' convention is the one documented in src/twin_fdtr.py, GKBoundary).
"""
import json
import mpmath as mp

mp.mp.dps = 60
OUT = {}


def Z(theta, w, bs=1, R=0, conv="local"):
    xi, b, tR, tl = theta
    p = mp.mpc(0, w)
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * mp.sqrt(p * g)
    B = b * mp.sqrt(p / g)
    S = bs * mp.sqrt(p)
    L = R + 1 / S
    if conv == "entropic":
        L = L / (1 + tl * p)
    E = mp.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    return (1 + r * E) / ((1 - r * E) * B)


def jac(theta, ws, bs=1, R=None, conv="local", h=mp.mpf("1e-25")):
    """Real Jacobian [Re; Im] of ln Z w.r.t. (xi, b, tauR, taul[, R]) by central differences."""
    pars = list(theta) + ([R] if R is not None else [])
    n = len(pars)
    rows = []
    cols = []
    for j in range(n):
        pp = list(pars); pm = list(pars)
        pp[j] += h; pm[j] -= h
        def f(q):
            th = q[:4]
            RR = q[4] if R is not None else 0
            return [mp.log(Z(th, w, bs, RR, conv)) for w in ws]
        fp, fm = f(pp), f(pm)
        cols.append([(a - c) / (2 * h) for a, c in zip(fp, fm)])
    M = mp.matrix(2 * len(ws), n)
    for j in range(n):
        for i in range(len(ws)):
            M[2 * i, j] = mp.re(cols[j][i])
            M[2 * i + 1, j] = mp.im(cols[j][i])
    return M


def svals(M):
    # column-normalise (scale-free rank test), then singular values relative to the largest
    n = M.cols
    Mn = M.copy()
    for j in range(n):
        nrm = mp.sqrt(sum(Mn[i, j] ** 2 for i in range(Mn.rows)))
        if nrm > 0:
            for i in range(Mn.rows):
                Mn[i, j] /= nrm
    s = mp.svd_r(Mn, compute_uv=False)
    s = sorted([abs(x) for x in s], reverse=True)
    return [float(x / s[0]) for x in s]


ws = [mp.mpf(10) ** (mp.mpf(k) / 5) for k in range(-15, 16)]  # 1e-3 .. 1e3, 31 points

# T1: certified two-frequency twin of Proposition TW
th_star = (mp.mpf(1), mp.mpf(1), mp.mpf("0.7"), mp.mpf("0.3"))
th_twin = (mp.mpf("0.38213304403876405621"), mp.mpf("0.99455364692717951823"),
           mp.mpf("1.58320216151666891086"), mp.mpf("0.77319048590574559013"))
t1 = {}
for w in (mp.mpf("0.5"), mp.mpf(4), mp.mpf("1.7"), mp.mpf("0.05"), mp.mpf(40)):
    z1, z2 = Z(th_star, w), Z(th_twin, w)
    t1[str(w)] = {"Z_star": [float(mp.re(z1)), float(mp.im(z1))], "abs_diff": float(abs(z1 - z2)),
                  "rel_diff": float(abs(z1 - z2) / abs(z1))}
OUT["T1_twin"] = t1

# T2: ranks under the two conventions, R = 0
pts = {
    "generic": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.5"), mp.mpf("0.1")),
    "resonance_b2": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.4"), mp.mpf("0.4")),
    "resonance_b=bs": (mp.mpf("0.8"), mp.mpf(1), mp.mpf("0.4"), mp.mpf("0.4")),
    "b=bs_offres": (mp.mpf("0.8"), mp.mpf(1), mp.mpf("0.5"), mp.mpf("0.1")),
    "fourier_point_b2": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("1e-30"), mp.mpf("1e-30")),
}
t2 = {}
for conv in ("local", "entropic"):
    for nm, th in pts.items():
        t2[f"{conv}:{nm}"] = svals(jac(th, ws, conv=conv))
OUT["T2_rank_4par"] = t2

# T3: thin-slab confound on the resonance at R = 0 (ratio dZ/dR / dZ/dxi = b bs^2/(bs^2 - b^2))
th = (mp.mpf("0.8"), mp.mpf("2.5"), mp.mpf("0.4"), mp.mpf("0.4"))
ratios = []
for w in (mp.mpf("0.01"), mp.mpf(1), mp.mpf(30)):
    hR = mp.mpf("1e-25")
    dR = (Z(th, w, R=hR) - Z(th, w, R=-hR)) / (2 * hR)
    dxi = (Z((th[0] + hR,) + th[1:], w) - Z((th[0] - hR,) + th[1:], w)) / (2 * hR)
    ratios.append([float(mp.re(dR / dxi)), float(mp.im(dR / dxi))])
OUT["T3_thin_slab_ratio"] = {"computed": ratios, "predicted": float(mp.mpf("2.5") / (1 - mp.mpf("6.25")))}
# same under the entropic convention
ratios_e = []
for w in (mp.mpf("0.01"), mp.mpf(1), mp.mpf(30)):
    hR = mp.mpf("1e-25")
    dR = (Z(th, w, R=hR, conv="entropic") - Z(th, w, R=-hR, conv="entropic")) / (2 * hR)
    dxi = (Z((th[0] + hR,) + th[1:], w, conv="entropic") - Z((th[0] - hR,) + th[1:], w, conv="entropic")) / (2 * hR)
    ratios_e.append([float(mp.re(dR / dxi)), float(mp.im(dR / dxi))])
OUT["T3_thin_slab_ratio_entropic"] = ratios_e

# T4: rank with resistance (5 params), local vs entropic, on the resonance with R = 0
t4 = {}
for conv in ("local", "entropic"):
    t4[f"{conv}:res_R0_b2"] = svals(jac(pts["resonance_b2"], ws, R=mp.mpf(0) + mp.mpf("1e-40"), conv=conv))
    t4[f"{conv}:res_R0.3_b2"] = svals(jac(pts["resonance_b2"], ws, R=mp.mpf("0.3"), conv=conv))
    t4[f"{conv}:generic_R0.3"] = svals(jac(pts["generic"], ws, R=mp.mpf("0.3"), conv=conv))
OUT["T4_rank_5par"] = t4

# T5: low-frequency Laurent coefficients c0..c3 by a trapezoidal Cauchy integral in s (p = s^2), R = 0
def F_s(theta, s, bs=1, R=0):
    xi, b, tR, tl = theta
    p = s * s
    g = (1 + tR * p) / (1 + tl * p)
    z = xi ** 2 * p * g
    rz = mp.sqrt(z)
    Cc = mp.cosh(rz)
    Ss = mp.sinh(rz) / rz if abs(rz) > 0 else mp.mpf(1)
    X = (xi / b) * g * Ss
    Y = xi * b * p * Ss
    L = R + 1 / (bs * s)
    return (Cc * L + X) / (Y * L + Cc)


th5 = (mp.mpf(7) / 10, mp.mpf(23) / 10, mp.mpf(3) / 10, mp.mpf(1) / 8)
bs5 = mp.mpf(11) / 10
rho = mp.mpf("0.15")
Nn = 256
ck = []
for k in range(0, 5):
    acc = mp.mpc(0)
    for j in range(Nn):
        s = rho * mp.e ** (2j * mp.pi * j / Nn)
        acc += (F_s(th5, s, bs5) - 1 / (bs5 * s)) * s ** (-k)
    ck.append(acc / Nn)
xi, b, tR, tl = th5
bs = bs5
c0 = xi * (bs ** 2 - b ** 2) / (b * bs ** 2)
c1 = xi ** 2 * (b ** 2 - bs ** 2) / bs ** 3
c2 = xi * (tR - tl) / b - xi ** 3 * (3 * b ** 4 - 4 * b ** 2 * bs ** 2 + bs ** 4) / (3 * b * bs ** 4)
c3 = xi ** 2 * (3 * bs ** 4 * (tl - tR) + (3 * b ** 4 - 5 * b ** 2 * bs ** 2 + 2 * bs ** 4) * xi ** 2) / (3 * bs ** 5)
OUT["T5_laurent"] = {
    "numeric": [[float(mp.re(c)), float(mp.im(c))] for c in ck],
    "closed_form_c0_c3": [float(c0), float(c1), float(c2), float(c3)],
    "rel_err_c0_c3": [float(abs(ck[i] - cf) / abs(cf)) for i, cf in enumerate((c0, c1, c2, c3))],
}
print(json.dumps(OUT, indent=1))
with open(__file__.replace(".py", ".json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
