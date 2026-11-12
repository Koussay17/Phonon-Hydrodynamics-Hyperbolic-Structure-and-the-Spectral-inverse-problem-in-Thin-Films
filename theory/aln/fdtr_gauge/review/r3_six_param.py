"""Rank of the 1D response with a KNOWN Fourier transducer on top and UNKNOWN top (R1) and bottom (R2) resistances:
parameters (xi, b, tauR, taul, R1, R2). Local interface temperature (Assumption M). Nondimensional: bs = 1.
Transducer: Fourier layer xi_A = 0.2, b_A = 3 (known). Z_tot = Mobius(T_A, R1 + Z_film(R2)).
"""
import json
import mpmath as mp

mp.mp.dps = 50
XA, BA = mp.mpf("0.2"), mp.mpf(3)


def zfilm(th, p, R2, bs=1):
    xi, b, tR, tl = th
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * mp.sqrt(p * g)
    B = b * mp.sqrt(p / g)
    S = bs * mp.sqrt(p)
    L = R2 + 1 / S
    E = mp.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    return (1 + r * E) / ((1 - r * E) * B)


def ztot(q, w):
    xi, b, tR, tl, R1, R2 = q
    p = mp.mpc(0, w)
    zb = R1 + zfilm((xi, b, tR, tl), p, R2)
    A = XA * mp.sqrt(p)
    B = BA * mp.sqrt(p)
    return (mp.cosh(A) * zb + mp.sinh(A) / B) / (B * mp.sinh(A) * zb + mp.cosh(A))


ws = [mp.mpf(10) ** (mp.mpf(k) / 5) for k in range(-20, 21)]


def svals(q, h=mp.mpf("1e-22")):
    n = len(q)
    M = mp.matrix(2 * len(ws), n)
    for j in range(n):
        qp = list(q); qm = list(q)
        qp[j] += h; qm[j] -= h
        for i, w in enumerate(ws):
            d = (mp.log(ztot(qp, w)) - mp.log(ztot(qm, w))) / (2 * h)
            M[2 * i, j] = mp.re(d)
            M[2 * i + 1, j] = mp.im(d)
    for j in range(n):
        nrm = mp.sqrt(sum(M[i, j] ** 2 for i in range(M.rows)))
        if nrm < mp.mpf("1e-40"):
            continue
        for i in range(M.rows):
            M[i, j] /= nrm
    s = sorted([abs(x) for x in mp.svd_r(M, compute_uv=False)], reverse=True)
    return [float(x / s[0]) for x in s]


tiny = mp.mpf("1e-35")
pts = {
    "generic": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.5"), mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.2")),
    "generic_R2=0": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.5"), mp.mpf("0.1"), mp.mpf("0.3"), tiny),
    "generic_b=bs": (mp.mpf("0.8"), mp.mpf(1), mp.mpf("0.5"), mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.2")),
    "cattaneo": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.5"), tiny, mp.mpf("0.3"), mp.mpf("0.2")),
    "resonance": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.4"), mp.mpf("0.4"), mp.mpf("0.3"), mp.mpf("0.2")),
    "resonance_R2=0": (mp.mpf("0.8"), mp.mpf(2), mp.mpf("0.4"), mp.mpf("0.4"), mp.mpf("0.3"), tiny),
    "resonance_b=bs_R2=0": (mp.mpf("0.8"), mp.mpf(1), mp.mpf("0.4"), mp.mpf("0.4"), mp.mpf("0.3"), tiny),
}
out = {k: svals(v) for k, v in pts.items()}
for k, v in out.items():
    print(k, ["%.2e" % x for x in v])
with open(__file__.replace(".py", ".json"), "w") as fh:
    json.dump(out, fh, indent=1)
