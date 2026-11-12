"""Independent re-certification of Theorem T4(ii) determinants (b2_05_certified.json), stack, N = 3.

Own jets over Arb complex balls; reflection form for the film block with the same coefficient at both faces,
    Phi = (N/B)(1 + r E)/(1 - r E),  r = (u - 1)/(u + 1),  u = B L/N,  N = 1 + eps tl p,  L = R2 + 1/sqrt(p),
tanh form for the Fourier transducer (xiA, bA).  Rows: Re, Im of d log Z_tot / d theta at each frequency, in the
order (omega_1 Re, omega_1 Im, omega_2 Re, ...), columns (xi, b, tR, tl, R1, R2) -- the same normalisation as
b2_05, so the certified determinants must agree.  256 bits.
Output: a07_certified_minors.json
"""
import os
import json
from flint import acb, arb, arb_mat, fmpq, ctx

ctx.prec = 256
HERE = os.path.dirname(os.path.abspath(__file__))
NP = 6


class J:
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = acb(v)
        self.d = d if d is not None else [acb(0)] * NP

    def __add__(self, o):
        o = o if isinstance(o, J) else J(o)
        return J(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
    __radd__ = __add__

    def __neg__(self):
        return J(-self.v, [-a for a in self.d])

    def __sub__(self, o):
        return self + (-(o if isinstance(o, J) else J(o)))

    def __rsub__(self, o):
        return J(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, J) else J(o)
        return J(self.v * o.v, [a * o.v + self.v * b for a, b in zip(self.d, o.d)])
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, J) else J(o)
        inv = 1 / o.v
        q = self.v * inv
        return J(q, [(a - q * b) * inv for a, b in zip(self.d, o.d)])

    def __rtruediv__(self, o):
        return J(o) / self


def jexp(x):
    e = x.v.exp()
    return J(e, [e * a for a in x.d])


def jsqrt(x):
    r = x.v.sqrt()
    h = 1 / (2 * r)
    return J(r, [h * a for a in x.d])


def q(sv):
    if "/" in sv:
        a, b = sv.split("/")
        return arb(fmpq(int(a), int(b)))
    return arb(sv)


def ztot(theta, w, eps, xiA, bA):
    V = []
    for i, t in enumerate(theta):
        d = [acb(0)] * NP
        d[i] = acb(1)
        V.append(J(acb(t), d))
    xi, b, tR, tl, R1, R2 = V
    p = J(acb(0, w))
    sq = jsqrt(p)
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * jsqrt(p * g)
    B = b * jsqrt(p / g)
    N = 1 + eps * tl * p
    u = B * (R2 + 1 / sq) / N
    r = (u - 1) / (u + 1)
    E = jexp(-2 * A)
    W = R1 + N * (1 + r * E) / (B * (1 - r * E))
    EA = jexp(-2 * xiA * sq)
    th = (1 - EA) / (1 + EA)
    BA = bA * sq
    return (W + th / BA) / (1 + BA * th * W)


def det_J3(theta_s, eps, ws_s, xiA="1/5", bA="3"):
    theta = [q(x) for x in theta_s]
    rows = []
    for w in ws_s:
        Z = ztot(theta, q(w), arb(eps), q(xiA), q(bA))
        gl = [dz / Z.v for dz in Z.d]
        rows.append([x.real for x in gl])
        rows.append([x.imag for x in gl])
    return arb_mat(rows).det()


W3 = ["1/2", "4", "17/10"]
CASES = [
    ("B4/local generic tR>tl, N=3", ["0.8", "2", "0.5", "0.1", "0.3", "0.2"], 0, W3, {}, "[1.584589472e-6 +/- 2.75e-16]"),
    ("B4/local R2=0, N=3", ["0.8", "2", "0.5", "0.1", "0.3", "0"], 0, W3, {}, "[-3.792198679e-6 +/- 2.23e-16]"),
    ("B4/entropic Delta R2=0, N=3", ["0.8", "2", "0.4", "0.4", "0.3", "0"], 1, W3, {}, "[4.323603453e-7 +/- 1.52e-17]"),
    ("B4/entropic Delta R2=0 b=bs, N=3", ["0.8", "1", "0.4", "0.4", "0.3", "0"], 1, W3, {}, "[3.477819722e-7 +/- 5.91e-19]"),
    ("B4/entropic Fourier corner R2>0, N=3", ["0.8", "2", "0", "0", "0.3", "0.2"], 1, W3, {}, "[-2.086057355e-8 +/- 2.34e-18]"),
    ("B4/local Au/AlN/sapphire Cattaneo (tl=0), N=3", ["1", "2.07", "0.00784", "0", "1.82", "3.65"], 0,
     ["1/10", "1", "10"], {"xiA": "0.183", "bA": "1.87"}, "[5.071334408e-10 +/- 1.29e-21]"),
    ("B4/entropic Au/AlN/sapphire GK (tl=10 ps), N=3", ["1", "2.07", "0.00784", "0.00314", "1.82", "3.65"], 1,
     ["1/10", "1", "10"], {"xiA": "0.183", "bA": "1.87"}, "[5.769660682e-10 +/- 3.13e-20]"),
]
OUT = []
for label, th, eps, ws, tr, rec in CASES:
    d = det_J3(th, eps, ws, **tr)
    OUT.append({"case": label, "det": d.str(10, radius=True), "nonzero": not d.contains(0), "record": rec})
    print(f"{label:52s} det {d.str(10, radius=True):32s} certified nonzero {not d.contains(0)}   record {rec}", flush=True)
with open(os.path.join(HERE, "a07_certified_minors.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
