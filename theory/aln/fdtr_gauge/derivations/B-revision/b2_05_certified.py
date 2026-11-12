"""B-2 / B-4, step 5: certified (Arb ball arithmetic, 256 bits) finite-frequency Jacobian determinants.

A nonzero n x n minor of the real finite-frequency Jacobian J_N(theta; omega) (rows Re/Im of d log Z/d theta at N
frequencies) certifies rank n of the function-level Jacobian at theta (the finite-frequency rank is a lower bound),
and, by real analyticity in (theta, omega) (Lemma F of proofs/B-analytic.md), generic sufficiency of those N
frequencies on the connected stratum containing theta. Derivatives: forward-mode automatic differentiation over
acb (exact rules for +, -, *, /, sqrt, exp, cosh, sinh), so no finite-difference error enters.
Models: models_mp.py conventions (b_s = 1). Output: b2_05_certified.json
"""
import os
import json
import itertools
from flint import acb, arb as _arb, arb_mat, ctx, fmpq


def arb(x):
    """Ball containing the decimal or rational number given as a string (or a number)."""
    if isinstance(x, str) and "/" in x:
        a, b = x.split("/")
        return _arb(fmpq(int(a), int(b)))
    return _arb(x)

ctx.prec = 256


class D:
    """acb value with gradient (list of acb)."""
    __slots__ = ("v", "d")

    def __init__(self, v, d):
        self.v, self.d = v, d

    @staticmethod
    def c(x, n):
        return D(acb(x), [acb(0)] * n)

    def _w(self, o):
        return o if isinstance(o, D) else D(acb(o), [acb(0)] * len(self.d))

    def __add__(self, o):
        o = self._w(o)
        return D(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
    __radd__ = __add__

    def __sub__(self, o):
        o = self._w(o)
        return D(self.v - o.v, [a - b for a, b in zip(self.d, o.d)])

    def __rsub__(self, o):
        return self._w(o) - self

    def __mul__(self, o):
        o = self._w(o)
        return D(self.v * o.v, [a * o.v + self.v * b for a, b in zip(self.d, o.d)])
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._w(o)
        q = self.v / o.v
        return D(q, [(a - q * b) / o.v for a, b in zip(self.d, o.d)])

    def __rtruediv__(self, o):
        return self._w(o) / self

    def __neg__(self):
        return D(-self.v, [-a for a in self.d])


def dsqrt(x):
    r = x.v.sqrt()
    return D(r, [a / (2 * r) for a in x.d])


def dexp(x):
    e = x.v.exp()
    return D(e, [e * a for a in x.d])


def dcosh(x):
    return D(x.v.cosh(), [x.v.sinh() * a for a in x.d])


def dsinh(x):
    return D(x.v.sinh(), [x.v.cosh() * a for a in x.d])


def mobius(M, x):
    (a, b_), (c, d) = M
    return (a * x + b_) / (c * x + d)


def response(par, names, spec, w):
    """Z(i w) with gradient w.r.t. `names`. par: dict of arb/strings."""
    n = len(names)
    V = {}
    for k, val in par.items():
        if k in names:
            d = [acb(0)] * n
            d[names.index(k)] = acb(1)
            V[k] = D(acb(arb(val)), d)
        else:
            V[k] = D.c(arb(val), n)
    p = D.c(acb(0, arb(w)), n)
    tl = V["tl"]
    iface, obs = spec.get("iface", "local"), spec.get("obs", "T")
    tn = D.c(0, n) if iface == "local" else (V["tn"] if iface == "free" else tl * arb(spec["eps"]))
    to = {"T": D.c(0, n), "theta_e": tl, "theta_q": tn}[obs]
    g = (1 + V["tR"] * p) / (1 + tl * p)
    A = V["xi"] * dsqrt(p * g)
    B = V["b"] * dsqrt(p / g)
    Mf = ((dcosh(A), dsinh(A) / B), (B * dsinh(A), dcosh(A)))
    S = dsqrt(p)
    if spec.get("transducer"):
        L = V["R2"] + 1 / S
        Win = V["R1"] + (1 + to * p) * mobius(Mf, L / (1 + tn * p))
        At = D.c(arb(spec["xiA"]), n) * S
        Bt = D.c(arb(spec["bA"]), n) * S
        return mobius(((dcosh(At), dsinh(At) / Bt), (Bt * dsinh(At), dcosh(At))), Win)
    L = V.get("R", D.c(0, n)) + 1 / S
    return (1 + to * p) * mobius(Mf, L / (1 + tn * p))


def real_jacobian(par, names, spec, ws):
    rows = []
    for w in ws:
        Z = response(par, names, spec, w)
        g = [dz / Z.v for dz in Z.d]          # d log Z
        rows.append([x.real for x in g])
        rows.append([x.imag for x in g])
    return rows


def certify(par, names, spec, ws):
    """Largest certified-nonzero n x n minor over row subsets (by midpoint magnitude)."""
    rows = real_jacobian(par, names, spec, ws)
    n = len(names)
    best = None
    for idx in itertools.combinations(range(len(rows)), n):
        det = arb_mat([rows[i] for i in idx]).det()
        mag = abs(float(det.mid()))
        if best is None or mag > best[0]:
            best = (mag, idx, det)
    mag, idx, det = best
    nonzero = not det.contains(0)
    return dict(rows=list(idx), det=det.str(10, radius=True), certified_nonzero=bool(nonzero))


T4 = ["xi", "b", "tR", "tl"]
T5 = T4 + ["R"]
T6 = T4 + ["R1", "R2"]
ent = dict(obs="T", iface="eps", eps="1")
kin = dict(obs="T", iface="eps", eps="0.084")
B4L = dict(transducer=True, xiA="0.2", bA="3", obs="theta_q", iface="local")
B4E = dict(transducer=True, xiA="0.2", bA="3", obs="theta_q", iface="eps", eps="1")
W2 = ["1/2", "4"]
W3 = ["1/2", "4", "17/10"]


def P(**kw):
    d = dict(xi="0.8", b="2", tR="0.5", tl="0.1")
    d.update({k: str(v) for k, v in kw.items()})
    return d


cases = [
    ("T/entropic Delta (t=0.4), N=2", P(tR="0.4", tl="0.4"), T4, ent, W2),
    ("T/entropic Delta b=bs, N=2", P(b="1", tR="0.4", tl="0.4"), T4, ent, W2),
    ("T/entropic Fourier corner b=2, N=2", P(tR="0", tl="0"), T4, ent, W2),
    ("T/entropic tR=0, N=2", P(tR="0", tl="0.3"), T4, ent, W2),
    ("T/entropic tl=0, N=2", P(tR="0.3", tl="0"), T4, ent, W2),
    ("T/kinetic eps=0.084 Delta, N=2", P(tR="0.4", tl="0.4"), T4, kin, W2),
    ("T/entropic +R: Delta R=0, N=3", P(tR="0.4", tl="0.4", R="0"), T5, ent, W3),
    ("T/entropic +R: Delta R=0.3, N=3", P(tR="0.4", tl="0.4", R="0.3"), T5, ent, W3),
    ("T/entropic +R: corner R=0.3, N=3", P(tR="0", tl="0", R="0.3"), T5, ent, W3),
    ("T/entropic +R: corner R=0 b=2.5 (rank 4: 4x4 minor of 5 cols)", P(b="2.5", tR="0", tl="0", R="0"),
     ["xi", "b", "tR", "tl"], ent, W3),
    ("B4/local generic tR>tl, N=3", P(R1="0.3", R2="0.2"), T6, B4L, W3),
    ("B4/local generic tR<tl, N=3", P(tR="0.1", tl="0.5", R1="0.3", R2="0.2"), T6, B4L, W3),
    ("B4/local b=bs, N=3", P(b="1", R1="0.3", R2="0.2"), T6, B4L, W3),
    ("B4/local Cattaneo tl=0, N=3", P(tl="0", R1="0.3", R2="0.2"), T6, B4L, W3),
    ("B4/local tR=0, N=3", P(tR="0", tl="0.3", R1="0.3", R2="0.2"), T6, B4L, W3),
    ("B4/local R1=0, N=3", P(R1="0", R2="0.2"), T6, B4L, W3),
    ("B4/local R2=0, N=3", P(R1="0.3", R2="0"), T6, B4L, W3),
    ("B4/local AlN-like (b=2.07, tR=0.03, tl=0), N=3", P(xi="1", b="2.07", tR="0.03", tl="0", R1="0.05", R2="0.1"), T6, B4L, W3),
    ("B4/entropic generic, N=3", P(R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic Delta R2>0, N=3", P(tR="0.4", tl="0.4", R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic Delta R2=0, N=3", P(tR="0.4", tl="0.4", R1="0.3", R2="0"), T6, B4E, W3),
    ("B4/entropic Delta R2=0 b=bs, N=3", P(b="1", tR="0.4", tl="0.4", R1="0.3", R2="0"), T6, B4E, W3),
    ("B4/entropic Fourier corner R2>0, N=3", P(tR="0", tl="0", R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic corner R2=0 b=2 (rank 5: cols xi,b,tR,tl,R1)", P(tR="0", tl="0", R1="0.3", R2="0"),
     ["xi", "b", "tR", "tl", "R1"], B4E, W3),
    ("B4/entropic b=bs, N=3", P(b="1", R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic Cattaneo tl=0, N=3", P(tl="0", R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic tR=0, N=3", P(tR="0", tl="0.3", R1="0.3", R2="0.2"), T6, B4E, W3),
    ("B4/entropic R1=0, N=3", P(R1="0", R2="0.2"), T6, B4E, W3),
    ("B4/entropic R2=0 (off Delta), N=3", P(R1="0.3", R2="0"), T6, B4E, W3),
    # Au/AlN/sapphire in units b_s = 1, time xi_1^2 = e^2/a = 3.19 ns: tau_R = 25 ps, R1 = 1e-8, R2 = 2e-8 m2K/W,
    # Au 80 nm: xi_A = 0.183, b_A = 1.87
    ("B4/local Au/AlN/sapphire Cattaneo (tl=0), N=3",
     P(xi="1", b="2.07", tR="0.00784", tl="0", R1="1.82", R2="3.65"), T6,
     dict(B4L, xiA="0.183", bA="1.87"), ["1/10", "1", "10"]),
    ("B4/entropic Au/AlN/sapphire GK (tl=10 ps), N=3",
     P(xi="1", b="2.07", tR="0.00784", tl="0.00314", R1="1.82", R2="3.65"), T6,
     dict(B4E, xiA="0.183", bA="1.87"), ["1/10", "1", "10"]),
]
out = []
for label, par, names, spec, ws in cases:
    r = certify(par, names, spec, ws)
    r["label"] = label
    out.append(r)
    print(f"{label:62s} certified={r['certified_nonzero']}  det={r['det']}")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_05_certified.json"), "w") as fh:
    json.dump(dict(prec_bits=256, cases=out, all_certified=all(c["certified_nonzero"] for c in out)), fh, indent=1)
print("all certified:", all(c["certified_nonzero"] for c in out))
