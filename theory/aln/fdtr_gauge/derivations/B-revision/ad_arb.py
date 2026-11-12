"""Forward-mode automatic differentiation over Arb complex balls, and the 1D responses of models_mp.py.

Used by b4_07_certify_twin.py (Krawczyk test). Same formulas as b2_05_certified.py, packaged without side effects.
"""
from flint import acb, arb as _arb, fmpq


def arb(x):
    if isinstance(x, str) and "/" in x:
        a, b = x.split("/")
        return _arb(fmpq(int(a), int(b)))
    if isinstance(x, _arb):
        return x
    return _arb(x)


class D:
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


def dsqrt(x):
    r = x.v.sqrt()
    return D(r, [a / (2 * r) for a in x.d])


def dcosh(x):
    return D(x.v.cosh(), [x.v.sinh() * a for a in x.d])


def dsinh(x):
    return D(x.v.sinh(), [x.v.cosh() * a for a in x.d])


def mobius(M, x):
    (a, b_), (c, d) = M
    return (a * x + b_) / (c * x + d)


def stack_response(vals, names, w, conv, xiA="0.2", bA="3"):
    """Z_tot(i w) of the transducer stack with gradient w.r.t. `names` (all six parameters).
    vals: dict name -> arb ball. conv: 'local' or 'entropic' (same coefficient at both film faces)."""
    n = len(names)
    V = {}
    for k, val in vals.items():
        d = [acb(0)] * n
        if k in names:
            d[names.index(k)] = acb(1)
        V[k] = D(acb(arb(val)), d)
    p = D.c(acb(0, arb(w)), n)
    tl = V["tl"]
    tn = tl if conv == "entropic" else D.c(0, n)
    g = (1 + V["tR"] * p) / (1 + tl * p)
    A = V["xi"] * dsqrt(p * g)
    B = V["b"] * dsqrt(p / g)
    Mf = ((dcosh(A), dsinh(A) / B), (B * dsinh(A), dcosh(A)))
    S = dsqrt(p)
    L = V["R2"] + 1 / S
    Win = V["R1"] + (1 + tn * p) * mobius(Mf, L / (1 + tn * p))
    At = D.c(arb(xiA), n) * S
    Bt = D.c(arb(bA), n) * S
    return mobius(((dcosh(At), dsinh(At) / Bt), (Bt * dsinh(At), dcosh(At))), Win)
