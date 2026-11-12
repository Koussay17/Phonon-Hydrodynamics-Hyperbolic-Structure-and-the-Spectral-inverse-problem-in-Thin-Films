"""Information measures for the comparison of a Gaussian measurement with a bounded static-data range (A-m1).

Prior: relative change x ~ U[-w, w] (the admissible half-width), measurement y = x + n, n ~ N(0, sigma^2) (the
response is linear in tau_R at these changes, A5). Everything depends on r = sigma / w only. With w = 1:
  p_y(y) = [Phi((y + 1)/r) - Phi((y - 1)/r)] / 2,
  I(x; y) = h(y) - (1/2) log2(2 pi e r^2)            [bits],
  E[Var(x | y)] = E_y[ variance of N(y, r^2) truncated to [-1, 1] ],   prior variance 1/3,
  prior variance removed = 1 - 3 E[Var(x | y)].
Method 'quad': composite Gauss-Legendre in y (64 nodes per panel, panels refined at y = +-1 on the scale r),
analytic truncated-normal variance in a cancellation-free form (erfcx). Method 'mc': Monte Carlo with a fixed seed,
independent of the quadrature (agreement within the MC standard error at r = 0.1 ... 10). Large-r limit:
I -> 1/(6 ln2 r^2), removed -> 1/(3 r^2), used for r > 100 (relative error < 4e-5 there).
"""
from __future__ import annotations

import numpy as np
from scipy import special

LN2 = np.log(2.0)
SQ2 = np.sqrt(2.0)
SQ2PI = np.sqrt(2.0 * np.pi)


def tn_var_std(a, b):
    """Variance of the standard normal truncated to [a, b] (vectorised, a < b), stable in both tails."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    flip = a > 0                                     # both positive -> reflect to both negative
    a, b = np.where(flip, -b, a), np.where(flip, -a, b)
    out = np.empty_like(a)
    neg = b <= 0                                     # both <= 0: erfcx form, factor exp(-b^2/2)
    if np.any(neg):
        aa, bb = a[neg], b[neg]
        ea = np.exp(-(aa * aa - bb * bb) / 2)        # <= 1
        den = special.erfcx(-bb / SQ2) - special.erfcx(-aa / SQ2) * ea          # = 2 Z exp(b^2/2)
        pa = 2.0 / SQ2PI * ea / den                  # phi(a)/Z
        pb = 2.0 / SQ2PI / den                       # phi(b)/Z
        out[neg] = 1 + aa * pa - bb * pb - (pa - pb) ** 2
    st = ~neg                                        # straddles 0: direct form
    if np.any(st):
        aa, bb = a[st], b[st]
        Z = special.ndtr(bb) - special.ndtr(aa)
        pa = np.exp(-aa * aa / 2) / SQ2PI / Z
        pb = np.exp(-bb * bb / 2) / SQ2PI / Z
        out[st] = 1 + aa * pa - bb * pb - (pa - pb) ** 2
    return np.maximum(out, 0.0)


def _py(y, r):
    return 0.5 * (special.ndtr((y + 1) / r) - special.ndtr((y - 1) / r))


def _nodes(r, nper=64):
    lim = 1 + 14 * r
    br = {0.0, lim}
    for k in (-12, -8, -6, -4, -3, -2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12):
        y = 1 + k * r
        if 0 < y < lim:
            br.add(y)
    br = np.array(sorted(br))
    # subdivide long panels
    edges = [br[0]]
    for lo, hi in zip(br[:-1], br[1:]):
        m = int(np.ceil((hi - lo) / max(min(r, 1.0), 1e-3) / 4)) if hi - lo > 0 else 1
        m = max(1, min(m, 200))
        edges.extend(np.linspace(lo, hi, m + 1)[1:])
    edges = np.array(edges)
    x, w = np.polynomial.legendre.leggauss(nper)
    a, b = edges[:-1, None], edges[1:, None]
    yy = (0.5 * (b - a) * x + 0.5 * (a + b)).ravel()
    ww = (0.5 * (b - a) * w).ravel()
    return yy, ww


def info_quad(r):
    r = float(r)
    yy, ww = _nodes(r)
    p = _py(yy, r)
    with np.errstate(divide="ignore", invalid="ignore"):
        hint = np.where(p > 0, -p * np.log2(p), 0.0)
    hy = 2 * float(np.sum(ww * hint))                # symmetric in y
    mi = hy - 0.5 * np.log2(2 * np.pi * np.e * r * r)
    v = r * r * tn_var_std((-1 - yy) / r, (1 - yy) / r)
    ev = 2 * float(np.sum(ww * p * v))
    norm = 2 * float(np.sum(ww * p))
    return {"bits": mi, "var_removed": 1 - 3 * ev, "norm_check": norm}


def info_mc(r, n=400000, seed=20261112):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1, 1, n)
    y = x + r * rng.standard_normal(n)
    lpyx = -0.5 * ((y - x) / r) ** 2 - np.log(r * SQ2PI)
    lpy = np.log(_py(y, r))
    d = (lpyx - lpy) / LN2
    ev = r * r * tn_var_std((-1 - y) / r, (1 - y) / r)
    return {"bits": float(np.mean(d)), "bits_se": float(np.std(d) / np.sqrt(n)),
            "var_removed": float(1 - 3 * np.mean(ev)), "var_removed_se": float(3 * np.std(ev) / np.sqrt(n))}


def info(r):
    """bits and prior variance removed for sigma/w = r (quadrature; r > 100 by the large-r limit, rel. error < 4e-5)."""
    r = float(r)
    if r > 100:
        return {"bits": 1.0 / (6 * LN2 * r * r), "var_removed": 1.0 / (3 * r * r)}
    q = info_quad(r)
    return {"bits": q["bits"], "var_removed": q["var_removed"]}
