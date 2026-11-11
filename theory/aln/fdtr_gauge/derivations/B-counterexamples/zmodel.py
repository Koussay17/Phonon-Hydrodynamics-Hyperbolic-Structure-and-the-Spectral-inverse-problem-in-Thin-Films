"""Front-face impedance of a Cattaneo/Guyer-Krumhansl film on a semi-infinite Fourier substrate.

Model (note 7, src/forward_model.response, note 16 eq. for Z_t):
    A(p) = xi1 sqrt(p (1 + tR p) / (1 + tl p))          (sigma e)
    B(p) = b   sqrt(p (1 + tl p) / (1 + tR p))          (lambda_eff sigma)
    Zs   = R + 1/(bs sqrt p)                            (interface resistance + substrate)
    Z    = (Zs + tanh(A)/B) / (1 + B Zs tanh(A))        variant "code"  (transmission-line form)
         = (1/B) (1 + rho E)/(1 - rho E),  E = exp(-2A), rho = (y - 1)/(y + 1), y = B Zs.

The brief of branch B-X writes Z = (1/B)(B tanh A + bs sqrt p)/(B + bs sqrt p tanh A); this is the
reciprocal ratio, i.e. y -> 1/y (rho -> -rho).  It fails the limit e -> 0 (Z -> bs sqrt(p)/B^2
instead of 1/(bs sqrt p)).  It is kept as variant "brief" only to show which conclusions depend on it.

Parameters: theta = (xi1, b, tR, tl) and optionally R.  bs > 0 is fixed and known.
Raw complex derivatives returned in the order
    dZ/dln(xi1), dZ/dln(b), dZ/dtR, dZ/dtl, dZ/dR
Structure used throughout:  with u = A dZ/dA, v = B dZ/dB (A fixed),
    dZ/dln xi1 = u,  dZ/dln b = v,  dZ/dtR = gR (u - v),  dZ/dtl = -gl (u - v),  g = p / (2 (1 + t p)).
"""
from __future__ import annotations

import numpy as np
import mpmath as mp

# --------------------------------------------------------------------------------------------
# numpy (double precision) implementation, vectorised over omega
# --------------------------------------------------------------------------------------------

def pieces_np(w, xi1, b, tR, tl, bs=1.0, R=0.0, variant="code"):
    w = np.asarray(w, dtype=float)
    p = 1j * w
    s = np.sqrt(p)
    num = 1.0 + tR * p
    den = 1.0 + tl * p
    r = np.sqrt(num / den)            # arg in (-pi/4, pi/4) for omega > 0
    A = xi1 * s * r
    B = b * s / r
    E = np.exp(-2.0 * A)
    one_m_E = -np.expm1(-2.0 * A)
    y = B * (R + 1.0 / (bs * s))
    if variant == "code":
        yy = y
        dy = yy                        # B d(yy)/dB
    elif variant == "brief":
        if R != 0.0:
            raise ValueError("variant 'brief' has no interface resistance")
        yy = 1.0 / y
        dy = -yy
    else:
        raise ValueError(variant)
    rho = (yy - 1.0) / (yy + 1.0)
    return dict(p=p, s=s, A=A, B=B, E=E, one_m_E=one_m_E, yy=yy, dy=dy, rho=rho)


def Z_np(w, xi1, b, tR, tl, bs=1.0, R=0.0, variant="code"):
    d = pieces_np(w, xi1, b, tR, tl, bs, R, variant)
    rE = d["rho"] * d["E"]
    return (1.0 + rE) / (d["B"] * (1.0 - rE))


def jac_np(w, xi1, b, tR, tl, bs=1.0, R=0.0, variant="code"):
    """Return Z (M,) and raw complex derivatives (M, 5): dln xi1, dln b, dtR, dtl, dR."""
    d = pieces_np(w, xi1, b, tR, tl, bs, R, variant)
    A, B, E, ome, yy, dy, rho, p = (d[k] for k in ("A", "B", "E", "one_m_E", "yy", "dy", "rho", "p"))
    D1 = 1.0 - rho * E
    Z = (1.0 + rho * E) / (B * D1)
    u = -4.0 * A * rho * E / (B * D1 ** 2)
    if variant == "code":
        # cancellation-free form of v = -Z + 4 E y/(B D1^2 (y+1)^2)
        v = -ome * ((1.0 + E) * (1.0 + yy ** 2) + 2.0 * yy * ome) / (B * D1 ** 2 * (yy + 1.0) ** 2)
    else:
        v = -Z + 4.0 * E * dy / (B * D1 ** 2 * (yy + 1.0) ** 2)
    gR = p / (2.0 * (1.0 + tR * p))
    gl = p / (2.0 * (1.0 + tl * p))
    J = np.empty(w.shape + (5,), dtype=complex)
    J[..., 0] = u
    J[..., 1] = v
    J[..., 2] = gR * (u - v)
    J[..., 3] = -gl * (u - v)
    if variant == "code":
        J[..., 4] = 4.0 * E / (D1 ** 2 * (yy + 1.0) ** 2)
    else:
        J[..., 4] = 0.0
    return Z, J


def col_scales(xi1, b, tR, tl, R=0.0, tscale="log", t0=None, R0=None):
    """Scale factors turning raw derivatives into the columns actually used.

    tscale = "log": tau d/dtau for tau > 0, t0 d/dtau at tau = 0 (t0 defaults to xi1^2).
    tscale = "lin": t0 d/dtau always.
    R: R d/dR for R > 0, R0 d/dR at R = 0 (R0 defaults to xi1/b, the film resistance e/lambda).
    """
    if t0 is None:
        t0 = xi1 ** 2
    if R0 is None:
        R0 = xi1 / b
    c = np.ones(5)
    for k, t in ((2, tR), (3, tl)):
        if tscale == "log" and t > 0:
            c[k] = t
        else:
            c[k] = t0
    c[4] = R if R > 0 else R0
    return c


def real_rel_jac_np(w, xi1, b, tR, tl, bs=1.0, R=0.0, variant="code", with_R=False,
                    tscale="log", weights=None, t0=None, R0=None):
    """Real Jacobian of (ln|Z|, arg Z) w.r.t. the scaled parameters; rows stacked (Re, Im) per omega."""
    Z, J = jac_np(w, xi1, b, tR, tl, bs, R, variant)
    c = col_scales(xi1, b, tR, tl, R, tscale, t0, R0)
    ncol = 5 if with_R else 4
    Jr = J[:, :ncol] * c[:ncol] / Z[:, None]
    if weights is not None:
        Jr = Jr * np.sqrt(weights)[:, None]
    M = np.empty((2 * len(w), ncol))
    M[0::2] = Jr.real
    M[1::2] = Jr.imag
    return M


def svals_np(M):
    return np.linalg.svd(M, compute_uv=False)


# --------------------------------------------------------------------------------------------
# mpmath (arbitrary precision) implementation
# --------------------------------------------------------------------------------------------

def pieces_mp(w, xi1, b, tR, tl, bs=1, R=0, variant="code"):
    p = mp.mpc(0, w)
    s = mp.sqrt(p)
    r = mp.sqrt((1 + tR * p) / (1 + tl * p))
    A = xi1 * s * r
    B = b * s / r
    E = mp.exp(-2 * A)
    ome = -mp.expm1(-2 * A)
    y = B * (R + 1 / (bs * s))
    if variant == "code":
        yy = y
        dy = yy
    else:
        yy = 1 / y
        dy = -yy
    rho = (yy - 1) / (yy + 1)
    return p, s, A, B, E, ome, yy, dy, rho


def Z_mp(w, xi1, b, tR, tl, bs=1, R=0, variant="code"):
    p, s, A, B, E, ome, yy, dy, rho = pieces_mp(w, xi1, b, tR, tl, bs, R, variant)
    return (1 + rho * E) / (B * (1 - rho * E))


def jac_mp(w, xi1, b, tR, tl, bs=1, R=0, variant="code"):
    p, s, A, B, E, ome, yy, dy, rho = pieces_mp(w, xi1, b, tR, tl, bs, R, variant)
    D1 = 1 - rho * E
    Z = (1 + rho * E) / (B * D1)
    u = -4 * A * rho * E / (B * D1 ** 2)
    if variant == "code":
        v = -ome * ((1 + E) * (1 + yy ** 2) + 2 * yy * ome) / (B * D1 ** 2 * (yy + 1) ** 2)
    else:
        v = -Z + 4 * E * dy / (B * D1 ** 2 * (yy + 1) ** 2)
    gR = p / (2 * (1 + tR * p))
    gl = p / (2 * (1 + tl * p))
    dR = 4 * E / (D1 ** 2 * (yy + 1) ** 2) if variant == "code" else mp.mpc(0)
    return Z, [u, v, gR * (u - v), -gl * (u - v), dR]


def real_rel_jac_mp(ws, xi1, b, tR, tl, bs=1, R=0, variant="code", with_R=False,
                    tscale="log", weights=None, t0=None, R0=None):
    if t0 is None:
        t0 = xi1 ** 2
    if R0 is None:
        R0 = xi1 / b
    c = [mp.mpf(1), mp.mpf(1),
         (tR if (tscale == "log" and tR > 0) else t0),
         (tl if (tscale == "log" and tl > 0) else t0),
         (R if R > 0 else R0)]
    ncol = 5 if with_R else 4
    rows = []
    for i, w in enumerate(ws):
        Z, J = jac_mp(w, xi1, b, tR, tl, bs, R, variant)
        sw = mp.sqrt(weights[i]) if weights is not None else 1
        jr = [J[k] * c[k] / Z * sw for k in range(ncol)]
        rows.append([x.real for x in jr])
        rows.append([x.imag for x in jr])
    return mp.matrix(rows)


def svals_mp(M):
    """Singular values of a tall real mp.matrix via Householder QR then SVD of the small factor."""
    m, n = M.rows, M.cols
    A = M.copy()
    # Householder QR, keep R only
    for k in range(n):
        x = [A[i, k] for i in range(k, m)]
        normx = mp.sqrt(mp.fsum(xi ** 2 for xi in x))
        if normx == 0:
            continue
        alpha = -normx if x[0] >= 0 else normx
        vvec = list(x)
        vvec[0] = vvec[0] - alpha
        vn2 = mp.fsum(vi ** 2 for vi in vvec)
        if vn2 == 0:
            continue
        for j in range(k, n):
            dot = mp.fsum(vvec[i - k] * A[i, j] for i in range(k, m))
            f = 2 * dot / vn2
            for i in range(k, m):
                A[i, j] -= f * vvec[i - k]
    Rm = mp.matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            Rm[i, j] = A[i, j]
    sv = mp.svd_r(Rm, compute_uv=False)
    return sorted([abs(x) for x in sv], reverse=True)


# --------------------------------------------------------------------------------------------
# frequency grids
# --------------------------------------------------------------------------------------------

def char_freqs(xi1, b, tR, tl, bs=1.0):
    """Characteristic frequencies of the problem (inverse times)."""
    beta = b / bs
    f = [1.0 / xi1 ** 2, beta ** 2 / xi1 ** 2, 1.0 / (beta ** 2 * xi1 ** 2)]
    for t in (tR, tl):
        if t > 0:
            f.append(1.0 / t)
    if tR > 0 and tl > 0:
        f.append(1.0 / (xi1 ** 2) * tl / tR)
        f.append(1.0 / (xi1 ** 2) * tR / tl)
    return np.array(f)


def char_freqs2(xi1, b, tR, tl, bs=1.0):
    """Physically relevant characteristic frequencies (inverse times):
    film diffusion 1/xi1^2, thin-film crossover min(beta^2, beta^-2)/xi1^2, 1/tR, 1/tl, and the
    GK thick-film crossover (tl/tR)/xi1^2 (A ~ 1 when omega tau >> 1)."""
    beta = b / bs
    f = [1.0 / xi1 ** 2, min(beta ** 2, beta ** -2) / xi1 ** 2]
    for t in (tR, tl):
        if t > 0:
            f.append(1.0 / t)
    if tR > 0 and tl > 0:
        f.append(tl / tR / xi1 ** 2)
    return np.array(f)


def resolved_grid(xi1, b, tR, tl, bs=1.0, R=0.0, pad=1.0, dlog=0.1, dphi=0.2, nmax=200000,
                  wmin=None, wmax=None, variant="code"):
    """Grid on [min f / 10^pad, max f * 10^pad] (f = char_freqs2) whose consecutive points differ by
    at most dlog in ln(omega) and at most min(dphi, (1 - |rho E|)/3) in Im(2A) (resolves the
    thermal-wave resonances of Cattaneo-like films).  Returns omega, d ln(omega) trapezoid weights,
    and a flag `ok` (False if more than nmax points would be needed)."""
    f = char_freqs2(xi1, b, tR, tl, bs)
    lo = np.log(f.min()) - pad * np.log(10) if wmin is None else np.log(wmin)
    hi = np.log(f.max()) + pad * np.log(10) if wmax is None else np.log(wmax)
    nf = int((hi - lo) / 0.01) + 2
    lw = np.linspace(lo, hi, nf)
    d = pieces_np(np.exp(lw), xi1, b, tR, tl, bs, R, variant)
    phi = np.imag(2.0 * d["A"])          # continuous: arg A in (0, pi/2) on the imaginary axis
    q = np.abs(d["rho"] * d["E"])
    step = np.minimum(dphi, np.maximum((1.0 - q) / 3.0, 1e-6))
    step = np.where(q > 1e-8, step, np.inf)   # oscillatory part below 1e-8 relative: no phase constraint
    dl = np.diff(lw); dp = np.abs(np.diff(phi))
    # local cost: number of steps needed between fine nodes
    cost = np.maximum(dl / dlog, dp / np.minimum(step[:-1], step[1:]))
    s = np.concatenate([[0.0], np.cumsum(cost)])
    n = int(np.ceil(s[-1])) + 1
    ok = n <= nmax
    n = min(n, nmax)
    lsel = np.interp(np.linspace(0, s[-1], n), s, lw)
    w = np.exp(lsel)
    wts = np.zeros(n)
    dls = np.diff(lsel)
    wts[:-1] += dls / 2; wts[1:] += dls / 2
    return w, wts, ok


def wide_grid(xi1, b, tR, tl, bs=1.0, pad=4.0, ppd=24, wmin=None, wmax=None):
    """Log grid covering every characteristic frequency with `pad` decades on each side."""
    f = char_freqs(xi1, b, tR, tl, bs)
    lo = np.log10(f.min()) - pad if wmin is None else np.log10(wmin)
    hi = np.log10(f.max()) + pad if wmax is None else np.log10(wmax)
    n = int(np.ceil((hi - lo) * ppd)) + 1
    w = np.logspace(lo, hi, n)
    wts = np.full(n, (hi - lo) * np.log(10) / (n - 1))   # d ln(omega) quadrature weights
    return w, wts
