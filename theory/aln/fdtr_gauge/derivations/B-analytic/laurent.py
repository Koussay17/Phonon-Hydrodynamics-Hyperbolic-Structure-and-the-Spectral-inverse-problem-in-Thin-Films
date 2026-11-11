"""Low-frequency Laurent coefficients of the front-face response on the s-plane (p = s^2).

F(s) = (C L + X)/(Y L + C),  L = R + 1/(b_s s),  C = cosh sqrt z, X = (xi/b) g S(z), Y = xi b s^2 S(z),
z = xi^2 s^2 g(s^2), S(w) = sinh(sqrt w)/sqrt w.  Multiplying numerator and denominator by b_s s:
    s F = Nh / D1,   Nh = C (1 + R b_s s) + X b_s s,   D1 = b_s C + xi b s S(z) (1 + R b_s s),
with D1(0) = b_s != 0, so the coefficients follow from truncated power-series division.
Returns [c_{-1}, c_0, c_1, ...] with F = c_{-1}/s + sum_k c_k s^k.
"""
import sympy as sp


def _trunc(expr, s, n):
    """Keep powers s^0..s^(n-1) of a polynomial in s."""
    poly = sp.Poly(sp.expand(expr), s)
    return sum(c*s**m[0] for m, c in zip(poly.monoms(), poly.coeffs()) if m[0] < n)


def laurent_coefficients(xi, b, bs, tR, tl, R, s, n):
    """First n coefficients q_0..q_{n-1} of s F(s); c_{k-1} = q_k."""
    # g(s^2) = (1 + tR s^2) * sum_m (-tl s^2)^m
    ginv = sum((-tl*s**2)**m for m in range(n//2 + 1))
    g = _trunc((1 + tR*s**2)*ginv, s, n)
    z = _trunc(xi**2*s**2*g, s, n)
    C = 0; Sz = 0; zp = 1
    for m in range(n//2 + 2):
        C += zp/sp.factorial(2*m)
        Sz += zp/sp.factorial(2*m + 1)
        zp = _trunc(zp*z, s, n)
    C = _trunc(C, s, n); Sz = _trunc(Sz, s, n)
    X = _trunc(xi/b*g*Sz, s, n)
    Nh = _trunc(C*(1 + R*bs*s) + X*bs*s, s, n)
    D1 = _trunc(bs*C + xi*b*s*Sz*(1 + R*bs*s), s, n)
    nk = [sp.expand(Nh).coeff(s, k) for k in range(n)]
    dk = [sp.expand(D1).coeff(s, k) for k in range(n)]
    q = []
    for k in range(n):
        acc = nk[k] - sum(dk[j]*q[k - j] for j in range(1, k + 1))
        q.append(sp.expand(sp.cancel(acc/dk[0])))
    return q
