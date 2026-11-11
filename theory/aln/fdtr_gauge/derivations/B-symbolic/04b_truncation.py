"""Task 4, truncation tests at full precision (replaces the double-precision coefficient conversion of the
numerical block of 04_expansions.py).

Test point: xi1 = 7/10, b = 23/10, bs = 11/10, tau_R = 3/10, tau_l = 1/8 (R = 0 or 1/5).
For each expansion, the error of the truncation divided by the first omitted term must tend to 1:
  small p     : |Z - sum_{n<=4} z_n q^n| / |z_5 q^5|                    (omega -> 0)
  generic     : |Z - sum_{k<=3} y_k p^(-1/2-k)| / |y_4 p^(-9/2)|        (omega -> infinity)
  Cattaneo    : |Z - sum_{k<=3} C_k(t) u^k| / |C_4(t) u^4|, u = p^(-1/2), t = tanh A exact
  tau_R = 0   : |Z - sum_{k=2..5} c_k u^k| / |c_6 u^6|
Coefficients are computed by SymPy series with the parameters substituted as exact rationals and evaluated at
60 digits; Z is evaluated with mpmath at 50 digits.  Output: results/04b_truncation.txt
"""
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

mp.mp.dps = 50
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


X1, BB, BS, TR, TL = (sp.Rational(7, 10), sp.Rational(23, 10), sp.Rational(11, 10), sp.Rational(3, 10),
                      sp.Rational(1, 8))
x1m, bbm, bsm, trm, tlm = (mp.mpf(7) / 10, mp.mpf(23) / 10, mp.mpf(11) / 10, mp.mpf(3) / 10, mp.mpf(1) / 8)
to_mp = lambda e: mp.mpmathify(str(sp.N(e, 60)).replace("I", "j").replace("*", "").replace(" ", "")) \
    if sp.im(sp.N(e, 60)) != 0 else mp.mpf(str(sp.N(e, 60)))
q, u, t = sp.symbols("q u t")


def Zexpr(A, B, Zs, tt=None):
    tt = sp.tanh(A) if tt is None else tt
    w = B * Zs
    return (w + tt) / (B * (1 + w * tt))


# ------------------------------------------------------------- small p
for Rv in (sp.Integer(0), sp.Rational(1, 5)):
    rho = sp.sqrt((1 + TR * q**2) / (1 + TL * q**2))
    Zq = Zexpr(X1 * q * rho, BB * q / rho, Rv + 1 / (BS * q))
    ser = sp.expand(sp.series(Zq * q, q, 0, 7).removeO())
    z = {n - 1: sp.nsimplify(ser.coeff(q, n)) for n in range(7)}
    log(f"small p, R = {Rv}: z_-1..z_5 = {[z[n] for n in range(-1, 6)]}")
    for om in ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6"):
        pv = mp.mpc(0, mp.mpf(om))
        qv = mp.sqrt(pv)
        exact = M.Z_mp(pv, x1m, bbm, trm, tlm, bsm, mp.mpf(Rv.p) / Rv.q)
        part = sum(mp.mpf(z[n].p) / z[n].q * qv**n for n in range(-1, 5))
        nxt = mp.mpf(z[5].p) / z[5].q * qv**5
        log(f"   omega={om}: |Z - sum_(n<=4)| / |z_5 q^5| = {mp.nstr(abs(exact - part) / abs(nxt), 10)}")

# ------------------------------------------------------------- large p, generic
invB = (1 / BB) * u * sp.sqrt((TR + u**2) / (TL + u**2))
serB = sp.expand(sp.series(invB, u, 0, 11).removeO())
y = [serB.coeff(u, 2 * k + 1) for k in range(5)]
ym = [mp.mpf(str(sp.N(c, 60))) for c in y]
log("large p, generic: y_0..y_4 =", [sp.nsimplify(c) for c in y])
for om in ("1e3", "1e4", "1e5", "1e6"):
    pv = mp.mpc(0, mp.mpf(om))
    exact = M.Z_mp(pv, x1m, bbm, trm, tlm, bsm)
    part = sum(ym[k] * pv ** (-mp.mpf(1) / 2 - k) for k in range(4))
    nxt = ym[4] * pv ** (-mp.mpf(9) / 2)
    A_, B_, _ = M.parts_mp(pv, x1m, bbm, trm, tlm)
    log(f"   omega={om}: |Z - sum_(k<=3)| / |y_4 p^-9/2| = {mp.nstr(abs(exact - part) / abs(nxt), 10)};"
        f"  |Z - 1/B| / |Z| = {mp.nstr(abs(exact - 1 / B_) / abs(exact), 3)};  exp(-2 Re A) = {mp.nstr(mp.exp(-2 * mp.re(A_)), 3)}")

# ------------------------------------------------------------- large p, Cattaneo (tau_l = 0)
BC = BB / sp.sqrt(TR + u**2)                 # b sqrt(p / (1 + tau_R p)), p = 1/u^2
serC = sp.expand(sp.series(Zexpr(None, BC, u / BS, tt=t), u, 0, 6).removeO())
Ck = [sp.factor(serC.coeff(u, k)) for k in range(5)]
log("large p, Cattaneo: C_0..C_4 (t = tanh A) =", Ck)
fC = [sp.lambdify(t, c, modules="mpmath") for c in Ck]
for om in ("1e3", "1e4", "1e5", "1e6", "1e7"):
    pv = mp.mpc(0, mp.mpf(om))
    uv = 1 / mp.sqrt(pv)
    A_, B_, _ = M.parts_mp(pv, x1m, bbm, trm, 0)
    tv = mp.tanh(A_)
    exact = M.Z_mp(pv, x1m, bbm, trm, 0, bsm)
    part = sum(fC[k](tv) * uv**k for k in range(4))
    nxt = fC[4](tv) * uv**4
    log(f"   omega={om}: |t| = {mp.nstr(abs(tv), 5)}, |Z - sum_(k<=3)| / |C_4(t) u^4| = {mp.nstr(abs(exact - part) / abs(nxt), 10)}")

# ------------------------------------------------------------- large p, tau_R = 0
BN = BB * sp.sqrt(1 + TL / u**2) / u
AN = X1 / sp.sqrt(TL + u**2)
serN = sp.expand(sp.series(Zexpr(AN, BN, u / BS), u, 0, 8).removeO())
cN = {k: serN.coeff(u, k) for k in range(2, 7)}
cNm = {k: mp.mpf(str(sp.N(v, 60))) for k, v in cN.items()}
log("large p, tau_R = 0: c_2..c_6 =", {k: mp.nstr(v, 20) for k, v in cNm.items()})
for om in ("1e3", "1e4", "1e5", "1e6"):
    pv = mp.mpc(0, mp.mpf(om))
    uv = 1 / mp.sqrt(pv)
    exact = M.Z_mp(pv, x1m, bbm, 0, tlm, bsm)
    part = sum(cNm[k] * uv**k for k in range(2, 6))
    nxt = cNm[6] * uv**6
    log(f"   omega={om}: |Z - sum_(k=2..5)| / |c_6 u^6| = {mp.nstr(abs(exact - part) / abs(nxt), 10)}")
(HERE / "results" / "04b_truncation.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
