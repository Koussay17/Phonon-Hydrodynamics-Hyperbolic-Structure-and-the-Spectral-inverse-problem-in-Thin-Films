"""Task 3. Gauge identity of the physical parametrisation.

Physical parameters q = (e, a, lambda, ell, tau_R), rho c = lambda / a eliminated, and
   xi1 = e / sqrt(a),  b = sqrt(lambda rho c) = lambda / sqrt(a),  tau_l = (1 + alpha) ell^2 / a
(historical GK: 1 + alpha = 3).
(G1) the Euler operator D = e d_e + 2a d_a + lambda d_lambda + ell d_ell annihilates each of
     xi1, b, tau_R, tau_l, hence Z (checked on the full composite expression, symbolically and at
     50 digits on p = i omega); the operator without the ell term does not; finite form mu.
(G2) Jacobian of phi: q -> theta (4 x 5): rank 4, kernel spanned by (e, 2a, lambda, ell, 0).
     The same with (e, lambda, rho c, ell, tau_R): kernel (e, lambda, -rho c, ell, 0).
(G3) kernel of the composite map q -> Z(i omega_j): dimension 1 off Sigma (= gauge), 2 on
     {tau_R = tau_l}, 3 on {tau_R = tau_l, b = bs}  (50-digit SVD, 8 frequencies).
(G4) dimensional analysis (Pint for the dimensions, SymPy for the linear algebra): the unit
     rescalings (m, kg, s, K) -> (L, M, T, Theta) that leave the data (omega, bs, Z) unchanged
     form a 2-parameter group {T = 1, Theta = M}; on the parameters, M acts trivially and L acts
     as the gauge.  No other dimensional direction acts on the parameters.
Output: results/03_gauge.txt
"""
import random
import sys
from pathlib import Path

import mpmath as mp
import pint
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


# ------------------------------------------------------------------ (G1)
e, a, lam, ell, tR, mu, k3 = sp.symbols("e a lambda ell tau_R mu kappa", positive=True)
Z, s, parts = M.sympy_model(with_R=False)
p, xi1, b, bs, tRs, tl = (s[k] for k in ("p", "xi1", "b", "bs", "tau_R", "tau_l"))
theta_phys = {xi1: e / sp.sqrt(a), b: lam / sp.sqrt(a), tRs: tR, tl: k3 * ell**2 / a}
Dop = lambda f: e * sp.diff(f, e) + 2 * a * sp.diff(f, a) + lam * sp.diff(f, lam) + ell * sp.diff(f, ell)
for name, f in [("xi1", e / sp.sqrt(a)), ("b", lam / sp.sqrt(a)), ("tau_R", tR), ("tau_l", k3 * ell**2 / a)]:
    log(f"(G1) D[{name}] =", sp.simplify(Dop(f)))
Zq = Z.subs(theta_phys)
DZ = Dop(Zq)
DZ_noell = DZ - ell * sp.diff(Zq, ell)
f_DZ = sp.lambdify((p, e, a, lam, ell, tR, k3, bs), [DZ, DZ_noell, e * sp.diff(Zq, e)], modules="mpmath")
mp.mp.dps = 50
rng = random.Random(3)
worst, worst_noell = mp.mpf(0), mp.mpf(10)
for _ in range(25):
    v = [mp.mpf(rng.uniform(0.1, 3)) for _ in range(5)]
    pv = mp.mpc(0, rng.uniform(0.01, 50))
    r, rno, ref = f_DZ(pv, *v, 3, mp.mpf(rng.uniform(0.2, 3)))
    worst = max(worst, abs(r) / abs(ref))
    worst_noell = min(worst_noell, abs(rno) / abs(ref))
log(f"(G1) |D Z| / |e dZ/de| on p = i omega, 25 random points, 50 digits: max {mp.nstr(worst, 3)}")
log(f"     operator without the ell term: min |residual| / |e dZ/de| = {mp.nstr(worst_noell, 3)} (nonzero)")
# symbolic: D Z expressed through the chain rule, and the finite transformation
Zmu = Zq.subs({e: mu * e, a: mu**2 * a, lam: mu * lam, ell: mu * ell}, simultaneous=True)
th_mu = [sp.simplify(x.subs({e: mu * e, a: mu**2 * a, lam: mu * lam, ell: mu * ell}, simultaneous=True) - x)
         for x in (e / sp.sqrt(a), lam / sp.sqrt(a), k3 * ell**2 / a)]
log("(G1) finite transformation (e, a, lambda, ell) -> (mu e, mu^2 a, mu lambda, mu ell): changes of (xi1, b, tau_l) =", th_mu)
f_mu = sp.lambdify((p, e, a, lam, ell, tR, k3, bs, mu), Zmu - Zq, modules="mpmath")
wmu = max(abs(f_mu(mp.mpc(0, rng.uniform(0.01, 50)), *[mp.mpf(rng.uniform(0.1, 3)) for _ in range(5)], 3,
                   mp.mpf(rng.uniform(0.2, 3)), mp.mpf(rng.uniform(0.1, 10)))) for _ in range(20))
log(f"(G1) max |Z(mu-transformed) - Z| at 20 random points (50 digits): {mp.nstr(wmu, 3)}")

# ------------------------------------------------------------------ (G2)
q = [e, a, lam, ell, tR]
phi = sp.Matrix([e / sp.sqrt(a), lam / sp.sqrt(a), tR, k3 * ell**2 / a])
Jphi = phi.jacobian(q)
log("(G2) rank of d phi / d(e, a, lambda, ell, tau_R):", Jphi.rank(simplify=True))
ns = Jphi.nullspace(simplify=True)
log("     kernel:", [sp.simplify(v.T * sp.diag(*[1 / x for x in q[:4]], 1)) for v in ns], "(rescaled by 1/q: must be proportional to (1,2,1,1,0))")
log("     kernel vector, raw:", [list(sp.simplify(v)) for v in ns])
rc = sp.Symbol("rho_c", positive=True)
phi2 = sp.Matrix([e * sp.sqrt(rc / lam), sp.sqrt(lam * rc), tR, k3 * ell**2 * rc / lam])
J2 = phi2.jacobian([e, lam, rc, ell, tR])
ns2 = J2.nullspace(simplify=True)
log("(G2) with (e, lambda, rho c, ell, tau_R): rank", J2.rank(simplify=True), ", kernel", [list(sp.simplify(v)) for v in ns2])

# ------------------------------------------------------------------ (G3)
jth = M.sympy_jacobian_mp(with_R=False)


def composite_sv(ev, av, lv, ellv, tRv, bsv, oms, kap=3):
    Jp = sp.lambdify(q, Jphi.subs(k3, kap), modules="mpmath")(ev, av, lv, ellv, tRv)
    Jp = mp.matrix(Jp.tolist()) if hasattr(Jp, "tolist") else mp.matrix(Jp)
    x1, bb, tlv = ev / mp.sqrt(av), lv / mp.sqrt(av), kap * ellv ** 2 / av
    rows = []
    for om in oms:
        g = jth(mp.mpc(0, om), x1, bb, tRv, tlv, bsv, 0)
        rows.append([sum(g[i] * Jp[i, j] for i in range(4)) for j in range(5)])
    Jr = M.real_jacobian(rows)
    # log-scaled columns of the physical parameters
    for j, val in enumerate((ev, av, lv, ellv, tRv)):
        for i in range(Jr.rows):
            Jr[i, j] *= val
    sv = M.singular_values(Jr)
    return [x / sv[0] for x in sv]


oms = [mp.mpf(10) ** (-1 + 3 * k / 7) for k in range(8)]
cases = []
for _ in range(3):
    cases.append(("off Sigma", [mp.mpf(rng.uniform(0.3, 2)) for _ in range(5)], mp.mpf(rng.uniform(0.3, 3))))
log("(G3) normalised singular values of the composite Jacobian d Z(i omega_j)/d ln(e, a, lambda, ell, tau_R), 8 frequencies, 50 digits")
for name, (ev, av, lv, ellv, tRv), bsv in cases:
    sv = composite_sv(ev, av, lv, ellv, tRv, bsv, oms)
    log(f"   {name}: [{', '.join(mp.nstr(x, 3) for x in sv)}]")
# on Sigma: tau_R = 3 ell^2 / a ; and b = bs additionally
for _ in range(2):
    ev, av, lv, ellv = [mp.mpf(rng.uniform(0.3, 2)) for _ in range(4)]
    tRv = 3 * ellv ** 2 / av
    bsv = mp.mpf(rng.uniform(0.3, 3))
    sv = composite_sv(ev, av, lv, ellv, tRv, bsv, oms)
    log(f"   on tau_R = tau_l: [{', '.join(mp.nstr(x, 3) for x in sv)}]")
    sv = composite_sv(ev, av, lv, ellv, tRv, lv / mp.sqrt(av), oms)
    log(f"   on tau_R = tau_l, b = bs: [{', '.join(mp.nstr(x, 3) for x in sv)}]")

# ------------------------------------------------------------------ (G4)
ureg = pint.UnitRegistry()
quant = {
    "e": ureg.meter, "a": ureg.meter ** 2 / ureg.second, "lambda": ureg.watt / (ureg.meter * ureg.kelvin),
    "rho_c": ureg.joule / (ureg.meter ** 3 * ureg.kelvin), "ell": ureg.meter, "tau_R": ureg.second,
    "b_s": ureg.watt * ureg.second ** 0.5 / (ureg.meter ** 2 * ureg.kelvin), "omega": 1 / ureg.second,
    "Z": ureg.kelvin * ureg.meter ** 2 / ureg.watt,
}
base = ["[length]", "[mass]", "[time]", "[temperature]"]
expo = {}
for k, u in quant.items():
    dim = (1 * u).to_base_units().dimensionality
    expo[k] = [sp.nsimplify(dim.get(d, 0)) for d in base]
    log(f"(G4) {k:6s} exponents (m, kg, s, K): {expo[k]}")
data = ["omega", "b_s", "Z"]
A_data = sp.Matrix([expo[k] for k in data])
ker = A_data.nullspace()
log("(G4) unit rescalings (ln L, ln M, ln T, ln Theta) leaving omega, b_s, Z unchanged: basis", [list(v) for v in ker])
params = ["e", "a", "lambda", "rho_c", "ell", "tau_R"]
A_par = sp.Matrix([expo[k] for k in params])
for v in ker:
    log(f"     action on ln(e, a, lambda, rho_c, ell, tau_R): {list(A_par * v)}")
log("     -> one direction acts as (1, 2, 1, -1, 1, 0) = the gauge; the other acts trivially.")

(HERE / "results" / "03_gauge.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
