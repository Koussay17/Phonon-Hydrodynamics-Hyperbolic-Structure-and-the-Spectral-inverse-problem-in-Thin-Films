"""Independent checks of Sec. IV (front-face model): Z(p) from a direct solution of the two-layer
boundary-value problem (not from the quadrupole), branch independence, gauge invariance (Prop. 6),
Corollary 2, the half-space phase (40) and the second-order entropy production.

Direct BVP (Laplace domain, rest state): film 0<z<e: q'' = kappa^2 q with
kappa^2 = (1+tau_R p)/(lambda/(rho c p) + L^2), theta = -q'/(rho c p); substrate: theta_s = T exp(-s (z-e)),
s = sqrt(p/a_s), q_s = lambda_s s theta_s. Conditions: q(0) = 1, theta and q continuous at z = e.
Usage: python -B audit_frontface.py > audit_frontface.out.txt
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


def Z_direct(p, e, lam, rc, tR, L2, lam_s, rc_s, branch=1):
    kap = branch * mp.sqrt((1 + tR * p) / (lam / (rc * p) + L2))
    s = mp.sqrt(p * rc_s / lam_s)  # Re s > 0 (principal root, Re p > 0)
    # q = c1 cosh(kap z) + c2 sinh(kap z); theta = -q'/(rho c p)
    # unknowns c1, c2, Tsub: q(0)=1 -> c1 = 1
    c1 = 1
    # continuity of theta: -(c1 kap sinh(kap e) + c2 kap cosh(kap e))/(rc p) = Tsub
    # continuity of q:      c1 cosh(kap e) + c2 sinh(kap e) = lam_s s Tsub
    A = mp.matrix([[-kap * mp.cosh(kap * e) / (rc * p), -1], [mp.sinh(kap * e), -lam_s * s]])
    b = mp.matrix([c1 * kap * mp.sinh(kap * e) / (rc * p), -c1 * mp.cosh(kap * e)])
    c2, Tsub = mp.lu_solve(A, b)
    theta0 = -(c2 * kap) / (rc * p)
    return theta0  # q(0) = 1


def Z_formula(p, e, lam, rc, tR, L2, lam_s, rc_s):
    a = lam / rc; b = mp.sqrt(lam * rc); bs = mp.sqrt(lam_s * rc_s)
    xi1 = e / mp.sqrt(a); tl = L2 / a
    sig_e = xi1 * mp.sqrt(p * (1 + tR * p) / (1 + tl * p))
    lsig = b * mp.sqrt(p * (1 + tl * p) / (1 + tR * p))
    rho = lsig / (bs * mp.sqrt(p))
    t = mp.tanh(sig_e)
    return (rho + t) / (lsig * (1 + rho * t))


pars = dict(e=mp.mpf("0.7"), lam=mp.mpf("2.3"), rc=mp.mpf("1.9"), tR=mp.mpf("0.4"), L2=mp.mpf("0.15"),
            lam_s=mp.mpf("1.1"), rc_s=mp.mpf("2.6"))
ps = [mp.mpc(0.3, 1.7), mp.mpc(2.0, 0.1), mp.mpc(0.01, 40), mp.mpc(5, -3), mp.mpc(0, 0.9) + mp.mpf(10) ** -30]
err = max(abs(Z_direct(p, **pars) - Z_formula(p, **pars)) / abs(Z_formula(p, **pars)) for p in ps)
check("Eq. (35)-(36): Z from the direct two-layer BVP equals the quotient form (5 complex p, Re p > 0)", err < 1e-30, f"max rel. diff {mp.nstr(err, 3)}")
errb = max(abs(Z_direct(p, branch=-1, **pars) - Z_direct(p, **pars)) / abs(Z_direct(p, **pars)) for p in ps)
check("branch independence: sigma -> -sigma leaves Z unchanged", errb < 1e-30, mp.nstr(errb, 3))
# gauge (Prop. 6): e -> mu e, lam -> mu lam, rho c -> rho c/mu, ell -> mu ell (L^2 -> mu^2 L^2), tau_R fixed
mu = mp.mpf("3.7")
g = dict(pars); g.update(e=mu * pars["e"], lam=mu * pars["lam"], rc=pars["rc"] / mu, L2=mu ** 2 * pars["L2"])
errg = max(abs(Z_direct(p, **g) - Z_direct(p, **pars)) / abs(Z_direct(p, **pars)) for p in ps)
check("Prop. 6(ii): Z invariant under the depth-scaling gauge G_mu (mu = 3.7), direct BVP", errg < 1e-30, mp.nstr(errg, 3))
# Corollary 2: (beta, ell) on (1+beta) ell^2 = const give the same Z; ell ratio 3/2 between beta = 1/3 and 2
L2c = pars["L2"]
ell_13 = mp.sqrt(L2c / (1 + mp.mpf(1) / 3)); ell_2 = mp.sqrt(L2c / 3)
check("Cor. 2: ell_{beta=1/3}/ell_{beta=2} = 3/2 at fixed tau_l", abs(ell_13 / ell_2 - mp.mpf(3) / 2) < 1e-35)
# half-space phase (40): thick film, branch Re sigma > 0
tl = pars["L2"] / (pars["lam"] / pars["rc"])
okph = True
for w in (mp.mpf("0.05"), mp.mpf(1), mp.mpf(30)):
    p = mp.mpc(0, w)
    # half-space film: decaying solution q = exp(-kappa z), Re kappa > 0, so theta(0) = kappa/(rho c p)
    kap = mp.sqrt((1 + pars["tR"] * p) / (pars["lam"] / (pars["rc"] * p) + pars["L2"]))
    if mp.re(kap) < 0:
        kap = -kap
    Zt = kap / (pars["rc"] * p)
    ph = mp.arg(Zt)
    pred = -mp.pi / 4 + (mp.atan(w * pars["tR"]) - mp.atan(w * tl)) / 2
    okph = okph and abs(ph - pred) < 1e-25
check("Eq. (40): thick-film phase = -pi/4 + [arctan(w tau_R) - arctan(w tau_l)]/2 (direct half-space solution, decaying branch)", okph)
# entropy production to second order (symbolic): sigma_s = d_t s + d_z J_s with the GK law
z, t = sp.symbols("z t")
T0, lamS, tRS, L2S, rcS = sp.symbols("T_0 lambda tau_R L2 rho_c", positive=True)
q = sp.Function("q")(z, t); th = sp.Function("theta")(z, t)
T = T0 + th
qt = (-q - lamS * sp.diff(th, z) + L2S * sp.diff(q, z, 2)) / tRS          # GK law solved for d_t q
et = -sp.diff(q, z)                                                        # energy balance d_t e = -d_z q
ds = et / T - tRS * q * qt / (lamS * T0 ** 2)                               # d_t s with s = s_eq(e) - tau_R q^2/(2 lam T0^2)
Js = q / T + L2S * q * sp.diff(q, z) / (lamS * T0 ** 2)
sigma = sp.simplify(ds + sp.diff(Js, z))
target = (q ** 2 + L2S * sp.diff(q, z) ** 2) / (lamS * T0 ** 2)
remainder = sp.simplify(sigma - target)
# remainder must be third order: q * d_z theta * (1/T0^2 - 1/T^2)
check("entropy production = (q^2 + L^2 q_z^2)/(lambda T0^2) + third-order terms",
      sp.simplify(remainder - q * sp.diff(th, z) * (1 / T0 ** 2 - 1 / T ** 2)) == 0, str(sp.factor(remainder)))

print()
print(f"{sum(res)} of {len(res)} checks passed")
