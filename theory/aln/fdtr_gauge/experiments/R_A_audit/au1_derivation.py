"""AU1: independent check of the first-order (k^2, p) expansion behind tau_l = |e|^2 (u.C^-1 u) / K^2, u = Q V x.

Three routes, none of which uses ra_kinetic.py:
 A. Exact rational arithmetic (SymPy) for the grey isotropic Callaway operator in an orthonormal Legendre basis
    psi_l = sqrt(2l+1) P_l(mu) (inner product <f,g> = (1/2) int f g dmu, so |e|^2 = 1):
      C = diag(0, 1/tau_R, 1/tau_c, 1/tau_c, ...),  V = v J (J the Legendre recurrence matrix),  b = V e.
    The exact static kernel kappa(k) = b.[C + i k Q V Q]^-1 b on e-perp (l >= 1, tridiagonal) is expanded in k
    symbolically; the GK time follows from L^2 = -[k^2] kappa / K and tau_l = L^2 / a, a = K / |e|^2.
    Also the "unprojected" variant (the energy component of V x relaxed at 1/tau_c as if it were not conserved),
    and the exact dispersion relation p + k^2 kappa(k, p) = 0 expanded to first order in p (lam_eff slope).
 B. Random non-grey parity-symmetric operators with an exact energy invariant (and optionally a momentum mode
    conserved by the normal part): the exact driven kernel from the bordered system
      [[p + i k V + C, -e], [e^T, 0]] [f; s] = [0; |e|^2],   kappa = (b.f) / (-i k),
    i.e. a source along e only, with the energy temperature fixed to 1 (no projection, no regularisation
    constant); first derivatives in p and k^2 by Richardson; compared with the closed formula evaluated through a
    KKT solve [[C, e], [e^T, 0]]. Then the diffusive root of p |e|^2 + k^2 kappa(k, p) = 0 by Newton and the
    slope of lam_eff = -p |e|^2 / k^2 against delta = tau_R - tau_l.
 C. Convention of the twin: src/twin_fdtr.py Material.tau_l = (1 + alpha) ell^2 C / lam, and the twin's k = 0
    kernel of a single GK slab on a Fourier substrate against an independent lam_eff quadrupole written here.
Output: out/au1_derivation.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
res = {}

# ============================================================== A. grey Callaway, exact rational arithmetic
tR, tN, k, p = sp.symbols("tau_R tau_N k p", positive=True)
v = sp.Integer(1)                # the speed scales out (K ~ v^2, L^2 ~ v^2)
tc = 1 / (1 / tR + 1 / tN)
LMAX = 4                         # e-perp block l = 1..4: exact through order k^6 for the tridiagonal chain


def J_entry(l):            # <psi_{l+1}, mu psi_l>
    return sp.Integer(l + 1) / sp.sqrt(sp.Integer((2 * l + 1) * (2 * l + 3)))


n = LMAX + 1
Jm = sp.zeros(n, n)
for l in range(n - 1):
    Jm[l, l + 1] = Jm[l + 1, l] = J_entry(l)
Vm = v * Jm
Cdiag = [0, 1 / tR] + [1 / tR + 1 / tN] * (n - 2)
Cm = sp.diag(*Cdiag)
e = sp.zeros(n, 1)
e[0] = 1
b = Vm * e                       # = v/sqrt(3) psi_1
# restriction to e-perp = indices 1..LMAX
Cr = Cm[1:, 1:]
Vr = Vm[1:, 1:]                  # Q V Q on e-perp (Q removes psi_0)
br = b[1:, :]
x = Cr.inv() * br
K = sp.simplify((br.T * x)[0])
tauR_kin = sp.simplify((x.T * x)[0] / K)
u = Vr * x                       # Q V x restricted (psi_0 component removed)
M = sp.simplify((u.T * Cr.inv() * u)[0])
a = K / 1                        # |e|^2 = 1
taul = sp.simplify(M / K / a)
res["A_grey"] = {"K": str(K), "tauR_kin": str(tauR_kin), "taul": str(taul),
                 "taul_over_tauc": str(sp.simplify(taul / tc))}
# unprojected variant: psi_0 component of V x kept and relaxed at 1/tau_c
Vx_full = Vm * sp.Matrix.vstack(sp.zeros(1, 1), x)
Chist = sp.diag(*([1 / tR + 1 / tN] + Cdiag[1:]))
M_hist = sp.simplify((Vx_full.T * Chist.inv() * Vx_full)[0])
res["A_grey"]["taul_unprojected_over_tauc"] = str(sp.simplify(M_hist / K / a / tc))
# exact static kernel (rational function of k from the exact LU solve of the tridiagonal block) and its expansion
Ak = Cr + sp.I * k * Vr
kap = sp.cancel(sp.together((br.T * Ak.LUsolve(br))[0]))
numer, denom = sp.fraction(kap)
ser = sp.expand(sp.series(numer / denom, k, 0, 6).removeO())
c0 = sp.cancel(ser.coeff(k, 0))
c2 = sp.cancel(-ser.coeff(k, 2))
c4 = sp.cancel(ser.coeff(k, 4))
res["A_grey"]["kernel_k0"] = str(c0)
res["A_grey"]["kernel_minus_k2_coef_equals_M"] = bool(sp.cancel(c2 - M) == 0)
res["A_grey"]["kernel_k4_coef"] = str(sp.factor(c4))
res["A_grey"]["k4_over_K_L4"] = str(sp.factor(sp.cancel(c4 / (K * (M / K) ** 2))))
# exact dispersion relation p + k^2 kappa(k, p) = 0, first order in p: lam_eff = -p / k^2
# coefficients of kappa(k, p) = K0 + Kp p + Kk2 k^2 + ... from derivatives of the exact rational kernel
Akp = Cr + p * sp.eye(n - 1) + sp.I * k * Vr
kapp = sp.cancel(sp.together((br.T * Akp.LUsolve(br))[0]))
K0 = sp.cancel(kapp.subs({k: 0, p: 0}))
Kp = sp.cancel(sp.diff(kapp, p).subs({k: 0, p: 0}))
Kk2 = sp.cancel(sp.diff(kapp, k, 2).subs({k: 0, p: 0}) / 2)
# root s = k^2 of p + s (K0 + Kp p + Kk2 s) = 0 to second order: s = -p/K0 (1 + eps p)
eps_ = sp.symbols("eps")
s_ = -p / K0 * (1 + eps_ * p)
eq = sp.expand(sp.series(p + s_ * (K0 + Kp * p + Kk2 * s_), p, 0, 3).removeO())
eps_sol = sp.solve(sp.Eq(eq.coeff(p, 2), 0), eps_)[0]
delta_disp = sp.simplify(eps_sol)          # lam_eff = -p/s = K0 (1 - eps p + ...), delta = eps
res["A_grey"]["delta_from_dispersion"] = str(sp.simplify(delta_disp))
res["A_grey"]["delta_equals_tauR_minus_taul"] = bool(sp.simplify(delta_disp - (tauR_kin - taul)) == 0)
res["A_grey"]["delta_equals_tauR_minus_4tauc_over_5"] = bool(sp.simplify(delta_disp - (tR - sp.Rational(4, 5) * tc)) == 0)
print("A grey:", json.dumps(res["A_grey"], indent=0), flush=True)

# ============================================================== B. random non-grey operators
rng = np.random.default_rng(424242)


def parity_operator(nh, conserve_momentum, seed):
    """Mode set of nh parity pairs; operator built in the parity-adapted basis (even/odd blocks), so it commutes
    with parity by construction (different construction from ra0_kinetic_tests.random_nongray)."""
    r = np.random.default_rng(seed)
    w = r.lognormal(0.0, 0.8, nh)                  # capacities (broad)
    vel = r.lognormal(0.0, 0.6, nh) * r.choice([-1, 1], nh)
    e_h = np.sqrt(w)
    # even/odd parity-adapted coordinates: f = [f+, f-]; even part (f+ + f-)/sqrt2, odd (f+ - f-)/sqrt2
    # energy invariant lives in the even block: e_even = sqrt2 * e_h ; current b = V e is odd: b_odd = sqrt2 * v e_h
    def spd(nn, scale, cond):
        Qm, _ = np.linalg.qr(r.normal(size=(nn, nn)))
        lam = scale * np.geomspace(1.0, cond, nn)
        r.shuffle(lam)
        return (Qm * lam) @ Qm.T
    ee = e_h / np.linalg.norm(e_h)
    Pe = np.eye(nh) - np.outer(ee, ee)
    Ceven = Pe @ spd(nh, 1.0, 300.0) @ Pe          # energy conserved exactly (even block)
    Codd = spd(nh, 0.3, 300.0)                       # odd block: momentum-like relaxation
    if conserve_momentum:
        # a normal-type part conserving an odd 'momentum' mode m (not parallel to the current): only part of
        # Codd; the resistive part keeps the odd block invertible
        mvec = e_h * vel * r.uniform(0.3, 1.7, nh)
        mm = mvec / np.linalg.norm(mvec)
        Pm = np.eye(nh) - np.outer(mm, mm)
        Codd = Pm @ spd(nh, 3.0, 50.0) @ Pm + spd(nh, 0.05, 30.0)
    # back to the mode basis
    T = np.block([[np.eye(nh), np.eye(nh)], [np.eye(nh), -np.eye(nh)]]) / np.sqrt(2.0)   # columns: even, odd
    Cblk = np.block([[Ceven, np.zeros((nh, nh))], [np.zeros((nh, nh)), Codd]])
    Cm_ = T @ Cblk @ T.T
    Cm_ = 0.5 * (Cm_ + Cm_.T)
    e_m = np.concatenate([e_h, e_h])
    v_m = np.concatenate([vel, -vel])
    return Cm_, e_m, v_m


def closed_form_kkt(Cm_, e_m, v_m):
    nn = e_m.size
    Kkt = np.zeros((nn + 1, nn + 1))
    Kkt[:nn, :nn] = Cm_
    Kkt[:nn, nn] = e_m
    Kkt[nn, :nn] = e_m
    bb = v_m * e_m
    sol = np.linalg.solve(Kkt, np.append(bb, 0.0))
    xx = sol[:nn]
    KK = bb @ xx
    tauR_ = xx @ xx / KK
    uu = v_m * xx
    uu = uu - e_m * (e_m @ uu) / (e_m @ e_m)
    y = np.linalg.solve(Kkt, np.append(uu, 0.0))[:nn]
    MM = uu @ y
    e2 = e_m @ e_m
    return {"K": KK, "tauR": tauR_, "taul": e2 * MM / KK ** 2, "M": MM, "e2": e2}


def kernel_bordered(Cm_, e_m, v_m, kk, pp):
    nn = e_m.size
    A = np.zeros((nn + 1, nn + 1), complex)
    A[:nn, :nn] = Cm_ + pp * np.eye(nn) + 1j * kk * np.diag(v_m)
    A[:nn, nn] = -e_m
    A[nn, :nn] = e_m
    rhs = np.zeros(nn + 1, complex)
    rhs[nn] = e_m @ e_m                 # e.f = |e|^2 theta with theta = 1
    sol = np.linalg.solve(A, rhs)
    f = sol[:nn]
    q = (v_m * e_m) @ f
    return q / (-1j * kk)


def coefs_numeric(Cm_, e_m, v_m, h_k, h_p):
    # K and tau_R from p-derivative at small k (even part in k removes the odd term), Richardson in p and k
    def kap_even(kk, pp):
        return 0.5 * (kernel_bordered(Cm_, e_m, v_m, kk, pp) + kernel_bordered(Cm_, e_m, v_m, -kk, pp))

    def dkdp(hk, hp):
        return (kap_even(hk, hp) - kap_even(hk, -hp)) / (2 * hp)

    def K0(hk):
        return kap_even(hk, 0.0)
    # K: Richardson in k^2 of kap_even(k, 0)
    K1, K2 = K0(h_k), K0(h_k / 2)
    Kx = (4 * K2 - K1) / 3
    # d/dk^2 at p = 0: (K(k) - K(0)) / k^2 with K(0) from extrapolation, Richardson
    c1 = (K1 - Kx) / h_k ** 2
    c2 = (K2 - Kx) / (h_k / 2) ** 2
    ck2 = (4 * c2 - c1) / 3
    # d/dp at k -> 0 (Richardson in p and in k)
    d11, d12 = dkdp(h_k, h_p), dkdp(h_k, h_p / 2)
    d21, d22 = dkdp(h_k / 2, h_p), dkdp(h_k / 2, h_p / 2)
    dA = (4 * d12 - d11) / 3
    dB = (4 * d22 - d21) / 3
    dp = (4 * dB - dA) / 3
    return {"K": Kx.real, "L2K": -ck2.real, "tauR": -dp.real / Kx.real, "imag_K": abs(Kx.imag)}


def dispersion_delta(Cm_, e_m, v_m, ps):
    """Diffusive root s = k^2 of p|e|^2 + s kappa(sqrt(s), p) = 0 by Newton (complex), lam_eff = -p|e|^2/s."""
    e2 = e_m @ e_m
    lam0 = closed_form_kkt(Cm_, e_m, v_m)["K"]
    out = []
    for pp in ps:
        s = -pp * e2 / lam0
        for _ in range(60):
            kk = np.sqrt(s)
            F = pp * e2 + s * kernel_bordered(Cm_, e_m, v_m, kk, pp)
            ds = 1e-7 * abs(s)
            kk2 = np.sqrt(s + ds)
            F2 = pp * e2 + (s + ds) * kernel_bordered(Cm_, e_m, v_m, kk2, pp)
            step = F / ((F2 - F) / ds)
            s = s - step
            if abs(step) < 1e-14 * abs(s):
                break
        lam_eff = -pp * e2 / s
        out.append(lam_eff / lam0)
    return np.array(out)


rows = []
for (nh, cm, seed) in ((30, False, 1), (60, False, 2), (60, True, 3), (120, True, 4)):
    Cm_, e_m, v_m = parity_operator(nh, cm, seed)
    cf = closed_form_kkt(Cm_, e_m, v_m)
    # scales for the finite differences: smallest relaxation rate and largest velocity
    ev_ = np.linalg.eigvalsh(Cm_)
    rmin = ev_[ev_ > 1e-10 * ev_.max()].min()
    hk = 2e-3 * rmin / np.abs(v_m).max()
    hp = 2e-3 * rmin
    nu = coefs_numeric(Cm_, e_m, v_m, hk, hp)
    # dispersion-relation slope: lam_eff/lam0 = 1 - delta p + O(p^2) -> Richardson in p
    ps = np.array([hp, hp / 2]) * 1j
    le = dispersion_delta(Cm_, e_m, v_m, ps)
    d1 = -(le[0] - 1) / ps[0]
    d2 = -(le[1] - 1) / ps[1]
    delta_disp = ((2 * d2 - d1)).real          # first-order Richardson (error O(p))
    a_ = cf["K"] / cf["e2"]
    row = {"n_modes": 2 * nh, "conserve_momentum_part": cm,
           "closed": {"K": cf["K"], "tauR": cf["tauR"], "taul": cf["taul"], "L2K": cf["M"]},
           "numeric": nu, "rel_err_K": abs(nu["K"] / cf["K"] - 1), "rel_err_tauR": abs(nu["tauR"] / cf["tauR"] - 1),
           "rel_err_L2K": abs(nu["L2K"] / cf["M"] - 1),
           "taul_numeric": nu["L2K"] / nu["K"] / a_, "delta_closed": cf["tauR"] - cf["taul"],
           "delta_dispersion": delta_disp,
           "rel_err_delta_dispersion": abs(delta_disp / (cf["tauR"] - cf["taul"]) - 1),
           "taul_over_tauR": cf["taul"] / cf["tauR"]}
    rows.append(row)
    print("B", nh, cm, f"relerr K {row['rel_err_K']:.1e} tauR {row['rel_err_tauR']:.1e} L2K {row['rel_err_L2K']:.1e}"
          f" delta(disp) {row['rel_err_delta_dispersion']:.1e} taul/tauR {row['taul_over_tauR']:.3f}", flush=True)
res["B_nongrey"] = rows

# ============================================================== C. twin convention
sys.path.insert(0, r"C:/Users/Koussay/these/src")
import twin_fdtr as tw  # noqa: E402

lam, Cv, tauR_, taul_, d = 150.0, 2.4e6, 25e-12, 220e-12, 500e-9
alpha = 1.0 / 3.0
ell = np.sqrt(taul_ * lam / ((1 + alpha) * Cv))
mat = tw.Material(lam, Cv, tau_r=tauR_, ell=ell, alpha=alpha)
res["C_twin"] = {"Material.tau_l": mat.tau_l, "expected": taul_, "rel": abs(mat.tau_l / taul_ - 1)}
# twin k = 0 kernel: Au transducer / R1 / GK film / R2 / sapphire (local temperature), compared with an
# independent 1D quadrupole product written here
fr = np.geomspace(1e4, 2e8, 25)
R1_, R2_ = 1e-8, 2e-8
au = tw.Material(150.0, 2.49e6)
sap = tw.Material(35.0, 3.03e6)
out_c = {}
for conv in ("local", "entropic"):
    st = tw.au_aln_sapphire(mat, aln_thickness=d, au_thickness=80e-9, au=au, sapphire=sap, r_au_aln=R1_,
                            r_aln_sapphire=R2_, boundary=tw.GKBoundary(slip=1.0, temperature=conv))
    Zt = tw.kernel(fr, 0.0, st)

    def quad(lamx, Cx, dx, pp):
        g = np.sqrt(Cx * pp / lamx)
        ch, sh = np.cosh(g * dx), np.sinh(g * dx)
        return np.array([[ch, sh / (lamx * g)], [lamx * g * sh, ch]])
    Zm = []
    for f_ in fr:
        pp = 2j * np.pi * f_
        lam_eff = lam * (1 + taul_ * pp) / (1 + tauR_ * pp)
        gfac = 1 + taul_ * pp if conv == "entropic" else 1.0
        # transfer matrices act on (theta, q); entropic: interface temperature theta = g T on both AlN faces
        Mau = quad(150.0, 2.49e6, 80e-9, pp)
        Mr1 = np.array([[1, R1_], [0, 1]])
        Gin = np.array([[gfac, 0], [0, 1]])         # (theta, q) at the film top = diag(g, 1) (T, q)
        Mf = quad(lam_eff, Cv, d, pp)
        Gout = np.array([[1 / gfac, 0], [0, 1]])    # (T, q) at the film bottom = diag(1/g, 1) (theta, q)
        Mr2 = np.array([[1, R2_], [0, 1]])
        Mtot = Mau @ Mr1 @ Gin @ Mf @ Gout @ Mr2
        Ys = np.sqrt(35.0 * 3.03e6 * pp)            # semi-infinite substrate: q = Ys T
        # (T0, q0) = Mtot (Tb, qb), qb = Ys Tb  ->  Z = T0/q0
        Zm.append((Mtot[0, 0] + Mtot[0, 1] * Ys) / (Mtot[1, 0] + Mtot[1, 1] * Ys))
    Zm = np.array(Zm)
    rat = Zt / Zm
    out_c[conv] = {"max_rel_diff": float(np.max(np.abs(rat - 1))), "ratio_spread": float(np.max(np.abs(rat / rat[0] - 1)))}
res["C_twin"]["kernel_vs_independent_quadrupole"] = out_c
print("C twin:", res["C_twin"], flush=True)
(OUT / "au1_derivation.json").write_text(json.dumps(res, indent=1, default=str))
