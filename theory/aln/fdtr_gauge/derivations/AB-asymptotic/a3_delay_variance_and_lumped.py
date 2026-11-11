"""A: (i) delay collinearity as a weighted variance of the film weight; (ii) lumped-transducer formula.

(i) Phase data only, equal noise, leading order in omega tau: dphi/dln tau_R = omega tau_R W(omega), delay column
    -omega. With mu_k = omega_k^2 / sum omega^2 on the sampled frequencies,
        cos^2(angle) = E_mu(W)^2 / E_mu(W^2),   1 - cos^2 = Var_mu(W) / E_mu(W^2)      (EXACT identity)
    and with the phase offset also free (columns 1 and -omega):
        1 - cos^2 = min_{a,c} E_mu[(W - a - c/omega)^2] / E_mu(W^2).
(ii) Au thermally thin (omega e_Au^2/a_Au <= 0.13 in band): Au = lumped heat capacity C_A = (rho c e)_Au, so
        dphi/w = omega tau_R Re[ zeta_f / (X (1 + i omega C_A X)) ],   X = R1 + Z_f,   zeta_f = -dZ_f/dln lambda,
    zeta_f -> R_f (1 + T_f + T_f^2)/3 (thin AlN), Z_f/2 (thick AlN). Compared with the exact 1D stack.
"""
import json
import os

import numpy as np

from ab_common import BASE, with_, dlnZ, Dual, stack_Z_d

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "a3_delay_variance_and_lumped.json")
f80 = np.logspace(4, np.log10(2e8), 80)
w80 = 2 * np.pi * f80
res = {}


def weight_stats(P, omega=w80):
    W = -dlnZ(omega, with_(P, tau_r=0.0), "lam").real          # film weight at leading order
    s = dlnZ(omega, P, "tau_r").imag                            # exact phase signature per unit w
    mu = omega ** 2 / np.sum(omega ** 2)
    EW, EW2 = np.sum(mu * W), np.sum(mu * W * W)
    one_minus_cos2_LO = (EW2 - EW ** 2) / EW2
    cos_exact = np.dot(s, omega) / (np.linalg.norm(s) * np.linalg.norm(omega))
    # with phase offset: project W onto span{1, 1/omega} in L2(mu)
    B = np.vstack([np.ones_like(omega), 1 / omega]).T * np.sqrt(mu)[:, None]
    coef, *_ = np.linalg.lstsq(B, W * np.sqrt(mu), rcond=None)
    r = W * np.sqrt(mu) - B @ coef
    # same with the exact signature and columns {1, omega}
    A = np.vstack([np.ones_like(omega), omega]).T
    c2, *_ = np.linalg.lstsq(A, s, rcond=None)
    rs = s - A @ c2
    return dict(E_mu_W=float(EW), rel_var=float(one_minus_cos2_LO), one_minus_cos2_exact=float(1 - cos_exact ** 2),
                with_offset_LO=float(np.dot(r, r) / EW2), with_offset_exact=float(np.dot(rs, rs) / np.dot(s, s)),
                W_top5=W[-5:].tolist())


for name, P in (("base 500nm", dict(BASE)), ("500nm R1=R2=0", with_(BASE, R1=1e-12, R2=1e-12)),
                ("100nm lam130 tau11", with_(BASE, e=100e-9, lam=130.0, tau_r=11e-12)),
                ("300nm lam180 tau20", with_(BASE, e=300e-9, lam=180.0, tau_r=20e-12)),
                ("bare AlN on sapphire", with_(BASE, e_au=1e-12, R1=1e-12, R2=1e-12))):
    st = weight_stats(P)
    res[name] = st
    print("%-22s E_mu W = %.4f  1-cos^2: LO %.4f exact %.4f | +phase offset: LO %.4f exact %.4f" %
          (name, st["E_mu_W"], st["rel_var"], st["one_minus_cos2_exact"], st["with_offset_LO"], st["with_offset_exact"]))

# ---------------------------------------------------------------- (ii) lumped transducer
P = dict(BASE)
ff = np.logspace(4, np.log10(2e8), 200)
om = 2 * np.pi * ff
p = 1j * om
bs = np.sqrt(P["lam_s"] * P["C_s"])
Z2 = 1 / (bs * np.sqrt(p)) + P["R2"]
xi1 = P["e"] * np.sqrt(P["C"] / P["lam"])
b = np.sqrt(P["lam"] * P["C"])
lamD = Dual(P["lam"], P["lam"])
from ab_common import dsqrt, _layer_d
xiD = P["e"] * dsqrt(P["C"] / lamD)
bD = dsqrt(lamD * P["C"])
Zf_d = _layer_d(Z2, xiD * np.sqrt(p), bD * np.sqrt(p))
Zf, zeta_f = Zf_d.v, -Zf_d.d
X = P["R1"] + Zf
C_A = P["C_au"] * P["e_au"]
lumped = om * P["tau_r"] * np.real(zeta_f / (X * (1 + 1j * om * C_A * X)))
exact = dlnZ(om, P, "tau_r").imag
# thin and thick film limits of zeta_f
Af = xi1 * np.sqrt(p)
Tf = (2 * np.exp(-Af) / (1 + np.exp(-2 * Af))) / (1 + b * np.sqrt(p) * Z2 * (1 - np.exp(-2 * Af)) / (1 + np.exp(-2 * Af)))
Rf = P["e"] / P["lam"]
zeta_thin = Rf * (1 + Tf + Tf ** 2) / 3
zeta_thick = 1 / (2 * b * np.sqrt(p))
thin = om * P["tau_r"] * np.real(zeta_thin / (X * (1 + 1j * om * C_A * X)))
thick = om * P["tau_r"] * np.real(zeta_thick / (X * (1 + 1j * om * C_A * X)))
i_max = np.argmax(np.abs(exact))
res["lumped"] = {
    "max_exact_deg_per_w": float(np.degrees(np.abs(exact[i_max]))), "f_max_Hz": float(ff[i_max]),
    "max_lumped_deg_per_w": float(np.degrees(np.max(np.abs(lumped)))),
    "max_rel_err_lumped_over_max": float(np.max(np.abs(lumped - exact)) / np.max(np.abs(exact))),
    "thin_over_exact_at": {str(f): float(np.interp(np.log(f), np.log(ff), thin / exact)) for f in (1e5, 1e6, 1e7, 3e7)},
    "thick_over_exact_at": {str(f): float(np.interp(np.log(f), np.log(ff), thick / exact)) for f in (1e8, 1.5e8, 2e8)},
}
print("lumped transducer: max exact %.4e deg/w at %.3g MHz; lumped max %.4e; max |lumped-exact|/max = %.3e" %
      (res["lumped"]["max_exact_deg_per_w"], ff[i_max] / 1e6, res["lumped"]["max_lumped_deg_per_w"], res["lumped"]["max_rel_err_lumped_over_max"]))
print("   thin-film zeta/exact:", {k: round(v, 3) for k, v in res["lumped"]["thin_over_exact_at"].items()})
print("   thick-film zeta/exact:", {k: round(v, 3) for k, v in res["lumped"]["thick_over_exact_at"].items()})
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("written", OUT)
