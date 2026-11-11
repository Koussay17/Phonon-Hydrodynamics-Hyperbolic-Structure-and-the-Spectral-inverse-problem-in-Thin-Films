"""A: phase sensitivity of the 1D Au/R1/AlN/R2/sapphire response to the AlN Cattaneo time tau_R.

Checks and outputs
 I1  EXACT identity: d ln Z/d ln tau_R = -(tau_R p/(1+tau_R p)) s_lambda(p; lambda_eff(p)),
     s_lambda = d ln Z/d ln lambda_AlN at fixed C, e (complexified conductivity).
 I2  EXACT reciprocity: dZ/dtau_R = p R_f <T^2>/(1+tau_l p), dZ/dR_i = T_i^2, T(z) = q(z)/q(0).
 I3  ASYMPTOTIC (thin AlN): <T^2> = (T1^2 + T1 T2 + T2^2)/3 + O(|A_f|^2); error slope vs omega.
 L1  Leading order in omega tau: dphi = omega tau_R w W0(omega), W0 = -Re s_lambda|_{tau=0} = -d ln|Z|/d ln lambda.
     Exact Delta phi vs first order in w (error O(w^2)) vs leading order in omega tau (error O(omega tau)).
 C1  Closed forms: substrate regime dphi ~ (tau w R_f b_s/sqrt2) omega^{3/2}; transducer-limited band
     dphi ~ K omega^3 (wa + wb)/((omega^2 + wa^2)(omega^2 + wb^2)), K = tau w R_f/(3 C_A R1^2),
     wa = 1/(C_f R1), wb = (1 + C_A/C_f)/(C_A R1); max ~ K/2 = tau w R_f/(6 C_A R1^2).
 S1  Scenario table: max |Delta phi| over 10 kHz - 200 MHz for lambda, tau, w, thickness, interfaces, transducer.
"""
import json
import os

import numpy as np

from ab_common import (BASE, stack_Z, stack_flux_transfer, with_, dlnZ, Dual, stack_Z_d)


def log_sensitivity(omega, P, name):
    return dlnZ(omega, P, name, log=True)


def abs_sensitivity(omega, P, name, h=None):
    return dlnZ(omega, P, name, log=False)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "a1_phase_sensitivity.json")
DEG = 180 / np.pi
f80 = np.logspace(np.log10(1e4), np.log10(2e8), 80)
w80 = 2 * np.pi * f80
ffine = np.logspace(3, 10, 701)
wfine = 2 * np.pi * ffine
res = {"base": BASE}


def dphi_exact(omega, P, w):
    return np.imag(np.log(stack_Z(omega, with_(P, tau_r=P["tau_r"] * (1 + w))) / stack_Z(omega, P)))


# ---------------------------------------------------------------- I1
P = dict(BASE)
s_tau = log_sensitivity(wfine, P, "tau_r")
# s_lambda at complex lambda_eff: d ln Z / d ln lambda evaluated with lambda -> lambda_eff, tau -> 0
pp = 1j * wfine
lam_eff = P["lam"] * (1 + P["tau_l"] * pp) / (1 + P["tau_r"] * pp)
# evaluate s_lambda(p; lambda_eff) frequency by frequency (complex conductivity, Fourier film)
def Z_complex_lambda(omega, lam_c):
    Q = with_(P, tau_r=0.0, tau_l=0.0)
    p = 1j * omega
    bs = np.sqrt(Q["lam_s"] * Q["C_s"])
    Z2 = 1 / (bs * np.sqrt(p)) + Q["R2"]
    A = Q["e"] * np.sqrt(p * Q["C"] / lam_c)
    B = np.sqrt(lam_c * Q["C"] * p)
    t = (1 - np.exp(-2 * A)) / (1 + np.exp(-2 * A))
    Zf = (Z2 + t / B) / (1 + B * Z2 * t)
    Z1 = Zf + Q["R1"]
    Aa = Q["e_au"] * np.sqrt(p * Q["C_au"] / Q["lam_au"])
    Ba = np.sqrt(Q["lam_au"] * Q["C_au"] * p)
    ta = (1 - np.exp(-2 * Aa)) / (1 + np.exp(-2 * Aa))
    return (Z1 + ta / Ba) / (1 + Ba * Z1 * ta)
def Z_complex_lambda_d(omega, lam_c):
    Q = with_(P, tau_r=0.0, tau_l=0.0)
    Q["lam"] = Dual(lam_c, lam_c)          # d/d ln lambda at the complex value lambda_eff(p)
    Z = stack_Z_d(omega, Q)
    return Z.d / Z.v
s_lam_eff = Z_complex_lambda_d(wfine, lam_eff)
pred = -(P["tau_r"] * pp / (1 + P["tau_r"] * pp)) * s_lam_eff
res["I1_max_rel_err"] = float(np.max(np.abs(pred - s_tau) / np.abs(s_tau)))
print("I1 exact identity, max rel. err:", res["I1_max_rel_err"])

# ---------------------------------------------------------------- I2, I3
ft = stack_flux_transfer(wfine, P, n_gauss=32)
Rf = P["e"] / P["lam"]
dZ_dtau = s_tau / P["tau_r"] * ft["parts"]["Z"]
pred_tau = pp * Rf * ft["meanT2"] / (1 + P["tau_l"] * pp)
dZ_dR1 = abs_sensitivity(wfine, P, "R1", 1e-11) * ft["parts"]["Z"]
dZ_dR2 = abs_sensitivity(wfine, P, "R2", 1e-11) * ft["parts"]["Z"]
res["I2_tau_rel_err"] = float(np.max(np.abs(pred_tau - dZ_dtau) / np.abs(dZ_dtau)))
res["I2_R1_rel_err"] = float(np.max(np.abs(ft["T1"] ** 2 - dZ_dR1) / np.abs(dZ_dR1)))
res["I2_R2_rel_err"] = float(np.max(np.abs(ft["T2"] ** 2 - dZ_dR2) / np.abs(dZ_dR2)))
print("I2 reciprocity rel. err: tau %.2e  R1 %.2e  R2 %.2e" % (res["I2_tau_rel_err"], res["I2_R1_rel_err"], res["I2_R2_rel_err"]))
T1, T2 = ft["T1"], ft["T2"]
# transfer relative to the AlN top (removes the Au factor, which is common)
lin = (T1 ** 2 + T1 * T2 + T2 ** 2) / 3
trap = (T1 ** 2 + T2 ** 2) / 2
err_lin = np.abs(lin - ft["meanT2"]) / np.abs(ft["meanT2"])
err_trap = np.abs(trap - ft["meanT2"]) / np.abs(ft["meanT2"])
Af = ft["parts"]["Af"]
lowf = ffine < 2e6
sl_lin = np.polyfit(np.log(np.abs(Af[lowf]) ** 2), np.log(err_lin[lowf]), 1)[0]
res["I3"] = {"slope_err_lin_vs_|A|^2": float(sl_lin),
             "err_lin_at_f": {str(f): float(np.interp(np.log(f), np.log(ffine), err_lin)) for f in (1e6, 1e7, 5e7, 1e8, 2e8)},
             "err_trap_at_f": {str(f): float(np.interp(np.log(f), np.log(ffine), err_trap)) for f in (1e6, 1e7, 5e7, 1e8, 2e8)}}
print("I3 linear-flux weight: err slope vs |A|^2 = %.3f; err at 1e7,1e8,2e8 Hz:" % sl_lin,
      ["%.2e" % res["I3"]["err_lin_at_f"][k] for k in ("10000000.0", "100000000.0", "200000000.0")])

# ---------------------------------------------------------------- L1
P0 = with_(P, tau_r=0.0)
s_lam0 = log_sensitivity(wfine, P0, "lam")
W0 = -np.real(s_lam0)
W0_80 = -np.real(log_sensitivity(w80, P0, "lam"))
lo = wfine * P["tau_r"] * W0                     # leading order, per unit w
first = np.imag(s_tau)                           # first order in w, exact in omega tau
rows = {}
for w in (0.03, 0.05, 0.08, -0.05):
    ex = dphi_exact(wfine, P, w)
    bnd = (ffine >= 1e4) & (ffine <= 2e8)
    rows[str(w)] = {"max_abs_deg": float(np.max(np.abs(ex[bnd])) * DEG),
                    "f_at_max_Hz": float(ffine[bnd][np.argmax(np.abs(ex[bnd]))]),
                    "rel_err_first_order_in_w_at_200MHz": float(abs(np.interp(np.log(2e8), np.log(ffine), ex / (w * first) - 1))),
                    "rel_err_leading_order_at_200MHz": float(abs(np.interp(np.log(2e8), np.log(ffine), ex / (w * lo) - 1)))}
res["L1"] = rows
band = (ffine >= 1e4) & (ffine <= 2e8)
# absolute error of the leading order, normalised by max|lo| on the band: predicted O((omega tau)^2)
absLO = np.abs(first[band] - lo[band]) / np.max(np.abs(lo[band]))
otau = (wfine * P["tau_r"])[band]
sel = (otau > 1e-4) & (otau < 1e-2)
sl_LO = np.polyfit(np.log(otau[sel]), np.log(absLO[sel]), 1)[0]
res["L1"]["slope_abs_err_LO_vs_omega_tau"] = float(sl_LO)
res["L1"]["abs_err_LO_over_max_at_200MHz"] = float(absLO[-1])
# tau -> 0 at fixed frequencies: (first - lo)/tau^2 -> const
tau_scan = []
for tau in (40e-12, 20e-12, 10e-12, 5e-12, 2.5e-12):
    Pt = with_(P, tau_r=tau)
    st = np.imag(dlnZ(w80[-5:], Pt, "tau_r"))
    lt = w80[-5:] * tau * W0_80[-5:]
    tau_scan.append([tau] + list((st - lt) / tau ** 2))
tau_scan = np.array(tau_scan)
res["L1"]["tau_scan_(first-lo)/tau^2_top5freqs"] = tau_scan.tolist()
res["L1"]["W0_at_f"] = {str(f): float(np.interp(np.log(f), np.log(ffine), W0)) for f in (1e5, 1e6, 1e7, 3e7, 1e8, 2e8)}
print("L1 max |Delta phi| (deg) for w = 0.03/0.05/0.08:", ["%.2e" % rows[k]["max_abs_deg"] for k in ("0.03", "0.05", "0.08")],
      " at f =", "%.3g Hz" % rows["0.05"]["f_at_max_Hz"])
print("L1 first-order-in-w error at 200 MHz (w=0.08): %.2e; |first - LO|/max|LO|: slope vs omega tau %.3f (expect 2), at 200 MHz %.2e" %
      (rows["0.08"]["rel_err_first_order_in_w_at_200MHz"], sl_LO, res["L1"]["abs_err_LO_over_max_at_200MHz"]))
print("L1 tau-scan (first - LO)/tau^2 at 200 MHz:", ["%.4e" % v for v in tau_scan[:, -1]])
print("L1 film weight W0(f):", {k: round(v, 5) for k, v in res["L1"]["W0_at_f"].items()})

# ---------------------------------------------------------------- C1 closed forms
bs = np.sqrt(P["lam_s"] * P["C_s"])
C_A = P["C_au"] * P["e_au"]
C_f = P["C"] * P["e"]
R1, R2 = P["R1"], P["R2"]
K = P["tau_r"] * Rf / (3 * C_A * R1 ** 2)
wa = 1 / (C_f * R1)
wb = (1 + C_A / C_f) / (C_A * R1)
closed_mid = K * wfine ** 3 * (wa + wb) / ((wfine ** 2 + wa ** 2) * (wfine ** 2 + wb ** 2))
closed_low = P["tau_r"] * Rf * bs / np.sqrt(2) * wfine ** 1.5
exact_lin = first                                  # per unit w
lowmask = ffine < 3e5
res["C1"] = {
    "K_rad_per_w": K, "f_a_Hz": wa / 2 / np.pi, "f_b_Hz": wb / 2 / np.pi,
    "f_s_Hz_(b_s R2)^-2": 1 / (bs * R2) ** 2 / 2 / np.pi, "f_2_Hz_1/(C_f R2)": 1 / (C_f * R2) / 2 / np.pi,
    "f_film_Hz_a/(2 pi e^2)": P["lam"] / P["C"] / P["e"] ** 2 / 2 / np.pi,
    "closed_max_K_over_2_deg_per_w": K / 2 * DEG,
    "exact_max_deg_per_w": float(np.max(np.abs(exact_lin[band])) * DEG),
    "f_exact_max_Hz": float(ffine[band][np.argmax(np.abs(exact_lin[band]))]),
    "closed_mid_max_deg_per_w": float(np.max(closed_mid[band]) * DEG),
    "f_closed_mid_max_Hz": float(ffine[band][np.argmax(closed_mid[band])]),
    "low_regime_ratio_exact_over_closed_at_1e3..3e5": [float(v) for v in (exact_lin[lowmask] / closed_low[lowmask])[::40]],
    "mid_ratio_exact_over_closed_at_f": {str(f): float(np.interp(np.log(f), np.log(ffine), exact_lin / closed_mid)) for f in (1e7, 3e7, 1e8, 2e8)},
}
print("C1 closed forms: K/2 = %.3e deg/w at f_b = %.3g MHz; exact max %.3e deg/w at %.3g MHz" %
      (res["C1"]["closed_max_K_over_2_deg_per_w"], res["C1"]["f_b_Hz"] / 1e6, res["C1"]["exact_max_deg_per_w"], res["C1"]["f_exact_max_Hz"] / 1e6))
print("   scales: f_s %.3g MHz, f_2 %.3g MHz, f_a %.3g MHz, f_film %.3g MHz, f_b %.3g MHz" %
      (res["C1"]["f_s_Hz_(b_s R2)^-2"] / 1e6, res["C1"]["f_2_Hz_1/(C_f R2)"] / 1e6, res["C1"]["f_a_Hz"] / 1e6, res["C1"]["f_film_Hz_a/(2 pi e^2)"] / 1e6, res["C1"]["f_b_Hz"] / 1e6))
print("   low-frequency ratio exact/closed:", ["%.3f" % v for v in res["C1"]["low_regime_ratio_exact_over_closed_at_1e3..3e5"]])
print("   mid ratio exact/closed:", {k: round(v, 3) for k, v in res["C1"]["mid_ratio_exact_over_closed_at_f"].items()})

# save curves
res["curves"] = {"f_Hz": ffine.tolist(), "dphi_per_w_rad_first_order": first.tolist(), "dphi_per_w_rad_leading": lo.tolist(),
                 "closed_mid": closed_mid.tolist(), "closed_low": closed_low.tolist(), "W0": W0.tolist(),
                 "dlnabs_per_w_first_order": np.real(s_tau).tolist()}

# ---------------------------------------------------------------- S1 scenario table (max over the 80-point band and fine band)
def scenario(P, label):
    s = log_sensitivity(wfine[band], P, "tau_r")
    s80 = log_sensitivity(w80, P, "tau_r")
    i = np.argmax(np.abs(s.imag))
    return {"label": label, "max_dphi_deg_per_w": float(np.abs(s.imag[i]) * DEG), "f_max_Hz": float(ffine[band][i]),
            "dphi_deg_per_w_at_200MHz": float(s80.imag[-1] * DEG),
            "rms80_dphi_deg_per_w": float(np.sqrt(np.mean(s80.imag ** 2)) * DEG),
            "sum80_sq_dphi_rad2_per_w2": float(np.sum(s80.imag ** 2)),
            "sum80_sq_dlnamp_per_w2": float(np.sum(s80.real ** 2))}

table = []
for lam in (150.0, 200.0, 250.0):
    for tau in (20e-12, 27.5e-12, 35e-12):
        table.append(scenario(with_(BASE, lam=lam, tau_r=tau), "500nm lam=%g tau=%gps" % (lam, tau * 1e12)))
table.append(scenario(with_(BASE, R1=2e-8, R2=2e-8), "500nm R1=R2=2e-8"))
table.append(scenario(with_(BASE, R1=1e-9, R2=2e-8), "500nm R1=1e-9"))
table.append(scenario(with_(BASE, R1=1e-12, R2=1e-12), "500nm R1=R2=0 (perfect contacts)"))
table.append(scenario(with_(BASE, e_au=1e-12, R1=1e-12), "500nm no transducer, R2 kept"))
table.append(scenario(with_(BASE, e_au=1e-12, R1=1e-12, R2=1e-12), "500nm bare AlN on sapphire"))
table.append(scenario(with_(BASE, e_au=40e-9), "500nm Au 40 nm"))
for e, lam, tau in ((100e-9, 130.0, 11e-12), (100e-9, 100.0, 11e-12), (300e-9, 180.0, 20e-12), (300e-9, 150.0, 20e-12)):
    table.append(scenario(with_(BASE, e=e, lam=lam, tau_r=tau), "%dnm lam=%g tau=%gps" % (e * 1e9, lam, tau * 1e12)))
res["S1"] = table
print("\n%-38s %12s %10s %12s" % ("scenario", "max deg/w", "f_max MHz", "200MHz deg/w"))
for r in table:
    print("%-38s %12.3e %10.3g %12.3e" % (r["label"], r["max_dphi_deg_per_w"], r["f_max_Hz"] / 1e6, r["dphi_deg_per_w_at_200MHz"]))

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("written", OUT)
