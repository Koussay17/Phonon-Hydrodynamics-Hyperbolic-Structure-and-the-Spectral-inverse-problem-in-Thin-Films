"""A: collinearity of the tau_R sensitivity with the instrument delay and the other nuisances (1D Fisher).

Model: H_obs = exp(g + i phi0 - i omega t_d) Z(omega; theta), 1D stack Au/R1/AlN/R2/sapphire (ab_common),
80 log-spaced frequencies 10 kHz - 200 MHz, Gaussian independent noise sigma_a on ln|H| and sigma_phi on arg H.
Nuisance sets (00-question.md): N0 none; N1 gain, phase offset, delay; N2 + R1, R2; N3 + lambda_z, C (AlN);
N4 + tau_l and Au (lambda_Au, C_Au). N5 (beam radii) has no 1D analogue.
Outputs: w*(1 sigma) = sqrt([F^-1]_{ln tau, ln tau}); information retention rho_N = 1/(F_tt [F^-1]_tt);
cosines of the tau_R phase signature with the delay (omega) and with span{1, omega}; the equivalent delay.
Also the leading-order structural identity dphi/dtau_R ~ omega R_f Re(<T^2>/Z) and its delay-collinear plateau.
"""
import json
import os

import numpy as np

from ab_common import BASE, with_, dlnZ, stack_flux_transfer

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "a2_nuisance_collinearity.json")
DEG = np.pi / 180
f80 = np.logspace(4, np.log10(2e8), 80)
w80 = 2 * np.pi * f80

THERMAL = {  # name: log-parameter?
    "tau_r": True, "R1": True, "R2": True, "lam": True, "C": True, "tau_l": False, "lam_au": True, "C_au": True,
}
SETS = {
    "N0": [],
    "N1": ["gain", "phi0", "delay"],
    "N2": ["gain", "phi0", "delay", "R1", "R2"],
    "N3": ["gain", "phi0", "delay", "R1", "R2", "lam", "C"],
    "N4": ["gain", "phi0", "delay", "R1", "R2", "lam", "C", "tau_l", "lam_au", "C_au"],
    # partial sets to locate the loss of information
    "N1_delay_only": ["delay"],
    "N1_phi0_only": ["phi0"],
    "N0+R1R2": ["R1", "R2"],
    "N0+tau_l": ["tau_l"],
    "N0+C": ["C"],
    "N0+lam": ["lam"],
    "N0+lam+C": ["lam", "C"],
    "N0+delay+C": ["delay", "C"],
    "N1+C": ["gain", "phi0", "delay", "C"],
    "N1+lam": ["gain", "phi0", "delay", "lam"],
    "N2+C": ["gain", "phi0", "delay", "R1", "R2", "C"],
    "N2+lam": ["gain", "phi0", "delay", "R1", "R2", "lam"],
    "N0+C_au": ["C_au"],
    "N1+C_au": ["gain", "phi0", "delay", "C_au"],
}


def columns(P, omega, names):
    """Complex d ln H / d theta for each name (instrument columns analytic)."""
    cols = {}
    for n in names:
        if n == "gain":
            cols[n] = np.ones_like(omega, dtype=complex)
        elif n == "phi0":
            cols[n] = 1j * np.ones_like(omega, dtype=complex)
        elif n == "delay":
            cols[n] = -1j * omega * 1e-12           # per ps
        elif n == "tau_l":
            cols[n] = dlnZ(omega, P, "tau_l", log=False) * 1e-12   # per ps (absolute, boundary-safe)
        else:
            cols[n] = dlnZ(omega, P, n, log=THERMAL[n])
    return cols


def fisher(P, sig_phi_deg, sig_a=0.01, use_amp=True, sets=SETS, omega=w80):
    names = ["tau_r", "gain", "phi0", "delay", "R1", "R2", "lam", "C", "tau_l", "lam_au", "C_au"]
    cols = columns(P, omega, names)
    sp = sig_phi_deg * DEG
    out = {}
    for sname, nuis in sets.items():
        use = ["tau_r"] + nuis
        rows = []
        for n in use:
            c = cols[n]
            rows.append(np.concatenate([c.real / sig_a if use_amp else np.zeros(0), c.imag / sp]))
        J = np.array(rows).T
        if not use_amp and "gain" in use:
            J = J[:, [i for i, n in enumerate(use) if n != "gain"]]
        F = J.T @ J
        # scale-invariant inverse via normalised columns
        nrm = np.sqrt(np.diag(F))
        Fn = F / np.outer(nrm, nrm)
        sv = np.linalg.svd(Fn, compute_uv=False)
        cond = sv[0] / sv[-1]
        Ci = np.linalg.inv(Fn) / np.outer(nrm, nrm)
        wstar = float(np.sqrt(Ci[0, 0]))
        retention = float(1.0 / (F[0, 0] * Ci[0, 0]))
        out[sname] = {"w_star_1sigma": wstar, "w_star_2sigma": 2 * wstar, "retention": retention,
                      "cond_normalised_fisher": float(cond)}
    return out


res = {}
P = dict(BASE)

# ------------------------------------------------------------------ structure of the tau_R phase signature
s_tau = dlnZ(w80, P, "tau_r")                 # d ln Z / d ln tau_R
phase_sig = s_tau.imag
amp_sig = s_tau.real
om = w80
c_delay = float(np.dot(phase_sig, om) / (np.linalg.norm(phase_sig) * np.linalg.norm(om)))
Aplane = np.vstack([np.ones_like(om), om]).T
coef, *_ = np.linalg.lstsq(Aplane, phase_sig, rcond=None)
resid = phase_sig - Aplane @ coef
sin_plane = float(np.linalg.norm(resid) / np.linalg.norm(phase_sig))
t_eq_per_w = float(-np.dot(phase_sig, om) / np.dot(om, om))     # s per unit w
res["structure"] = {"cos_phase_signature_vs_delay": c_delay,
                    "sin_angle_to_span{1,omega}": sin_plane,
                    "equivalent_delay_ps_per_w": t_eq_per_w * 1e12,
                    "equivalent_delay_ps_for_w=0.05": 0.05 * t_eq_per_w * 1e12,
                    "max_phase_sig_deg_per_w": float(np.max(np.abs(phase_sig)) / DEG),
                    "max_amp_sig_per_w": float(np.max(np.abs(amp_sig)))}
print("tau_R phase signature: cos(angle with delay) = %.5f; sin(angle to span{1,omega}) = %.4f" % (c_delay, sin_plane))
print("equivalent delay: %.4f ps per unit w (%.4f ps at w = 0.05)" % (t_eq_per_w * 1e12, 0.05 * t_eq_per_w * 1e12))

# band-restricted cosines (where is the non-delay information?)
for lo_, hi_ in ((1e4, 2e8), (1e6, 2e8), (1e7, 2e8), (3e7, 2e8), (1e7, 1e8)):
    m = (f80 >= lo_) & (f80 <= hi_)
    c = np.dot(phase_sig[m], om[m]) / (np.linalg.norm(phase_sig[m]) * np.linalg.norm(om[m]))
    res["structure"]["cos_delay_band_%g_%g" % (lo_, hi_)] = float(c)

# R sensitivities in the phase and their relation to the tau_R signature
sR1 = dlnZ(w80, P, "R1"); sR2 = dlnZ(w80, P, "R2"); slam = dlnZ(w80, P, "lam")
def cosv(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
res["structure"]["cos_phase_tau_vs_phase_R1"] = cosv(phase_sig, sR1.imag)
res["structure"]["cos_phase_tau_vs_phase_R2"] = cosv(phase_sig, sR2.imag)
res["structure"]["cos_phase_tau_vs_phase_lambda"] = cosv(phase_sig, slam.imag)
res["structure"]["cos_phase_tau_vs_omega_x_amp_R1"] = cosv(phase_sig, om * sR1.real)
# leading-order identity: phase_sig ~ omega tau_R W0, W0 = -Re s_lambda(tau=0); and ~ omega tau R_f Re(<T^2>/Z)
P0 = with_(P, tau_r=0.0)
W0 = -dlnZ(w80, P0, "lam").real
ft = stack_flux_transfer(w80, P0, n_gauss=32)
Rf = P["e"] / P["lam"]
Wrec = Rf * np.real(ft["meanT2"] / ft["parts"]["Z"])
res["structure"]["max_rel_dev_W0_vs_reciprocity"] = float(np.max(np.abs(W0 - Wrec)) / np.max(np.abs(W0)))
res["structure"]["max_rel_dev_phase_sig_vs_LO"] = float(np.max(np.abs(phase_sig - om * P["tau_r"] * W0)) / np.max(np.abs(phase_sig)))
print("cos(phase tau, phase R1) %.4f, (phase tau, phase R2) %.4f, (phase tau, phase lambda) %.4f, (phase tau, omega*amp R1) %.4f" %
      (res["structure"]["cos_phase_tau_vs_phase_R1"], res["structure"]["cos_phase_tau_vs_phase_R2"],
       res["structure"]["cos_phase_tau_vs_phase_lambda"], res["structure"]["cos_phase_tau_vs_omega_x_amp_R1"]))
print("band cosines with delay:", {k[15:]: round(v, 4) for k, v in res["structure"].items() if k.startswith("cos_delay_band")})

# ------------------------------------------------------------------ Fisher: base and scenarios
scen = {
    "500nm lam200 tau27.5": dict(BASE),
    "500nm lam150 tau35": with_(BASE, lam=150.0, tau_r=35e-12),
    "500nm lam250 tau20": with_(BASE, lam=250.0, tau_r=20e-12),
    "300nm lam180 tau20": with_(BASE, e=300e-9, lam=180.0, tau_r=20e-12),
    "100nm lam130 tau11": with_(BASE, e=100e-9, lam=130.0, tau_r=11e-12),
    "500nm no R (R1=R2=0)": with_(BASE, R1=1e-12, R2=1e-12),
}
res["fisher"] = {}
for name, Ps in scen.items():
    res["fisher"][name] = {}
    for sig in (0.01, 0.1):
        for amp in (True, False):
            key = "sigphi=%g_%s" % (sig, "amp+phase" if amp else "phase")
            res["fisher"][name][key] = fisher(Ps, sig, use_amp=amp)
print("\nw*(1 sigma) [relative change of tau_R], 1D, 80 freq. 10 kHz-200 MHz, sigma_a = 1 %")
hdr = ["N0", "N1_phi0_only", "N1_delay_only", "N1", "N0+R1R2", "N2", "N3", "N0+tau_l", "N4"]
print("%-24s %-22s" % ("scenario", "noise") + "".join("%13s" % h for h in hdr))
for name in scen:
    for key in ("sigphi=0.01_amp+phase", "sigphi=0.1_amp+phase", "sigphi=0.01_phase"):
        r = res["fisher"][name][key]
        print("%-24s %-22s" % (name, key) + "".join("%13.3g" % r[h]["w_star_1sigma"] for h in hdr))
print("\ninformation retention rho_N (base, sigma_phi = 0.01 deg, amp+phase):")
r = res["fisher"]["500nm lam200 tau27.5"]["sigphi=0.01_amp+phase"]
print({h: "%.3e" % r[h]["retention"] for h in hdr})

print("structural identity checks: max dev W0 vs reciprocity %.2e; phase signature vs LO %.2e" %
      (res["structure"]["max_rel_dev_W0_vs_reciprocity"], res["structure"]["max_rel_dev_phase_sig_vs_LO"]))
hdr2 = ["N0+C", "N0+lam", "N0+lam+C", "N0+delay+C", "N1+C", "N1+lam", "N2+C", "N2+lam", "N0+C_au", "N1+C_au"]
print("partial sets, base, sigma_phi = 0.01 deg, amp+phase: retention")
r = res["fisher"]["500nm lam200 tau27.5"]["sigphi=0.01_amp+phase"]
print({h: "%.3e" % r[h]["retention"] for h in hdr2})
# whitened single-column cosines with the tau_R column (amp + phase, sigma_phi = 0.01 deg)
names = ["gain", "phi0", "delay", "R1", "R2", "lam", "C", "tau_l", "lam_au", "C_au"]
cols = columns(P, w80, ["tau_r"] + names)
def wvec(c, sp=0.01 * DEG, sa=0.01):
    return np.concatenate([c.real / sa, c.imag / sp])
vt = wvec(cols["tau_r"])
res["single_cosines"] = {n: float(np.dot(vt, wvec(cols[n])) / (np.linalg.norm(vt) * np.linalg.norm(wvec(cols[n])))) for n in names}
print("whitened cosines with the tau_R column:", {k: round(v, 4) for k, v in res["single_cosines"].items()})

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("written", OUT)
