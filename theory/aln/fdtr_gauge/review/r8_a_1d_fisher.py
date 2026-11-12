"""Minimal independent 1D check of workstream A (large-beam limit; my own code, no import from the project).
Stack: Au (Fourier, 80 nm, 150 W/mK, 2.49e6) / R1 / AlN (Cattaneo, L, lam_z, C = 2.41e6, tau_R) / R2 / sapphire (35, 3.03e6).
Observable ln H = ln Z(i omega) + g + i phi0 - i omega t_d. Noise 1 % on ln|H|, sigma_phi on arg H, 80 log-spaced 10 kHz-200 MHz.
Fisher: sigma(ln tau_R) with nuisance sets N0 {}, N1 {g, phi0, t_d}, N2 N1 + {ln R1, ln R2}, N3 N2 + {ln lam_z, ln C}.
Also: max |d phi / d ln tau_R| and sensitivity of sigma to (R1, R2)."""
import json
import numpy as np


def zstack(f, P):
    p = 2j * np.pi * f
    # substrate
    Z = 1 / (np.sqrt(35.0 * 3.03e6) * np.sqrt(p))
    Z = Z + P["R2"]
    # AlN Cattaneo layer: lam_eff = lam/(1 + tau p)
    lam_eff = P["lam"] / (1 + P["tau"] * p)
    s = np.sqrt(P["C"] * p / lam_eff)
    Y = lam_eff * s
    t = np.tanh(s * P["L"])
    Z = (Z + t / Y) / (1 + Y * Z * t)
    Z = Z + P["R1"]
    # Au layer
    s = np.sqrt(2.49e6 * p / 150.0)
    Y = 150.0 * s
    t = np.tanh(s * 80e-9)
    Z = (Z + t / Y) / (1 + Y * Z * t)
    return Z


def sigma_lntau(P, sphi_deg, nuis, f, sA=0.01):
    sphi = np.radians(sphi_deg)
    names = ["tau"] + [n for n in nuis if n in ("R1", "R2", "lam", "C")]
    cols = []
    for n in names:
        h = 1e-5
        Pp = dict(P); Pm = dict(P)
        Pp[n] = P[n] * np.exp(h); Pm[n] = P[n] * np.exp(-h)
        d = (np.log(zstack(f, Pp)) - np.log(zstack(f, Pm))) / (2 * h)
        cols.append(d)
    if "inst" in nuis:
        cols.append(np.ones_like(f, dtype=complex))            # gain
        cols.append(1j * np.ones_like(f, dtype=complex))       # phase offset (rad)
        cols.append(-1j * 2 * np.pi * f * 1e-9)                # delay (ns)
    J = np.array([np.concatenate([c.real / sA, c.imag / sphi]) for c in cols]).T
    Fm = J.T @ J
    cov = np.linalg.inv(Fm)
    return float(np.sqrt(cov[0, 0]))


f = np.logspace(4, np.log10(2e8), 80)
base = {"500": dict(L=500e-9, lam=188.6, C=2.41e6, tau=25.0e-12, R1=1e-8, R2=2e-8),
        "300": dict(L=300e-9, lam=167.25, C=2.41e6, tau=19.70e-12, R1=1e-8, R2=2e-8),
        "100": dict(L=100e-9, lam=115.62, C=2.41e6, tau=11.23e-12, R1=1e-8, R2=2e-8)}
out = {}
for Lnm, P in base.items():
    row = {}
    # max phase sensitivity
    h = 1e-5
    Pp = dict(P); Pm = dict(P); Pp["tau"] *= np.exp(h); Pm["tau"] *= np.exp(-h)
    dphi = np.degrees((np.angle(zstack(f, Pp) / zstack(f, Pm))) / (2 * h))
    row["max_dphi_dlntau_deg"] = float(np.max(np.abs(dphi)))
    row["f_at_max_MHz"] = float(f[np.argmax(np.abs(dphi))] / 1e6)
    for sp in (0.01, 0.1):
        row[f"N0_{sp}"] = sigma_lntau(P, sp, [], f)
        row[f"N1_{sp}"] = sigma_lntau(P, sp, ["inst"], f)
        row[f"N2_{sp}"] = sigma_lntau(P, sp, ["inst", "R1", "R2"], f)
        row[f"N3_{sp}"] = sigma_lntau(P, sp, ["inst", "R1", "R2", "lam", "C"], f)
    out[Lnm] = row
# resistance map at 500 nm, N2
rmap = {}
for R1 in (2e-9, 5e-9, 1e-8, 2e-8):
    for R2 in (2e-9, 5e-9, 1e-8, 2e-8):
        P = dict(base["500"]); P["R1"] = R1; P["R2"] = R2
        rmap[f"R1={R1:.0e},R2={R2:.0e}"] = {"N2_0.01": sigma_lntau(P, 0.01, ["inst", "R1", "R2"], f),
                                           "N2_0.1": sigma_lntau(P, 0.1, ["inst", "R1", "R2"], f),
                                           "N3_0.1": sigma_lntau(P, 0.1, ["inst", "R1", "R2", "lam", "C"], f)}
out["R_map_500nm"] = rmap
print(json.dumps(out, indent=1))
with open(__file__.replace(".py", ".json"), "w") as fh:
    json.dump(out, fh, indent=1)
