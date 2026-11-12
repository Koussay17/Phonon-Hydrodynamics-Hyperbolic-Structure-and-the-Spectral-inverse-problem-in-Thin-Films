"""B-4, step 8: practical conditioning of the transducer configuration in the FDTR band (1D, k = 0).

Stack of note 16 / Workstream A: Au 80 nm (150 W/mK, 2.49e6 J/m3K) / R1 / AlN 500 nm (lam = 188.6 W/mK,
C = 2.41e6 J/m3K, tau_R = 25 ps, tau_l in {0, 10, 25} ps) / R2 / sapphire (35 W/mK, 3.03e6 J/m3K);
R1 = 1e-8, R2 = 2e-8 m2K/W; 80 log-spaced frequencies 10 kHz - 200 MHz; noise 1 % on ln|H| and sigma_phi on arg H.
Conventions: local (tau_n = 0) and entropic (tau_n = tau_l at both film faces, the twin's "entropic").
Parameter sets (all in log except the two times, which enter in seconds because tau_l can be 0):
  S6     = (ln lam, ln C, tau_R, tau_l, ln R1, ln R2)            [structural problem of B-4, e known]
  S6+I   = S6 + instrument (log gain, phase offset, delay)        [A's N3 + tau_l]
  S5+I   = S5 (tau_l = 0 known, Cattaneo) + instrument            [A's N3, 1D; review r8 reproduces it]
Reports: whitened singular values (relative), sigma(tau_R)/tau_R, sigma(tau_l) [ps], and the standard deviation of
the leading low-frequency time combination delta_eps = tau_R - tau_l + eps (1 + b^2/b_s^2) tau_l (Sec. 6 of the
proof file). Derivatives: central differences, relative step 1e-5 (log parameters) / 1e-15 s (times).
Output: b4_08_practical.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import numpy as np

f = np.logspace(4, np.log10(2e8), 80)
w = 2 * np.pi * f
AU = (150.0, 2.49e6, 80e-9)
SAP = (35.0, 3.03e6)
E_ALN = 500e-9


def quad_fourier(lam, C, d, p):
    s = np.sqrt(C * p / lam)
    B = lam * s
    return np.cosh(s * d), np.sinh(s * d) / B, B * np.sinh(s * d)


def ztot(P, conv):
    p = 1j * w
    lam, C, tR, tl, R1, R2 = P["lam"], P["C"], P["tR"], P["tl"], P["R1"], P["R2"]
    a = lam / C
    xi = E_ALN / np.sqrt(a)
    b = np.sqrt(lam * C)
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(SAP[0] * SAP[1]) * np.sqrt(p)
    N = 1 + (tl * p if conv == "entropic" else 0)
    L = (R2 + 1 / S) / N
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    W = R1 + N * (1 + r * E) / ((1 - r * E) * B)
    c, x, y = quad_fourier(*AU, p)
    return (c * W + x) / (y * W + c)


def columns(P, conv, names):
    cols = []
    for nm in names:
        if nm in ("tR", "tl"):
            h = 1e-15
            Pp, Pm = dict(P), dict(P)
            Pp[nm] += h
            Pm[nm] -= h
            cols.append((np.log(ztot(Pp, conv)) - np.log(ztot(Pm, conv))) / (2 * h) * 1e-12)  # per ps
        elif nm == "gain":
            cols.append(np.ones_like(w, dtype=complex))
        elif nm == "phase":
            cols.append(1j * np.ones_like(w, dtype=complex))
        elif nm == "delay":
            cols.append(-1j * w * 1e-9)                                   # per ns
        else:
            h = 1e-5
            Pp, Pm = dict(P), dict(P)
            Pp[nm] *= np.exp(h)
            Pm[nm] *= np.exp(-h)
            cols.append((np.log(ztot(Pp, conv)) - np.log(ztot(Pm, conv))) / (2 * h))
    return cols


def analyse(P, conv, names, sphi_deg, sA=0.01):
    sphi = np.radians(sphi_deg)
    cols = columns(P, conv, names)
    J = np.array([np.concatenate([c.real / sA, c.imag / sphi]) for c in cols]).T
    s = np.linalg.svd(J / np.linalg.norm(J, axis=0), compute_uv=False)
    cov = np.linalg.inv(J.T @ J)
    out = {"rel_svals": [float(x / s[0]) for x in s]}
    if "tR" in names:
        i = names.index("tR")
        out["sigma_tauR_over_tauR"] = float(np.sqrt(cov[i, i]) / (P["tR"] * 1e12))
    if "tl" in names:
        j = names.index("tl")
        out["sigma_taul_ps"] = float(np.sqrt(cov[j, j]))
        i = names.index("tR")
        eps = 1.0 if conv == "entropic" else 0.0
        bb = np.sqrt(P["lam"] * P["C"]) ** 2 / (SAP[0] * SAP[1])
        v = np.zeros(len(names))
        v[i], v[j] = 1.0, -1.0 + eps * (1 + bb)
        out["delta_eps_weights"] = [1.0, float(v[j])]
        out["sigma_delta_eps_ps"] = float(np.sqrt(v @ cov @ v))
        u = np.zeros(len(names))
        u[i], u[j] = 1.0, -1.0
        out["sigma_(tauR-taul)_ps"] = float(np.sqrt(u @ cov @ u))
        # best-determined combination tau_R + w tau_l (smallest variance per unit tau_R weight)
        C2 = cov[np.ix_([i, j], [i, j])]
        evals, evecs = np.linalg.eigh(C2)
        v0 = evecs[:, 0] / evecs[0, 0]
        out["best_combination_tauR_plus_w_taul"] = {"w": float(v0[1]),
                                                    "sigma_ps": float(np.sqrt(v0 @ C2 @ v0))}
    return out


base = dict(lam=188.6, C=2.41e6, tR=25e-12, R1=1e-8, R2=2e-8)
S6 = ["lam", "C", "tR", "tl", "R1", "R2"]
S6I = S6 + ["gain", "phase", "delay"]
S5I = ["lam", "C", "tR", "R1", "R2", "gain", "phase", "delay"]
OUT = {}
for tl in (0.0, 10e-12, 25e-12):
    P = dict(base, tl=tl)
    for conv in ("local", "entropic"):
        for sp in (0.01, 0.1):
            key = f"tl={tl*1e12:.0f}ps,{conv},sigma_phi={sp}"
            OUT[key] = {"S6": analyse(P, conv, S6, sp), "S6+I": analyse(P, conv, S6I, sp)}
            if tl == 0.0:
                OUT[key]["S5+I (A's N3, 1D)"] = analyse(P, conv, S5I, sp)
            print(key, "S6 s_min/s_max=%.2e" % OUT[key]["S6"]["rel_svals"][-1],
                  "sig(tR)/tR S6=%.3g S6+I=%.3g" % (OUT[key]["S6"]["sigma_tauR_over_tauR"],
                                                   OUT[key]["S6+I"]["sigma_tauR_over_tauR"]),
                  "sig(tl) S6+I=%.3g ps" % OUT[key]["S6+I"]["sigma_taul_ps"],
                  "sig(delta_eps) S6+I=%.3g ps" % OUT[key]["S6+I"]["sigma_delta_eps_ps"],
                  "best: tR%+.3g tl, sig=%.3g ps" % (OUT[key]["S6+I"]["best_combination_tauR_plus_w_taul"]["w"],
                                                    OUT[key]["S6+I"]["best_combination_tauR_plus_w_taul"]["sigma_ps"]),
                  ("S5+I sig(tR)/tR=%.3g" % OUT[key]["S5+I (A's N3, 1D)"]["sigma_tauR_over_tauR"]) if tl == 0.0 else "")
# phase difference produced by the convention itself at fixed parameters (GK film, transducer stack)
dphi = {}
for tl in (10e-12, 25e-12):
    P = dict(base, tl=tl)
    d = np.degrees(np.angle(ztot(P, "entropic") / ztot(P, "local")))
    dphi[f"tl={tl*1e12:.0f}ps"] = {"max_abs_deg": float(np.max(np.abs(d))), "f_at_max_MHz": float(f[np.argmax(np.abs(d))] / 1e6)}
    Pf = dict(base, tl=tl, tR=tl)
    PF = dict(base, tl=0.0, tR=0.0)
    d2 = np.degrees(np.angle(ztot(Pf, "entropic") / ztot(PF, "entropic")))
    dphi[f"resonant tR=tl={tl*1e12:.0f}ps vs Fourier, entropic"] = float(np.max(np.abs(d2)))
    d3 = np.degrees(np.angle(ztot(Pf, "local") / ztot(PF, "local")))
    dphi[f"resonant tR=tl={tl*1e12:.0f}ps vs Fourier, local"] = float(np.max(np.abs(d3)))
OUT["phase_effect_of_convention"] = dphi
print(json.dumps(dphi, indent=1))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_08_practical.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
