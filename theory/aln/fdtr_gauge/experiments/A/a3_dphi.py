"""A3: phase change Delta phi(f) = phi(f; tau_R (1 + w)) - phi(f; tau_R), everything else fixed (exact, nonlinear).

Thicknesses 30, 100, 300, 500 nm with their own base values (A1); w = +-2.8, 3.7, 5.2, 8 %.
Models: Cattaneo with finite beams (10/8 um) and in the 1D large-beam limit; GK (tau_l = 0.5 and 2 tau_R) with
finite beams; 500 nm also with the RTA base values. Also d ln|H| and the linear prediction ln(1+w) dphi/dln tau.
Usage: python a3_dphi.py -> out/a3_dphi.json, out/a3_dphi_curves.csv
"""
from __future__ import annotations

import json
import time

import numpy as np

from common import FREQ, BANDS, OUT, W_ADMISSIBLE, model_for, band_mask, jacobian, load_base, Model, jsonable

T0 = time.perf_counter()
base = load_base()
WS = [s * w for w in W_ADMISSIBLE for s in (-1, 1)]
variants = []
for L in (30, 100, 300, 500):
    variants += [(L, "cat_beam", dict(law="cattaneo", geometry="beam")),
                 (L, "cat_1d", dict(law="cattaneo", geometry="1d")),
                 (L, "gk0.5_beam", dict(law="gk", taul_factor=0.5, geometry="beam")),
                 (L, "gk2_beam", dict(law="gk", taul_factor=2.0, geometry="beam"))]
res = {"w": WS, "freq_Hz": FREQ.tolist(), "cases": {}}
rows = []
for L, name, kw in variants:
    m = model_for(L, **kw)
    l0 = m.lnH(FREQ)
    J = jacobian(m, FREQ, ["ln_tauR"], h=2e-3, scheme="richardson")[:, 0]
    case = {"tauR_ps": m.tauR * 1e12, "lamz": m.lamz, "lamr": m.phys()["lamr"], "taul_ps": m.taul * 1e12,
            "phase_deg": np.degrees(l0.imag).tolist(),
            "dphi_dlntau_deg": np.degrees(J.imag).tolist(), "dlnA_dlntau": J.real.tolist(), "w": {}}
    for w in WS:
        lw = m.lnH(FREQ, {"ln_tauR": np.log1p(w)})
        dphi = np.degrees((lw - l0).imag)
        dlnA = (lw - l0).real
        lin = np.degrees(np.log1p(w) * J.imag)
        e = {"dphi_deg": dphi.tolist(), "dlnA": dlnA.tolist(),
             "nonlinearity_max_abs_deg": float(np.max(np.abs(dphi - lin))),
             "max_abs_dphi_deg": float(np.max(np.abs(dphi))), "f_at_max_Hz": float(FREQ[np.argmax(np.abs(dphi))]),
             "max_abs_dlnA": float(np.max(np.abs(dlnA))), "f_at_max_dlnA_Hz": float(FREQ[np.argmax(np.abs(dlnA))]),
             "rms_dphi_deg": float(np.sqrt(np.mean(dphi ** 2)))}
        for b in BANDS:
            mk = band_mask(FREQ, b)
            e[f"max_abs_dphi_deg_{b}"] = float(np.max(np.abs(dphi[mk])))
        case["w"][f"{w:+.3f}"] = e
        rows.append((L, name, w, dphi, dlnA))
    res["cases"][f"L{L}_{name}"] = case
    mx = {k: round(v["max_abs_dphi_deg"], 6) for k, v in case["w"].items()}
    print(L, name, f"tau {m.tauR*1e12:.3f} ps", "max|dphi| deg", mx, f"{time.perf_counter()-T0:.0f}s", flush=True)

# 500 nm with the RTA base values (sensitivity of the base point)
b = base["500_rta"]
m = Model(L_nm=500.0, tauR=b["tau_R_ps"] * 1e-12, lamz=b["lam_z"], lamr=b["lam_r"])
l0 = m.lnH(FREQ)
case = {"tauR_ps": m.tauR * 1e12, "lamz": m.lamz, "lamr": m.lamr, "w": {}}
for w in WS:
    dphi = np.degrees((m.lnH(FREQ, {"ln_tauR": np.log1p(w)}) - l0).imag)
    case["w"][f"{w:+.3f}"] = {"max_abs_dphi_deg": float(np.max(np.abs(dphi))),
                             "f_at_max_Hz": float(FREQ[np.argmax(np.abs(dphi))]), "dphi_deg": dphi.tolist()}
res["cases"]["L500_cat_beam_RTAbase"] = case
print("500 RTA base", {k: round(v["max_abs_dphi_deg"], 6) for k, v in case["w"].items()})

(OUT / "a3_dphi.json").write_text(json.dumps(jsonable(res), indent=1))
with open(OUT / "a3_dphi_curves.csv", "w") as fh:
    fh.write("L_nm,variant,w," + ",".join(f"dphi_deg@{f:.6g}" for f in FREQ) + "\n")
    for L, name, w, dphi, _ in rows:
        fh.write(f"{L},{name},{w:+.3f}," + ",".join(f"{x:.6e}" for x in dphi) + "\n")
print("done", time.perf_counter() - T0)
