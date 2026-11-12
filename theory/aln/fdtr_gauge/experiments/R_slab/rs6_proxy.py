"""RS6: does a cheap regularised quantity predict the slab-visible delta of the witnesses?

Proxy: the first-order kinetic delta of the regularised operator C + QDQ(|v_z|/L_eff) (and RTA V1 for the RTA slab) at an
effective length L_eff = L/c, with c calibrated so that REF reproduces its slab delta_fit (N3 and N2lam, 0.01 deg,
note-16, wall port). Compared with the witnesses' slab delta_fit excursions (rs2b): slope, correlation, max error.
A proxy that works would make delta-optimised witnesses for the slab affordable (one KKT solve instead of a slab).
usage: python -B rs6_proxy.py [freqset]       Output: out/rs6_proxy.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs_core as rc  # noqa: E402
from rs2_slab import rta_film  # noqa: E402

FS = sys.argv[1] if len(sys.argv) > 1 else "sub4"
R2B = json.loads((OUT / f"rs2b_m775_film2_{FS}.json").read_text())
G = rc.geometry("m775")
names = ["REF"] + rc.WIT12


def dreg_C(C, Leff, op):
    if op == "full":
        return rc.kinetic_times(G, rc.regulated(G, C, Leff, "QDQ"))["delta"] * 1e12
    return rta_film(G["cap"], G["vz"], 1.0 / np.diag(C), Leff)["delta"] * 1e12


# pass 1: calibrate c on REF for every (operator family, thickness, nuisance set)
CREF = rc.operator(G, "film2", "REF")
combos = [(op, L, N) for op in ("full", "rta") for L in ("300", "500") for N in ("N3", "N2lam")]
cal = {}
for op, L, N in combos:
    d = float(L) * 1e-9
    target = R2B[op]["REF"][f"fit_{L}|wall|note16"][f"{N}|0.01"]["delta_fit_ps"]
    c = brentq(lambda c_: dreg_C(CREF, d / c_, op) - target, 1.0, 60.0, xtol=1e-6)
    cal[(op, L, N)] = c
# pass 2: proxies, one operator in memory at a time
prox = {k: {} for k in combos}
for nm in names:
    C = CREF if nm == "REF" else rc.operator(G, "film2", nm)
    for (op, L, N), c in cal.items():
        prox[(op, L, N)][nm] = dreg_C(C, float(L) * 1e-9 / c, op)
    del C
res = {}
for (op, L, N), c in cal.items():
    d = float(L) * 1e-9
    key = f"{N}|0.01"
    fit = {nm: R2B[op][nm][f"fit_{L}|wall|note16"][key]["delta_fit_ps"] for nm in names}
    pr = prox[(op, L, N)]
    ex_f = np.array([fit[w] - fit["REF"] for w in rc.WIT12])
    ex_p = np.array([pr[w] - pr["REF"] for w in rc.WIT12])
    res[f"{op}|{L}|{N}"] = {"c": c, "L_eff_nm": d / c * 1e9, "corr": float(np.corrcoef(ex_p, ex_f)[0, 1]),
                            "slope_fit_vs_proxy": float(np.polyfit(ex_p, ex_f, 1)[0]),
                            "max_abs_err_ps": float(np.max(np.abs(ex_p - ex_f))),
                            "rms_excursion_fit_ps": float(np.sqrt(np.mean(ex_f ** 2))),
                            "per_witness": {w: {"fit": float(a), "proxy": float(b)} for w, a, b in zip(rc.WIT12, ex_f, ex_p)}}
    print(op, L, N, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res[f"{op}|{L}|{N}"].items()
                     if k != "per_witness"}, flush=True)
(OUT / "rs6_proxy.json").write_text(json.dumps(res, indent=1))
