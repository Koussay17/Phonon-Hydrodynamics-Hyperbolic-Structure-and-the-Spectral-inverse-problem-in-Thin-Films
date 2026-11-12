"""RS2c: accuracy of the AAA interpolation of the slab two-ports (REF, 7x7x5, full and RTA, all wall models).

From the exact 'all' set (80 FDTR frequencies + 100 Hz, 1 kHz, 2 kHz) take every k-th FDTR frequency (k = 3, 4, 5, 6;
plus the last one and the three low frequencies) as AAA support, evaluate at the 80 FDTR frequencies and compare:
max relative error of each two-port entry; max |Delta ln H| through the note-16 stack (wall port); delta_fit (N2lam,
N3, 0.01 deg) against the exact data.  Output: out/rs2c_interp_check.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import AAA

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs2b_fits as fb  # noqa: E402

a5, a6 = fb.a5, fb.a6
F = a5.FREQ
res = {}
stat = json.loads((OUT / "rs2_m775_film2_static_REF.json").read_text())
for op in ("full", "rta"):
    walls = fb.walls_of("m775", "film2", op, "all", "REF")
    for iw, w in enumerate(walls):
        f, d, Y, Th, _ = fb.load("m775", "film2", op, "all", "REF", iw)
        Yx, Tx, _ = fb.at_fdtr(f, Y, Th, "all")
        low = np.where(f < 9e3)[0]
        fi = np.where(f >= 9e3)[0]
        for k in (3, 4, 5, 6):
            sub = np.unique(np.r_[np.arange(0, 80, k), 79])
            idx = np.r_[low, fi[sub]]
            pz = 2j * np.pi * f[idx]
            ent = {"n_support": int(idx.size)}
            for i, dd in enumerate(d):
                L = f"{dd * 1e9:.0f}"
                Yi = np.zeros((F.size, 2, 2), complex)
                Ti = np.zeros_like(Yi)
                for a in range(2):
                    for b in range(2):
                        Yi[:, a, b] = AAA(pz, Y[i, idx, a, b], rtol=1e-14)(2j * np.pi * F)
                        Ti[:, a, b] = AAA(pz, Th[i, idx, a, b], rtol=1e-14)(2j * np.pi * F)
                eY = float(np.max(np.abs(Yi - Yx[i])) / np.max(np.abs(Yx[i])))
                eT = float(np.max(np.abs(Ti - Tx[i])) / np.max(np.abs(Tx[i])))
                p = 2j * np.pi * F
                y_ex = np.log(a6.Z_stack_slab(list(Yx[i]), list(Tx[i]), p, 1e-8, 2e-8, port="wall"))
                y_in = np.log(a6.Z_stack_slab(list(Yi), list(Ti), p, 1e-8, 2e-8, port="wall"))
                kr = stat[f"reg_full_QDQ_{L}"] if op == "full" else stat[f"reg_rta_V1_{L}"]
                base = {"tauR": kr["tauR"], "taul": kr["taul"], "lamz": kr["K"] if op == "full" else kr["lam"],
                        "C": stat["C_tot"], "R1": 1e-8, "R2": 2e-8, "d": dd}
                dfit = {}
                for N in ("N2lam", "N3"):
                    c1, t1 = a5.gn_refit(y_ex, base, fb.SETS[N], 0.01, free=(1,))
                    c2, t2 = a5.gn_refit(y_in, base, fb.SETS[N], 0.01, free=(1,))
                    dfit[N] = {"exact_ps": (kr["tauR"] - kr["taul"] - t1[1] * 1e-12) * 1e12,
                               "interp_ps": (kr["tauR"] - kr["taul"] - t2[1] * 1e-12) * 1e12}
                ent[L] = {"max_rel_err_Y": eY, "max_rel_err_Th": eT,
                          "max_abs_dlnH_real": float(np.max(np.abs((y_in - y_ex).real))),
                          "max_abs_dphase_deg": float(np.rad2deg(np.max(np.abs((y_in - y_ex).imag)))), "delta_fit": dfit}
            res[f"{op}|{fb.wall_tag(w)}|every{k}"] = ent
            print(op, fb.wall_tag(w), k, json.dumps(ent), flush=True)
(OUT / "rs2c_interp_check.json").write_text(json.dumps(res, indent=1))
