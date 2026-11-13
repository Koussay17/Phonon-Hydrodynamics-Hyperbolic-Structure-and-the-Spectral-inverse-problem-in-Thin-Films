"""SA0b: which part of the z mesh controls the error (5x5x3, REF, 500 nm, 200 MHz and 10 kHz)? Variation of
hmax (interior), hmin (wall), gamma (grading); Richardson extrapolation over bisection levels.
Output: out/sa0b_meshstudy553.json"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sa_core as sc  # noqa: E402
import numpy as np  # noqa: E402

RS = HERE.parent / "R_slab" / "out"
T0 = time.perf_counter()
G = sc.load_modes("m553")
C = sc.operator(G, "filmw", "REF")
e, v = G["e"], G["vz"]
Zs = np.load(RS / "rs2_m553_filmw_full_sub4_REF.npz")
fs, ds, Ys = Zs["f"], Zs["d"], Zs["Y"]
Ys = Ys[0] if Ys.ndim == 5 else Ys
res = {}
for fi in (3, fs.size - 1):
    p = 2j * np.pi * fs[fi]
    di = 1
    d = ds[di]
    Yref = Ys[di, fi]
    nrm = np.max(np.abs(Yref))
    for tag, kw in (("base", {}), ("hmax4", {"hmax": 4e-9}), ("hmin0.005", {"hmin": 0.005e-9}),
                    ("gamma0.05", {"gamma": 0.05}), ("hmax4_gamma0.08", {"hmax": 4e-9, "gamma": 0.08})):
        Yl = []
        rows = []
        for lev in (0, 1, 2):
            z = sc.make_mesh(d, level=lev, **kw)
            Y, Th, dg, _ = sc.twoport(C, v, e, p, z, rtol=1e-13, restart=60, maxiter=40)
            Yl.append(Y)
            rows.append({"level": lev, "N": int(z.size - 1), "err": float(np.max(np.abs(Y - Yref)) / nrm),
                         "iters": dg["gmres"][0]["iters"], "t": dg["t_total"]})
        R12 = (4 * Yl[2] - Yl[1]) / 3
        R01 = (4 * Yl[1] - Yl[0]) / 3
        RR = (16 * R12 - R01) / 15
        p_obs = float(np.log2(np.max(np.abs(Yl[1] - Yl[0])) / np.max(np.abs(Yl[2] - Yl[1]))))
        ent = {"rows": rows, "p_obs": p_obs, "err_rich12": float(np.max(np.abs(R12 - Yref)) / nrm),
               "err_rich01": float(np.max(np.abs(R01 - Yref)) / nrm), "err_rich012": float(np.max(np.abs(RR - Yref)) / nrm),
               "err_rich12_ImY11": float(abs((R12[0, 0] - Yref[0, 0]).imag) / abs(Yref[0, 0].imag)),
               "Y_levels": [[[complex(x) for x in r] for r in Y_] for Y_ in Yl]}
        res[f"{fs[fi]:.4g}|{tag}"] = ent
        print(f"{fs[fi]:.4g}", tag, [(r["N"], f"{r['err']:.2e}") for r in rows], f"p_obs {p_obs:.2f}",
              f"rich12 {ent['err_rich12']:.2e} rich012 {ent['err_rich012']:.2e} ImY11(rich12) {ent['err_rich12_ImY11']:.2e}",
              f"[{time.perf_counter() - T0:.0f}s]", flush=True)
        (HERE / "out" / "sa0b_meshstudy553.json").write_text(json.dumps(res, indent=1, default=str))
