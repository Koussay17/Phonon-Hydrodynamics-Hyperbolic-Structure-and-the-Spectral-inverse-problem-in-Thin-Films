"""SA0c: fourth-order (cubic-source) variant of the independent solver on 5x5x3 (REF, 300/500 nm, 10 kHz and
200 MHz) against the stored modal two-ports; order-2 regression at level 0. Output: out/sa0c_order4_553.json"""
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
fs, ds, Ys, Ts = Zs["f"], Zs["d"], Zs["Y"], Zs["Th"]
Ys, Ts = (Ys[0], Ts[0]) if Ys.ndim == 5 else (Ys, Ts)
MESH = {"hmin": 0.02e-9, "gamma": 0.08, "hmax": 4e-9}
res = {}
for fi in (3, fs.size - 1):
    p = 2j * np.pi * fs[fi]
    for di, d in enumerate(ds):
        Yref, Tref = Ys[di, fi], Ts[di, fi]
        nrm = np.max(np.abs(Yref))
        for order, levels in ((2, (0,)), (4, (0, 1, 2))):
            Yl, rows = [], []
            for lev in levels:
                z = sc.make_mesh(d, level=lev, **MESH)
                Y, Th, dg, _ = sc.twoport(C, v, e, p, z, rtol=2e-12, restart=60, maxiter=4, order=order)
                Yl.append(Y)
                rows.append({"level": lev, "N": int(z.size - 1), "err": float(np.max(np.abs(Y - Yref)) / nrm),
                             "err_Th": float(np.max(np.abs(Th - Tref)) / np.max(np.abs(Tref))),
                             "err_ImY11_rel": float(abs((Y[0, 0] - Yref[0, 0]).imag) / abs(Yref[0, 0].imag)),
                             "ebal": dg["ebal"], "iters": [g["iters"] for g in dg["gmres"]], "t": dg["t_total"]})
            ent = {"rows": rows}
            if len(Yl) == 3:
                ent["p_obs"] = float(np.log2(np.max(np.abs(Yl[1] - Yl[0])) / np.max(np.abs(Yl[2] - Yl[1]))))
                R = (16 * Yl[2] - Yl[1]) / 15
                ent["err_rich12_o4"] = float(np.max(np.abs(R - Yref)) / nrm)
            key = f"{fs[fi]:.4g}|{d * 1e9:.0f}|order{order}"
            res[key] = ent
            print(key, [(r["N"], f"{r['err']:.2e}", f"Im {r['err_ImY11_rel']:.1e}", f"eb {r['ebal']:.1e}", r["iters"][0]) for r in rows],
                  {k: f"{v:.3g}" for k, v in ent.items() if k != "rows"}, f"[{time.perf_counter() - T0:.0f}s]", flush=True)
            (HERE / "out" / "sa0c_order4_553.json").write_text(json.dumps(res, indent=1, default=str))
