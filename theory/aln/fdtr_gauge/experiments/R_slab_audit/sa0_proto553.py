"""SA0: prototype of the independent slab solver on the 5x5x3 event set (REF): GMRES convergence, z-mesh
refinement, comparison with the stored modal two-ports (rs2_m553_filmw_full_sub4_REF.npz) and, for the RTA
reduction, with the audit's RTA slab code au6_slab.slab_twoport. Output: out/sa0_proto553.json"""
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


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


G = sc.load_modes("m553")
log("modes", G["n"], "orbits", G["n_orb"], "missing images", G["missing_images"])
C = sc.operator(G, "filmw", "REF")
e, v = G["e"], G["vz"]
log("energy residual", np.linalg.norm(C @ e) / np.linalg.norm(np.diag(C) * e), "sym", np.abs(C - C.T).max() / np.abs(C).max())
kb = sc.kinetic_bulk(C, e, v)
log("bulk", kb)
Zs = np.load(RS / "rs2_m553_filmw_full_sub4_REF.npz")
fs, ds, Ys, Ths = Zs["f"], Zs["d"], Zs["Y"], Zs["Th"]
if Ys.ndim == 5:
    Ys, Ths = Ys[0], Ths[0]
res = {"bulk": kb, "rows": []}
for fi in (3, 10, fs.size - 1):
    f = fs[fi]
    p = 2j * np.pi * f
    for di, d in enumerate(ds):
        Yref, Tref = Ys[di, fi], Ths[di, fi]
        prev = None
        for lev in (0, 1, 2):
            z = sc.make_mesh(d, level=lev)
            Y, Th, dg, _ = sc.twoport(C, v, e, p, z, rtol=1e-12, restart=60, maxiter=40)
            row = {"f": f, "d_nm": d * 1e9, "level": lev, "N": int(z.size - 1),
                   "rel_Y": float(np.max(np.abs(Y - Yref)) / np.max(np.abs(Yref))),
                   "rel_ImY11": float(abs((Y[0, 0] - Yref[0, 0]).imag) / abs(Yref[0, 0].imag)),
                   "rel_Th": float(np.max(np.abs(Th - Tref)) / np.max(np.abs(Tref))),
                   "ebal": dg["ebal"], "iters": [g["iters"] for g in dg["gmres"]],
                   "relres": [g["true_relres"] for g in dg["gmres"]], "t": dg["t_total"]}
            if prev is not None:
                row["rel_change_Y"] = float(np.max(np.abs(Y - prev)) / np.max(np.abs(Yref)))
            prev = Y
            res["rows"].append(row)
            log(row)
            (HERE / "out" / "sa0_proto553.json").write_text(json.dumps(res, indent=1, default=float))
