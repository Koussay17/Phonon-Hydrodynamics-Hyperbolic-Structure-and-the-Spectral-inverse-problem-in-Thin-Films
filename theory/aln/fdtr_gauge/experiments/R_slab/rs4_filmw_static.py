"""RS4a: regularised (first-order) film times of REF and the ten film witnesses (filmw_pp) that exist on both 5x5x3
and 7x7x5, for the mesh-trend check of the slab transfer. Full operator C + QDQ (|v_z|/L) and RTA V1, 300/500 nm.
Writes out/rs2_<mesh>_filmw_static_<name>.json (same format as rs2_slab.py) and out/rs4_filmw_static.json (summary with
the delta-extreme witnesses per mesh and thickness).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs_core as rc  # noqa: E402
from rs2_slab import DS, rta_film  # noqa: E402

res = {}
for mesh in ("m553", "m775"):
    G = rc.geometry(mesh)
    rows = {}
    for nm in ["REF"] + rc.WITW:
        C = rc.operator(G, "filmw", nm)
        meta = {"name": nm, "mesh": mesh, "group": "filmw", "C_tot": G["C_tot"], "bulk_full": rc.kinetic_times(G, C)}
        r = np.diag(C)
        for L in DS:
            for var in ("QDQ", "phys"):
                meta[f"reg_full_{var}_{L * 1e9:.0f}"] = rc.kinetic_times(G, rc.regulated(G, C, L, var))
            meta[f"reg_rta_V1_{L * 1e9:.0f}"] = rta_film(G["cap"], G["vz"], 1.0 / r, L)
        meta["bulk_rta"] = rc.kinetic_times(G, rc.rta_operator(G, C))
        (OUT / f"rs2_{mesh}_filmw_static_{nm}.json").write_text(json.dumps(meta, indent=1))
        rows[nm] = {L: {"delta_reg_full_ps": meta[f"reg_full_QDQ_{L}"]["delta"] * 1e12,
                        "delta_reg_rta_ps": meta[f"reg_rta_V1_{L}"]["delta"] * 1e12,
                        "tauR_full_ps": meta[f"reg_full_QDQ_{L}"]["tauR"] * 1e12} for L in ("300", "500")}
        print(mesh, nm, {L: {k: round(v, 2) for k, v in rows[nm][L].items()} for L in rows[nm]}, flush=True)
    ext = {}
    for L in ("300", "500"):
        ref = rows["REF"][L]["delta_reg_full_ps"]
        exc = {w: rows[w][L]["delta_reg_full_ps"] - ref for w in rc.WITW}
        ext[L] = {"max": max(exc, key=exc.get), "min": min(exc, key=exc.get), "excursions_ps": exc}
    res[mesh] = {"rows": rows, "extremes": ext, "n_modes": G["n"], "pairs": int(G["ip"].size), "local": int(G["i0"].size)}
    print(mesh, "extremes", {L: (v["max"], v["min"]) for L, v in ext.items()}, flush=True)
(OUT / "rs4_filmw_static.json").write_text(json.dumps(res, indent=1))
