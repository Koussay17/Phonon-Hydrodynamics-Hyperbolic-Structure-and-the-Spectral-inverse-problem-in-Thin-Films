"""RS7: band dependence of the FDTR-visible delta of the REF slabs (exact 80-frequency data, black walls, wall port,
note-16 stack). The GK-film fit (tau_l free; N2lam, N3; 0.01 and 0.1 deg) is repeated on sub-bands by giving the
out-of-band points zero weight in the whitened residual and Jacobian (au5_fisher1d.whiten masked); same multi-start
and warm-start fitting as rs2b. A delta that drifts with the band would mean that the GK form is only an effective fit.
Output: out/rs7_band.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs2b_fits as fb  # noqa: E402

a5, a6 = fb.a5, fb.a6
F = a5.FREQ
BANDS = {"full_10k_200M": (1e4, 2e8), "low_10k_50M": (1e4, 5e7), "high_1M_200M": (1e6, 2e8), "top_10M_200M": (1e7, 2e8)}
_orig_whiten = a5.whiten
res = {}
for op in ("full", "rta"):
    f, d, Y, Th, _ = fb.load("m775", "film2", op, "all", "REF", 0)
    Yx, Tx, _ = fb.at_fdtr(f, Y, Th, "all")
    stat = json.loads((OUT / "rs2_m775_film2_static_REF.json").read_text())
    for i, dd in enumerate(d):
        L = f"{dd * 1e9:.0f}"
        kr = stat[f"reg_full_QDQ_{L}"] if op == "full" else stat[f"reg_rta_V1_{L}"]
        base = {"tauR": kr["tauR"], "taul": kr["taul"], "lamz": kr["K"] if op == "full" else kr["lam"], "C": stat["C_tot"]}
        for bname, (lo, hi) in BANDS.items():
            m = ((F >= lo * 0.999) & (F <= hi * 1.001)).astype(float)
            w = np.concatenate([m, m])
            a5.whiten = lambda Jc, sp, w=w: _orig_whiten(Jc, sp) * w[:, None]
            try:
                r_, _ = fb.fits_one(Yx[i], Tx[i], dd, base, "wall", "note16", extra=False)
            finally:
                a5.whiten = _orig_whiten
            res[f"{op}|{L}|{bname}"] = {k: {"delta_fit_ps": v["delta_fit_ps"], "dchi2": v["dchi2"], "lamz_fit": v["lamz_fit"],
                                            "n_points": int(m.sum())} for k, v in r_.items()}
            print(op, L, bname, {k: (round(v["delta_fit_ps"], 2), round(v["dchi2"], 4)) for k, v in r_.items()}, flush=True)
(OUT / "rs7_band.json").write_text(json.dumps(res, indent=1))
