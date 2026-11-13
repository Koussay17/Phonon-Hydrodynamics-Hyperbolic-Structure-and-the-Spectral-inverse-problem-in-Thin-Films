"""SA5: accuracy of the AAA interpolation for the witnesses (the margins use AAA-interpolated witness two-ports).

(1) Exact 80-frequency two-ports of the range-setting witness cal3w_m3L100min (black walls, 300/500 nm) from the
    modal solver rs_core.slab_paired (verified independently to <= 1.7e-9 by sa1 at two frequencies of this tensor);
    AAA from the stored sub4 frequencies vs exact: two-port error, |Delta ln H|, delta_fit (N2lam, N3; 0.01 deg).
(2) Leave-out cross-validation on all 13 tensors: AAA rebuilt from the sub4 set without every second in-band support
    point, evaluated at the removed points (relative two-port error), and delta_fit with the thinned AAA.
Output: out/sa5_aaa_witness.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "R_slab"))
import numpy as np  # noqa: E402
import rs_core as rc  # noqa: E402
import sa_fdtr as fd  # noqa: E402
from sa2_fits import WIT, aaa_interp, base_par, data, load_sub  # noqa: E402
from scipy.interpolate import AAA  # noqa: E402

T0 = time.perf_counter()
RES = HERE / "out" / "sa5_aaa_witness.json"
F = fd.FREQ
NM = "cal3w_m3L100min"


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


res = json.loads(RES.read_text()) if RES.exists() else {}
exact_path = HERE / "out" / f"sa5_exact80_{NM}.npz"
if not exact_path.exists():
    G = rc.geometry("m775")
    C = rc.operator(G, "film2", NM)
    Bk = rc.blocks(G, C)
    del C
    Y = np.zeros((2, F.size, 2, 2), complex)
    for k, f in enumerate(F):
        outs, info = rc.slab_paired(Bk, 2j * np.pi * f, (300e-9, 500e-9))
        for i, (Y_, Th_, dg) in enumerate(outs):
            Y[i, k] = Y_
        if k % 10 == 0:
            log("exact", k, f"{f:.4g}")
    np.savez(exact_path, f=F, Y=Y)
Ye = np.load(exact_path)["Y"]
fs, ds, Ys, _ = load_sub(NM)
out = {}
for i, L in enumerate(("300", "500")):
    dd = float(L) * 1e-9
    b = base_par(L)
    Ya = aaa_interp(fs, Ys[i])
    # the stored sub4 values must equal the exact values at the shared frequencies
    idx = [int(np.argmin(np.abs(F - x))) for x in fs if x >= 9e3]
    ent = {"sub4_vs_exact_shared": float(np.max(np.abs(Ys[i][fs >= 9e3] - Ye[i][idx])) / np.max(np.abs(Ye[i]))),
           "aaa_max_rel_err_Y": float(np.max(np.abs(Ya - Ye[i])) / np.max(np.abs(Ye[i]))),
           "aaa_max_abs_dlnH": float(np.max(np.abs(data(Ya) - data(Ye[i]))))}
    for N in ("N2lam", "N3"):
        r1 = fd.fit_delta(data(Ye[i]), dd, b, 0.01, N)["delta_ps"]
        r2 = fd.fit_delta(data(Ya), dd, b, 0.01, N)["delta_ps"]
        ent[f"delta|{N}"] = {"exact": r1, "aaa": r2, "diff": r2 - r1}
    out[L] = ent
    log(L, ent)
res["exact_witness"] = out
RES.write_text(json.dumps(res, indent=1, default=float))
# (2) leave-out cross-validation
cv = {}
for nm in ["REF"] + WIT:
    fs_, ds_, Ys_, _ = load_sub(nm)
    inb = np.where(fs_ >= 9e3)[0]
    drop = inb[1:-1:2]
    keep = np.setdiff1d(np.arange(fs_.size), drop)
    row = {}
    for i, L in enumerate(("300", "500")):
        dd = float(L) * 1e-9
        b = base_par(L)
        err = 0.0
        Yt = np.zeros((F.size, 2, 2), complex)
        for a in range(2):
            for bb in range(2):
                r = AAA(2j * np.pi * fs_[keep], Ys_[i][keep, a, bb], rtol=1e-14)
                err = max(err, float(np.max(np.abs(r(2j * np.pi * fs_[drop]) - Ys_[i][drop, a, bb])) / np.max(np.abs(Ys_[i][:, a, bb]))))
                Yt[:, a, bb] = r(2j * np.pi * F)
        Yf = aaa_interp(fs_, Ys_[i])
        row[L] = {"leaveout_max_rel_err": err}
        for N in ("N2lam", "N3"):
            d_full = fd.fit_delta(data(Yf), dd, b, 0.01, N)["delta_ps"]
            d_thin = fd.fit_delta(data(Yt), dd, b, 0.01, N)["delta_ps"]
            row[L][f"delta_thin_minus_full|{N}"] = d_thin - d_full
    cv[nm] = row
    log("cv", nm, row)
    res["leave_out_cv"] = cv
    RES.write_text(json.dumps(res, indent=1, default=float))
res["max_leaveout_err"] = max(cv[n][L]["leaveout_max_rel_err"] for n in cv for L in cv[n])
res["max_abs_delta_thin_minus_full"] = max(abs(cv[n][L][k]) for n in cv for L in cv[n] for k in cv[n][L] if k.startswith("delta"))
RES.write_text(json.dumps(res, indent=1, default=float))
log("max leave-out err", res["max_leaveout_err"], "max |delta change|", res["max_abs_delta_thin_minus_full"])
