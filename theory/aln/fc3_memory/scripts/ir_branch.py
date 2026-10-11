"""E23b: per-branch low-frequency structure of the RTA memory moment (phono3py per-mode data at 300 K).
For each of the three acoustic branches (band index 0,1,2 sorted by frequency at each q): exponent n of gamma ~ f^n
over frequency windows; cumulative share of N_x, N_z, K_x, K_z from modes below f.
usage: python -B ir_branch.py <kappa.hdf5> <out.json>
"""
import json
import sys

import h5py
import numpy as np

f = h5py.File(sys.argv[1], "r")
it = int(np.argmin(np.abs(f["temperature"][:] - 300)))
g, mk, fr, w = f["gamma"][it], f["mode_kappa"][it], f["frequency"][:], f["weight"][:]
q = f["qpoint"][:]
ok = (g > 0) & (fr > 1e-3)
tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, g, 1)), 0)
out = {"mesh": f["mesh"][:].tolist(), "branches": {}, "cumulative": {}}
for b in range(3):
    m = ok[:, b]
    fb, gb = fr[m, b], g[m, b]
    rec = {}
    for lo, hi in ((0, 1.0), (0, 1.5), (1.0, 2.0), (1.5, 3.0), (2.0, 4.0)):
        s = (fb >= lo) & (fb < hi)
        if s.sum() >= 4:
            rec[f"{lo}-{hi}"] = {"n": float(np.polyfit(np.log(fb[s]), np.log(gb[s]), 1)[0]), "count": int(s.sum())}
    out["branches"][b] = rec
# direction-resolved: modes on the c axis (qx=qy=0) and in the basal plane (qz=0)
for name, sel in (("c_axis", (np.abs(q[:, 0]) < 1e-9) & (np.abs(q[:, 1]) < 1e-9)), ("basal", np.abs(q[:, 2]) < 1e-9)):
    for b in range(3):
        m = sel & ok[:, b]
        if m.sum() >= 4:
            fb, gb = fr[m, b], g[m, b]
            s = fb < 3.0
            if s.sum() >= 3:
                out["branches"].setdefault(f"{name}_{b}", {})["<3THz"] = {
                    "n": float(np.polyfit(np.log(fb[s]), np.log(gb[s]), 1)[0]), "count": int(s.sum()),
                    "f": fb[s].tolist(), "gamma": gb[s].tolist()}
for c, nm in ((0, "x"), (2, "z")):
    k = mk[..., c]
    N, K = (k * tau)[ok].sum(), k[ok].sum()
    o = np.argsort(fr[ok])
    fs = fr[ok][o]
    cn, ck = np.cumsum((k * tau)[ok][o]) / N, np.cumsum(k[ok][o]) / K
    out["cumulative"][nm] = {str(x): [float(cn[max(np.searchsorted(fs, x) - 1, 0)]), float(ck[max(np.searchsorted(fs, x) - 1, 0)])]
                             for x in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0)}
json.dump(out, open(sys.argv[2], "w"), indent=1)
print(json.dumps({"branches": {k: {kk: (vv["n"], vv["count"]) for kk, vv in v.items()} for k, v in out["branches"].items()},
                  "cumulative [N share, K share]": out["cumulative"]}, indent=0))
