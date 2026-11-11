"""E23d: branch decomposition of the infrared part (f < fc) of the RTA memory moment; band index 0,1 = TA, 2 = LA.
usage: python -B ir_branch_share.py <out.json> label=kappa.hdf5 ...
"""
import json
import sys

import h5py
import numpy as np

res = {}
for s in sys.argv[2:]:
    lab, p = s.split("=", 1)
    f = h5py.File(p, "r")
    it = int(np.argmin(np.abs(f["temperature"][:] - 300)))
    g, mk, fr = f["gamma"][it], f["mode_kappa"][it], f["frequency"][:]
    ok = (g > 0) & (fr > 1e-3)
    tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, g, 1)), 0)
    r = {}
    for c, nm in ((0, "x"), (2, "z")):
        k = mk[..., c]
        K = k[ok].sum()
        for fc in (2.0, 3.0):
            low = ok & (fr < fc)
            r[f"{nm}_below{fc}"] = {b: float((k * tau)[:, b][low[:, b]].sum() / K) for b in range(4)}
    res[lab] = r
    print(lab, json.dumps(r))
json.dump(res, open(sys.argv[1], "w"), indent=1)
