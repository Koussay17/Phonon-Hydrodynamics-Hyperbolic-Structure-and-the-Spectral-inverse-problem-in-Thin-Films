"""E23c: memory moment above a fixed frequency cutoff, tau_>(fc) = sum_{f>fc} kappa_mode tau / sum_all kappa_mode,
compared across meshes; if the mesh dependence of tau_mem is infrared, tau_>(fc) converges for fixed fc.
usage: python -B ir_cutoff.py <out.json> label=kappa.hdf5 ...
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
        r[nm] = {str(fc): float((k * tau)[ok & (fr > fc)].sum() / K) for fc in (0.0, 1.5, 2.0, 2.5, 3.0, 4.0)}
    res[lab] = r
    print(lab, json.dumps(r))
json.dump(res, open(sys.argv[1], "w"), indent=1)
