"""E23 (T1): infrared behaviour of the memory moment at RTA level from phono3py per-mode data.
RTA: C diagonal, x = b/r, K = sum b^2/r, N = sum b^2/r^2, so tau_mem(RTA) = sum kappa_mode tau / sum kappa_mode
with tau = 1/(4 pi gamma) ps (gamma = phono3py half width in THz). Also: share of N from modes below f, low-f
exponent of gamma, and f_min of each mesh.
usage: python -B ir_rta_moment.py <out.json> label=kappa.hdf5 ...
"""
import json
import sys

import h5py
import numpy as np

out, specs = sys.argv[1], sys.argv[2:]
res = {}
for s in specs:
    lab, path = s.split("=", 1)
    f = h5py.File(path, "r")
    T = f["temperature"][:]
    it = int(np.argmin(np.abs(T - 300)))
    g = f["gamma"][it]
    mk = f["mode_kappa"][it]
    fr = f["frequency"][:]
    ok = g > 0
    tau = np.where(ok, 1.0 / (4 * np.pi * np.where(ok, g, 1)), 0.0)
    r = {"mesh": f["mesh"][:].tolist(), "kappa": f["kappa"][it][[0, 2]].tolist()}
    for c, nm in ((0, "x"), (2, "z")):
        k = mk[..., c]
        K, N = float(k[ok].sum()), float((k * tau)[ok].sum())
        r[f"kappa_sum_{nm}"] = K
        r[f"tau_mem_{nm}"] = N / K
        order = np.argsort(fr[ok])
        cumN = np.cumsum((k * tau)[ok][order]) / N
        frs = fr[ok][order]
        r[f"N_share_below_{nm}"] = {str(x): float(cumN[np.searchsorted(frs, x) - 1]) if np.searchsorted(frs, x) > 0 else 0.0
                                    for x in (1.0, 2.0, 3.0, 5.0)}
    pos = fr > 1e-3
    r["f_min_THz"] = float(fr[pos].min())
    lo = pos & (fr < 3.0) & ok
    if lo.sum() > 5:
        p = np.polyfit(np.log(fr[lo]), np.log(g[lo]), 1)
        r["gamma_lowf_exponent"] = float(p[0])
        r["n_lowf"] = int(lo.sum())
    res[lab] = r
    print(lab, json.dumps({k: v for k, v in r.items() if not isinstance(v, dict)}), flush=True)
json.dump(res, open(out, "w"), indent=1)
