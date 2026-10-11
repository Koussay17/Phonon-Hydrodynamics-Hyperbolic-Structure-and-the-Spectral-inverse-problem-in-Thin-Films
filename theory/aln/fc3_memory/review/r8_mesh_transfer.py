"""Review R8: mesh transfer of the witnesses. The fc3 tensors found on the 5x5x3 event set are evaluated on the
7x7x5 event set (memory-gap export, sigma 0.1 THz, 4 sigma, 300 K) with the independent event operator of r5
(conserving convention, orbit-averaged rates) and, optionally, the phono3py-form detuning convention of r6d.
Reported per witness: kappa_xx, kappa_zz and tau_x, tau_z ratios to the 7x7x5 reference; mode lifetimes
(r = diag C) after degeneracy averaging: min, max, rms relative change and fraction within 10 %.
usage: python -B r8_mesh_transfer.py <events.npz> <pp_ref.npz> <outjson> <conv: cons|p3form> name=pp.npz ...
"""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

evp, ppref, outp, conv = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4]
specs = sys.argv[5:]
saved = list(sys.argv)
t0 = time.time()
mod = "r6d_detuning_convention.py" if conv == "p3form" else "r5_event_operator.py"
sys.argv = ["x", evp, ppref, str(outp.with_name(outp.stem + "_dummy.json"))]
sp_ = importlib.util.spec_from_file_location("m", str(Path(__file__).with_name(mod)))
M = importlib.util.module_from_spec(sp_)
sp_.loader.exec_module(M)
sys.argv = saved
r5 = M.r5 if conv == "p3form" else M
print("setup", time.time() - t0, flush=True)
n, nb = r5.n, r5.Z["freqs"].shape[1]
qv = np.nonzero(r5.valid)[0] // nb
f = r5.eps
grp = -np.ones(n, int)
g = 0
for q in np.unique(qv):
    ii = np.nonzero(qv == q)[0]
    ii = ii[np.argsort(f[ii])]
    st = 0
    for k in range(1, len(ii) + 1):
        if k == len(ii) or f[ii[k]] - f[ii[k - 1]] > 1e-6:
            grp[ii[st:k]] = g
            g += 1
            st = k
da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]


def evaluate(ppfile):
    if conv == "p3form":
        C = M.op_p3form(M.amp(ppfile))
        kx, tx, kz, tz = M.kt(C)
    else:
        C = r5.operator(r5.rates(ppfile), "cons")
        w, V = np.linalg.eigh(C)
        eh = r5.e_vec / np.linalg.norm(r5.e_vec)
        keep = np.ones(n, bool)
        keep[np.argmax(np.abs(V.T @ eh))] = False
        out = []
        for col in (0, 2):
            bc = V[:, keep].T @ r5.b_all[:, col]
            K = float(np.sum(bc ** 2 / w[keep]))
            N = float(np.sum(bc ** 2 / w[keep] ** 2))
            out += [r5.kf * K, N / K / (4 * np.pi)]
        kx, tx, kz, tz = out
    return {"kappa_x": kx, "tau_x": tx, "kappa_z": kz, "tau_z": tz, "r": np.diag(C).copy()}


ref = evaluate(ppref)
res = {"conv": conv, "events": evp, "REF": {k: v for k, v in ref.items() if k != "r"}}
print("REF", res["REF"], flush=True)
for spec in specs:
    name, path = spec.split("=", 1)
    e = evaluate(path)
    lr = da(e["r"]) / da(ref["r"]) - 1
    res[name] = {"kappa_ratio": [e["kappa_x"] / ref["kappa_x"], e["kappa_z"] / ref["kappa_z"]],
                 "tau_ratio": [e["tau_x"] / ref["tau_x"], e["tau_z"] / ref["tau_z"]],
                 "tau_ps": [e["tau_x"], e["tau_z"]],
                 "lifetime_degavg": {"min": float(lr.min()), "max": float(lr.max()), "rms": float(np.sqrt(np.mean(lr ** 2))),
                                     "frac_within_10pct": float(np.mean(np.abs(lr) <= 0.1)),
                                     "frac_within_20pct": float(np.mean(np.abs(lr) <= 0.2))}}
    print(name, json.dumps(res[name]), flush=True)
outp.write_text(json.dumps(res, indent=1))
