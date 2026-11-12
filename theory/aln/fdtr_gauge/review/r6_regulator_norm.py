"""Sensitivity of the film memory time to the normalisation of the heuristic boundary rate c |v_z| / L (c = 1 in E26/E30).
RTA 64x64x35 mode data (read-only). Reports tau_mem,z(L) and kappa_z(L) for c = 0.5, 0.75, 1, 4/3, 2."""
import json
from pathlib import Path
import numpy as np
D = Path(r"D:\ResearchLab\orchestration\campaigns\20261005-fc3-memory\results_c\rta\m646435")
meta = json.loads((D / "meta.json").read_text())
fs = sorted(f for f in D.glob("chunk_*.npz") if ".tmp" not in f.name)
dat = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk")}
ok = (dat["gamma"] > 0) & (dat["f"] > 1e-3)
tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, dat["gamma"], 1)), 0.0)
vz = np.abs(dat["gv"][..., 2]) * 100.0
kz = np.where(ok, dat["mk"][..., 2], 0.0)
wgt = np.where(ok, kz / np.where(ok, tau, 1), 0.0)
out = {}
for c in (0.5, 0.75, 1.0, 4/3, 2.0):
    row = {}
    for L in (30, 100, 300, 500):
        tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + c * vz / (L * 1e-9) * 1e-12), 0.0)
        row[str(L)] = {"tau_mem_ps": float((wgt * tl**2).sum() / (wgt * tl).sum()),
                       "kappa_z": float((wgt * tl).sum() / meta["sum_weights"])}
    out[f"c={c:.3f}"] = row
print(json.dumps(out, indent=1))
Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1))
