"""Is the film-regularised cross-plane kernel a single pole in the FDTR band?

RTA mode data (phono3py 4.5.0, 300 K, 64x64x35; chunks of fc3-memory E27, read-only). Same regulator as rta_chunks.py:
1/tau_L = 1/tau + |v_z|/L. Weights w_i = mode_kappa_zz / tau (proportional to c v_z^2).
kappa_L(omega)/kappa_L(0) = sum w tau_L/(1 + i omega tau_L) / sum w tau_L   (exact multi-pole RTA kernel)
single pole: 1/(1 + i omega tau1), tau1 = <tau^2>/<tau> = tau_mem,z(L)   (what the twin uses as tau_R)
"""
import json
from pathlib import Path
import numpy as np

D = Path(r"D:\ResearchLab\orchestration\campaigns\20261005-fc3-memory\results_c\rta\m646435")
meta = json.loads((D / "meta.json").read_text())
fs = sorted(f for f in D.glob("chunk_*.npz") if ".tmp" not in f.name)
assert len(fs) == meta["n_chunks"]
dat = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk")}
ok = (dat["gamma"] > 0) & (dat["f"] > 1e-3)
tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, dat["gamma"], 1)), 0.0)          # ps
vz = np.abs(dat["gv"][..., 2]) * 100.0                                              # m/s
kz = np.where(ok, dat["mk"][..., 2], 0.0)
wgt = np.where(ok, kz / np.where(ok, tau, 1), 0.0)

out = {}
for L in (30, 100, 300, 500, 1000):
    tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vz / (L * 1e-9) * 1e-12), 0.0)
    m1, m2, m3, m4 = [(wgt * tl ** n).sum() for n in (1, 2, 3, 4)]
    tau1 = m2 / m1
    ratio = m3 * m1 / m2 ** 2          # = 1 for a single relaxation time
    row = {"tau1_ps": tau1, "moment_ratio_<t3><t>/<t2>^2": ratio, "max_tauL_ps_weighted99": None}
    # weighted quantiles of tau_L (weight w tau_L, i.e. share of kappa_z(L))
    sel = (wgt * tl) > 0
    order = np.argsort(tl[sel])
    cw = np.cumsum((wgt * tl)[sel][order]); cw /= cw[-1]
    row["tauL_quantiles_ps_(50,90,99,99.9)"] = [float(tl[sel][order][np.searchsorted(cw, q)]) for q in (0.5, 0.9, 0.99, 0.999)]
    for fMHz in (65.0, 200.0):
        om = 2 * np.pi * fMHz * 1e6 * 1e-12   # rad/ps
        K = (wgt * tl / (1 + 1j * om * tl)).sum() / m1
        Ks = 1 / (1 + 1j * om * tau1)
        row[f"{int(fMHz)}MHz"] = {
            "omega_tau1": om * tau1,
            "exact_ReK-1": float(K.real - 1), "singlepole_ReK-1": float(Ks.real - 1),
            "exact_ImK": float(K.imag), "singlepole_ImK": float(Ks.imag),
            "ratio_Re_dev_exact/single": float((K.real - 1) / (Ks.real - 1)),
            "Re_dev_exact/|ImK|": float(abs(K.real - 1) / abs(K.imag)),
        }
    row.pop("max_tauL_ps_weighted99")
    out[str(L)] = row
print(json.dumps(out, indent=1))
Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1))
