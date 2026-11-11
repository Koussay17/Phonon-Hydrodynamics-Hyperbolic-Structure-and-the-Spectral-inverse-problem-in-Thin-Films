"""E26: boundary-regularised RTA memory time tau_mem(L) = sum kappa_mode tau_L / sum kappa_mode tau_L ... see below.
RTA with Matthiessen boundary rate: 1/tau_L = 1/tau_3ph + |v_dir| / L (grey diffusive boundary, HEURISTIC regulator;
dir = z for a c-oriented film, or |v| for an isotropic Casimir length). With x = b/r, K_L = sum c v_a^2 tau_L and
N_L = sum c v_a^2 tau_L^2, so tau_mem(L) = N_L / K_L and kappa_L = kappa-scale * K_L. Mode weights from mode_kappa / tau_3ph.
usage: python -B ir_film.py <out.json> label=kappa.hdf5 ...
"""
import json
import sys

import h5py
import numpy as np

Ls_nm = [30, 100, 300, 1000, 3000, 10000, 30000]
res = {"L_nm": Ls_nm}
for s in sys.argv[2:]:
    lab, p = s.split("=", 1)
    f = h5py.File(p, "r")
    it = int(np.argmin(np.abs(f["temperature"][:] - 300)))
    g, mk, gv = f["gamma"][it], f["mode_kappa"][it], f["group_velocity"][:]
    ok = g > 0
    tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, g, 1)), 0.0)           # ps
    vabs = np.linalg.norm(gv, axis=-1) * 100.0                                # THz A -> m/s
    vz = np.abs(gv[..., 2]) * 100.0
    r = {}
    for reg, v in (("vz", vz), ("vabs", vabs)):
        for c, nm in ((0, "x"), (2, "z")):
            w = np.where(ok, mk[..., c] / np.where(ok, tau, 1), 0.0)          # c v_a^2 (per ps)
            row = {"inf": [float((w * tau ** 2).sum() / (w * tau).sum()), 1.0]}
            K0 = (w * tau).sum()
            for L in Ls_nm:
                rb = v / (L * 1e-9) * 1e-12                                   # 1/ps
                tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + rb), 0.0)
                row[str(L)] = [float((w * tl ** 2).sum() / (w * tl).sum()), float((w * tl).sum() / K0)]
            r[f"{reg}_{nm}"] = row                                            # [tau_mem ps, kappa_L/kappa_bulk]
    res[lab] = r
    print(lab, json.dumps({k: {L: round(v[0], 2) for L, v in row.items()} for k, row in r.items() if k.startswith("vz")}))
json.dump(res, open(sys.argv[1], "w"), indent=1)
