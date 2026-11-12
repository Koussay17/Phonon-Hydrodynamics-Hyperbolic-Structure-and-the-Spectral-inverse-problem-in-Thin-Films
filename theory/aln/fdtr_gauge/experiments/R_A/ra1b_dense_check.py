"""R-A1b: independent check of the A-2 quantities of ra1_film_event_ext.py on 7x7x5 (REF and one witness).

Route: dense operator 4 pi [C(g) + Q D Q] (C assembled by FastOps.C, not by the padded bincount buffer), reduced to
e-perp by a Householder basis and solved with LU (ra_kinetic.kinetic_times_dense): no regularisation constant, no
Cholesky. Also the kinetic times of the RTA operator with energy conservation on the same diagonal,
R - (R e)(R e)^T / (e^T R e), R = diag(4 pi (r + d)), for which M = sum u_i^2 / R_i exactly.
Output: out/ra1b_dense_check.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
sys.path.insert(0, str(FC3 / "scripts_c"))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ra_kinetic import kinetic_times_dense  # noqa: E402

ev = FC3 / "results_c/m775x/events_m775_s0.1.npz"
ppd = FC3 / "results_c/m775x/film2_pp"
M = FC3Class(str(ev), ppd / "pp_r01_REF.npz", allow_missing=True)
geom = M.geom
n = geom.n
vz = geom.vel[:, 2]
e = geom.e
ehat = geom.ehat
b = geom.b[:, 2]
ra1 = json.loads((HERE / "out" / "ra1_film_m775.json").read_text())
res = {}
for nm in ("REF", "orb0.1_m3L300max"):
    g = M.rates(np.load(ppd / f"pp_r01_{nm}.npz")["pp"])
    C = M.ops.C(g).copy()
    res[nm + "_energy_residual"] = float(np.linalg.norm(C @ e) / (np.linalg.norm(np.diag(C) * e)))
    r = np.diag(C).copy()
    for L in (30.0, 100.0, 300.0, 500.0, 1000.0):
        d = np.abs(vz) / (FOURPI * L * 10.0)
        Q = np.eye(n) - np.outer(ehat, ehat)
        Ch = FOURPI * (C + Q @ np.diag(d) @ Q)
        kt = kinetic_times_dense(Ch, e, vz)
        # energy-conserving RTA on the same diagonal
        R = FOURPI * (r + d)
        x = b / R
        K = float(b @ x)
        u = vz * x
        u = u - ehat * (ehat @ u)
        taul_rta = float(e @ e) * float(np.sum(u * u / R)) / K ** 2
        tau_rta = float(x @ x) / K
        key = f"{L:.6g}"
        ref = ra1[nm][key]
        res[f"{nm}|{key}"] = {"tau_dense": kt["tauR"], "taul_dense": kt["taul"],
                              "tau_ra1": ref["tau_z"], "taul_ra1": ref["taul_kin"],
                              "rel_diff_tau": abs(kt["tauR"] / ref["tau_z"] - 1),
                              "rel_diff_taul": abs(kt["taul"] / ref["taul_kin"] - 1),
                              "gamma_g": kt["gamma_g"], "tau_RTA_ec": tau_rta, "taul_RTA_ec": taul_rta}
        print(nm, key, res[f"{nm}|{key}"], flush=True)
res["max_rel_diff"] = max(v["rel_diff_taul"] for k, v in res.items() if "|" in k)
(HERE / "out" / "ra1b_dense_check.json").write_text(json.dumps(res, indent=1))
print("max rel diff taul", res["max_rel_diff"])
