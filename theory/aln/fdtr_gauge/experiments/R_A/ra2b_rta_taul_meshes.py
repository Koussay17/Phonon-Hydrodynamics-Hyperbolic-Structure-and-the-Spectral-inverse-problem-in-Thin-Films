"""R-A2b: mesh convergence and mode decomposition of the RTA kinetic nonlocal time tau_l(L) (objection A-2).

RTA with energy conservation, Matthiessen film rate |v_z|/L (c = 1), on the three converged phono3py meshes of
fc3-memory E27 (48x48x26, 56x56x30, 64x64x35; chunks read-only):
  tau_l = <tau (v_z^2 tau - s)^2> / s^2,   s = <v_z^2 tau>,   <.> capacity-weighted average,
  tau_mem = <v_z^2 tau^2> / s  (A1 value, checked).
Decomposition of tau_l = sum_i c_i tau_i (v_z,i^2 tau_i - s)^2 / (s^2 sum c):
  (i) by the plane q_z = 0 of the uniform mesh (v_z = 0 exactly, never boundary-regularised, weight 1/N_z),
  (ii) by |v_z| / v_max bins and by the bulk lifetime tau_3ph,
  (iii) the 'projection' part (the -s term alone: sum c tau s^2) against the 'direct' part.
If tau_l is carried by modes the regulator does not limit (v_z ~ 0), it inherits the infrared growth of the bulk
memory time and is not a converged film quantity. Variant: isotropic regulator |v|/L (all modes capped), to measure
how much of tau_l the grazing modes carry.
Output: out/ra2b_rta_taul_meshes.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
RTA = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/rta")
H_OVER_KB = 47.99243073366221
LS = (30, 100, 300, 500, 1000, 3000)
res = {}
for mesh in ("m484826", "m565630", "m646435"):
    D0 = RTA / mesh
    meta = json.loads((D0 / "meta.json").read_text())
    fs = sorted(f for f in D0.glob("chunk_*.npz") if ".tmp" not in f.name)
    assert len(fs) == meta["n_chunks"]
    D = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk", "q")}
    ok = (D["gamma"] > 0) & (D["f"] > 1e-3)
    tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, D["gamma"], 1)), 0.0)
    vz = np.abs(D["gv"][..., 2]) * 100.0
    mult = np.broadcast_to(D["w"][:, None], tau.shape).astype(float)
    xq = H_OVER_KB * np.where(ok, D["f"], 1.0) / 300.0
    cap = np.where(ok, mult * xq ** 2 * np.exp(xq) / np.expm1(xq) ** 2, 0.0)
    qz = np.abs(D["q"][:, 2])
    qz = np.minimum(qz, 1 - qz)
    plane = np.broadcast_to((qz < 1e-9)[:, None], tau.shape)
    vmax = vz.max()
    vabs = np.linalg.norm(D["gv"], axis=-1) * 100.0              # |v| for the isotropic (Casimir-type) variant
    row = {"mesh": meta["mesh"], "f_min_THz": float(D["f"][D["f"] > 1e-3].min()),
           "tau3ph_max_ps": float(tau.max()), "plane_weight": float(mult[plane].sum() / mult.sum() / 1.0 * 1.0)}
    for L in LS:
        tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vz / (L * 1e-9) * 1e-12), 0.0)
        S = cap.sum()
        s = (cap * vz ** 2 * tl).sum() / S
        term = cap * tl * (vz ** 2 * tl - s) ** 2 / (S * s * s)
        proj = cap * tl * s * s / (S * s * s)
        taul = float(term.sum())
        bins = {}
        edges = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.01)
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (vz / vmax >= lo) & (vz / vmax < hi) & ok
            bins[f"{lo}-{hi}"] = float(term[m].sum() / taul)
        tb = {}
        for lo, hi in ((0, 10), (10, 100), (100, 1000), (1000, 1e9)):
            m = (tau >= lo) & (tau < hi) & ok
            tb[f"{lo}-{hi:g}ps"] = float(term[m].sum() / taul)
        # regulator variant |v|/L (caps the grazing modes too; not the c-axis film rate of A1)
        ta = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vabs / (L * 1e-9) * 1e-12), 0.0)
        sa = (cap * vz ** 2 * ta).sum() / S
        taul_vabs = float((cap * ta * (vz ** 2 * ta - sa) ** 2).sum() / (S * sa * sa))
        tmem_vabs = float((cap * vz ** 2 * ta ** 2).sum() / (cap * vz ** 2 * ta).sum())
        row[str(L)] = {"tau_mem": float((cap * vz ** 2 * tl ** 2).sum() / (cap * vz ** 2 * tl).sum()),
                       "tau_mem_vabs_regulator": tmem_vabs, "taul_RTA_vabs_regulator": taul_vabs,
                       "taul_RTA": taul, "share_qz0_plane": float(term[plane].sum() / taul),
                       "share_projection_term": float(proj.sum() / taul), "share_by_vz_over_vmax": bins,
                       "share_by_tau3ph": tb}
    res[mesh] = row
    print(mesh, "fmin", round(row["f_min_THz"], 3), {L: (round(row[str(L)]["tau_mem"], 3), round(row[str(L)]["taul_RTA"], 2),
                                                        round(row[str(L)]["share_qz0_plane"], 3))
                                                    for L in LS}, flush=True)
(OUT / "ra2b_rta_taul_meshes.json").write_text(json.dumps(res, indent=1))
for L in LS:
    v = [res[m][str(L)]["taul_RTA"] for m in res]
    va = [res[m][str(L)]["taul_RTA_vabs_regulator"] for m in res]
    print(L, "taul_RTA |v|/L regulator", [round(x, 2) for x in va], "tau_mem |v|/L", round(res["m646435"][str(L)]["tau_mem_vabs_regulator"], 3))
    print(L, "taul_RTA meshes", [round(x, 2) for x in v], "steps %", [round(100 * (v[i + 1] / v[i] - 1), 2) for i in range(2)],
          "plane share 64:", round(res["m646435"][str(L)]["share_qz0_plane"], 3),
          "by vz:", {k: round(x, 3) for k, x in res["m646435"][str(L)]["share_by_vz_over_vmax"].items()},
          "by tau:", {k: round(x, 3) for k, x in res["m646435"][str(L)]["share_by_tau3ph"].items()})
