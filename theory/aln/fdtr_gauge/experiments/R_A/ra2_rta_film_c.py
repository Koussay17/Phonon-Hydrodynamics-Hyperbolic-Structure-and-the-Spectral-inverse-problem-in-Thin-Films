"""R-A2: RTA film route on the converged phono3py mesh 64x64x35 (fc3-memory E27 chunks, read-only) for the regulator
normalisation c in {1/2, 2/3, 3/4, 1, 4/3, 3/2, 2} (rate c |v_z| / L), and the RTA kinetic nonlocal time.

Per mode (irreducible points, multiplicity w): 1/tau_L = 1/tau_3ph + c |v_z| / L.
  tau_mem,z(L; c) = sum W tau_L^2 / sum W tau_L,   W = mode_kappa_zz / tau_3ph          (A1, a1_base_values.py),
  kappa_a(L; c)   = kappa_a,bulk sum W_a tau_L / sum W_a tau_3ph                          (a = x, z),
so that tau(L; c) = tau(L/c; 1) identically (checked).
RTA kinetic nonlocal time with energy conservation (ra_kinetic.py; operator R - (R e)(R e)^T / (e^T R e)):
  tau_l = <tau (v_z^2 tau - <v_z^2 tau>)^2> / <v_z^2 tau>^2 ,   <.> = capacity-weighted mode average,
and the kinematic factor gamma_g = <v_z^4> / <v_z^2>^2 - 1. |v_z| and the capacity are invariant on the star of each
irreducible point (point group 6mm with time reversal), so the irreducible sums are exact for these z quantities.
Base values for the FDTR box (event operator at c = 1, A1; c-dependence from this RTA route):
  tau_R(L; c) = tau_R,event(L) tau_RTA(L/c) / tau_RTA(L),  lam_z(L; c) = lam_z,event(L) kz_RTA(L/c) / kz_RTA(L),
  lam_r(L; c) = kx_RTA(L/c)  (A1 adopted the RTA in-plane value).
Output: out/ra2_rta_film_c.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
RTA = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/rta/m646435")
A1 = Path("D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/experiments/A/out/a1_base_values.json")
CS = (0.5, 2 / 3, 0.75, 1.0, 4 / 3, 1.5, 2.0)
LB = (30, 100, 300, 500)
H_OVER_KB = 47.99243073366221        # h / k_B in K per THz

meta = json.loads((RTA / "meta.json").read_text())
fs = sorted(f for f in RTA.glob("chunk_*.npz") if ".tmp" not in f.name)
assert len(fs) == meta["n_chunks"]
D = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk")}
assert int(D["w"].sum()) == meta["sum_weights"]
T = 300.0
ok = (D["gamma"] > 0) & (D["f"] > 1e-3)
tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, D["gamma"], 1)), 0.0)          # ps
vz = np.abs(D["gv"][..., 2]) * 100.0                                            # THz*A -> m/s
mult = np.broadcast_to(D["w"][:, None], tau.shape).astype(float)
xq = H_OVER_KB * np.where(ok, D["f"], 1.0) / T
cap = np.where(ok, mult * xq ** 2 * np.exp(xq) / np.expm1(xq) ** 2, 0.0)       # capacity weight x multiplicity
kz = np.where(ok, D["mk"][..., 2], 0.0)
kx = np.where(ok, D["mk"][..., 0], 0.0)
Wz = np.where(ok, kz / np.where(ok, tau, 1), 0.0)
Wx = np.where(ok, kx / np.where(ok, tau, 1), 0.0)
kz_bulk, kx_bulk = kz.sum() / meta["sum_weights"], kx.sum() / meta["sum_weights"]


def film(Lnm, c):
    tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + c * vz / (Lnm * 1e-9) * 1e-12), 0.0)
    out = {"tau_z": float((Wz * tl ** 2).sum() / (Wz * tl).sum()),
           "kappa_z": float(kz_bulk * (Wz * tl).sum() / (Wz * tau).sum()),
           "kappa_x": float(kx_bulk * (Wx * tl).sum() / (Wx * tau).sum())}
    # kinetic times from capacity weights and v_z (independent of mode_kappa); same tau_mem as a check
    s = (cap * vz ** 2 * tl).sum() / cap.sum()                      # <v^2 tau>
    out["tau_z_capweights"] = float((cap * vz ** 2 * tl ** 2).sum() / (cap * vz ** 2 * tl).sum())
    out["taul_RTA"] = float((cap * tl * (vz ** 2 * tl - s) ** 2).sum() / cap.sum() / s ** 2)
    out["a_m2s"] = float(s * 1e-12)                                 # <v^2 tau> in m^2/s (diffusivity)
    return out


gam_g = float((cap * vz ** 4).sum() * cap.sum() / (cap * vz ** 2).sum() ** 2 - 1.0)
a1 = json.loads(A1.read_text())
res = {"mesh": meta["mesh"], "c_values": list(CS), "gamma_g_RTA64": gam_g, "film": {}, "base_box": {}}
Lall = sorted({round(L / c, 6) for L in LB for c in CS} | {1000.0})
for L in Lall:
    res["film"][f"{L:.6g}"] = film(L, 1.0)
# identity tau(L; c) = tau(L/c; 1)
idmax = 0.0
for L in LB:
    for c in CS:
        a = film(L, c)
        bq = res["film"][f"{round(L / c, 6):.6g}"]
        idmax = max(idmax, max(abs(a[k] / bq[k] - 1) for k in ("tau_z", "kappa_z", "kappa_x", "taul_RTA")))
res["identity_c_vs_L_over_c_max_rel"] = idmax
# reproduction of A1 RTA values at c = 1
res["check_A1_rta_max_rel"] = max(abs(res["film"][f"{L:.6g}"]["tau_z"] / a1["base"][f"{L}_rta"]["tau_R_ps"] - 1)
                                  for L in LB)
res["check_capweights_vs_modekappa_tau_max_rel"] = max(abs(v["tau_z_capweights"] / v["tau_z"] - 1)
                                                       for v in res["film"].values())
for L in LB:
    ev = a1["base"][str(L)]
    f1 = res["film"][f"{float(L):.6g}"]
    for c in CS:
        fc = res["film"][f"{round(L / c, 6):.6g}"]
        res["base_box"][f"{L}|{c:.6g}"] = {
            "tau_R_ps": ev["tau_R_ps"] * fc["tau_z"] / f1["tau_z"],
            "lam_z": ev["lam_z"] * fc["kappa_z"] / f1["kappa_z"],
            "lam_r": fc["kappa_x"],
            "tau_R_ps_rta": fc["tau_z"], "lam_z_rta": fc["kappa_z"],
            "L_eff_nm": L / c}
print(json.dumps({k: v for k, v in res.items() if k not in ("film", "base_box")}, indent=1))
for L in LB:
    print(L, {f"{c:.3g}": round(res["base_box"][f"{L}|{c:.6g}"]["tau_R_ps"], 3) for c in CS})
print("taul_RTA(64x64x35):", {k: round(v["taul_RTA"], 2) for k, v in res["film"].items()
                              if k in ("30", "100", "300", "500", "1000")})
(OUT / "ra2_rta_film_c.json").write_text(json.dumps(res, indent=1))
