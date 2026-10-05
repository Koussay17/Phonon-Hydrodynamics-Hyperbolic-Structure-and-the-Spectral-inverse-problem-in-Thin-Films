"""Attack 4 (part 1): observables the campaign does not fix, evaluated on the witnesses; and the
information content of multi-temperature lifetime data with temperature-independent |Phi|^2.

(i) For reference and witnesses: K(z)/K0 at z_phys = 1/(10 ns), 1/ns, 1/(100 ps), 1/(10 ps)
    (z = z_phys/(4 pi) in operator units), second moment S = b^T C b, accumulation
    sum_{f_mu < fc} b_mu x_mu / K at fc = 5, 10, 15 THz, and the implied lifetimes and K at
    T = 100, 200, 500, 800 K if |Phi|^2 delta is temperature independent (w = g/Bf(300)).
(ii) Rank / singular values of the stacked lifetime-constraint matrix over event orbits
    [W_gamma(T_k) diag(Bf(T_k))]/r(T_k) for growing temperature sets.
Output: review/cx_attack4_observables.json
"""
import json
import sys

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI, their_ev_orb, load_witness  # noqa: E402

M = AlN()
o = M.orbits()
gref = M.g_ref()
forb = gref < 1e-12 * gref.max()
ev_orb_theirs, _ = their_ev_orb()
temps_other = [100.0, 200.0, 500.0, 800.0]
MT = {T: AlN(T=T) for T in temps_other}
zs_phys = {"1/10ns": 1e-4, "1/ns": 1e-3, "1/100ps": 1e-2, "1/10ps": 1e-1}


def observables(g):
    C = M.C(g)
    gam = np.trace(C) / M.n
    cf = sla.cho_factor(C + gam * np.outer(M.ehat, M.ehat), lower=True)
    rec = {}
    for c in (0, 2):
        b = M.b[:, c]
        x = sla.cho_solve(cf, b)
        K = float(b @ x)
        d = {"K": K, "tau_ps": float(x @ x) / K / FOURPI, "S2": float(b @ C @ b)}
        for lab, zp in zs_phys.items():
            z = zp / FOURPI
            xz = np.linalg.solve(C + z * np.eye(M.n), b)
            d[f"K({lab})/K"] = float(b @ xz) / K
        for fc in (5.0, 10.0, 15.0):
            mask = M.eps < fc
            d[f"acc<{fc:g}THz"] = float(b[mask] @ x[mask]) / K
        rec[c] = d
    # implied data at other temperatures (T-independent |Phi|^2 delta)
    w = g / M.Bf
    for T, MM in MT.items():
        gT = w * MM.Bf
        mo = MM.moments(gT, (0, 2))
        rec[f"T{int(T)}"] = {"r": MM.diag(gT), "K0": mo[0]["K"], "K2": mo[2]["K"], "tau0": mo[0]["tau_ps"], "tau2": mo[2]["tau_ps"]}
    return rec


ref = observables(gref)
out = {"reference": {c: ref[c] for c in (0, 2)}}
out["reference_other_T"] = {k: {"K0": v["K0"], "K2": v["K2"], "tau0_ps": v["tau0"], "tau2_ps": v["tau2"]}
                            for k, v in ref.items() if str(k).startswith("T")}
for nm in ["full_c0_max_08_vertex", "full_c0_max_04_logn3.0", "full_c2_max_01_logn0.5", "full_c0_min_02_logn1.0"]:
    g = load_witness(nm, ev_orb_theirs)
    ob = observables(g)
    rec = {c: {k: (v / ref[c][k] if k in ("K", "S2") else v) for k, v in ob[c].items()} for c in (0, 2)}
    for c in (0, 2):
        rec[c]["S2_over_ref"] = rec[c].pop("S2")
        rec[c]["K_over_ref"] = rec[c].pop("K")
    for T in temps_other:
        key = f"T{int(T)}"
        rr = ob[key]["r"] / ref[key]["r"]
        rec[key] = {"lifetime_ratio_min": float(rr.min()), "lifetime_ratio_max": float(rr.max()),
                    "lifetime_rel_dev_median": float(np.median(np.abs(rr - 1))),
                    "K_x_ratio": ob[key]["K0"] / ref[key]["K0"], "K_z_ratio": ob[key]["K2"] / ref[key]["K2"],
                    "tau_x_ps": ob[key]["tau0"], "tau_z_ps": ob[key]["tau2"]}
    out[nm] = rec
    print(nm, json.dumps(rec, default=float), flush=True)

# (ii) rank of stacked lifetime constraints over event orbits (and over allowed orbits only)
Porb = sp.csr_matrix((np.ones(M.m), (o["ev_orb"], np.arange(M.m))), shape=(o["n_ev_orb"], M.m))
w0 = gref / M.Bf
allowed_orb = np.bincount(o["ev_orb"], weights=(~forb).astype(float), minlength=o["n_ev_orb"]) > 0


def rowblock(T):
    MM = M if T == 300.0 else AlN(T=T)
    gT = w0 * MM.Bf
    r = MM.diag(gT)
    Wr = MM.W[o["mode_reps"], :].multiply(gT[None, :]) @ Porb.T     # derivative wrt log w per orbit
    return (Wr.toarray() / r[o["mode_reps"]][:, None])


sets = {"1T": [300.0], "2T": [300.0, 800.0], "3T": [100.0, 300.0, 800.0], "5T": [100.0, 200.0, 300.0, 500.0, 800.0],
        "10T": list(np.geomspace(50, 1500, 10)), "25T": list(np.geomspace(30, 3000, 25)), "40T": list(np.geomspace(20, 5000, 40))}
rk = {}
cache = {}
for name, Ts in sets.items():
    blocks = []
    for T in Ts:
        T = float(T)
        if T not in cache:
            cache[T] = rowblock(T)
        blocks.append(cache[T])
    A = np.vstack(blocks)[:, allowed_orb]
    s = np.linalg.svd(A, compute_uv=False)
    rk[name] = {"rows": int(A.shape[0]), "cols": int(A.shape[1]),
                "rank_rel_1e-3": int(np.sum(s > 1e-3 * s[0])), "rank_rel_1e-6": int(np.sum(s > 1e-6 * s[0])),
                "rank_rel_1e-10": int(np.sum(s > 1e-10 * s[0])), "rank_rel_1e-13": int(np.sum(s > 1e-13 * s[0]))}
    print(name, rk[name], flush=True)
out["multiT_rank_allowed_orbits"] = rk
out["n_allowed_orbits"] = int(allowed_orb.sum())
(CAMP / "review" / "cx_attack4_observables.json").write_text(json.dumps(out, indent=1, default=float))
print("done")
