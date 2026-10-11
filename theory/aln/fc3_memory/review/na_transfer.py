"""Numerical audit, part 2: witness tensors on 5x5x3 / 7x7x5 / 9x9x6 (S5, E19-E21).

For every pp file listed below (phono3py 4.5.0 Rust interaction strengths of a witness tensor on the mesh's event
list), the independent operator of na_common.py gives tau and kappa by full eigendecomposition; ratios to the
mesh's reference are compared with the campaign's stored ratios. Diagnostics per tensor: lambda_min on H, cond on H,
share of N carried by the lowest 5 / 20 eigenpairs and by the largest single eigenpair, degeneracy-averaged
lifetime deviations (own grouping: same q, |df| <= 1e-6 THz), and the sensitivity of the tau ratio to the event cut
(3 / 3.5 / 4 sigma) and to the smearing width (sigma' = 0.08 / 0.09 / 0.10 THz, exact Gaussian re-weighting).
usage: python -B na_transfer.py m553|m775|m996 [out.json]
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import na_common as nc  # noqa: E402

RES = nc.RES
CASES = {
    "m553": {"dir": RES + "/m553/joint_pp", "ref": RES + "/m553/joint_pp/pp_r01_REF.npz",
             "stored": [RES + "/joint/verify_m553.json"],
             "names": ["cal1nz_J_c0_max", "cal1nz_J_c2_max", "orb0.1_J_c0_max", "orb0.1_J_c2_max"]},
    "m775": {"dir": RES + "/m775x", "ref": RES + "/m775x/pp_r01_REF.npz",
             "stored": [RES + "/m775x/mesh_transfer.json", RES + "/joint/verify_m775.json"],
             "names": ["cal1_c0_max", "cal1_c2_max", "cal1_c0_min", "cal1_c2_min", "cal2_c0_max", "cal2_c2_max",
                       "orb0.1_c0_max", "orb0.1_c2_max", "cal1nz_J_c0_max", "cal1nz_J_c2_max", "orb0.1_J_c0_max",
                       "orb0.1_J_c2_max"]},
    "m996": {"dir": RES + "/m996", "ref": RES + "/m996/pp_r01_REF.npz",
             "stored": [RES + "/m996/transfer_controls.json", RES + "/joint/verify_m996.json"],
             "names": ["alamode", "cal1_c0_max", "joint/cal1nz_J_c0_max", "joint/cal1nz_J_c2_max",
                       "joint/orb0.1_J_c0_max", "joint/orb0.1_J_c2_max"]},
}


def degeneracy_groups(op, tol=1e-6):
    grp = -np.ones(op.n, dtype=np.int64)
    key = np.lexsort((op.f, op.q_of))
    gid, prev_q, prev_f = -1, -1, None
    for i in key:
        if op.q_of[i] != prev_q or op.f[i] - prev_f > tol:
            gid += 1
        grp[i] = gid
        prev_q, prev_f = op.q_of[i], op.f[i]
    return grp


def analyse(op, g, want_eig=True):
    C = op.C(g)
    out = {}
    if want_eig:
        E = nc.eig_response(op, C)
        out = {"tau": [E["cols"]["x"]["tau_ps"], E["cols"]["z"]["tau_ps"]],
               "kappa": [E["cols"]["x"]["kappa"], E["cols"]["z"]["kappa"]],
               "lam_min_H": E["lam_min_H"], "cond_H": E["cond_H"], "null_overlap": E["null_overlap"],
               "N_frac_lowest5": [E["cols"][a]["N_frac_lowest"]["5"] for a in "xz"],
               "N_frac_lowest20": [E["cols"][a]["N_frac_lowest"]["20"] for a in "xz"],
               "N_top_pair": [E["cols"][a]["top_N_contributors"][0] for a in "xz"],
               "n_pairs_N50": [E["cols"][a]["n_pairs_for_N_50_90_99"][0] for a in "xz"]}
    out["r"] = op.diagC(g)
    del C
    return out


def variants(op, g):
    res = {}
    for ns in (3.0, 3.5):
        gg = np.where(np.abs(op.dl) <= ns * op.sigma, g, 0.0)
        Cc = op.C(gg)
        rr = nc.chol_response(op, Cc, np.trace(Cc) / op.n)
        res[f"cut{ns:g}"] = [rr["x"]["tau_ps"], rr["z"]["tau_ps"], rr["x"]["kappa"], rr["z"]["kappa"]]
        del Cc
    for s2 in (0.08, 0.09):
        wgt = (op.sigma / s2) * np.exp(-op.dl ** 2 / 2 * (1 / s2 ** 2 - 1 / op.sigma ** 2))
        gg = np.where(np.abs(op.dl) <= 4 * s2, g * wgt, 0.0)
        Cc = op.C(gg)
        rr = nc.chol_response(op, Cc, np.trace(Cc) / op.n)
        res[f"sigma{s2:g}"] = [rr["x"]["tau_ps"], rr["z"]["tau_ps"], rr["x"]["kappa"], rr["z"]["kappa"]]
        del Cc
    return res


def main():
    mk = sys.argv[1]
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"transfer_{mk}.json"
    cs = CASES[mk]
    t0 = time.perf_counter()
    op = nc.Operator(nc.MESHES[mk][0])
    stored = {}
    for s in cs["stored"]:
        stored.update(json.load(open(s)))
    grp = degeneracy_groups(op)
    da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
    P0 = np.load(cs["ref"])
    g0 = op.group_average(op.rates_from_pp(P0["pp"], P0["g0w"], float(P0["conv"])))
    A0 = analyse(op, g0)
    V0 = variants(op, g0)
    R = {"mesh": op.mesh, "REF": {k: v for k, v in A0.items() if k != "r"}, "REF_variants": V0, "witness": {}}
    print("REF", R["REF"]["tau"], R["REF"]["kappa"], time.perf_counter() - t0, flush=True)
    for nm in cs["names"]:
        path = cs["dir"] + "/" + (nm if "/" not in nm else nm.split("/")[0] + "/pp_r01_" + nm.split("/")[1])
        path = path if "/" in nm else cs["dir"] + f"/pp_r01_{nm}"
        path += ".npz"
        P = np.load(path)
        dg0w = float(np.max(np.abs(P["g0w"] - P0["g0w"])))
        g = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], float(P["conv"])))
        A = analyse(op, g)
        V = variants(op, g)
        key = nm.split("/")[-1]
        rec = {"g0w_maxdiff_vs_REF": dg0w,
               "tau_ratio": [A["tau"][i] / A0["tau"][i] for i in range(2)],
               "kappa_ratio": [A["kappa"][i] / A0["kappa"][i] for i in range(2)],
               "lifetime_degavg_minmax": [float((da(A["r"]) / da(A0["r"]) - 1).min()), float((da(A["r"]) / da(A0["r"]) - 1).max())],
               "lam_min_H_ratio": A["lam_min_H"] / A0["lam_min_H"], "cond_H": A["cond_H"],
               "N_frac_lowest5": A["N_frac_lowest5"], "N_frac_lowest20": A["N_frac_lowest20"],
               "N_top_pair": A["N_top_pair"], "n_pairs_N50": A["n_pairs_N50"],
               "variants_tau_ratio": {k: [V[k][0] / V0[k][0], V[k][1] / V0[k][1]] for k in V},
               "variants_kappa_ratio": {k: [V[k][2] / V0[k][2], V[k][3] / V0[k][3]] for k in V}}
        if key in stored:
            st = stored[key]
            rec["stored_tau_ratio"] = st["tau_ratio"]
            rec["stored_kappa_ratio"] = st["kappa_ratio"]
            rec["diff_tau_ratio"] = [rec["tau_ratio"][i] - st["tau_ratio"][i] for i in range(2)]
            rec["diff_kappa_ratio"] = [rec["kappa_ratio"][i] - st["kappa_ratio"][i] for i in range(2)]
            rec["stored_lifetime_minmax"] = st["lifetime_degavg_minmax"]
        R["witness"][key] = rec
        print(key, json.dumps({k: rec[k] for k in ("tau_ratio", "diff_tau_ratio", "kappa_ratio", "lifetime_degavg_minmax",
                                                   "variants_tau_ratio") if k in rec}), round(time.perf_counter() - t0), flush=True)
    R["t_total"] = time.perf_counter() - t0
    out_path.write_text(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
