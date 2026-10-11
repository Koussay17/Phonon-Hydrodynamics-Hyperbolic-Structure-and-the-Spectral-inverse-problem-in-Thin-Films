"""Numerical audit, part 3 (S1): the exact linearisation of the event amplitudes at 7x7x5 (and 5x5x3).

Checks, reading the stored scaled matrices G_s = G M^-1/2 (memory-mapped, row chunks):
  (a) G_s (M^-1/2)^-1 c0 = Phi0 for each scaling tag (orb, cal1, cal1nz): the stored scaled matrices, i.e. the ones
      the joint optimisation (E20) actually used, reproduce the reference amplitudes;
  (b) at the witness coordinates u: amplitude-level comparison of the linear model, uc |Phi_row(u)|^2, with
      phono3py's Rust interaction strength pp of the rebuilt witness tensor at the same (representative) event,
      for all orbits with a non-degenerate block (blk = 1);
  (c) tau / kappa ratios of the linear model (block-average convention, own operator from na_common) versus the
      phono3py route (orbit average) and versus the stored model values of the optimiser history;
  (d) the campaign's event-orbit partition (symmetry.orbits, union-find) versus the group-average partition of
      na_common (orbit representative = lowest event index).
usage: python -B na_linear.py m775|m553 [out.json]
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

RES = Path(nc.RES)
TAGPRE = {"orb": "scaled_", "cal1": "scaled_cal1_", "cal1nz": "scaled_cal1nz_"}
WIT = {  # tag -> list of (u file, pp name on the mesh, history json for the stored model values or None)
    "cal1": [("m553/t3_cal/cal1_c0_max_y.npy", "cal1_c0_max", None), ("m553/t3_cal/cal1_c2_max_y.npy", "cal1_c2_max", None),
             ("m553/t3_cal/cal1_c0_min_y.npy", "cal1_c0_min", None), ("m553/t3_cal/cal1_c2_min_y.npy", "cal1_c2_min", None)],
    "cal1nz": [("joint/cal1nz_x_max_u.npy", "cal1nz_J_c0_max", "joint/cal1nz_x_max.json"),
               ("joint/cal1nz_z_max_u.npy", "cal1nz_J_c2_max", "joint/cal1nz_z_max.json")],
    "orb": [("joint/orb01_x_max_u.npy", "orb0.1_J_c0_max", "joint/orb01_x_max.json"),
            ("joint/orb01_z_max_u.npy", "orb0.1_J_c2_max", "joint/orb01_z_max.json")],
}


def main():
    mk = sys.argv[1]
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"linear_{mk}.json"
    t0 = time.perf_counter()
    rdir = RES / ("m775x" if mk == "m775" else "m553")
    ppdir = rdir if mk == "m775" else rdir / "joint_pp"
    meta = np.load(rdir / ("lin_meta.npz" if mk == "m775" else "linearization.npz"))
    Phi0, row_orb, blk, W_orb = meta["Phi0"], meta["row_orb"], meta["blk"], meta["W_orb"]
    allowed, reps, uc, c0 = meta["allowed"], meta["reps"], float(meta["uc"]), meta["c0"]
    g_ref_orb = meta["g_ref_orb"]
    n_orb = len(blk)
    op = nc.Operator(nc.MESHES[mk][0])
    img = op.images()
    my_id = img.min(axis=0)                                   # lowest event index of the orbit
    lookup = -np.ones(op.m, dtype=np.int64)
    lookup[reps] = np.arange(n_orb)
    camp_orb = lookup[my_id]
    R = {"mesh": op.mesh, "n_rows": int(len(Phi0)), "n_orbits": int(n_orb),
         "orbit_partition_agrees": bool((camp_orb >= 0).all() and len(np.unique(my_id)) == n_orb)}
    # campaign partition check: every event's orbit label equals the label of its lowest-index orbit member
    del img
    w_row = W_orb[row_orb] * uc / blk[row_orb]
    g_fixed = np.where(allowed, 0.0, g_ref_orb)

    def model_metrics(Phi):
        g_orb = g_fixed + np.bincount(row_orb, weights=w_row * (Phi.real ** 2 + Phi.imag ** 2), minlength=n_orb)
        C = op.C(g_orb[camp_orb])
        rr = nc.chol_response(op, C, np.trace(C) / op.n)
        return g_orb, [rr["x"]["tau_ps"], rr["z"]["tau_ps"]], [rr["x"]["kappa"], rr["z"]["kappa"]]

    _, tau_b0, kap_b0 = model_metrics(Phi0)
    R["model_base_block_avg"] = {"tau": tau_b0, "kappa": kap_b0}
    P0 = np.load(ppdir / "pp_r01_REF.npz")
    g0 = op.group_average(op.rates_from_pp(P0["pp"], P0["g0w"], float(P0["conv"])))
    C0 = op.C(g0)
    r0 = nc.chol_response(op, C0, np.trace(C0) / op.n)
    del C0
    tau_o0, kap_o0 = [r0["x"]["tau_ps"], r0["z"]["tau_ps"]], [r0["x"]["kappa"], r0["z"]["kappa"]]
    R["p3_base_orbit_avg"] = {"tau": tau_o0, "kappa": kap_o0}
    nondeg = allowed & (blk == 1)
    rows_nd = np.nonzero(nondeg[row_orb])[0]
    R["tags"] = {}
    for tag, wl in WIT.items():
        if mk == "m553" and tag == "cal1":
            continue
        pre = TAGPRE[tag]
        Mmh = np.load(RES / "m553" / f"{pre}Mmh.npy")
        u0 = np.linalg.solve(Mmh, c0)
        res_solve = float(np.linalg.norm(Mmh @ u0 - c0) / np.linalg.norm(c0))
        del Mmh
        U = np.column_stack([u0] + [np.load(RES / f) for f, _, _ in wl])
        Gr = np.load(rdir / f"{pre}Gr.npy", mmap_mode="r")
        Gi = np.load(rdir / f"{pre}Gi.npy", mmap_mode="r")
        GU = np.empty((Gr.shape[0], U.shape[1]), dtype=complex)
        for s in range(0, Gr.shape[0], 1500):
            e = min(s + 1500, Gr.shape[0])
            GU[s:e] = np.asarray(Gr[s:e]) @ U + 1j * (np.asarray(Gi[s:e]) @ U)
        del Gr, Gi
        d0 = np.abs(GU[:, 0] - Phi0)
        T = {"solve_residual": res_solve, "Gs_u0_vs_Phi0_relmax_rowwise": float(np.max(d0 / np.abs(Phi0))),
             "Gs_u0_vs_Phi0_rel_to_maxPhi0": float(d0.max() / np.abs(Phi0).max()),
             "Gs_u0_vs_Phi0_rowwise_q999": float(np.quantile(d0 / np.abs(Phi0), 0.999)), "witnesses": {}}
        for k, (f, nm, hist) in enumerate(wl, start=1):
            Phi = Phi0 + GU[:, k]
            g_orb, tau_b, kap_b = model_metrics(Phi)
            P = np.load(ppdir / f"pp_r01_{nm}.npz")
            # (b) amplitude level, non-degenerate blocks: row of orbit o is its representative event reps[o]
            ev = reps[row_orb[rows_nd]]
            same_tp = np.all(P["tp"][ev] == meta["rows_tp"][rows_nd], axis=1) & np.all(P["bands"][ev] == meta["rows_b"][rows_nd], axis=1)
            pp_lin = uc * np.abs(Phi[rows_nd]) ** 2
            pp_p3 = P["pp"][ev]
            big = pp_p3 > 1e-6 * pp_p3.max()
            rel = np.abs(pp_lin - pp_p3) / pp_p3
            gw = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], float(P["conv"])))
            Cw = op.C(gw)
            rw = nc.chol_response(op, Cw, np.trace(Cw) / op.n)
            del Cw
            tau_o, kap_o = [rw["x"]["tau_ps"], rw["z"]["tau_ps"]], [rw["x"]["kappa"], rw["z"]["kappa"]]
            rec = {"rows_nondegenerate": int(len(rows_nd)), "rep_triplet_and_bands_match": float(same_tp.mean()),
                   "pp_rel_max_all_nondeg": float(rel.max()), "pp_rel_max_nondeg_above_1e-6max": float(rel[big].max()),
                   "pp_rel_median": float(np.median(rel)),
                   "pp_abs_max_over_maxpp": float(np.max(np.abs(pp_lin - pp_p3)) / pp_p3.max()),
                   "model_tau_ratio_block": [tau_b[i] / tau_b0[i] for i in range(2)],
                   "model_kappa_ratio_block": [kap_b[i] / kap_b0[i] for i in range(2)],
                   "p3_tau_ratio_orbit": [tau_o[i] / tau_o0[i] for i in range(2)],
                   "p3_kappa_ratio_orbit": [kap_o[i] / kap_o0[i] for i in range(2)]}
            rec["model_minus_p3_tau_ratio"] = [rec["model_tau_ratio_block"][i] - rec["p3_tau_ratio_orbit"][i] for i in range(2)]
            rec["model_minus_p3_kappa_ratio"] = [rec["model_kappa_ratio_block"][i] - rec["p3_kappa_ratio_orbit"][i] for i in range(2)]
            rec["model_minus_p3_tau_ratio_rel"] = [rec["model_minus_p3_tau_ratio"][i] / rec["p3_tau_ratio_orbit"][i] for i in range(2)]
            if hist:
                H = json.load(open(RES / hist))
                mkey = "m775" if mk == "m775" else "m553"
                rec["history_best_model_tau_ratio"] = H["best"]["mesh"][mkey]["tau_ratio"]
                rec["history_best_model_kappa_ratio"] = H["best"]["mesh"][mkey]["kappa_ratio"]
                rec["own_model_minus_history"] = [rec["model_tau_ratio_block"][i] - H["best"]["mesh"][mkey]["tau_ratio"][i] for i in range(2)]
            T["witnesses"][nm] = rec
            print(tag, nm, json.dumps({k2: rec[k2] for k2 in ("pp_rel_max_nondeg_above_1e-6max", "pp_rel_median",
                                                              "model_minus_p3_tau_ratio", "model_minus_p3_kappa_ratio")}),
                  flush=True)
        R["tags"][tag] = T
        print(tag, {k2: T[k2] for k2 in T if k2 != "witnesses"}, round(time.perf_counter() - t0), flush=True)
        del GU
    R["t_total"] = time.perf_counter() - t0
    out_path.write_text(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
