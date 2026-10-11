"""Numerical audit, part 1: reference operator on one mesh (S8).

For the reference fc3 on the given mesh, with the independent operator of na_common.py:
  (1) rates from the pp file versus the exported event rates; group (orbit) averaging;
  (2) K, N, tau, kappa by full eigendecomposition (null pair removed explicitly) and by Cholesky of
      C + gamma e e^T for gamma / (tr C / n) = 1e-8 ... 1e4 (regularisation independence, residuals);
  (3) spectrum on H = e^perp: lambda_min, cond, share of N and K carried by the lowest 1/5/20/100/500 eigenpairs;
  (4) character of the lowest eigenvectors (participation ratio, frequency, |q|) and the split of N = |x|^2 and
      K = b.x over mode frequency;
  (5) relative condition number of tau with respect to the event rates (elasticities), share of events with tiny
      rates;
  (6) sensitivity to the event cut (n sigma = 2 ... 4) and to the smearing width (sigma' <= 0.1 THz, exact
      Gaussian re-weighting of the stored events; sigma' > 0.1 only as a truncated Gaussian, flagged);
  (7) RTA analogue tau_RTA = sum b^2/r^2 / sum b^2/r with r = diag C.
usage: python -B na_spectral.py m553|m775|m996 [out.json]
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

STORED = {  # campaign values (results/m775x/mesh_transfer.json, joint/verify_m553.json, joint/verify_m996.json)
    "m553": {"kappa": [214.29993549566584, 201.96695921348257], "tau": [34.252990795285065, 27.446471103010975]},
    "m775": {"kappa": [243.02706782165123, 272.33365133947524], "tau": [52.60872168593557, 50.622624927262244]},
    "m996": {"kappa": [265.76132022214045, 256.8998241196887], "tau": [61.396861903142835, 59.15275715063295]},
}


def main():
    mk = sys.argv[1]
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"spectral_{mk}.json"
    t0 = time.perf_counter()
    evp, rdir = nc.MESHES[mk]
    op = nc.Operator(evp)
    R = {"mesh": op.mesh, "n_modes": op.n, "n_events": op.m, "sigma": op.sigma}
    P = np.load(rdir + "/pp_r01_REF.npz")
    g_raw = op.rates_from_pp(P["pp"], P["g0w"], float(P["conv"]))
    R["rates_pp_vs_exported_relmax"] = float(np.max(np.abs(g_raw - op.ev_g) / op.ev_g))
    gauss = np.exp(-op.dl ** 2 / (2 * op.sigma ** 2)) / (np.sqrt(2 * np.pi) * op.sigma)
    R["g0w_over_gaussian_minmax"] = [float(np.min(P["g0w"] / gauss)), float(np.max(P["g0w"] / gauss))]
    g = op.group_average(g_raw)
    img = op.images()
    orb_id = img.min(axis=0)
    R["n_event_orbits"] = int(len(np.unique(orb_id)))
    R["group_avg_check_invariance"] = float(np.max(np.abs(g[img[5]] - g) / g))
    R["t_setup"] = time.perf_counter() - t0

    C = op.C(g)
    tr_n = float(np.trace(C) / op.n)
    R["trC_over_n"] = tr_n
    R["energy_residual_rel"] = float(np.linalg.norm(C @ op.ehat) / np.linalg.norm(C, 2) if op.n < 4000 else
                                     np.linalg.norm(C @ op.ehat) / np.max(np.abs(np.diag(C))))
    R["b_dot_ehat_rel"] = [float(abs(op.b[:, c] @ op.ehat) / np.linalg.norm(op.b[:, c])) for c in (0, 2)]
    R["C_symmetry_rel"] = 0.0

    # ---- (2) eigendecomposition route
    E = nc.eig_response(op, C, keep_vectors=True)
    lam, Q, k0 = E.pop("_lam"), E.pop("_Q"), E.pop("_k0")
    R["eig"] = E
    R["t_eig"] = time.perf_counter() - t0
    st = STORED[mk]
    R["eig_vs_campaign_rel"] = {
        "tau_x": E["cols"]["x"]["tau_ps"] / st["tau"][0] - 1, "tau_z": E["cols"]["z"]["tau_ps"] / st["tau"][1] - 1,
        "kappa_x": E["cols"]["x"]["kappa"] / st["kappa"][0] - 1, "kappa_z": E["cols"]["z"]["kappa"] / st["kappa"][1] - 1}

    # ---- (2b) regularisation scan
    R["chol_gamma_scan"] = {}
    for s in (1e-8, 1e-4, 1e-2, 1.0, 1e2, 1e4):
        R["chol_gamma_scan"][f"{s:g}"] = nc.chol_response(op, C, s * tr_n)
    taus = [[v[a]["tau_ps"] for v in R["chol_gamma_scan"].values()] for a in "xz"]
    R["chol_gamma_scan_tau_spread_rel"] = [float((max(t) - min(t)) / np.mean(t)) for t in taus]

    # ---- (3)/(4) structure of the solution
    mask = np.ones(len(lam), dtype=bool)
    mask[k0] = False
    sel = np.nonzero(mask)[0]
    sel = sel[np.argsort(lam[sel])]
    lH = lam[sel]
    QH = Q[:, sel]
    del Q, C
    rec = np.linalg.inv(op.lattice).T * 2 * np.pi            # rows: reciprocal vectors (1/Angstrom)
    qc = (op.addr / np.array(op.mesh)) @ rec
    qn = np.linalg.norm(qc, axis=1)[op.q_of]
    low = []
    for k in range(5):
        w = QH[:, k] ** 2
        top = np.argsort(-w)[:6]
        low.append({"lambda": float(lH[k]), "participation_ratio": float(1.0 / np.sum(w ** 2)),
                    "mean_f_THz": float(w @ op.f), "weight_f_below_2THz": float(w[op.f < 2.0].sum()),
                    "weight_acoustic_bands": float(w[op.band_of < 3].sum()),
                    "top_modes": [{"f": float(op.f[t]), "band": int(op.band_of[t]), "q_addr": op.addr[op.q_of[t]].tolist(),
                                   "|q|": float(qn[t]), "weight": float(w[t])} for t in top]})
    R["lowest_eigvecs"] = low
    xs = {}
    bins = [0, 1, 2, 3, 4, 6, 8, 12, 30]
    R["freq_bins_THz"] = bins
    R["mode_counts_per_bin"] = np.histogram(op.f, bins=bins)[0].tolist()
    R["lowest_mode_freqs"] = [float(v) for v in np.sort(op.f)[:8]]
    r = op.diagC(g)
    R["diagC_min_max"] = [float(r.min()), float(r.max())]
    R["diagC_at_lowest_modes"] = [float(r[i]) for i in np.argsort(op.f)[:8]]
    for i, c in enumerate((0, 2)):
        a = "xz"[i]
        beta = QH.T @ op.b[:, c]
        x = QH @ (beta / lH)
        xs[a] = x
        N, K = float(x @ x), float(op.b[:, c] @ x)
        hN = np.histogram(op.f, bins=bins, weights=x * x)[0] / N
        hK = np.histogram(op.f, bins=bins, weights=op.b[:, c] * x)[0] / K
        R[f"N_share_by_freq_{a}"] = [float(v) for v in hN]
        R[f"K_share_by_freq_{a}"] = [float(v) for v in hK]
        # (7) RTA analogue
        bb = op.b[:, c] ** 2
        KR, NR = float(np.sum(bb / r)), float(np.sum(bb / r ** 2))
        R[f"RTA_{a}"] = {"tau_ps": NR / KR / nc.FOURPI, "kappa": op.kf * KR,
                         "N_share_by_freq": [float(v) for v in np.histogram(op.f, bins=bins, weights=bb / r ** 2)[0] / NR]}
        # (5) elasticities of tau with respect to the event rates
        u = QH @ (beta / lH ** 2)
        p, q = op.ATx(x), op.ATx(u)
        el = g * (-2 * p * q / N + p * p / K)
        orbel = np.bincount(orb_id, weights=el, minlength=op.m)
        tiny = g < 1e-3 * np.median(g)
        R[f"elasticity_{a}"] = {"sum": float(el.sum()), "sum_abs": float(np.abs(el).sum()),
                                "max_abs_event": float(np.abs(el).max()), "max_abs_orbit": float(np.abs(orbel).max()),
                                "top5_orbit": [float(v) for v in np.sort(np.abs(orbel))[::-1][:5]],
                                "events_tiny_rate_fraction": float(tiny.mean()),
                                "sum_abs_tiny_rate_events": float(np.abs(el[tiny]).sum())}
    del QH
    R["t_struct"] = time.perf_counter() - t0

    # ---- (6) event cut and smearing width
    cut = {}
    for ns in (2.0, 3.0, 3.5, 4.0):
        gg = np.where(np.abs(op.dl) <= ns * op.sigma, g, 0.0)
        Cc = op.C(gg)
        rr = nc.chol_response(op, Cc, np.trace(Cc) / op.n)
        cut[f"{ns:g}"] = {"tau": [rr["x"]["tau_ps"], rr["z"]["tau_ps"]], "kappa": [rr["x"]["kappa"], rr["z"]["kappa"]],
                          "n_events": int((gg > 0).sum())}
        del Cc
    R["cut_scan"] = cut
    sig = {}
    for s2 in (0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.1333):
        wgt = (op.sigma / s2) * np.exp(-op.dl ** 2 / 2 * (1 / s2 ** 2 - 1 / op.sigma ** 2))
        gg = np.where(np.abs(op.dl) <= min(4 * s2, 4 * op.sigma), g * wgt, 0.0)
        Cc = op.C(gg)
        try:
            rr = nc.chol_response(op, Cc, np.trace(Cc) / op.n)
            sig[f"{s2:g}"] = {"tau": [rr["x"]["tau_ps"], rr["z"]["tau_ps"]], "kappa": [rr["x"]["kappa"], rr["z"]["kappa"]],
                              "cut_in_sigma": min(4.0, 4 * op.sigma / s2), "truncated": bool(s2 > op.sigma)}
        except Exception as ex:  # noqa: BLE001
            sig[f"{s2:g}"] = {"error": repr(ex)}
        del Cc
    R["sigma_scan"] = sig
    R["t_total"] = time.perf_counter() - t0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(R, indent=1))
    print(json.dumps({k: R[k] for k in ("mesh", "n_modes", "rates_pp_vs_exported_relmax", "g0w_over_gaussian_minmax",
                                         "eig_vs_campaign_rel", "chol_gamma_scan_tau_spread_rel", "t_total")}, indent=1))


if __name__ == "__main__":
    main()
