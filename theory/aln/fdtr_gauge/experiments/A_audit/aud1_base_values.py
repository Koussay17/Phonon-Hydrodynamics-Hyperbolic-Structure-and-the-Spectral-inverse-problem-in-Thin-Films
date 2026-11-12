"""Audit of A1 (base values of the AlN layer): independent re-assembly from the stored phono3py RTA chunks and
from the REF blocks of the event-operator result files. Written without reference to a1_base_values.py code.

Definitions (fc3-memory E26/E27/E30):
  tau_alpha = 1/(4 pi Gamma_alpha) [ps], Gamma the phono3py linewidth in THz;
  boundary rate (heuristic Matthiessen, diffuse grey film, c-oriented): r_alpha = |v_z,alpha| / L;
  film lifetime 1/tau_L = 1/tau + r;
  memory time = first moment of the normalised RTA kernel in the Laplace variable,
      tau_mem(L) = sum c_alpha tau_L^2 / sum c_alpha tau_L,  c_alpha = C_alpha v_z,alpha^2  (= mode_kappa_zz / tau);
  kappa_z(L) = sum c tau_L / N_q (weights of the irreducible points are already inside mode_kappa).
Checks: the archived assembled.json of each mesh, the REF blocks of all *_m12128.json files, the 500 nm ratio
method and its alternatives, the admissible film tau_z ratios (E33/E37 tables) from the result files.
Usage: python aud1_base_values.py -> out/aud1_base_values.json
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
RTA = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/rta")
RTA_REPO = Path("C:/Users/Koussay/these/theory/aln/fc3_memory/results_c/rta")
FILM3 = Path("C:/Users/Koussay/these/theory/aln/fc3_memory/results_c/film3")
A1 = Path("D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/experiments/A/out/a1_base_values.json")
LENGTHS_NM = (30.0, 100.0, 300.0, 500.0, 1000.0, 3000.0, 10000.0)


def load_mesh(mdir: Path):
    meta = json.loads((mdir / "meta.json").read_text())
    files = sorted(mdir.glob("chunk_*.npz"))
    files = [f for f in files if f.name.count(".") == 1]
    assert len(files) == meta["n_chunks"], (mdir, len(files))
    cols = {k: [] for k in ("w", "f", "gamma", "gv", "mk", "q")}
    for f in files:
        with np.load(f) as z:
            for k in cols:
                cols[k].append(z[k])
    D = {k: np.concatenate(v) for k, v in cols.items()}
    assert int(D["w"].sum()) == meta["sum_weights"] == int(np.prod(meta["mesh"]))
    return meta, D


def film_values(meta, D):
    nq = float(meta["sum_weights"])
    good = (D["gamma"] > 0.0) & (D["f"] > 1e-3)
    tau = np.zeros_like(D["gamma"])
    tau[good] = 1.0 / (4.0 * np.pi * D["gamma"][good])                       # ps
    vz = np.abs(D["gv"][:, :, 2]) * 1e2                                      # THz*Angstrom -> m/s
    out = {"mesh": meta["mesh"], "n_ir": int(D["w"].size), "n_good_modes_ir": int(good.sum())}
    for comp, j in (("x", 0), ("z", 2)):
        mk = np.where(good, D["mk"][:, :, j], 0.0)
        c = np.zeros_like(mk)
        c[good] = mk[good] / tau[good]                                       # proportional to C v_j^2 (with weights)
        out[f"kappa_{comp}_bulk"] = float(mk.sum() / nq)
        out[f"tau_{comp}_bulk"] = float((c * tau ** 2).sum() / (c * tau).sum())
        tl_all = {}
        for L in LENGTHS_NM:
            rate = vz / (L * 1e-9) * 1e-12                                     # 1/ps
            tl = np.zeros_like(tau)
            tl[good] = tau[good] / (1.0 + tau[good] * rate[good])
            tl_all[L] = tl
            out.setdefault(f"tau_{comp}", {})[f"{L:g}"] = float((c * tl ** 2).sum() / (c * tl).sum())
            out.setdefault(f"kappa_{comp}", {})[f"{L:g}"] = float((c * tl).sum() / nq)
    # the q_z = 0 plane (v_z = 0 by C2z x time reversal)
    qz = np.abs(D["q"][:, 2]) % 1.0
    plane = np.minimum(qz, 1.0 - qz) < 1e-9
    out["qz0_weight_fraction"] = float(D["w"][plane].sum() / nq)
    out["qz0_max_abs_vz_m_s"] = float(vz[plane].max())
    return out


def main():
    res = {"rta": {}}
    for m in ("m484826", "m565630", "m646435"):
        meta, D = load_mesh(RTA / m)
        r = film_values(meta, D)
        arch = json.loads((RTA / m / "assembled.json").read_text())
        arch_repo = json.loads((RTA_REPO / m / "assembled.json").read_text()) if (RTA_REPO / m / "assembled.json").exists() else None
        devs = []
        for comp in ("x", "z"):
            for L, v in arch[f"film_vz_{comp}"].items():
                devs.append(abs(r[f"tau_{comp}"][f"{float(L):g}"] / v - 1.0))
        r["max_rel_dev_vs_archived_film_tau"] = float(max(devs))
        r["rel_dev_kappa_bulk_vs_archived"] = float(max(abs(r["kappa_x_bulk"] / arch["kappa_x"] - 1),
                                                       abs(r["kappa_z_bulk"] / arch["kappa_z"] - 1)))
        r["repo_copy_identical"] = (arch_repo == arch) if arch_repo is not None else None
        res["rta"][m] = r
        print(m, "max dev vs archive", r["max_rel_dev_vs_archived_film_tau"], "kappa dev", r["rel_dev_kappa_bulk_vs_archived"],
              "tau_z(500) %.4f kappa_z(500) %.3f kappa_x(500) %.3f" % (r["tau_z"]["500"], r["kappa_z"]["500"],
                                                                      r["kappa_x"]["500"]), flush=True)
    # mesh sequence N = 48, 56, 64 (in-plane); observed convergence of tau_z(500) and a Richardson-type estimate
    seq = [res["rta"][m]["tau_z"]["500"] for m in ("m484826", "m565630", "m646435")]
    nx = np.array([48.0, 56.0, 64.0])
    d1, d2 = seq[1] - seq[0], seq[2] - seq[1]
    res["rta_tau_z500_sequence"] = seq
    res["rta_tau_z500_increments"] = [d1, d2]
    # fit tau(N) = tau_inf + a N^-p through three points
    from scipy.optimize import brentq
    def g(p):
        return (seq[1] - seq[0]) / (nx[1] ** -p - nx[0] ** -p) - (seq[2] - seq[1]) / (nx[2] ** -p - nx[1] ** -p)
    try:
        p = brentq(g, 0.05, 8.0)
        a = (seq[2] - seq[1]) / (nx[2] ** -p - nx[1] ** -p)
        res["rta_tau_z500_fit"] = {"p": p, "tau_inf": seq[2] - a * nx[2] ** -p}
    except ValueError:
        res["rta_tau_z500_fit"] = "no power-law fit through three points (non-monotone ratio)"
    for p in (1.0, 2.0):
        a = (seq[2] - seq[1]) / (nx[2] ** -p - nx[1] ** -p)
        res[f"rta_tau_z500_extrap_p{p:g}"] = seq[2] - a * nx[2] ** -p

    # event operator REF blocks
    files = sorted(glob.glob(str(FILM3 / "*_m12128.json")))
    blocks = [json.loads(Path(f).read_text())["REF"] for f in files]
    spread = 0.0
    for b in blocks[1:]:
        for L in blocks[0]:
            for q in blocks[0][L]:
                spread = max(spread, abs(b[L][q] / blocks[0][L][q] - 1.0))
    ref = blocks[0]
    res["event_REF_files"] = [Path(f).name for f in files]
    res["event_REF_spread"] = spread
    res["event_REF"] = ref
    conv = res["rta"]["m646435"]

    def interp_loglog(L0, y0, L1, y1, L):
        t = np.log(L / L0) / np.log(L1 / L0)
        return float(np.exp((1 - t) * np.log(y0) + t * np.log(y1)))

    est = {}
    for q, qr in (("tau_z", "tau_z"), ("kappa_z", "kappa_z"), ("kappa_x", "kappa_x")):
        r300 = ref["300"][q] / conv[qr]["300"]
        r1000 = ref["1000"][q] / conv[qr]["1000"]
        r500 = interp_loglog(300, r300, 1000, r1000, 500)
        est[q] = {"ratio_event12_over_rta64_300": r300, "ratio_event12_over_rta64_1000": r1000,
                  "ratio_interp_500": r500, "ratio_method_500": r500 * conv[qr]["500"],
                  "event_loglog_300_1000": interp_loglog(300, ref["300"][q], 1000, ref["1000"][q], 500),
                  "rta64_direct_500": conv[qr]["500"]}
    # alternative: event/RTA ratio on the SAME 12x12x8 mesh (REF block carries tauRTA_z on the same diagonal), applied to the
    # converged RTA value: an estimate of the event value corrected for the 12x12x8 mesh bias (tau only)
    s300 = ref["300"]["tau_z"] / ref["300"]["tauRTA_z"]
    s1000 = ref["1000"]["tau_z"] / ref["1000"]["tauRTA_z"]
    s500 = interp_loglog(300, s300, 1000, s1000, 500)
    est["tau_z"]["same_mesh_ratio_300"] = s300
    est["tau_z"]["same_mesh_ratio_1000"] = s1000
    est["tau_z"]["same_mesh_ratio_500"] = s500
    est["tau_z"]["mesh_corrected_event_500"] = s500 * conv["tau_z"]["500"]
    # mesh bias of RTA on 12x12x8 against 64x64x35 at the tabulated lengths
    est["tau_z"]["rta12_over_rta64"] = {L: ref[L]["tauRTA_z"] / conv["tau_z"][f"{float(L):g}"]
                                        for L in ("30", "100", "300", "1000", "3000", "10000")}
    res["L500"] = est

    # compare with the branch values
    a1 = json.loads(A1.read_text())["base"]
    cmp_ = {}
    for L in ("30", "100", "300"):
        cmp_[L] = {"tau_R_ps": (a1[L]["tau_R_ps"], ref[L]["tau_z"]), "lam_z": (a1[L]["lam_z"], ref[L]["kappa_z"]),
                   "lam_r_rta": (a1[L + "_rta"]["lam_r"], conv["kappa_x"][f"{float(L):g}"])}
    cmp_["500"] = {"tau_R_ps": (a1["500"]["tau_R_ps"], est["tau_z"]["ratio_method_500"]),
                   "lam_z": (a1["500"]["lam_z"], est["kappa_z"]["ratio_method_500"]),
                   "lam_r_rta": (a1["500_rta"]["lam_r"], conv["kappa_x"]["500"]),
                   "tau_R_rta_ps": (a1["500_rta"]["tau_R_ps"], conv["tau_z"]["500"]),
                   "lam_z_rta": (a1["500_rta"]["lam_z"], conv["kappa_z"]["500"])}
    maxdev = 0.0
    for L, d in cmp_.items():
        for k, (x, y) in d.items():
            maxdev = max(maxdev, abs(x / y - 1.0))
    res["branch_vs_audit"] = cmp_
    res["branch_vs_audit_max_rel_dev"] = maxdev

    # admissible film tau_z ratios from the result files (orb 0.1, three data meshes + held-out 12x12x8)
    adm = {}
    for f in sorted(glob.glob(str(FILM3 / "orb0.1_m3L*_m12128.json"))) + sorted(glob.glob(str(FILM3 / "cal*_m3L100*_m12128.json"))):
        d = json.loads(Path(f).read_text())
        name = [k for k in d if k.endswith("_ratio")][0]
        Lsel = name.split("m3L")[1].split("m")[0].rstrip("ax").rstrip("i")
        adm[name] = {L: d[name][L]["tau_z"] for L in ("30", "100", "300", "1000")}
    res["heldout_12x12x8_tau_z_ratios"] = adm
    # the three-mesh (5x5x3, 7x7x5, 9x9x6) values are in the companion files of the other meshes
    tri = {}
    for f in sorted(glob.glob(str(FILM3 / "orb0.1_m3L*_m*.json"))):
        if "transfer" in f:
            continue
        d = json.loads(Path(f).read_text())
        name = [k for k in d if k.endswith("_ratio")]
        if not name:
            continue
        tri[Path(f).stem] = {L: d[name[0]][L]["tau_z"] for L in ("30", "100", "300")}
    res["orb01_ratio_files"] = tri
    (OUT / "aud1_base_values.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(est, indent=1))
    print("branch vs audit max rel dev", maxdev)
    print(json.dumps(adm, indent=1))
    print(json.dumps(tri, indent=1))


if __name__ == "__main__":
    main()
