"""AU8: witness-specific transfer of the delta range from the regularised film to a kinetic slab (RTA diagonal).

For REF and the twelve witnesses on the event meshes 7x7x5 and 9x9x6:
  * single-mode rates r_i = diag C(g) = sum_alpha g_alpha a_alpha,i^2 with orbit-tied g (orbits = minimum event index
    over the 24 stored grid maps, one pass per map: the maps form a group, so the orbit of alpha is {R alpha});
    tau_i = 1/(4 pi r_i) ps; |v_z|, capacities from the event geometry (aln_geometry conventions re-derived, au3);
  * RTA (V1) regularised film times tau_R,M, tau_l,M, delta_M (|v_z|/L, unbinned modes);
  * kinetic slab with diffuse black walls (au6 model) on mode classes whose membership is fixed by the REF lifetimes and
    velocities (so that witness differences carry no binning noise); the class tau and |v| are recomputed from each
    tensor's own lifetimes; FDTR data of each slab fitted by the GK film of the REF classes (tau_l free; N2lam, N3).
Outputs per (mesh, L): excursions of delta_M (RTA diagonal), of the full-operator delta (ra1_film_*.json) and of the
slab-visible delta_fit; witness-specific transfer T_W = Delta delta_fit / Delta delta_M; class half-ranges.
Output: out/au8_witness_slab.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au5_fisher1d as a5  # noqa: E402
import au6_slab as a6  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
RA = HERE.parent / "R_A" / "out"
MESH = {"m775": ("results_c/m775x/events_m775_s0.1.npz", "results_c/m775x/film2_pp"),
        "m996": ("results_c/m996/events_m996_s0.1.npz", "results_c/m996/film2_pp")}
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
KB = 1.380649e-23
NB = 24


def geometry(mesh):
    Z = np.load(FC3 / MESH[mesh][0])
    f = Z["freqs"]
    Nq, nbnd = f.shape
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nbnd, dtype=np.int64)
    newidx[valid] = np.arange(valid.sum())
    eps = f.ravel()[valid]
    kT = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    vz = np.abs(Z["gv"].reshape(Nq * nbnd, 3)[valid][:, 2]) * 100.0          # m/s
    pf, af, bf = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    P, A, B = newidx[pf], newidx[af], newidx[bf]
    x = eps / kT
    nocc = 1.0 / np.expm1(x)
    dsc = np.sqrt(nocc * (1 + nocc))
    rep = A == B
    lp = np.sqrt((eps[A] + eps[B]) / eps[P])
    ld = 1.0 / lp
    w2P = (lp / dsc[P]) ** 2
    w2A = (np.where(rep, 2 * ld, ld) / dsc[A]) ** 2
    w2B = np.where(rep, 0.0, (ld / dsc[B]) ** 2)
    s = np.sinh(x / 2)
    cap = KB * x * x * nocc * (1 + nocc) / (Nq * float(Z["volume"]) * 1e-30)       # J/(m^3 K) per mode
    # orbit labels: minimum event index over the images under the stored maps (group action, one pass)
    NF = Nq * nbnd
    key = (pf * NF + np.minimum(af, bf)) * NF + np.maximum(af, bf)
    order = np.argsort(key)
    skey = key[order]
    m = key.size
    lab = np.arange(m)
    missing = 0
    for R in Z["rot_maps"]:
        def img(full):
            return R[full // nbnd] * nbnd + full % nbnd
        ip, ia, ib = img(pf), img(af), img(bf)
        k2 = (ip * NF + np.minimum(ia, ib)) * NF + np.maximum(ia, ib)
        pos = np.clip(np.searchsorted(skey, k2), 0, m - 1)
        found = skey[pos] == k2
        missing += int((~found).sum())
        cand = np.where(found, order[pos], lab)
        lab = np.minimum(lab, cand)
    _, orb = np.unique(lab, return_inverse=True)
    return {"eps": eps, "vz": vz, "cap": cap, "P": P, "A": A, "B": B, "rep": rep, "w2P": w2P, "w2A": w2A, "w2B": w2B,
            "s": s, "orb": orb, "n_orb": int(orb.max() + 1), "missing": missing, "n": int(eps.size), "Nq": Nq}


def diag_rates(G, pp, conv, g0w, tie=True):
    P, A, B, rep, s = G["P"], G["A"], G["B"], G["rep"], G["s"]
    W = np.where(rep, conv * g0w / (8 * s[P] * s[A] * s[A]), conv * g0w / (4 * s[P] * s[A] * s[np.where(rep, A, B)]))
    g = W * pp
    if tie:
        g = (np.bincount(G["orb"], weights=g, minlength=G["n_orb"]) / np.bincount(G["orb"], minlength=G["n_orb"]))[G["orb"]]
    n = G["n"]
    r = np.bincount(P, weights=g * G["w2P"], minlength=n) + np.bincount(A, weights=g * G["w2A"], minlength=n) \
        + np.bincount(B, weights=g * G["w2B"], minlength=n)
    return r


def rta_film(cap, vz, tau, L):
    te = 1.0 / (1.0 / tau + vz / (L * 1e-9))
    S = cap.sum()
    lam = (cap * vz * vz * te).sum()
    s = lam / S
    tauR = (cap * vz * vz * te ** 2).sum() / lam
    taul = (cap * te * (vz * vz * te - s) ** 2).sum() / (S * s * s)
    return {"lam": lam, "tauR": tauR, "taul": taul, "delta": tauR - taul}


def fixed_classes(vz, tau_ref, cap, nb, vfloor=1.0):
    mov = vz > vfloor
    lt = np.log10(tau_ref)
    lv = np.log10(np.maximum(vz, vfloor))
    et = np.linspace(lt.min() - 1e-9, lt.max() + 1e-9, nb + 1)
    ev = np.linspace(lv[mov].min() - 1e-9, lv[mov].max() + 1e-9, nb + 1)
    it = np.clip(np.searchsorted(et, lt) - 1, 0, nb - 1)
    iv = np.clip(np.searchsorted(ev, lv) - 1, 0, nb - 1)
    kb = np.where(mov, it * nb + iv, -1)
    k0 = np.where(~mov, it, -1)
    return kb, k0


def classes_from(kb, k0, cap, vz, tau, nb):
    mov = kb >= 0
    C = np.bincount(kb[mov], weights=cap[mov], minlength=nb * nb)
    Ct = np.bincount(kb[mov], weights=(cap * tau)[mov], minlength=nb * nb)
    Cv = np.bincount(kb[mov], weights=(cap * vz * vz * tau)[mov], minlength=nb * nb)
    ok = C > 0
    loc = k0 >= 0
    C0 = np.bincount(k0[loc], weights=cap[loc], minlength=nb)
    Ct0 = np.bincount(k0[loc], weights=(cap * tau)[loc], minlength=nb)
    ok0 = C0 > 0
    return {"c": C[ok] / 2, "tau": Ct[ok] / C[ok], "v": np.sqrt(Cv[ok] / Ct[ok]), "c0": C0[ok0], "tau0": Ct0[ok0] / C0[ok0]}


def main():
    t0 = time.perf_counter()
    res = {}
    p = 2j * np.pi * a5.FREQ
    for mesh in ("m775", "m996"):
        G = geometry(mesh)
        ref = np.load(FC3 / MESH[mesh][1] / "pp_r01_REF.npz")
        conv, g0w = float(ref["conv"]), ref["g0w"]
        ra1 = json.loads((RA / f"ra1_film_{mesh}.json").read_text())
        rc = json.loads((RA / "ra1c_ref_rta.json").read_text())[mesh]
        taus = {}
        for nm in ["REF"] + WIT:
            pp = ref["pp"] if nm == "REF" else np.load(FC3 / MESH[mesh][1] / f"pp_r01_{nm}.npz")["pp"]
            r = diag_rates(G, pp, conv, g0w)
            taus[nm] = 1.0 / (4 * np.pi * r) * 1e-12                         # s
        r_untied = diag_rates(G, ref["pp"], conv, g0w, tie=False)
        out = {"n_modes": G["n"], "n_orbits": G["n_orb"], "missing_images": G["missing"],
               "REF_tie_effect_max_rel": float(np.max(np.abs(taus["REF"] * 4 * np.pi * 1e12 * r_untied - 1)))}
        print(mesh, "modes", G["n"], "orbits", G["n_orb"], "missing", G["missing"], "REF tie effect",
              out["REF_tie_effect_max_rel"], f"{time.perf_counter() - t0:.0f}s", flush=True)
        for L in (300, 500):
            d = L * 1e-9
            rows = {}
            # check against the campaign's RTA kinetic times on the same diagonal (REF, ra1c)
            fm = rta_film(G["cap"], G["vz"], taus["REF"], L)
            out[f"{L}|REF_check"] = {"tauR_ps": fm["tauR"] * 1e12, "taul_ps": fm["taul"] * 1e12,
                                     "campaign_tauR_ps": rc[f"{float(L):.6g}"]["tau_RTA"],
                                     "campaign_taul_ps": rc[f"{float(L):.6g}"]["taul_RTA"]}
            kb, k0 = fixed_classes(G["vz"], taus["REF"], G["cap"], NB)
            cl_ref = classes_from(kb, k0, G["cap"], G["vz"], taus["REF"], NB)
            Mref = a6.matthiessen(cl_ref, d)
            for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8)}.items():
                base = {"tauR": Mref["tauR"], "taul": Mref["taul"], "lamz": Mref["lam"], "C": Mref["C"], "R1": R1,
                        "R2": R2, "d": d}
                for nm in ["REF"] + WIT:
                    fmw = rta_film(G["cap"], G["vz"], taus[nm], L)
                    clw = classes_from(kb, k0, G["cap"], G["vz"], taus[nm], NB)
                    Mw = a6.matthiessen(clw, d)
                    Ys, Ths = [], []
                    for pp_ in p:
                        Y, Th, _ = a6.slab_twoport(clw, d, pp_)
                        Ys.append(Y), Ths.append(Th)
                    y = np.log(a6.Z_stack_slab(Ys, Ths, p, R1, R2, port="wall"))
                    ent = {"delta_M_unbinned_ps": fmw["delta"] * 1e12, "delta_M_classes_ps": Mw["delta"] * 1e12,
                           "lam_M": Mw["lam"], "delta_full_operator_ps": ra1[nm][f"{float(L):.6g}"]["delta_kin"]}
                    for N, nuis in {"N2lam": a5.SETS["N2"] + [2], "N3": a5.SETS["N3"]}.items():
                        chi, th = a5.gn_refit(y, base, nuis, 0.01, free=(1,))
                        ent[f"delta_fit_{N}_ps"] = (Mref["tauR"] - Mref["taul"] - th[1] * 1e-12) * 1e12
                        ent[f"chi_{N}"] = chi
                    rows[nm] = ent
                # excursions and transfers
                R0 = rows["REF"]
                for nm in WIT:
                    e_ = rows[nm]
                    dM = e_["delta_M_classes_ps"] - R0["delta_M_classes_ps"]
                    e_["exc_delta_M_unbinned"] = e_["delta_M_unbinned_ps"] / R0["delta_M_unbinned_ps"] - 1
                    e_["exc_delta_full_operator"] = e_["delta_full_operator_ps"] / R0["delta_full_operator_ps"] - 1
                    for N in ("N2lam", "N3"):
                        dS = e_[f"delta_fit_{N}_ps"] - R0[f"delta_fit_{N}_ps"]
                        e_[f"transfer_{N}"] = dS / dM if dM != 0 else None
                        e_[f"exc_delta_fit_{N}"] = e_[f"delta_fit_{N}_ps"] / R0[f"delta_fit_{N}_ps"] - 1
                summ = {}
                for q in ("delta_M_classes_ps", "delta_fit_N2lam_ps", "delta_fit_N3_ps", "delta_full_operator_ps",
                          "delta_M_unbinned_ps"):
                    dv = [rows[nm][q] - R0[q] for nm in WIT]
                    summ[q] = {"up": max(dv + [0.0]), "low": min(dv + [0.0]),
                               "half_range_ps": 0.5 * (max(dv + [0.0]) - min(dv + [0.0]))}
                summ["range_transfer_N2lam"] = summ["delta_fit_N2lam_ps"]["half_range_ps"] / summ["delta_M_classes_ps"]["half_range_ps"]
                summ["range_transfer_N3"] = summ["delta_fit_N3_ps"]["half_range_ps"] / summ["delta_M_classes_ps"]["half_range_ps"]
                # correlation of the RTA-diagonal excursions with the full-operator ones
                xa = np.array([rows[nm]["exc_delta_M_unbinned"] for nm in WIT])
                xb = np.array([rows[nm]["exc_delta_full_operator"] for nm in WIT])
                summ["corr_RTAdiag_vs_full_excursions"] = float(np.corrcoef(xa, xb)[0, 1])
                summ["slope_RTAdiag_vs_full"] = float(np.polyfit(xb, xa, 1)[0])
                out[f"{L}|{Rtag}"] = {"rows": rows, "summary": summ}
                print(mesh, L, Rtag, "REF check", {k: round(v, 3) for k, v in out[f"{L}|REF_check"].items()},
                      "| range transfer N2lam %.3f N3 %.3f | corr(RTA-diag, full) %.2f slope %.2f" % (
                          summ["range_transfer_N2lam"], summ["range_transfer_N3"], summ["corr_RTAdiag_vs_full_excursions"],
                          summ["slope_RTAdiag_vs_full"]),
                      "| T_W N3:", {nm[:18]: round(rows[nm]["transfer_N3"], 2) for nm in WIT if rows[nm]["transfer_N3"] is not None},
                      f"{time.perf_counter() - t0:.0f}s", flush=True)
            res[mesh] = out
            (OUT / "au8_witness_slab.json").write_text(json.dumps(res, indent=1))
    res["seconds"] = time.perf_counter() - t0
    (OUT / "au8_witness_slab.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
