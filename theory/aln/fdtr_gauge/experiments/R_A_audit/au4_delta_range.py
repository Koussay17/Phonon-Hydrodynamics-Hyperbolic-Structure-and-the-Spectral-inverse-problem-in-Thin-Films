"""AU4: audit of the static-data range of delta = tau_R - tau_l (R-A6 section 4), from the stored per-witness values.

Recomputed from ra1_film_m775.json / ra1_film_m996.json (not from ra5_delta.json):
  * robust excursion per witness (campaign rule: smaller magnitude of the 7x7x5 and 9x9x6 excursions if they have the
    same sign, else 0), class ends and half-width w = (up + low)/2 for tau_mem, tau_l, delta, kappa_z;
  * the absolute half-range of delta in units of tau_mem,0 = w(delta) |delta_0| / tau_0 (9x9x6 REF), its ratio to
    w(tau_mem), compared with the R-A6 table;
  * alternative readings: (a) 9x9x6 only; (b) 7x7x5 only; (c) a shrink-extrapolation of the excursions to finer
    meshes using the observed 7x7x5 -> 9x9x6 ratio (geometric, i.e. the next refinement shrinks by the same
    factor), as a pessimistic estimate for the held-out mesh where no tau_l exists; (d) the absolute half-range with
    the delta_0 of the FDTR base point (A1 tau_R, mesh-corrected tau_l) instead of the 9x9x6 REF ratio;
  * mesh trend: excursion ratio 9x9x6/7x7x5 per witness for tau_mem and tau_l;
  * film kappa_z excursions of the witnesses (the FDTR side holds lam_z at the base for N1/N2);
  * normal/umklapp labels: which quantities depend on them.
Output: out/au4_delta_range.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
RA = HERE.parent / "R_A" / "out"
A1 = HERE.parent / "A" / "out" / "a1_base_values.json"
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
CLASSES = {"orb0.1": WIT[:6], "all": WIT}
LS = (30, 100, 300, 500, 1000)
key = lambda L: f"{float(L):.6g}"  # noqa: E731
ra = {m: json.loads((RA / f"ra1_film_{m}.json").read_text()) for m in ("m775", "m996")}
r5 = json.loads((RA / "ra5_delta.json").read_text())
a1 = json.loads(A1.read_text())["base"]


def val(m, W, L, q):
    rec = ra[m][W][key(L)]
    return {"tau": rec["tau_z"], "taul": rec["taul_kin"], "delta": rec["delta_kin"], "kappa": rec["kappa_z"],
            "taul_check": rec["tau_z"] - rec["delta_kin"]}[q]


def excursion(m, W, L, q):
    return val(m, W, L, q) / val(m, "REF", L, q) - 1


def class_range(xs):
    up = max([x for x in xs] + [0.0])
    lo = max([-x for x in xs] + [0.0])
    return up, lo, 0.5 * (up + lo)


res = {"checks": {}, "table": {}}
# internal consistency of the stored records: delta_kin = tau_z - taul_kin
res["checks"]["max_abs_delta_identity_ps"] = max(abs(val(m, W, L, "taul_check") - val(m, W, L, "taul"))
                                                 for m in ra for W in ["REF"] + WIT for L in LS)
for L in LS:
    row = {}
    exc = {}
    for W in WIT:
        exc[W] = {}
        for q in ("tau", "taul", "delta", "kappa"):
            x7, x9 = excursion("m775", W, L, q), excursion("m996", W, L, q)
            rob = (x7 if abs(x7) < abs(x9) else x9) if x7 * x9 > 0 else 0.0
            shrink = x9 / x7 if x7 != 0 else np.nan
            # pessimistic held-out estimate: next refinement shrinks by the same factor (only if |shrink| < 1, same sign)
            x12 = x9 * shrink if (0 < shrink < 1) else x9
            exc[W][q] = {"m775": x7, "m996": x9, "robust": rob, "ratio_996_over_775": shrink, "x12_geometric": x12}
    row["excursions"] = exc
    R9, R7 = ra["m996"]["REF"][key(L)], ra["m775"]["REF"][key(L)]
    for cl, ws in CLASSES.items():
        c = {}
        for q in ("tau", "taul", "delta", "kappa"):
            for rule in ("robust", "m996", "m775", "x12_geometric"):
                up, lo, w = class_range([exc[W][q][rule] for W in ws])
                c[f"{q}|{rule}"] = {"up": up, "low": lo, "w": w}
        d0_9 = abs(R9["delta_kin"]) / R9["tau_z"]
        out = {}
        for rule in ("robust", "m996", "m775", "x12_geometric"):
            wd, wt = c[f"delta|{rule}"]["w"], c[f"tau|{rule}"]["w"]
            out[rule] = {"w_tau": wt, "w_taul": c[f"taul|{rule}"]["w"], "w_delta": wd,
                         "halfrange_delta_over_tau0": wd * d0_9, "ratio_to_w_tau": wd * d0_9 / wt if wt > 0 else None}
        # FDTR base point normalisation (A1 tau_R, mesh-corrected tau_l) where available
        if str(L) in a1:
            est = r5["table"][key(L)]["taul0_estimates_ps"]["m996_meshcorr_rta64"]
            tauR0 = a1[str(L)]["tau_R_ps"]
            d0_base = abs(tauR0 - est)
            out["base_point"] = {"tauR0_ps": tauR0, "taul0_ps": est, "delta0_ps": tauR0 - est,
                                 "halfrange_delta_ps_campaign_rule": c["delta|robust"]["w"] * d0_9 * tauR0,
                                 "halfrange_delta_ps_base_delta0": c["delta|robust"]["w"] * d0_base,
                                 "halfrange_tau_ps": c["tau|robust"]["w"] * tauR0}
        out["kappa_w_robust"] = c["kappa|robust"]["w"]
        out["kappa_ends_robust"] = [c["kappa|robust"]["up"], -c["kappa|robust"]["low"]]
        # comparison with the stored ra5 class values
        st = r5["table"][key(L)]["classes"][cl]
        out["campaign"] = {"w_tau": st["tau"]["w_half"], "w_taul": st["taul"]["w_half"], "w_delta": st["delta"]["w_half"],
                           "halfrange_delta_over_tau0": st["delta"]["abs_over_tauR0_m996"],
                           "ratio_to_w_tau": st["rho_kin_over_wtau"]}
        out["agrees_with_campaign"] = bool(abs(out["robust"]["ratio_to_w_tau"] - st["rho_kin_over_wtau"]) < 1e-9)
        row[cl] = out
    # mesh trend summary
    rt = [exc[W]["tau"]["ratio_996_over_775"] for W in WIT if abs(exc[W]["tau"]["m775"]) > 0.005]
    rl = [exc[W]["taul"]["ratio_996_over_775"] for W in WIT if abs(exc[W]["taul"]["m775"]) > 0.005]
    row["mesh_trend"] = {"tau_ratio_996_775_median": float(np.median(rt)) if rt else None,
                         "tau_ratio_range": [float(min(rt)), float(max(rt))] if rt else None,
                         "taul_ratio_996_775_median": float(np.median(rl)) if rl else None,
                         "taul_ratio_range": [float(min(rl)), float(max(rl))] if rl else None,
                         "REF_taul_775_996": [R7["taul_kin"], R9["taul_kin"]],
                         "REF_tau_775_996": [R7["tau_z"], R9["tau_z"]]}
    res["table"][key(L)] = row
    a = row["all"]
    print(f"L={L}: all-class robust w_tau {a['robust']['w_tau']:.4f} w_taul {a['robust']['w_taul']:.4f} "
          f"w_delta {a['robust']['w_delta']:.4f} ratio {a['robust']['ratio_to_w_tau']:.2f} (campaign {a['campaign']['ratio_to_w_tau']:.2f}, "
          f"agree {a['agrees_with_campaign']}) | m996-only {a['m996']['ratio_to_w_tau']:.2f} m775-only {a['m775']['ratio_to_w_tau']:.2f} "
          f"geom-x12 {a['x12_geometric']['ratio_to_w_tau']:.2f} | kappa_z w {a['kappa_w_robust']:.4f} | trend taul "
          f"{row['mesh_trend']['taul_ratio_996_775_median']:.3f} tau {row['mesh_trend']['tau_ratio_996_775_median']:.3f}",
          flush=True)
    if "base_point" in a:
        print("    base point:", {k_: round(v_, 3) for k_, v_ in a["base_point"].items()})
res["labels"] = {"kinetic_tau_l_uses_NU_labels": False,
                 "note": "taul_kin, delta_kin are computed from the full operator (no normal/umklapp split); the labels "
                         "enter only the gray-closure normal times tauN_G1/G2/RQ and their delta_gray."}
(OUT / "au4_delta_range.json").write_text(json.dumps(res, indent=1))
