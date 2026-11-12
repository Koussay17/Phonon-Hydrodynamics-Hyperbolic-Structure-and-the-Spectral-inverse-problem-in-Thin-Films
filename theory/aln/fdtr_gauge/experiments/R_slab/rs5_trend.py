"""RS5: mesh trend 5x5x3 -> 7x7x5 of the slab/regularised transfer (REF and the ten film witnesses of filmw_pp).

Inputs: rs4_filmw_static.json (regularised first-order delta, both meshes), rs2b_<mesh>_filmw_<fs>.json (slab fits,
black walls; 5x5x3 all witnesses, 7x7x5 the subset run), rs2b_m775_film2_<fs>.json (REF cross-check: the REF tensor is
the same in both amplitude sets). Ranges by the campaign rule on the witness subset available on both meshes.
usage: python -B rs5_trend.py [freqset]       Output: out/rs5_trend.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
FS = sys.argv[1] if len(sys.argv) > 1 else "sub4"
ST = json.loads((OUT / "rs4_filmw_static.json").read_text())
FIT = {m: json.loads((OUT / f"rs2b_{m}_filmw_{FS}.json").read_text()) for m in ("m553", "m775")
       if (OUT / f"rs2b_{m}_filmw_{FS}.json").exists()}
F2 = json.loads((OUT / f"rs2b_m775_film2_{FS}.json").read_text()) if (OUT / f"rs2b_m775_film2_{FS}.json").exists() else None


def hr(dv):
    up, low = max(list(dv) + [0.0]), min(list(dv) + [0.0])
    return {"up": up, "low": low, "half_range_ps": 0.5 * (up - low)}


res = {"freqset": FS}
for L in ("300", "500"):
    rl = {}
    common = None
    for m in FIT:
        have = {nm for nm in FIT[m].get("full", {}) if nm != "REF"} & {nm for nm in FIT[m].get("rta", {}) if nm != "REF"}
        common = have if common is None else common & have
    common = sorted(common or [])
    for m in ("m553", "m775"):
        if m not in FIT:
            continue
        rows = ST[m]["rows"]
        ent = {"REF": {}}
        for op in ("full", "rta"):
            fr = FIT[m][op]
            for N in ("N2lam", "N3"):
                key = f"{N}|0.01"
                ref_fit = fr["REF"][f"fit_{L}|wall|note16"][key]["delta_fit_ps"]
                ref_reg = fr["REF"][f"delta_reg_{L}_ps"]
                ent["REF"][f"{op}|{N}"] = {"delta_fit_ps": ref_fit, "delta_reg_ps": ref_reg, "ratio": ref_fit / ref_reg,
                                           "lamz_fit": fr["REF"][f"fit_{L}|wall|note16"][key]["lamz_fit"]}
                for sel, ws in (("common", common), ("all_available", sorted(nm for nm in fr if nm != "REF"))):
                    if not ws:
                        continue
                    dfit = [fr[w][f"fit_{L}|wall|note16"][key]["delta_fit_ps"] - ref_fit for w in ws]
                    dreg = [fr[w][f"delta_reg_{L}_ps"] - ref_reg for w in ws]
                    a, b = hr(dfit), hr(dreg)
                    ent[f"{op}|{N}|{sel}"] = {"witnesses": ws, "slab": a, "reg": b,
                                              "transfer": a["half_range_ps"] / b["half_range_ps"] if b["half_range_ps"] else None,
                                              "per_witness": {w: {"d_fit": x, "d_reg": y} for w, x, y in zip(ws, dfit, dreg)}}
            ent["REF"][f"{op}|G0"] = fr["REF"]["static"][L]["G0"]
            ent["REF"][f"{op}|lam_eff"] = fr["REF"]["static"][L]["lam_eff"]
        # regularised-only ranges over all ten witnesses (cheap, both meshes)
        exc = {w: rows[w][L]["delta_reg_full_ps"] - rows["REF"][L]["delta_reg_full_ps"] for w in rows if w != "REF"}
        ent["reg_full_all10"] = dict(hr(exc.values()), ref_ps=rows["REF"][L]["delta_reg_full_ps"])
        excr = {w: rows[w][L]["delta_reg_rta_ps"] - rows["REF"][L]["delta_reg_rta_ps"] for w in rows if w != "REF"}
        ent["reg_rta_all10"] = dict(hr(excr.values()), ref_ps=rows["REF"][L]["delta_reg_rta_ps"])
        rl[m] = ent
    if F2 is not None and "m775" in FIT:
        rl["REF_crosscheck_film2_vs_filmw"] = {
            op: abs(F2[op]["REF"][f"fit_{L}|wall|note16"]["N3|0.01"]["delta_fit_ps"]
                    - FIT["m775"][op]["REF"][f"fit_{L}|wall|note16"]["N3|0.01"]["delta_fit_ps"]) for op in ("full", "rta")}
    res[L] = rl
    print(L, json.dumps({m: {"REF": {k: (round(v["delta_fit_ps"], 2), round(v["ratio"], 3)) for k, v in rl[m]["REF"].items()
                                     if isinstance(v, dict)},
                             "transfer_common": {k: (round(v["transfer"], 3) if v["transfer"] else None,
                                                     round(v["slab"]["half_range_ps"], 2), round(v["reg"]["half_range_ps"], 2))
                                                 for k, v in rl[m].items() if k.endswith("common")}}
                         for m in rl if m.startswith("m")}), flush=True)
(OUT / "rs5_trend.json").write_text(json.dumps(res, indent=1))
