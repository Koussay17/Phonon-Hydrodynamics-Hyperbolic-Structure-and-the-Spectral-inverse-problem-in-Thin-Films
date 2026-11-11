"""Summary of the three-mesh film optima (E33) and their held-out 12x12x8 checks (E37).

Reads results_c/film_opt/m3_<tag>_z_L<L>_<sense>.json and results_c/film3/<tag>_m3L<L><sense>_m*.json
(+ _transfer.json); writes results_c/film3/e33_e37_summary.json and prints a Markdown table.
Run from the campaign root: python scripts_c/summarize_e33_e37.py
"""
import json
import os

RUNS = [("orb0.1", 30), ("orb0.1", 100), ("orb0.1", 300),
        ("cal1nz", 100), ("cal3wnz", 100), ("cal3w", 100)]
MESHES = ["m775", "m996", "m12128"]


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


out = {}
for tag, L in RUNS:
    for sense in ("max", "min"):
        opt = load(f"results_c/film_opt/m3_{tag}_z_L{L}_{sense}.json")
        if opt is None:
            continue
        b = opt["best"]
        key = f"{tag}_m3L{L}{sense}"
        rec = {"tag": tag, "L_nm": L, "sense": sense, "outer": b["outer"], "t_h": b["t"] / 3600,
               "viol": b["viol"], "feasible": b["feasible"], "robust_ratio": b["robust_ratio"],
               "fitted_film_z": {m: v["film_tau_ratio"][1] for m, v in b["mesh"].items()},
               "check": {}}
        for m in MESHES:
            ev = load(f"results_c/film3/{key}_{m}.json")
            tr = load(f"results_c/film3/{key}_{m}_transfer.json")
            if ev is None:
                continue
            ref, w = ev["REF"], ev[key]
            c = {"film_z": {l: w[l]["tau_z"] / ref[l]["tau_z"] for l in ("30", "100", "300")},
                 "bulk_z": w["bulk"]["tau_z"] / ref["bulk"]["tau_z"]}
            if tr is not None:
                t = tr[key]
                # the quantile fields were added to mesh_transfer_check.py after the first checks (E37)
                c.update(kappa_ratio=t["kappa_ratio"], lifetime_minmax=t["lifetime_degavg_minmax"],
                         frac_outside_10pct=t.get("lifetime_frac_outside_10pct", float("nan")),
                         absdev_q50_q90_q99=t.get("lifetime_absdev_q50_q90_q99", [float("nan")] * 3))
            rec["check"][m] = c
        out[key] = rec

json.dump(out, open("results_c/film3/e33_e37_summary.json", "w"), indent=1)

print("| run | outer / h | viol | robust ratio | film z, fitted 9x9x6 | film z, held-out 12x12x8 | kept | "
      "kappa 12x12x8 (xx, zz) | lifetimes min/max | outside ±10 % | q99 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for key, r in out.items():
    h = r["check"].get("m12128")
    L = str(r["L_nm"])
    f9 = r["fitted_film_z"].get("m996")
    if h is None:
        print(f"| {key} | {r['outer']} / {r['t_h']:.1f} | {r['viol']:.1e} | {r['robust_ratio']:.4f} | {f9:.4f} | pending | | | | | |")
        continue
    fz = h["film_z"][L]
    kept = (fz - 1) / (r["robust_ratio"] - 1)
    k = h.get("kappa_ratio", [float("nan")] * 2)
    lm = h.get("lifetime_minmax", [float("nan")] * 2)
    print(f"| {key} | {r['outer']} / {r['t_h']:.1f} | {r['viol']:.1e} | {r['robust_ratio']:.4f} | {f9:.4f} | "
          f"{fz:.4f} | {100 * kept:.0f} % | {k[0]:.4f} / {k[1]:.4f} | {100 * lm[0]:+.1f} % / {100 * lm[1]:+.1f} % | "
          f"{100 * h.get('frac_outside_10pct', float('nan')):.2f} % | {100 * h.get('absdev_q50_q90_q99', [0, 0, float('nan')])[2]:.1f} % |")
