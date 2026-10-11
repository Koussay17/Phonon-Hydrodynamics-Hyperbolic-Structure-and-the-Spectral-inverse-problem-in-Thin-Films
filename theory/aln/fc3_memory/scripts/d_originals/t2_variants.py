"""Task 2 analysis: kappa, tau_mem, tau_CS, lambda_min for every fc3 variant (phono3py-route pp, orbit-averaged rates).

usage: python -B t2_variants.py <events.npz> <resultdir> <label> [families ...]
Writes <resultdir>/t2_variants_<label>.json (per variant) and .csv.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, summary  # noqa: E402


def main():
    ev_path, rdir, label = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
    fams = sys.argv[4:] or ["cut", "iso", "mul", "scale"]
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    pp_ref = np.load(rdir / "pp_r01_REF.npz")["pp"]
    g_ref = M.rates(pp_ref)
    mref = M.metrics(g_ref)
    allowed = g_ref > 1e-9 * g_ref.max()
    rows = []
    out = {"reference": summary(mref), "variants": {}}
    for fam in fams:
        fn = rdir / f"variants_{fam}.npz"
        if not fn.exists():
            continue
        V = np.load(fn)
        for name, pp, eps in zip(V["names"], V["pp"], V["eps"]):
            g = M.rates(pp)
            m = M.metrics(g)
            dr = m["r"] / mref["r"] - 1
            ratio = g[allowed] / g_ref[allowed]
            rec = summary(m)
            rec.update({"family": fam, "eps": float(eps),
                        "tau_x_ratio": m["tau_x"] / mref["tau_x"], "tau_z_ratio": m["tau_z"] / mref["tau_z"],
                        "kappa_x_ratio": m["kappa_x"] / mref["kappa_x"], "kappa_z_ratio": m["kappa_z"] / mref["kappa_z"],
                        "lifetime_relchange_maxabs": float(np.abs(dr).max()), "lifetime_relchange_rms": float(np.sqrt(np.mean(dr ** 2))),
                        "event_rate_ratio_q(0.001,0.5,0.999)": np.quantile(ratio, [0.001, 0.5, 0.999]).tolist(),
                        "event_rate_ratio_min": float(ratio.min()), "event_rate_ratio_max": float(ratio.max())})
            out["variants"][str(name)] = rec
            rows.append({"name": str(name), **{k: v for k, v in rec.items() if not isinstance(v, list)}})
            print(f"{name:18s} eps={eps:.3f} kx={m['kappa_x']:.2f} kz={m['kappa_z']:.2f} tx={m['tau_x']:.3f} "
                  f"tz={m['tau_z']:.3f} tCSx={m['tauCS_x']:.3f} lmin={m['lam_min_H']:.4e} dr_max={np.abs(dr).max():.3f}",
                  flush=True)
    # aggregates for random families
    agg = {}
    for name, rec in out["variants"].items():
        if rec["family"] in ("iso", "mul"):
            key = f"{rec['family']}_{rec['eps']:.2f}"
            agg.setdefault(key, []).append(rec)
    out["aggregates"] = {}
    for key, recs in agg.items():
        d = {}
        for q in ("tau_x_ratio", "tau_z_ratio", "kappa_x_ratio", "kappa_z_ratio", "lifetime_relchange_maxabs",
                  "lifetime_relchange_rms", "lam_min_H"):
            v = np.array([r[q] for r in recs])
            d[q] = {"mean": float(v.mean()), "std": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                    "min": float(v.min()), "max": float(v.max())}
        d["n"] = len(recs)
        out["aggregates"][key] = d
    (rdir / f"t2_variants_{label}.json").write_text(json.dumps(out, indent=1))
    with open(rdir / f"t2_variants_{label}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(out["aggregates"], indent=1))


if __name__ == "__main__":
    main()
