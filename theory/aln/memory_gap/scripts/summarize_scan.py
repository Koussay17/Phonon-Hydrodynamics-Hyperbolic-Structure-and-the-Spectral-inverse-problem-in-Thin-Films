"""Merge scan task outputs: per configuration, best witnesses over seeds and summary table."""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

SCAN = Path(__file__).resolve().parent.parent / "results" / "scan"


def collect(pattern="*.json"):
    groups = defaultdict(list)
    for p in sorted(SCAN.glob(pattern)):
        rec = json.loads(p.read_text())
        tag = p.name.split("__s")[0]
        groups[tag].append((p, rec))
    rows = []
    for tag, recs in groups.items():
        good = [r for _, r in recs if "error" not in r and r.get("max") and r.get("min")]
        if not good:
            rows.append({"tag": tag, "error": recs[0][1].get("error", "no result"), "n_seeds": len(recs)})
            continue
        info = good[0]["info"]
        cfg = good[0]["cfg"]
        tmax = max(r["max"]["tau"] for r in good)
        tmin = min(r["min"]["tau"] for r in good)
        smax = [r["max"]["tau"] for r in good]
        smin = [r["min"]["tau"] for r in good]
        # count of starts within 1 % of the best (search reliability indicator)
        allmax = [l[2] for r in good for l in r["log"] if l[1] == 1 and isinstance(l[2], float)]
        allmin = [l[2] for r in good for l in r["log"] if l[1] == -1 and isinstance(l[2], float)]
        hit_max = sum(1 for v in allmax if v >= 0.99 * tmax)
        hit_min = sum(1 for v in allmin if v <= 1.01 * tmin)
        rows.append({"tag": tag, "d": cfg["d"], "N": cfg["N"], "tie": cfg["tie"], "alpha": cfg.get("alpha"),
                     "n": info["n"], "m": info["m"], "m_var": info["m_var"], "K0": info["K0"],
                     "tau_ref": info["tau_ref"], "tau_RTA": info["tau_RTA"], "CS_lower": info["CS_lower"],
                     "K_RTA_over_K0": info["K_RTA"] / info["K0"], "alpha_fit": info.get("alpha_fit"),
                     "tau_min": tmin, "tau_max": tmax, "ratio": tmax / tmin,
                     "ref_over_min": info["tau_ref"] / tmin, "max_over_ref": tmax / info["tau_ref"],
                     "seed_max": smax, "seed_min": smin, "starts_hit_max": hit_max, "starts_total_max": len(allmax),
                     "starts_hit_min": hit_min, "starts_total_min": len(allmin), "n_seeds": len(good),
                     "time_s": sum(r["time"] for r in good)})
    return rows


def main():
    rows = collect()
    rows.sort(key=lambda r: (r.get("d", 9), str(r.get("alpha")), r.get("tie", ""), r.get("N", 0)))
    out = SCAN.parent / "scan_summary.json"
    out.write_text(json.dumps(rows, indent=1))
    hdr = f"{'tag':28s} {'n':>5s} {'m':>6s} {'tau_RTA':>10s} {'tau_ref':>10s} {'tau_min':>10s} {'tau_max':>11s} {'max/min':>10s} {'hits':>9s}"
    print(hdr)
    for r in rows:
        if "error" in r:
            print(f"{r['tag']:28s} ERROR {r['error']}")
            continue
        print(f"{r['tag']:28s} {r['n']:5d} {r['m']:6d} {r['tau_RTA']:10.4g} {r['tau_ref']:10.4g} {r['tau_min']:10.4g} "
              f"{r['tau_max']:11.4g} {r['ratio']:10.4g} {r['starts_hit_max']:>3d}/{r['starts_total_max']:<3d}")


if __name__ == "__main__":
    main()
