"""Collect task-3 optimisation results into tables (json + csv + text).

usage: python -B t3_summary.py <resultdir> <subdir1> [<subdir2> ...]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def main():
    rdir = Path(sys.argv[1])
    rows = []
    for sd in sys.argv[2:]:
        for fn in sorted((rdir / sd).glob("*.json")):
            d = json.loads(fn.read_text())
            if "history" not in d:
                continue
            b = d.get("best")
            last = d["history"][-1] if d["history"] else {}
            rows.append({"set": sd, "col": "xyz"[d["col"]], "sense": d["sense"], "rho": d["rho"], "dr": d["dr"],
                         "dK": d["dK"], "seed": d["seed"], "init": d["init_scale"], "tau_base": d["tau_base_ps"],
                         "tau_best": b["tau_ps"] if b else None,
                         "ratio": (b["tau_ps"] / d["tau_base_ps"]) if b else None,
                         "tau_x": b["tau_x_ps"] if b else None, "tau_z": b["tau_z_ps"] if b else None,
                         "rho_used": b["rho"] if b else None, "max_dr": b["max_lifetime_dev"] if b else None,
                         "kx": b["kappa_x_ratio"] if b else None, "kz": b["kappa_z_ratio"] if b else None,
                         "rate_ratio_min": b["orbit_rate_ratio_min"] if b else None,
                         "rate_ratio_max": b["orbit_rate_ratio_max"] if b else None,
                         "n_supp_10x": b["n_orbits_suppressed_10x"] if b else None,
                         "n_supp_100x": b["n_orbits_suppressed_100x"] if b else None,
                         "n_outer": len(d["history"]), "last_viol": last.get("viol"),
                         "seconds": d.get("seconds")})
    rows.sort(key=lambda r: (r["set"], r["col"], r["sense"], r["dr"], r["rho"]))
    (rdir / "t3_summary.json").write_text(json.dumps(rows, indent=1))
    if rows:
        with open(rdir / "t3_summary.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    for r in rows:
        tb = f"{r['tau_best']:.3f}" if r["tau_best"] is not None else "none"
        rt = f"{r['ratio']:.3f}" if r["ratio"] is not None else "-"
        print(f"{r['set']:10s} {r['col']} {r['sense']} rho={r['rho']:<6g} dr={r['dr']:<5g} tau={tb:>12s} "
              f"ratio={rt:>8s} |y|={r['rho_used'] if r['rho_used'] is None else round(r['rho_used'], 4)} "
              f"maxdr={r['max_dr'] if r['max_dr'] is None else round(r['max_dr'], 4)} "
              f"supp10={r['n_supp_10x']} supp100={r['n_supp_100x']} outer={r['n_outer']}")


if __name__ == "__main__":
    main()
