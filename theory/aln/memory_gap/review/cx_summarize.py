"""Collect the review's optimisation results (review/cx_runs/*.json) into one table.
Output: review/cx_summary.json and a printed table.
"""
import glob
import json
import os

base = r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review\cx_runs"
rows = []
for f in sorted(glob.glob(os.path.join(base, "*.json"))):
    name = os.path.basename(f)[:-5]
    if name.startswith("events_"):
        continue
    d = json.load(open(f))
    if "job" in d:
        r = d.get("result")
        j = d["job"]
        rows.append({"id": name, "col": j.get("col"), "sign": j.get("sign"), "lo": j.get("lo"), "hi": j.get("hi"),
                     "start": j.get("start"), "forbid_zero": j.get("forbid_zero", False), "temps": j.get("temps", [300.0]),
                     "lifetimes_T": j.get("lifetimes_T"), "extras": j.get("extras"),
                     "tau_ps": None if not r else r["tau_ps"], "tau_x_ps": None if not r else r.get("tau_x_ps"),
                     "tau_z_ps": None if not r else r.get("tau_z_ps"), "c_max": None if not r else r["c_max"],
                     "iters": None if not r else r["iters"], "seconds": d.get("seconds"), "error": d.get("error"),
                     "frac_at_lo": None if not r else r.get("frac_at_lo"), "hist_tail": None if not r else r.get("hist_tail")})
    elif name.startswith("smooth_"):
        rows.append({"id": name, **{k: d[k] for k in ("degree", "n_coef", "F", "col", "sign", "tau_ps", "tau_x_ps", "tau_z_ps",
                                                      "c_max", "u_min", "u_max", "nit", "message")}})
json.dump(rows, open(os.path.join(os.path.dirname(base), "cx_summary.json"), "w"), indent=1)
for r in rows:
    t = r.get("tau_ps")
    print(f"{r['id']:28s} tau={t if t is None else round(t, 3)!s:>14} c_max={r.get('c_max')!s:>10.10} "
          f"iters={r.get('iters', r.get('nit'))!s:>5} err={r.get('error')}")
