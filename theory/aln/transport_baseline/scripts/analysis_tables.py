"""Collect campaign results into comparison tables (results/tables.json + results/tables.md).

Re-runnable: missing runs are reported as absent.  Reads only run.json / kappa-*.hdf5 /
solve_T*.json files written by run_kappa.py and lbte_ooc.py, and
results/olympics_benchmark.json.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

C = common.CAMPAIGN
R = C / "runs"
bench = json.loads((C / "results" / "olympics_benchmark.json").read_text())
md = []
tables = {}


def kappa_file(run):
    fs = sorted((R / run).glob("kappa-m*.hdf5"))
    fs = [f for f in fs if "-g" not in f.name]
    return fs[0] if fs else None


def read_kappa(run, key="kappa"):
    f = kappa_file(run)
    if f is None:
        return None
    with h5py.File(f, "r") as h:
        if key not in h:
            return None
        return np.array(h["temperature"][()]), np.array(h[key][()])


def at_T(run, T, key="kappa"):
    d = read_kappa(run, key)
    if d is None:
        return None
    t, k = d
    i = np.where(np.isclose(t, T))[0]
    return None if len(i) == 0 else (float(k[i[0], 0]), float(k[i[0], 2]))


def ooc(run, k=0, tag=""):
    p = R / run / f"solve_T{k}{tag}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d


def rel(a, b):
    return (a - b) / b


# ---------------------------------------------------------------- reproduction by mesh
olymp_mesh = {tuple(r["mesh"]): r for r in bench["phono3py_mesh_series_300K"]}
rep_rows = []
for mesh, rta_run, lbte_run, ooc_run in [
        ((9, 9, 5), "t1-rta-m995-v2compat-C", "t1-lbte-m995-v2compat-C", "ooc-m995-v2C"),
        ((15, 15, 9), None, "rep-lbte-m15159-v2C", "ooc-m15159-v2C"),
        ((19, 19, 11), None, None, "ooc-m191911-v2C"),
        ((23, 23, 13), None, None, "ooc-m232313-v2C"),
        ((27, 27, 15), None, None, "ooc-m272715-v2C"),
        ((31, 31, 17), "rep-rta-m313117-v2C", None, "ooc-m313117-v2C")]:
    o = olymp_mesh.get(mesh)
    row = {"mesh": list(mesh), "olympics": o}
    if lbte_run and (R / lbte_run / "run.json").exists():
        j = json.loads((R / lbte_run / "run.json").read_text())["result"]
        kk = j["kappa_WmK[sigma][T][xx,yy,zz,yz,xz,xy]"][0][0]
        kr = j["kappa_RTA_WmK[sigma][T][xx,yy,zz,yz,xz,xy]"][0][0]
        row["phono3py_dense"] = {"lbte": [kk[0], kk[2]], "rta": [kr[0], kr[2]]}
    if rta_run and at_T(rta_run, 300.0):
        row["phono3py_rta"] = list(at_T(rta_run, 300.0))
    s = ooc(ooc_run)
    if s:
        row["ooc"] = {"lbte": [s["kappa_xx"], s["kappa_zz"]], "rta": [s["kappa_RTA_xx"], s["kappa_RTA_zz"]],
                      "pcg_iterations": s["pcg_iterations"], "residual": s["true_rel_residual"]}
    rep_rows.append(row)
tables["reproduction_by_mesh_300K"] = rep_rows
md.append("### Reproduction by mesh, 300 K (phono3py-2.1.0-like settings), W/(m K)\n")
md.append("| mesh | Olympics LBTE xx/zz | ours LBTE xx/zz | dLBTE % | Olympics RTA xx/zz | ours RTA xx/zz | dRTA % | route |")
md.append("|---|---|---|---|---|---|---|---|")
for row in rep_rows:
    o = row["olympics"]
    if "phono3py_dense" in row:
        lb, rt, route = row["phono3py_dense"]["lbte"], row["phono3py_dense"]["rta"], "phono3py dense"
    elif "ooc" in row:
        lb, rt, route = row["ooc"]["lbte"], row["ooc"]["rta"], "rows+PCG"
    else:
        lb, rt, route = None, row.get("phono3py_rta"), "phono3py RTA"
    if "phono3py_rta" in row:
        rt = row["phono3py_rta"]
    def f2(v):
        return "-" if v is None else f"{v[0]:.3f} / {v[1]:.3f}"
    def d2(v, ref):
        if v is None or ref is None:
            return "-"
        return f"{100*rel(v[0], ref[0]):+.3f} / {100*rel(v[1], ref[1]):+.3f}"
    ol = None if o is None else [o["lbte_xx"], o["lbte_zz"]]
    orr = None if o is None else [o["rta_xx"], o["rta_zz"]]
    md.append(f"| {'x'.join(map(str,row['mesh']))} | {f2(ol)} | {f2(lb)} | {d2(lb, ol)} | {f2(orr)} | {f2(rt)} | {d2(rt, orr)} | {route} |")

# ---------------------------------------------------------------- RTA kappa(T) reproduction
olyT = {r["T"]: r for r in bench["phono3py_kappa_T_31x31x17"]}
d = read_kappa("rep-rta-m313117-v2C")
if d is not None:
    t, k = d
    rows = []
    for i, T in enumerate(t):
        o = olyT.get(float(T))
        rows.append({"T": float(T), "ours_rta": [float(k[i, 0]), float(k[i, 2])],
                     "olympics_rta": None if o is None else [o["rta_xx"], o["rta_zz"]]})
    tables["reproduction_rta_T_31x31x17"] = rows
    md.append("\n### Reproduction of the phono3py-team RTA kappa(T), 31x31x17, W/(m K)\n")
    md.append("| T (K) | Olympics xx / zz | ours xx / zz | rel. diff xx / zz (%) |")
    md.append("|---|---|---|---|")
    for r in rows:
        o, u = r["olympics_rta"], r["ours_rta"]
        if o is None or o[0] is None:
            continue
        md.append(f"| {r['T']:.0f} | {o[0]:.3f} / {o[1]:.3f} | {u[0]:.3f} / {u[1]:.3f} | "
                  f"{100*rel(u[0], o[0]):+.4f} / {100*rel(u[1], o[1]):+.4f} |")

(C / "results" / "tables.json").write_text(json.dumps(tables, indent=1))
(C / "results" / "tables.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
