"""Convergence, kappa(T), isotope and cost tables for the production setting.

Reads runs/ooc-*-prod*/solve_T*.json, runs/prod-rta-*/run.json and kappa files,
results/olympics_benchmark.json.  Writes results/convergence.json and appends markdown
to results/tables_production.md.  Missing runs are skipped.
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
out, md = {}, []


def solve(run, k=0, tag=""):
    p = R / run / f"solve_T{k}{tag}.json"
    return json.loads(p.read_text()) if p.exists() else None


def rows_seconds(run):
    p = R / run / "progress.json"
    return json.loads(p.read_text()).get("rows_seconds") if p.exists() else None


# ------------------------------------------------------------ mesh series (tetrahedron)
meshes = [(11, 11, 7), (15, 15, 9), (19, 19, 11), (23, 23, 13), (27, 27, 15), (31, 31, 17), (35, 35, 19)]
series = []
for m in meshes:
    tag = "".join(map(str, m))
    run = f"ooc-m{tag}-prod"
    k = 2 if m == (31, 31, 17) else 0          # 31x31x17 run holds 5 temperatures, 300 K is index 2
    s = solve(run, k)
    if s is None:
        continue
    N = int(np.prod(m))
    series.append({"mesh": list(m), "N": N, "lbte": [s["kappa_xx"], s["kappa_zz"]],
                   "rta": [s["kappa_RTA_xx"], s["kappa_RTA_zz"]], "pcg_iterations": s["pcg_iterations"],
                   "rows_seconds": rows_seconds(run), "solve_seconds": s["seconds"]})
out["mesh_series_prod_300K"] = series


def extrapolate(xs, ys, p):
    """least-squares kappa = k_inf + b N^-p over given points; returns k_inf, b, rms residual"""
    A = np.vstack([np.ones(len(xs)), np.asarray(xs, float) ** (-p)]).T
    coef, *_ = np.linalg.lstsq(A, np.asarray(ys, float), rcond=None)
    res = np.asarray(ys) - A @ coef
    return float(coef[0]), float(coef[1]), float(np.sqrt(np.mean(res**2)))


if len(series) >= 3:
    fin = series[-4:]
    ex = {}
    for key in ("lbte", "rta"):
        for ia, ax in enumerate(("xx", "zz")):
            xs = [r["N"] for r in fin]
            ys = [r[key][ia] for r in fin]
            ex[f"{key}_{ax}"] = {f"p={p}": extrapolate(xs, ys, p) for p in (1.0 / 3, 2.0 / 3, 1.0)}
            ex[f"{key}_{ax}"]["last_two_rel_diff"] = (ys[-1] - ys[-2]) / ys[-1]
            ex[f"{key}_{ax}"]["spread_finest4_rel"] = (max(ys) - min(ys)) / ys[-1]
    out["extrapolation_finest_meshes"] = ex
md.append("### Production mesh series, 300 K, tetrahedron (W/(m K))\n")
md.append("| mesh | N | LBTE xx / zz | RTA xx / zz | LBTE/RTA xx / zz | PCG it | rows s |")
md.append("|---|---|---|---|---|---|---|")
for r in series:
    md.append(f"| {'x'.join(map(str, r['mesh']))} | {r['N']} | {r['lbte'][0]:.2f} / {r['lbte'][1]:.2f} | "
              f"{r['rta'][0]:.2f} / {r['rta'][1]:.2f} | {r['lbte'][0]/r['rta'][0]:.4f} / {r['lbte'][1]/r['rta'][1]:.4f} | "
              f"{r['pcg_iterations']} | {r['rows_seconds'] and round(r['rows_seconds'])} |")

# ------------------------------------------------------------ Gaussian vs tetrahedron (RTA)
gauss = []
for m in [(15, 15, 9), (19, 19, 11), (23, 23, 13)]:
    tag = "".join(map(str, m))
    p = R / f"prod-rta-m{tag}-gauss" / "run.json"
    if not p.exists():
        continue
    j = json.loads(p.read_text())
    if j.get("status") != "finished":
        continue
    sig = j["result"]["sigmas"]
    kk = j["result"]["kappa_WmK[sigma][T][xx,yy,zz,yz,xz,xy]"]
    tet = next((r for r in series if r["mesh"] == list(m)), None)
    gauss.append({"mesh": list(m), "sigmas": sig, "rta": [[k[0][0], k[0][2]] for k in kk],
                  "tetra_rta": None if tet is None else tet["rta"]})
for s in [0.1, 0.2]:
    so = solve(f"ooc-m191911-prod-s{s}")
    if so:
        out[f"gauss_lbte_m191911_s{s}"] = {"lbte": [so["kappa_xx"], so["kappa_zz"]],
                                           "rta": [so["kappa_RTA_xx"], so["kappa_RTA_zz"]]}
out["gaussian_rta"] = gauss
if gauss:
    md.append("\n### Gaussian smearing vs tetrahedron, RTA, 300 K (W/(m K), xx / zz)\n")
    md.append("| mesh | tetrahedron | sigma=0.05 THz | 0.1 | 0.2 |")
    md.append("|---|---|---|---|---|")
    for g in gauss:
        t = g["tetra_rta"]
        cells = " | ".join(f"{v[0]:.2f} / {v[1]:.2f}" for v in g["rta"])
        md.append(f"| {'x'.join(map(str, g['mesh']))} | {'-' if t is None else f'{t[0]:.2f} / {t[1]:.2f}'} | {cells} |")

# ------------------------------------------------------------ kappa(T)
olyT = {r["T"]: r for r in bench["phono3py_kappa_T_31x31x17"]}
shT = {r["T"]: r for r in bench["shengbte_kappa_T"]}
fk = R / "prod-rta-m313117-NU" / "kappa-m313117.hdf5"
lbteT = {}
meta = R / "ooc-m313117-prod" / "meta.npz"
if meta.exists():
    temps = np.load(meta)["temps"]
    for k, T in enumerate(temps):
        s = solve("ooc-m313117-prod", k)
        if s:
            lbteT[float(T)] = s
kT = []
if fk.exists():
    with h5py.File(fk, "r") as h:
        t, kap = h["temperature"][()], h["kappa"][()]
    fi = R / "prod-rta-m313117-iso" / "kappa-m313117.hdf5"
    kap_iso = None
    if fi.exists():
        with h5py.File(fi, "r") as h:
            kap_iso = h["kappa"][()]
    for i, T in enumerate(t):
        row = {"T": float(T), "prod_rta": [float(kap[i, 0]), float(kap[i, 2])]}
        if kap_iso is not None:
            row["prod_rta_iso"] = [float(kap_iso[i, 0]), float(kap_iso[i, 2])]
        if float(T) in lbteT:
            s = lbteT[float(T)]
            row["prod_lbte"] = [s["kappa_xx"], s["kappa_zz"]]
        o = olyT.get(float(T))
        if o:
            row["olympics_phono3py_rta"] = [o["rta_xx"], o["rta_zz"]]
            row["olympics_phono3py_lbte"] = [o["lbte_xx"], o["lbte_zz"]]
        sh = shT.get(float(T))
        if sh:
            row["shengbte_rta"] = [sh["rta_a"], sh["rta_c"]]
            row["shengbte_iter"] = [sh["iter_a"], sh["iter_c"]]
        kT.append(row)
out["kappa_T"] = kT
if kT:
    md.append("\n### kappa(T), 31x31x17, W/(m K) (xx / zz)\n")
    md.append("| T | prod RTA | prod LBTE | Olympics phono3py RTA | Olympics phono3py LBTE | ShengBTE RTA | ShengBTE iterative |")
    md.append("|---|---|---|---|---|---|---|")
    def f2(v):
        return "-" if v is None else f"{v[0]:.1f} / {v[1]:.1f}"
    for r in kT:
        md.append(f"| {r['T']:.0f} | {f2(r.get('prod_rta'))} | {f2(r.get('prod_lbte'))} | {f2(r.get('olympics_phono3py_rta'))} | "
                  f"{f2(r.get('olympics_phono3py_lbte'))} | {f2(r.get('shengbte_rta'))} | {f2(r.get('shengbte_iter'))} |")

# ------------------------------------------------------------ isotope
iso = {}
s0, s1 = solve("ooc-m313117-prod", 2), solve("ooc-m313117-prod", 2, "_iso")
if s0 and s1:
    iso["lbte_300K"] = {"no_iso": [s0["kappa_xx"], s0["kappa_zz"]], "iso": [s1["kappa_xx"], s1["kappa_zz"]],
                        "rel_change": [(s1["kappa_xx"] - s0["kappa_xx"]) / s0["kappa_xx"],
                                       (s1["kappa_zz"] - s0["kappa_zz"]) / s0["kappa_zz"]],
                        "rta_no_iso": [s0["kappa_RTA_xx"], s0["kappa_RTA_zz"]],
                        "rta_iso_from_same_assembly": [s1["kappa_RTA_xx"], s1["kappa_RTA_zz"]]}
for r in kT:
    if "prod_rta_iso" in r:
        iso.setdefault("rta_rel_change_by_T", []).append(
            {"T": r["T"], "xx": (r["prod_rta_iso"][0] - r["prod_rta"][0]) / r["prod_rta"][0],
             "zz": (r["prod_rta_iso"][1] - r["prod_rta"][1]) / r["prod_rta"][1]})
out["isotope"] = iso
common.write_json(C / "results" / "convergence.json", out)
(C / "results" / "tables_production.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
print(json.dumps({k: out[k] for k in ("extrapolation_finest_meshes", "isotope") if k in out}, indent=1)[:6000])
