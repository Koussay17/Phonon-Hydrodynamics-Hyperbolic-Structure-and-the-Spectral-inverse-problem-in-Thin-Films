"""Audit, part 4: comparisons.

Inputs: out/aud_julia_core.json (independent model, dual-number Jacobian, adaptive Gauss-Kronrod),
out/aud2_twin_checks.json (twin through an own wrapper), out/aud3_mp_check.json (50-digit reference),
experiments/A/out/a4_fisher.json (branch table).
Fisher algebra here is a fourth implementation: normal matrix with Jacobi scaling, Cholesky in float64, and the
marginal variance from the mpmath (40 digits) inverse; no SVD, no code of common.py.
Usage: python aud4_compare.py -> out/aud4_compare.json (and a printed summary)
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
JL = json.loads((OUT / "aud_julia_core.json").read_text())
TW = json.loads((OUT / "aud2_twin_checks.json").read_text())
MPR = json.loads((OUT / "aud3_mp_check.json").read_text()) if (OUT / "aud3_mp_check.json").exists() else None
BR = json.loads((HERE.parent / "A" / "out" / "a4_fisher.json").read_text())
F = np.array(JL["freq80"])
PN = ["ln_tauR", "ln_lamz", "ln_lamr", "ln_C", "ln_R1", "ln_R2", "ln_lamAu", "ln_CAu", "ln_dAu", "ln_dAlN",
      "ln_lamS", "ln_CS", "ln_wscale"]
SETS = {"N0": ["ln_tauR"], "N1": ["ln_tauR", "g", "phi0", "td_ns"]}
SETS["N2"] = SETS["N1"] + ["ln_R1", "ln_R2"]
SETS["N3"] = SETS["N2"] + ["ln_lamz", "ln_C"]
SETS["N4"] = SETS["N3"] + ["ln_lamAu", "ln_CAu"]
SETS["N5"] = SETS["N4"] + ["ln_wscale"]
BANDS = {"10k-200M": (1e4, 2e8), "1M-200M": (1e6, 2e8), "10k-50M": (1e4, 5e7)}


def cvec(a):
    a = np.asarray(a, float)
    return a[..., 0] + 1j * a[..., 1]


def design(Jth, names, f):
    cols = []
    for n in names:
        if n == "g":
            cols.append(np.ones(f.size, complex))
        elif n == "phi0":
            cols.append(np.full(f.size, 1j))
        elif n == "td_ns":
            cols.append(-2j * np.pi * f * 1e-9)
        else:
            cols.append(Jth[:, PN.index(n)])
    return np.stack(cols, axis=1)


def sigma_normal(D, sa, sp_deg, mpcheck=False):
    A = np.vstack([D.real / sa, D.imag / np.radians(sp_deg)])
    Fm = A.T @ A
    d = np.sqrt(np.diag(Fm))
    Fs = Fm / np.outer(d, d)
    L = np.linalg.cholesky(Fs)
    e0 = np.zeros(Fs.shape[0]); e0[0] = 1.0
    y = np.linalg.solve(L, e0)
    s = float(np.sqrt(y @ y) / d[0])
    if not mpcheck:
        return s, None
    with mp.workdps(40):
        M = mp.matrix(Fs.tolist())
        Mi = M ** -1
        s_mp = float(mp.sqrt(Mi[0, 0]) / d[0])
    return s, s_mp


def relmax(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


res = {}
# ------------------------------------------------------------- 1. audit table vs branch table (Julia sigmas)
cmp_ = {}
worst = 0.0
for Lk, key in (("30", "L30_cat_aniso_beam"), ("100", "L100_cat_aniso_beam"), ("300", "L300_cat_aniso_beam"),
                ("500", "L500_cat_aniso_beam"), ("500_rta", "L500_cat_aniso_beam_RTAbase")):
    c = JL["cases"][Lk]
    for band, tk in (("10k-200M", "table"), ("1M-200M", "table_1M_200M"), ("10k-50M", "table_10k_50M")):
        for N in SETS:
            for s in ("0.01", "0.03", "0.1"):
                a = c[tk][N][s]["sigma"]
                b = BR["cases"][key]["sets"][N][band][s]["sigma"]
                r = a / b - 1
                cmp_[f"{Lk}|{band}|{N}|{s}"] = [a, b, r]
                worst = max(worst, abs(r))
res["julia_vs_branch_sigma"] = cmp_
res["julia_vs_branch_sigma_max_rel"] = worst
# 1D limit
worst1d = 0.0
for Lk, key in (("30", "L30_cat_aniso_1d"), ("100", "L100_cat_aniso_1d"), ("300", "L300_cat_aniso_1d"),
                ("500", "L500_cat_aniso_1d")):
    for N in ("N0", "N1", "N2", "N3", "N4"):
        for s in ("0.01", "0.1"):
            a = JL["cases"][Lk]["table_1d"][N][s]["sigma"]
            b = BR["cases"][key]["sets"][N]["10k-200M"][s]["sigma"]
            worst1d = max(worst1d, abs(a / b - 1))
res["julia_vs_branch_sigma_1d_max_rel"] = worst1d

# ------------------------------------------------------------- 2. own Fisher algebra on both Jacobians
fa = {}
for Lk in ("30", "100", "300", "500"):
    Jj = np.stack([cvec(col) for col in JL["cases"][Lk]["J_all"]], axis=1)
    Jt = np.stack([cvec(col) for col in TW["jacobian"][Lk]["J"]], axis=1)
    e = {"col_rel_err_twinFD_vs_AD": {n: relmax(Jt[:, j], Jj[:, j]) for j, n in enumerate(PN)}}
    for N in SETS:
        for s in (0.01, 0.1):
            sj, sj_mp = sigma_normal(design(Jj, SETS[N], F), 0.01, s, mpcheck=True)
            st, _ = sigma_normal(design(Jt, SETS[N], F), 0.01, s)
            e[f"{N}|{s}"] = {"sigma_AD_cholesky": sj, "sigma_AD_mp40": sj_mp, "sigma_twinFD_cholesky": st,
                             "julia_qr": JL["cases"][Lk]["table"][N][str(s)]["sigma"],
                             "rel_twinFD_vs_AD": st / sj - 1, "rel_chol_vs_qr": sj / JL["cases"][Lk]["table"][N][str(s)]["sigma"] - 1}
    fa[Lk] = e
res["fisher_crosscheck"] = fa

# ------------------------------------------------------------- 3. quadrature of the twin vs the independent reference
qd = {}
for Lk in ("30", "100", "300", "500"):
    ref = cvec(JL["cases"][Lk]["lnH"])
    q = {}
    for k, v in TW["quadrature"][Lk].items():
        q[k] = float(np.max(np.abs(cvec(v) - ref)))
    qd[Lk] = q
res["twin_lnH_vs_independent_max_abs"] = qd

# ------------------------------------------------------------- 4. FD error curve of the twin tau_R column
fdc = {}
for Lk in ("500", "30"):
    ref = cvec(JL["cases"][Lk]["J_tau"])
    fdc[Lk] = {k: relmax(cvec(v), ref) for k, v in TW["fd_curve"][Lk].items()}
res["twin_fd_curve_tau_col_rel_err"] = fdc

# ------------------------------------------------------------- 5. 50-digit reference
if MPR:
    mpc = {}
    for Lk in ("500", "30"):
        ref_lnH = cvec(JL["cases"][Lk]["lnH"])
        ref_J = cvec(JL["cases"][Lk]["J_tau"])
        tJ = cvec(TW["jacobian"][Lk]["J"][0])
        tl = cvec(TW["quadrature"][Lk]["step0.04_xlow-22"])
        rows = {}
        for i, d in MPR[Lk].items():
            i = int(i)
            lnH_mp = complex(float(d["lnH"][0]), float(d["lnH"][1]))
            J_mp = complex(float(d["Jtau"][0]), float(d["Jtau"][1]))
            rows[str(i)] = {"f": d["f"], "julia_lnH_abs_err": abs(ref_lnH[i] - lnH_mp), "twin_lnH_abs_err": abs(tl[i] - lnH_mp),
                            "julia_Jtau_rel_err": abs(ref_J[i] / J_mp - 1), "twin_Jtau_rel_err": abs(tJ[i] / J_mp - 1),
                            "Jtau_mp": [J_mp.real, J_mp.imag], "quad_self_check": d["quad_self_check_rel"]}
        mpc[Lk] = rows
    res["mp50_reference"] = mpc

# ------------------------------------------------------------- 6. conditioning (Julia, column-normalised whitened)
cond = {}
for Lk in ("100", "500", "30", "300"):
    cond[Lk] = {N: {s: JL["cases"][Lk]["table"][N][s]["cond"] for s in ("0.01", "0.1")} for N in SETS}
res["cond"] = cond

(OUT / "aud4_compare.json").write_text(json.dumps(res, indent=1))

print("Julia vs branch sigma, all L/sets/bands/noise (finite beams): max rel", worst, "; 1D:", worst1d)
for Lk in ("30", "100", "300", "500"):
    e = fa[Lk]
    print(Lk, "twin FD col err (tau, R1, lamz, CAu, wscale):",
          ["%.1e" % e["col_rel_err_twinFD_vs_AD"][n] for n in ("ln_tauR", "ln_R1", "ln_lamz", "ln_CAu", "ln_wscale")],
          "| max |rel| twinFD vs AD sigma:", "%.1e" % max(abs(v["rel_twinFD_vs_AD"]) for k, v in e.items() if "|" in k),
          "| chol vs qr:", "%.1e" % max(abs(v["rel_chol_vs_qr"]) for k, v in e.items() if "|" in k),
          "| mp40 vs chol:", "%.1e" % max(abs(v["sigma_AD_mp40"] / v["sigma_AD_cholesky"] - 1) for k, v in e.items() if "|" in k))
for Lk in qd:
    print(Lk, "twin lnH vs independent:", {k: "%.1e" % v for k, v in qd[Lk].items()})
print("FD curve:", json.dumps(fdc, indent=0))
if MPR:
    for Lk in res["mp50_reference"]:
        r = res["mp50_reference"][Lk]
        print(Lk, "mp50: julia lnH err max %.1e, twin lnH err max %.1e, julia J err max %.1e, twin J err max %.1e" % (
            max(v["julia_lnH_abs_err"] for v in r.values()), max(v["twin_lnH_abs_err"] for v in r.values()),
            max(v["julia_Jtau_rel_err"] for v in r.values()), max(v["twin_Jtau_rel_err"] for v in r.values())))
print("cond:", json.dumps(cond))
