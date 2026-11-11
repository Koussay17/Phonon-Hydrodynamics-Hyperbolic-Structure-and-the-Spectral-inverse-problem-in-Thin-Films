"""A4: local Fisher analysis of ln tau_R with the nuisance sets N0..N5 (declared in 00-question.md).

Noise: sigma(ln|H|) = 1 %, sigma_phi in {0.01, 0.03, 0.1} deg, independent Gaussian; 80 log-spaced frequencies
10 kHz - 200 MHz and the sub-bands 1-200 MHz, 10 kHz-50 MHz (subsets of the same 80 points).
Jacobian: twin, Richardson central differences in the log parameters (h = 2e-3), fixed trapezoid rule
(step 0.04, x_low = -22); instrument columns exact. Marginal sigma(ln tau_R) from the SVD of the column-normalised
whitened Jacobian; sets with condition number > 1e10 are flagged and not inverted.
w*_1sigma = sigma(ln tau_R), w*_2sigma = 2 sigma(ln tau_R) (changes of ln tau_R; relative change exp(w) - 1).
Every Cattaneo row is recomputed with exact autodiff of the independent model (indep.py) and the relative
difference recorded. Also: the phase noise needed to resolve the admissible half-width w(L) at 2 sigma.
Usage: python a4_fisher.py -> out/a4_fisher.json, out/a4_wstar_table.csv
"""
from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import brentq

from common import (FREQ, BANDS, SIG_PHI_DEG, SIG_A, OUT, SETS_CAT, SETS_GK, EXTRA, W_BY_L, Model, model_for,
                    jacobian, weighted, fisher, band_mask, load_base, jsonable)
from indep import Indep, base_from_model

T0 = time.perf_counter()
base = load_base()

ALLP_CAT = SETS_CAT["N5"] + ["ln_lamr", "ln_dAu", "ln_dAlN", "ln_wp", "ln_wr"]
ALLP_GK = SETS_GK["N5"] + ["ln_dAu", "ln_dAlN", "ln_wp", "ln_wr"]


def variants():
    for L in (30, 100, 300, 500):
        for geo in ("beam", "1d"):
            yield L, f"cat_aniso_{geo}", dict(law="cattaneo", geometry=geo, lamr_mode="rta")
            yield L, f"cat_iso_{geo}", dict(law="cattaneo", geometry=geo, lamr_mode="iso")
            yield L, f"gk0.5_{geo}", dict(law="gk", taul_factor=0.5, geometry=geo)
            yield L, f"gk2_{geo}", dict(law="gk", taul_factor=2.0, geometry=geo)
        yield L, "cat_xextrap_beam", dict(law="cattaneo", geometry="beam", lamr_mode="xextrap")
    b = base["500_rta"]
    for geo in ("beam", "1d"):
        yield 500, f"cat_aniso_{geo}_RTAbase", dict(model=Model(L_nm=500.0, tauR=b["tau_R_ps"] * 1e-12,
                                                                 lamz=b["lam_z"], lamr=b["lam_r"], geometry=geo))


def analyse(J, names, sets, extra_sets=None, mask_full=None):
    out = {}
    allsets = dict(sets)
    if extra_sets:
        allsets.update(extra_sets)
    for N, nm in allsets.items():
        if any(n not in names for n in nm):
            continue
        idx = [names.index(n) for n in nm]
        out[N] = {}
        for band in BANDS:
            mk = band_mask(FREQ, band)
            out[N][band] = {}
            for s in SIG_PHI_DEG:
                r = fisher(weighted(J[:, idx], sig_phi_deg=s, mask=mk))
                r["params"] = nm
                out[N][band][str(s)] = r
            # phase only and amplitude only (full band information split), sigma_phi = 0.1 deg
            if band == "10k-200M":
                Jw = weighted(J[:, idx], sig_phi_deg=0.1, mask=mk)
                nf = int(mk.sum())
                out[N][band]["phase_only_0.1"] = fisher(Jw[nf:])
                out[N][band]["amp_only"] = fisher(Jw[:nf])
    return out


def sigma_phi_needed(J, names, nm, w_target, band="10k-200M"):
    """Phase noise (deg) at which 2 sigma(ln tau_R) = ln(1 + w_target), sigma_A = 1 % fixed."""
    idx = [names.index(n) for n in nm]
    mk = band_mask(FREQ, band)

    def g(ls):
        r = fisher(weighted(J[:, idx], sig_phi_deg=10 ** ls, mask=mk))
        return np.log(2 * r["sigma"]) - np.log(np.log1p(w_target)) if r["sigma"] else np.inf
    try:
        return float(10 ** brentq(g, -9, 1, xtol=1e-6))
    except ValueError:
        return None


res = {"settings": {"freq": "80 log 1e4..2e8", "bands": BANDS, "sigma_A": SIG_A, "sigma_phi_deg": SIG_PHI_DEG,
                    "fd": "richardson h=2e-3", "quadrature": "trapezoid step 0.04, x_low -22",
                    "cond_max": 1e10}, "cases": {}}
csv = ["L_nm,variant,set,band,sigma_phi_deg,sigma_ln_tauR,wstar_1sigma_rel,wstar_2sigma_rel,cond,flag,indep_rel_diff"]
for L, name, kw in variants():
    t1 = time.perf_counter()
    m = kw.pop("model") if "model" in kw else model_for(L, **kw)
    gk = m.law == "gk"
    sets = SETS_GK if gk else SETS_CAT
    names = ALLP_GK if gk else ALLP_CAT
    J = jacobian(m, FREQ, names, h=2e-3, scheme="richardson")
    extra = None
    if not gk and name.startswith("cat_aniso"):
        extra = EXTRA
    elif gk:
        extra = {"N4-taul": SETS_CAT["N4"], "N3+taul": SETS_GK["N3"] + ["ln_taul"]}
    case = {"L_nm": L, "tauR_ps": m.tauR * 1e12, "lamz": m.lamz, "lamr": m.phys()["lamr"], "taul_ps": m.taul * 1e12,
            "law": m.law, "geometry": m.geometry, "names": names, "warnings": getattr(m, "_warnings", []),
            "sets": analyse(J, names, sets, extra)}
    # independent check (Cattaneo): exact autodiff of the quadrupole-recursion model
    if not gk:
        nm_i = [n for n in names if n not in ("ln_wp", "ln_wr")]
        ind = Indep(base_from_model(m), FREQ, geometry=m.geometry)
        Ji = ind.jacobian(nm_i)
        Jt = J[:, [names.index(n) for n in nm_i]]
        chk = analyse(Ji, nm_i, sets)
        case["indep_rel_diff"] = {}
        for N in chk:
            for band in BANDS:
                for s in SIG_PHI_DEG:
                    a = case["sets"][N][band][str(s)]["sigma"]
                    b_ = chk[N][band][str(s)]["sigma"]
                    case["indep_rel_diff"][f"{N}|{band}|{s}"] = (abs(a / b_ - 1) if (a and b_) else None)
        case["indep_lnH_max_abs_diff"] = float(np.max(np.abs(ind.lnH() - m.lnH(FREQ))))
    # phase noise needed for 2 sigma at the admissible half-width
    wt = W_BY_L.get(L, 0.052)
    case["w_admissible_used"] = wt
    case["sigma_phi_needed_deg_2sigma"] = {N: sigma_phi_needed(J, names, nm, wt) for N, nm in sets.items()
                                           if all(n in names for n in nm) and not (m.geometry == "1d" and N == "N5")}
    if L == 100:
        case["sigma_phi_needed_deg_2sigma_w0.08"] = {N: sigma_phi_needed(J, names, nm, 0.08)
                                                     for N, nm in sets.items()
                                                     if not (m.geometry == "1d" and N == "N5")}
    res["cases"][f"L{L}_{name}"] = case
    for N in case["sets"]:
        for band in BANDS:
            for s in SIG_PHI_DEG:
                r = case["sets"][N][band][str(s)]
                sg = r["sigma"]
                d = case.get("indep_rel_diff", {}).get(f"{N}|{band}|{s}")
                # H is linear in tau_R (A5, A7): the resolvable relative change of tau_R is sigma itself
                csv.append(f"{L},{name},{N},{band},{s},{sg if sg is not None else 'NA'},"
                           f"{sg if sg is not None else 'NA'},{2 * sg if sg is not None else 'NA'},"
                           f"{r['cond']:.4g},{r['flag'].replace(',', ';')},{d if d is not None else ''}")
    s01 = {N: case["sets"][N]["10k-200M"]["0.1"]["sigma"] for N in sets}
    s001 = {N: case["sets"][N]["10k-200M"]["0.01"]["sigma"] for N in sets}
    print(L, name, "sigma(ln tau) 0.01deg:", {k: (round(v, 4) if v else None) for k, v in s001.items()},
          "| 0.1deg:", {k: (round(v, 3) if v else None) for k, v in s01.items()},
          f"{time.perf_counter() - t1:.0f}s", flush=True)
    (OUT / "a4_fisher.json").write_text(json.dumps(jsonable(res), indent=1))
    (OUT / "a4_wstar_table.csv").write_text("\n".join(csv) + "\n")
res["seconds"] = time.perf_counter() - T0
(OUT / "a4_fisher.json").write_text(json.dumps(jsonable(res), indent=1))
print("done", res["seconds"])
