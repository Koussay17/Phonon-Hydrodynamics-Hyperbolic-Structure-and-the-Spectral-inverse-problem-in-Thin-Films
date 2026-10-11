"""Numerical audit, part 9: mesh-refinement table and observed-order analysis for the reference tau_mem and kappa.

Collects the reference values of na_spectral / na_bigmesh outputs (orbit convention as reported by the campaign, and
the gauge-invariant block convention), and for the similar-aspect family 6x6x4 / 9x9x6 / 12x12x8 (h = 1/N, ratios 1.5
and 4/3) solves Q(h) = Q_inf + C h^p for p from the ratio of successive differences
    F(p) = (Q_6 - Q_9) / (Q_9 - Q_12) = (6^-p - 9^-p) / (9^-p - 12^-p),
F increasing in p, F(0) = ln(3/2)/ln(4/3) = 1.4094. A sequence with F < F(0) is not in any convergent power-law
regime (p <= 0: growth at least logarithmic); Richardson extrapolation is then not applicable.
usage: python -B na_refine.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent / "na_out"


def load(name):
    return json.load(open(OUT / name))


def F(p):
    if abs(p) < 1e-12:
        return np.log(1.5) / np.log(4 / 3)
    return (6.0 ** -p - 9.0 ** -p) / (9.0 ** -p - 12.0 ** -p)


def order(q6, q9, q12):
    r = (q6 - q9) / (q9 - q12)
    try:
        p = brentq(lambda t: F(t) - r, -3.0, 6.0)
    except ValueError:
        p = float("nan")
    out = {"ratio_of_differences": r, "p": p}
    if p > 0:
        out["richardson_Q_inf"] = q12 + (q12 - q9) * 12.0 ** -p / (9.0 ** -p - 12.0 ** -p)
    return out


def main():
    rows = {}
    s553, s775, s996 = load("spectral_m553.json"), load("spectral_m775.json"), load("spectral_m996.json")
    for nm, S in (("5x5x3", s553), ("7x7x5", s775), ("9x9x6", s996)):
        rows[nm] = {"orbit": {"tau": [S["eig"]["cols"]["x"]["tau_ps"], S["eig"]["cols"]["z"]["tau_ps"]],
                              "kappa": [S["eig"]["cols"]["x"]["kappa"], S["eig"]["cols"]["z"]["kappa"]]},
                    "lam_min_H": S["eig"]["lam_min_H"], "cond_H": S["eig"]["cond_H"], "f_min": S["lowest_mode_freqs"][0],
                    "N_over_Nq_x": S["eig"]["cols"]["x"]["N"] / (np.prod(S["mesh"])),
                    "K_over_Nq_x": S["eig"]["cols"]["x"]["K"] / (np.prod(S["mesh"]))}
    for nm, f in (("6x6x4", "bigmesh_m664.json"), ("12x12x8", "bigmesh_m12128.json")):
        B = load(f)
        rows[nm] = {"orbit": {"tau": [B["x"]["tau_ps"], B["z"]["tau_ps"]], "kappa": [B["x"]["kappa"], B["z"]["kappa"]]},
                    "lam_min_H": B.get("lam_min_H"), "f_min": B["lowest_mode_freqs"][0],
                    "N_over_Nq_x": B["x"]["N"] / np.prod(B["mesh"]), "K_over_Nq_x": B["x"]["K"] / np.prod(B["mesh"])}
    for nm, f in (("5x5x3", "bigmesh_m553_block.json"), ("6x6x4", "bigmesh_m664_block.json"), ("7x7x5", "bigmesh_m775_block.json"),
                  ("9x9x6", "bigmesh_m996_block.json"), ("12x12x8", "bigmesh_m12128_block.json")):
        B = load(f)
        rows[nm]["block"] = {"tau": [B["x"]["tau_ps"], B["z"]["tau_ps"]], "kappa": [B["x"]["kappa"], B["z"]["kappa"]]}
    for nm, f in (("12x12x8", "bigmesh_m12128_sig0.07.json"), ("12x12x8", "bigmesh_m12128_sig0.05.json")):
        B = load(f)
        rows[nm][f"orbit_sigma{B['sigma_prime']}"] = {"tau": [B["x"]["tau_ps"], B["z"]["tau_ps"]], "kappa": [B["x"]["kappa"], B["z"]["kappa"]]}
    res = {"table": rows, "family_6_9_12": {}}
    for conv in ("orbit", "block"):
        for q in ("tau", "kappa"):
            for i, a in enumerate("xz"):
                v = [rows[m][conv][q][i] for m in ("6x6x4", "9x9x6", "12x12x8")]
                res["family_6_9_12"][f"{conv}_{q}_{a}"] = {"values": v, **order(*v)}
    (OUT / "refine_summary.json").write_text(json.dumps(res, indent=1))
    for m, r in rows.items():
        print(m, "tau", [round(t, 3) for t in r["orbit"]["tau"]], "kappa", [round(t, 2) for t in r["orbit"]["kappa"]],
              "| block tau", [round(t, 3) for t in r["block"]["tau"]], "kappa", [round(t, 2) for t in r["block"]["kappa"]],
              "| lam_min", r["lam_min_H"], "f_min", round(r["f_min"], 3), "N/Nq_x %.3e K/Nq_x %.3e" % (r["N_over_Nq_x"], r["K_over_Nq_x"]))
    for k, v in res["family_6_9_12"].items():
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else [round(t, 3) for t in vv]) for kk, vv in v.items()})


if __name__ == "__main__":
    main()
