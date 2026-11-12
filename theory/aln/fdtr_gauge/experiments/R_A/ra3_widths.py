"""R-A3: admissible half-width of the film memory time at 500 nm and at the regulator-normalisation lengths L/c,
from the existing three-mesh optima (E33) without new optimisation. INNER BOUNDS: the witnesses were optimised at
30/100/300 nm (orb 0.1) or 100 nm (calibrated classes), not at the lengths evaluated here.

Inputs: out/ra1_film_m775.json, out/ra1_film_m996.json (direct evaluation, R-A1), the stored E37 12x12x8 film metrics
(results_c/film3/*_m12128.json, L = 30, 100, 300, 1000 nm) and the E33 optimiser summary (e33_e37_summary.json).
Per witness W, mesh m, length L: ratio rho = tau_z,W,m(L) / tau_z,REF,m(L).
  robust excursion over the data meshes evaluated (7x7x5, 9x9x6): the one of smaller magnitude if both have the same
  sign, else 0 (the 5x5x3 film values at these L are not available; they were never binding at the optimised L);
  upper end U = max_W robust(+), lower end Lo = max_W robust(-), half-width w = (U + Lo) / 2 (campaign convention,
  e.g. 300 nm [0.9491, 1.0537] -> 5.2 %), also w_max = max(U, Lo).
12x12x8 (held out) at L not stored: (a) log-linear interpolation of rho in ln L between stored neighbours,
(b) the same with the shape of the 9x9x6 curve: rho12(L) = rho12(L0) + [rho12(L1) - rho12(L0)] s9(L),
s9 = [rho9(L) - rho9(L0)] / [rho9(L1) - rho9(L0)] (ill-conditioned when the excursion is not monotone in L),
(c) kept fraction: kappa = (rho12 - 1)/(rho9 - 1) at the stored lengths, interpolated in ln L and applied to the
direct 9x9x6 value (primary; validated by predicting 9x9x6 from 7x7x5). Below 30 nm (no stored 12x12x8 value): (b).
Output: out/ra3_widths.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
FILM3 = Path("C:/Users/Koussay/these/theory/aln/fc3_memory/results_c/film3")
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
CLASSES = {"orb0.1": WIT[:6], "cal1nz": WIT[6:8], "cal3wnz": WIT[8:10], "cal3w": WIT[10:12], "all": WIT}
STORED = ("30", "100", "300", "1000")

ra = {m: json.loads((OUT / f"ra1_film_{m}.json").read_text()) for m in ("m775", "m996")}
Ls = [float(x) for x in ra["m996"]["L_nm"]]
key = lambda L: f"{L:.6g}"  # noqa: E731


def rho(m, W, L, q="tau_z"):
    return ra[m][W][key(L)][q] / ra[m]["REF"][key(L)][q]


st12 = {}
for W in WIT:
    d = json.loads((FILM3 / f"{W}_m12128.json").read_text())
    st12[W] = {L: d[W][L]["tau_z"] / d["REF"][L]["tau_z"] for L in STORED}


def interp_loglin(r0, r1, L0, L1, L):
    t = (np.log(L) - np.log(L0)) / (np.log(L1) - np.log(L0))
    return r0 + t * (r1 - r0)


def bracket(L):
    s = [float(x) for x in STORED]
    if L <= s[0]:
        return s[0], s[1]
    for a, b in zip(s[:-1], s[1:]):
        if a <= L <= b:
            return a, b
    return s[-2], s[-1]


def kept_pred(r_fine0, r_fine1, r_coarse0, r_coarse1, r_coarseL, L0, L1, L):
    """Kept-fraction interpolation: kappa = (rho_fine - 1)/(rho_coarse - 1) at the stored lengths, interpolated in
    ln L and applied to the coarse value computed at L; falls back to log-linear in rho_fine if ill-conditioned."""
    x0, x1 = r_coarse0 - 1, r_coarse1 - 1
    if abs(x0) > 3e-3 and abs(x1) > 3e-3:
        k0, k1 = (r_fine0 - 1) / x0, (r_fine1 - 1) / x1
        if 0 < k0 < 2 and 0 < k1 < 2:
            return 1 + interp_loglin(k0, k1, L0, L1, L) * (r_coarseL - 1)
    return interp_loglin(r_fine0, r_fine1, L0, L1, L)


def rho12(W, L, method):
    if key(L) in st12[W]:
        return st12[W][key(L)]
    L0, L1 = bracket(L)
    r0, r1 = st12[W][key(L0)], st12[W][key(L1)]
    if method == "loglin":
        return interp_loglin(r0, r1, L0, L1, L)
    if method == "kept":
        return kept_pred(r0, r1, rho("m996", W, L0), rho("m996", W, L1), rho("m996", W, L), L0, L1, L)
    a0, a1, aL = rho("m996", W, L0), rho("m996", W, L1), rho("m996", W, L)
    s9 = (aL - a0) / (a1 - a0) if abs(a1 - a0) > 1e-12 else 0.0
    return r0 + s9 * (r1 - r0)


res = {"note": "inner bounds; witnesses optimised at 30/100/300 nm (orb 0.1) and 100 nm (calibrated)", "per_L": {},
       "interp_validation": {}, "w500": {}}
# validation of the 12x12x8 interpolation where the direct value exists
for m, other in (("m996", "m775"), ("m775", "m996")):
    errs_ll, errs_sh = [], []
    for W in WIT:
        for L in Ls:
            if key(L) in STORED:
                continue
            L0, L1 = bracket(L)
            r0, r1, rL = rho(m, W, L0), rho(m, W, L1), rho(m, W, L)
            ll = interp_loglin(r0, r1, L0, L1, L)
            b0, b1, bL = rho(other, W, L0), rho(other, W, L1), rho(other, W, L)
            s = (bL - b0) / (b1 - b0) if abs(b1 - b0) > 1e-12 else 0.0
            sh = r0 + s * (r1 - r0)
            if L >= 30:
                errs_ll.append((abs(ll - rL), W, L))
                errs_sh.append((abs(sh - rL), W, L))
    errs_kp = []
    if m == "m996":
        for W in WIT:
            for L in Ls:
                if key(L) in STORED or L < 30:
                    continue
                L0, L1 = bracket(L)
                pr = kept_pred(rho("m996", W, L0), rho("m996", W, L1), rho("m775", W, L0), rho("m775", W, L1),
                               rho("m775", W, L), L0, L1, L)
                errs_kp.append((abs(pr - rho("m996", W, L)), W, L))
    res["interp_validation"].setdefault("kept_m775_to_m996", {})
    if errs_kp:
        res["interp_validation"]["kept_m775_to_m996"] = {
            "max_abs": max(errs_kp)[0], "worst": list(max(errs_kp)[1:]),
            "500nm_max_abs": max(e[0] for e in errs_kp if abs(e[2] - 500) < 1e-9)}
    res["interp_validation"][m] = {
        "loglin_max_abs": max(errs_ll)[0], "loglin_worst": list(max(errs_ll)[1:]),
        "shape_from_" + other + "_max_abs": max(errs_sh)[0], "shape_worst": list(max(errs_sh)[1:]),
        "loglin_500nm_max_abs": max(e[0] for e in errs_ll if abs(e[2] - 500) < 1e-9),
        "shape_500nm_max_abs": max(e[0] for e in errs_sh if abs(e[2] - 500) < 1e-9)}

for L in Ls:
    row = {"witness": {}}
    for W in WIT:
        r7, r9 = rho("m775", W, L), rho("m996", W, L)
        x7, x9 = r7 - 1, r9 - 1
        rob = (x7 if abs(x7) < abs(x9) else x9) if x7 * x9 > 0 else 0.0
        row["witness"][W] = {"rho_m775": r7, "rho_m996": r9, "robust_exc": rob,
                             "rho_m12128_loglin": rho12(W, L, "loglin") if L >= 30 else None,
                             "rho_m12128_shape": rho12(W, L, "shape"),
                             "rho_m12128_kept": rho12(W, L, "kept") if L >= 30 else rho12(W, L, "shape")}
    for cname, ws in CLASSES.items():
        up = max([row["witness"][W]["robust_exc"] for W in ws] + [0.0])
        lo = max([-row["witness"][W]["robust_exc"] for W in ws] + [0.0])
        u12 = max([row["witness"][W]["rho_m12128_kept"] - 1 for W in ws] + [0.0])
        l12 = max([1 - row["witness"][W]["rho_m12128_kept"] for W in ws] + [0.0])
        us = max([row["witness"][W]["rho_m12128_shape"] - 1 for W in ws] + [0.0])
        ls_ = max([1 - row["witness"][W]["rho_m12128_shape"] for W in ws] + [0.0])
        row[cname] = {"upper": up, "lower": lo, "w_half": 0.5 * (up + lo), "w_max": max(up, lo),
                      "upper_m12128": u12, "lower_m12128": l12, "w_half_m12128": 0.5 * (u12 + l12),
                      "w_half_m12128_shape": 0.5 * (us + ls_)}
    res["per_L"][key(L)] = row

# 500 nm summary
r5 = res["per_L"][key(500.0)]
summ = {}
for cname in CLASSES:
    summ[cname] = dict(r5[cname])
    ws = CLASSES[cname]
    summ[cname]["best_upper_witness"] = max(ws, key=lambda W: r5["witness"][W]["robust_exc"])
    summ[cname]["best_lower_witness"] = min(ws, key=lambda W: r5["witness"][W]["robust_exc"])
    ll_u = max([r5["witness"][W]["rho_m12128_loglin"] - 1 for W in ws] + [0.0])
    ll_l = max([1 - r5["witness"][W]["rho_m12128_loglin"] for W in ws] + [0.0])
    summ[cname]["w_half_m12128_loglin"] = 0.5 * (ll_u + ll_l)
res["w500"] = summ
# optimised half-widths at 30/100/300 nm (E33 robust ratios) and a power-law trend (heuristic, not a bound)
e33 = json.loads((FILM3 / "e33_e37_summary.json").read_text())
opt = {}
for L in (30, 100, 300):
    mx, mn = e33[f"orb0.1_m3L{L}max"]["robust_ratio"], e33[f"orb0.1_m3L{L}min"]["robust_ratio"]
    opt[str(L)] = 0.5 * ((mx - 1) + (1 - mn))
lx = np.log([30.0, 100.0, 300.0])
ly = np.log([opt["30"], opt["100"], opt["300"]])
p_all = np.polyfit(lx, ly, 1)
p_hi = np.polyfit(lx[1:], ly[1:], 1)
res["trend_heuristic"] = {"opt_half_width_orb0.1": opt, "slope_30_300": p_all[0], "slope_100_300": p_hi[0],
                          "w500_powerlaw_30_300": float(np.exp(np.polyval(p_all, np.log(500.0)))),
                          "w500_powerlaw_100_300": float(np.exp(np.polyval(p_hi, np.log(500.0))))}
# how much of its optimised excursion each witness keeps at its own L in this evaluation (reproduction check)
rep = {}
for W in WIT:
    Lw = float(W.split("_m3L")[1].rstrip("maxin"))
    rep[W] = {"e33_robust": e33[W]["robust_ratio"], "e33_m996": e33[W]["fitted_film_z"]["m996"],
              "ra1_m996": rho("m996", W, Lw), "ra1_m775": rho("m775", W, Lw),
              "e33_m775": e33[W]["fitted_film_z"]["m775"]}
res["reproduction_at_optimised_L"] = rep
(OUT / "ra3_widths.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res["interp_validation"], indent=1))
for c, v in summ.items():
    print(c, {k: (round(x, 4) if isinstance(x, float) else x) for k, v2 in [(0, 0)] for k, x in v.items()})
print(json.dumps(res["trend_heuristic"], indent=1))
for L in Ls:
    o, a = res["per_L"][key(L)]["orb0.1"], res["per_L"][key(L)]["all"]
    print(f"{L:8.2f} orb w {o['w_half']:.4f} [{-o['lower']:+.4f},{o['upper']:+.4f}] 12: {o['w_half_m12128']:.4f} | "
          f"all w {a['w_half']:.4f} [{-a['lower']:+.4f},{a['upper']:+.4f}] 12: {a['w_half_m12128']:.4f}")
