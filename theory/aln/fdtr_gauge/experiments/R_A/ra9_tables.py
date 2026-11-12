"""R-A9: markdown tables for 07-experiments.md (entries R-A1 ... R-A6) printed from the raw outputs in out/.
usage: python -B ra9_tables.py > out/ra9_tables.txt
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
J = lambda name: json.loads((OUT / name).read_text())  # noqa: E731
key = lambda L: f"{float(L):.6g}"  # noqa: E731
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]

r7, r9 = J("ra1_film_m775.json"), J("ra1_film_m996.json")
print("## R-A1 checks")
for m, r in (("7x7x5", r7), ("9x9x6", r9)):
    print(f"* {m}: n = {r['n_modes']}, events {r['n_events']} (normal {r['n_normal']}), setup {r['setup_seconds']:.0f} s, "
          f"max rel. deviation from E37 (tau_z, kappa_z, tauRTA_z at 30/100/300/1000 nm, all tensors) "
          f"{r['check_vs_E37_max_rel']:.1e}, peak working set {max(v for k, v in r.items() if k.endswith('_peak_wset_gb')):.2f} GB, "
          f"total {r['seconds']:.0f} s; mixed N/U orbits {r['n_mixed_NU_orbits']}/{r['n_event_orbits']} "
          f"(REF rate share {r['REF_rate_share_mixed']:.3f}), normal share {r['REF_rate_share_normal']:.3f}, "
          f"gamma_g {r['gamma_g']:.3f}")

w = J("ra3_widths.json")
print("\n## R-A3 500 nm, witness ratios tau_z(500 nm)/REF")
print("| witness | 7x7x5 | 9x9x6 | robust | 12x12x8 (shape) | 12x12x8 (log-lin) |")
print("|---|---|---|---|---|---|")
r5 = w["per_L"][key(500)]["witness"]
for W in WIT:
    x = r5[W]
    print(f"| {W} | {x['rho_m775']:.4f} | {x['rho_m996']:.4f} | {x['robust_exc']:+.4f} | {x['rho_m12128_shape']:.4f} | "
          f"{x['rho_m12128_loglin']:.4f} |")
print("\n| class | upper | lower | w_half (robust) | w_half 12x12x8 (shape) | w_half 12x12x8 (log-lin) |")
print("|---|---|---|---|---|---|")
for c, v in w["w500"].items():
    print(f"| {c} | {v['upper']:+.4f} | {-v['lower']:+.4f} | {v['w_half']:.4f} | {v['w_half_m12128']:.4f} | "
          f"{v['w_half_m12128_loglin']:.4f} |")
print("\ninterpolation validation:", json.dumps(w["interp_validation"]))
print("trend:", json.dumps(w["trend_heuristic"]))
print("reproduction at optimised L:", {k: (round(v["ra1_m996"], 4), round(v["e33_m996"], 4)) for k, v in
                                       w["reproduction_at_optimised_L"].items()})
print("\n## R-A3 half-width vs L (inner bounds): L, orb w_half [lo, up], 12x12x8; all w_half [lo, up], 12x12x8")
for Lk in sorted(w["per_L"], key=float):
    o, a = w["per_L"][Lk]["orb0.1"], w["per_L"][Lk]["all"]
    print(f"| {float(Lk):.4g} | {o['w_half']:.4f} | [{-o['lower']:+.4f}, {o['upper']:+.4f}] | {o['w_half_m12128']:.4f} | "
          f"{a['w_half']:.4f} | [{-a['lower']:+.4f}, {a['upper']:+.4f}] | {a['w_half_m12128']:.4f} |")

r2 = J("ra2_rta_film_c.json")
print("\n## R-A2 base values (L, c): tau_R ps / lam_z / lam_r ; RTA identity", r2["identity_c_vs_L_over_c_max_rel"],
      "A1 check", r2["check_A1_rta_max_rel"], "capweights", r2["check_capweights_vs_modekappa_tau_max_rel"])
CS = r2["c_values"]
print("| L | " + " | ".join(f"c = {c:.3g}" for c in CS) + " |")
print("|---|" + "---|" * len(CS))
for L in (30, 100, 300, 500):
    print(f"| {L} | " + " | ".join(
        f"{r2['base_box'][f'{L}|{c:.6g}']['tau_R_ps']:.2f} / {r2['base_box'][f'{L}|{c:.6g}']['lam_z']:.0f}"
        for c in CS) + " |")
print("\nwidths w(L/c), orb robust / all robust:")
for L in (30, 100, 300, 500):
    print(f"| {L} | " + " | ".join(
        f"{w['per_L'][key(round(L / c, 6))]['orb0.1']['w_half']:.4f} / {w['per_L'][key(round(L / c, 6))]['all']['w_half']:.4f}"
        for c in CS) + " |")

b = J("ra4_box.json")
print("\n## R-A4 checks:", json.dumps(b["checks"]["A4_reproduction_rel"]), "twin corners max rel",
      b["checks"]["twin_vs_indep_corners_max_rel"])
for cl in ("orb0.1", "all"):
    for basis in ("robust", "m12128"):
        print(f"\n### box summary, width class {cl}, basis {basis}: base m (bits, % var removed) | box min m [at] (bits) | "
              f"box max m | fraction m <= 1")
        print("| L | set | 0.01 deg | 0.03 deg | 0.1 deg |")
        print("|---|---|---|---|---|")
        for L in (30, 100, 300, 500):
            for N in ("N1", "N2", "N3"):
                cells = []
                for sp in (0.01, 0.03, 0.1):
                    v = b["summary"][f"{cl}|{basis}|{L}|{N}|{sp}"]
                    at = v["at"]
                    cells.append(f"{v['base_ratio']:.3g} ({v['bits_at_base']:.2g} b, {100 * v['var_removed_at_base']:.2g} %) | "
                                 f"min {v['min_ratio']:.3g} [c={at['c']:.3g}, {at['var']}, R1={at['R1']:.0e}, R2={at['R2']:.0e}] "
                                 f"({v['bits_at_min']:.2g} b) | max {v['max_ratio']:.3g} | flip {v['frac_flip_1sigma']:.2f}")
                print(f"| {L} | {N} | " + " || ".join(cells) + " |")
if "au_subbox" in b:
    print("\nAu sub-box:", json.dumps({k: (round(v["min_sigma"], 4), v["at"], round(v["min_ratio_all"], 3))
                                       for k, v in b["au_subbox"].items()}))

print("\n## R-A2 check: c-scaling of tau_z, event operator (REF) against the RTA 64x64x35 route, tau(L/c)/tau(L)")
for L in (30, 100, 300, 500):
    row = []
    for c in r2["c_values"]:
        Lc = key(round(L / c, 6))
        ev9 = r9["REF"][Lc]["tau_z"] / r9["REF"][key(L)]["tau_z"]
        ev7 = r7["REF"][Lc]["tau_z"] / r7["REF"][key(L)]["tau_z"]
        rt = r2["film"][Lc]["tau_z"] / r2["film"][key(L)]["tau_z"]
        row.append(f"{rt:.3f}/{ev7:.3f}/{ev9:.3f}")
    print(f"| {L} | " + " | ".join(row) + " |")
