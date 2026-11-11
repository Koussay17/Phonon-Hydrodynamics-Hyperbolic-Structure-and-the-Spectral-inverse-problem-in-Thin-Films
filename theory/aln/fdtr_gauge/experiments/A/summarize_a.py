"""Summary tables of workstream A from the saved JSON (A2, A3, A4, A5, A6, A7). Prints Markdown.

Usage: python summarize_a.py > out/summary_tables.md
"""
from __future__ import annotations

import json

import numpy as np

from common import OUT

SETS = ["N0", "N1", "N2", "N3", "N4", "N5"]
S = ["0.01", "0.03", "0.1"]
WL = {30: 0.028, 100: 0.037, 300: 0.052, 500: None}


def fmt(x, d=3):
    if x is None:
        return "flag"
    if x >= 100:
        return f"{x:.3g}"
    if x >= 1:
        return f"{x:.2f}"
    return f"{x:.{d}f}"


def load(name):
    p = OUT / name
    return json.loads(p.read_text()) if p.exists() else None


a4 = load("a4_fisher.json")
if a4:
    c = a4["cases"]
    for band in ("10k-200M", "1M-200M", "10k-50M"):
        print(f"\n### w*_1sigma = sigma(ln tau_R), Cattaneo (lam_r = RTA kappa_x), beams 10/8 um, band {band}\n")
        print("| L (nm) | tau_R (ps) | sigma_phi (deg) | " + " | ".join(SETS) + " | w(L) |")
        print("|---|---|---|" + "---|" * len(SETS) + "---|")
        for L in (30, 100, 300, 500):
            k = f"L{L}_cat_aniso_beam"
            if k not in c:
                continue
            for s in S:
                row = [fmt(c[k]["sets"][N][band][s]["sigma"]) for N in SETS]
                wl = WL[L]
                print(f"| {L} | {c[k]['tauR_ps']:.2f} | {s} | " + " | ".join(row) +
                      f" | {('%.3f' % wl) if wl else '5.5-8 % (trend, see text)'} |")
    print("\n### Variants, full band, sigma_phi = 0.1 deg and 0.01 deg (sigma(ln tau_R))\n")
    print("| case | " + " | ".join(f"{N} (0.01/0.1)" for N in SETS) + " |")
    print("|---|" + "---|" * len(SETS))
    for k in c:
        r = []
        for N in SETS:
            if N not in c[k]["sets"]:
                r.append("—")
                continue
            a = c[k]["sets"][N]["10k-200M"]["0.01"]["sigma"]
            b = c[k]["sets"][N]["10k-200M"]["0.1"]["sigma"]
            r.append(f"{fmt(a)} / {fmt(b)}")
        print(f"| {k} | " + " | ".join(r) + " |")
    print("\n### Extra sets (Cattaneo aniso, beams, full band): sigma(ln tau_R) at 0.01 / 0.1 deg\n")
    for L in (30, 100, 300, 500):
        k = f"L{L}_cat_aniso_beam"
        if k not in c:
            continue
        e = {N: (c[k]["sets"][N]["10k-200M"]["0.01"]["sigma"], c[k]["sets"][N]["10k-200M"]["0.1"]["sigma"],
                 c[k]["sets"][N]["10k-200M"]["0.1"]["cond"], c[k]["sets"][N]["10k-200M"]["0.1"]["flag"])
             for N in c[k]["sets"] if N not in SETS}
        print(f"* {L} nm: " + "; ".join(f"{N}: {fmt(a)} / {fmt(b)} (cond {cd:.3g}{', ' + fl if fl else ''})"
                                         for N, (a, b, cd, fl) in e.items()))
    print("\n### Phase-only and amplitude-only (sigma_phi = 0.1 deg), Cattaneo aniso, beams, full band\n")
    for L in (30, 100, 300, 500):
        k = f"L{L}_cat_aniso_beam"
        if k not in c:
            continue
        print(f"* {L} nm: " + "; ".join(
            f"{N}: both {fmt(c[k]['sets'][N]['10k-200M']['0.1']['sigma'])}, phase {fmt(c[k]['sets'][N]['10k-200M']['phase_only_0.1']['sigma'])}, amp {fmt(c[k]['sets'][N]['10k-200M']['amp_only']['sigma'])}"
            for N in SETS))
    print("\n### Main correlation partner of ln tau_R (full band, 0.1 deg)\n")
    for L in (100, 500):
        k = f"L{L}_cat_aniso_beam"
        if k not in c:
            continue
        for N in SETS[1:]:
            r = c[k]["sets"][N]["10k-200M"]["0.1"]
            if r.get("corr_target") is None:
                continue
            ct = np.array(r["corr_target"])
            ct[0] = 0
            j = int(np.argmax(np.abs(ct)))
            print(f"* {L} nm {N}: cond {r['cond']:.3g}, sigma_cond {r['sigma_conditional']:.3g}, "
                  f"inflation {r['sigma'] / r['sigma_conditional']:.3g}, strongest |corr| with {r['params'][j]} "
                  f"({ct[j]:+.4f}); v_min {np.round(r['v_min'], 3).tolist()}")
    print("\n### Phase noise needed for 2 sigma(ln tau_R) = ln(1 + w(L)) (sigma_A = 1 %, full band)\n")
    for k in c:
        if "cat_aniso" in k or "gk" in k:
            sp = c[k]["sigma_phi_needed_deg_2sigma"]
            print(f"* {k} (w = {c[k]['w_admissible_used']}): " +
                  ", ".join(f"{N} {('%.2e' % v) if v else 'none'}" for N, v in sp.items()))
            if "sigma_phi_needed_deg_2sigma_w0.08" in c[k]:
                sp = c[k]["sigma_phi_needed_deg_2sigma_w0.08"]
                print(f"  (w = 0.08): " + ", ".join(f"{N} {('%.2e' % v) if v else 'none'}" for N, v in sp.items()))
    print("\n### Independent model (exact autodiff) vs twin (Richardson FD): max relative difference of sigma(ln tau_R)\n")
    for k in c:
        if "indep_rel_diff" in c[k]:
            v = [x for x in c[k]["indep_rel_diff"].values() if x is not None]
            print(f"* {k}: max {max(v):.2e} over {len(v)} (set, band, noise); max |d lnH| {c[k]['indep_lnH_max_abs_diff']:.2e}")

a3 = load("a3_dphi.json")
if a3:
    print("\n### Max |Delta phi| (deg) and frequency, all other parameters fixed\n")
    print("| case | w | max abs dphi (deg) | f at max (Hz) | max abs d ln|H| | nonlinearity (deg) |")
    print("|---|---|---|---|---|---|")
    for k, v in a3["cases"].items():
        for w, e in v["w"].items():
            if w.startswith("-") and "RTAbase" not in k:
                continue
            print(f"| {k} | {w} | {e['max_abs_dphi_deg']:.3e} | {e['f_at_max_Hz']:.3g} | "
                  f"{e.get('max_abs_dlnA', float('nan')):.2e} | {e.get('nonlinearity_max_abs_deg', float('nan')):.1e} |")

a2 = load("a2_numerics.json")
if a2:
    print("\n### A2 numerics\n")
    for k, v in a2["quadrature"].items():
        print(f"* quadrature {k}: lnH err vs adaptive {v['lnH_err_vs_adaptive']}; x_low -30: {v['lnH_err_xlow30_vs_adaptive']:.1e}")
        for st, d in v["sigma_rel_err_vs_step0.01"].items():
            m = max(max(x.values()) for x in d.values())
            print(f"    step {st}: tau column err {v['tau_col_err_vs_step0.01'][st]:.1e}, max sigma rel err {m:.1e}")
        if "sigma_adaptive_fd_h1e-2_rel_err_vs_step0.01" in v:
            d = v["sigma_adaptive_fd_h1e-2_rel_err_vs_step0.01"]
            print(f"    adaptive rule + FD h=1e-2: max sigma rel err {max(max(x.values()) for x in d.values()):.1e} "
                  f"({ {N: {s: f'{e:.1e}' for s, e in x.items()} for N, x in d.items()} })")
    for k, v in a2["fd_step"].items():
        print(f"* FD {k} (reference {v['reference']}):")
        for sch, d in v["sigma_rel_err"].items():
            m = max(max(x.values()) for x in d.values())
            print(f"    {sch}: tau column err {v['col_err'][sch]['ln_tauR']:.1e}, max sigma rel err {m:.1e} "
                  f"(N4: {max(d['N4'].values()):.1e})")
    for k, v in a2["one_d"].items():
        if k == "large_beam_convergence":
            for kk, vv in v.items():
                print(f"* large beam {kk}: max |d lnH| {np.array(vv['max_abs_dlnH'])}, order {np.round(vv['observed_order'], 4)}")
            continue
        m = max(max(x.values()) for x in v["sigma_rel_diff"].values())
        print(f"* 1D {k}: max |lnH twin - indep| (10 kHz-10 GHz) {v['lnH_max_abs_diff_10k_10G']:.1e}; "
              f"max sigma rel diff {m:.1e}")
    for k, v in a2["indep_beam"].items():
        m = max(max(x.values()) for x in v["sigma_rel_diff"].values())
        print(f"* finite beam {k}: max |d lnH| {v['lnH_max_abs_diff']:.1e}, tau col err {v['col_err']['ln_tauR']:.1e}, "
              f"max sigma rel diff {m:.1e}")
    for k, v in a2["frequency_grid_sigma_sqrtN"].items():
        print(f"* grid {k}: sigma*sqrt(N), N2 at 0.1 deg: " + ", ".join(f"N={N}: {r['N2']['0.1']:.3f}" for N, r in v.items())
              + "; N0 at 0.01: " + ", ".join(f"{N}: {r['N0']['0.01']:.4f}" for N, r in v.items()))
    for k, v in a2["inplane_tau"].items():
        print(f"* in-plane tau {k}: tau_x {v['tau_x_ps']:.1f} ps, max |d lnH| {v['max_abs_dlnH']:.1e}; "
              f"N2 0.1: scalar {v['sigma_scalar_tau']['N2']['0.1']:.4g}, z-only {v['sigma_tau_ip_eq_taux_dlnTauZ']['N2']['0.1']:.4g}; "
              f"N0 0.01: {v['sigma_scalar_tau']['N0']['0.01']:.4g} / {v['sigma_tau_ip_eq_taux_dlnTauZ']['N0']['0.01']:.4g}")

a5 = load("a5_profile.json")
if a5:
    print("\n### A5 profile likelihood (500 nm)\n")
    for k, v in a5["asimov"].items():
        print(f"* Asimov {k}: sigma_F {v['sigma_F']:.4g}; " + "; ".join(
            f"w={p['w']:+.3f}: {p['dchi2']:.4g} (Fisher {p['fisher_pred']:.4g})" if p["fisher_pred"] is not None
            else f"Fourier: {p['dchi2']:.4g}" for p in v["profile"]))
    for k, v in a5["seeds"].items():
        print(f"* seeds {k}: sigma_F {v['sigma_F']:.4g}, std {v['std_hat']:.4g} (ratio {v['std_over_sigmaF']:.3f}), "
              f"robust {v['robust_sigma_half_16_84']:.4g}, mean {v['mean_hat']:+.4f}, cov1 {v['coverage_1sigmaF']:.3f}, "
              f"cov2 {v['coverage_2sigmaF']:.3f}, at bounds {v['frac_at_bounds']:.3f}, chi2/dof {v['mean_chi2'] / v['dof']:.3f}")
    if "twin_refit_N2_0.01" in a5:
        for r in a5["twin_refit_N2_0.01"]:
            print(f"* twin refit seed {r['seed']}: twin {r['twin_hat']:+.6f} indep {r['indep_hat']:+.6f}")

a7 = load("a7_fourier_limit.json")
if a7:
    print("\n### A7 Fourier film and large changes (noiseless Delta chi2_min, nuisances refitted)\n")
    print("| L | set | sigma_phi | sigma_F | Fourier (tau_R=0) | x2 | x10 | lin. pred. Fourier |")
    print("|---|---|---|---|---|---|---|---|")
    for k, v in a7.items():
        L, N, s = k.split("|")
        print(f"| {L} | {N} | {s} | {v['sigma_F']:.4g} | {v['fourier']['dchi2']:.4g} | {v['x2']['dchi2']:.4g} | "
              f"{v['x10']['dchi2']:.4g} | {v['fourier']['linear_pred_(w/sigmaF)^2']:.4g} |")

a6 = load("a6_design.json")
if a6:
    print("\n### A6 large-beam Fisher convergence (Cattaneo 500 nm), sigma/sigma_1D - 1 at 0.1 deg\n")
    for k, v in a6["large_beam"].items():
        if k == "1d":
            continue
        print(f"* {k}: " + ", ".join(f"{N} {v['rel_diff_to_1d'][N]['0.1']:+.2e}" for N in v["rel_diff_to_1d"]))
    print("\n### A6 design sweep (independent model): sigma(ln tau_R) at 0.1 / 0.01 deg\n")
    print("| case | n_f | max dphi/dlntau (deg) | N0 | N1 | N2 | N3 |")
    print("|---|---|---|---|---|---|---|")
    for k, v in a6["design"].items():
        print(f"| {k} | {v['n_freq']} | {v['max_abs_dphi_dlntau_deg']:.3g} | " +
              " | ".join(f"{fmt(v[f'{N}|0.1'])} / {fmt(v[f'{N}|0.01'])}" for N in ("N0", "N1", "N2", "N3")) + " |")
    for k, v in a6.get("design_twin_check", {}).items():
        print(f"* twin check {k}: " + ", ".join(f"{N} {d['rel_diff']:.1e}" for N, d in v.items()))
