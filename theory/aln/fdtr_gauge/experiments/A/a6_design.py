"""A6: (1) large-beam convergence of the Fisher numbers to the 1D limit (twin), (2) model-internal design sweep:
what band edge / transducer thickness would make the admissible gap resolvable.

(1) Cattaneo, 500 nm: sigma(ln tau_R) for N0..N4 with Gaussian beams wp = 10 ... 3000 um (wr = 0.8 wp) against the
    1D limit; twin, Richardson FD (h = 2e-3), fixed rule.
(2) Independent model (exact autodiff), Cattaneo, L = 300 and 500 nm: f_max in {0.2, 1, 5, 20} GHz at constant
    point density (80 points over 1e4..2e8 Hz, i.e. 18.6 per decade), Au thickness 80 / 40 / 20 nm, sets N0..N3,
    sigma_phi 0.01 and 0.1 deg. Spot checks against the twin. MODEL-INTERNAL ONLY: surface absorption, no
    electron-phonon coupling, no optical penetration, Fourier Au and sapphire; above a few hundred MHz these
    omissions are not small, so (2) only indicates the scale of the required bandwidth.
Usage: python a6_design.py -> out/a6_design.json
"""
from __future__ import annotations

import json
import time

import numpy as np

from common import FREQ, OUT, SETS_CAT, W_BY_L, model_for, jacobian, weighted, fisher, jsonable
from indep import Indep, base_from_model

T0 = time.perf_counter()
res = {"large_beam": {}, "design": {}}
names = SETS_CAT["N4"]

# (1) large-beam Fisher convergence
m1 = model_for(500, geometry="1d")
J1 = jacobian(m1, FREQ, names, h=2e-3, scheme="richardson")
s1d = {N: {s: fisher(weighted(J1[:, [names.index(n) for n in SETS_CAT[N]]], sig_phi_deg=s))["sigma"]
           for s in (0.01, 0.1)} for N in ("N0", "N1", "N2", "N3", "N4")}
res["large_beam"]["1d"] = s1d
for w in (10e-6, 30e-6, 100e-6, 300e-6, 1e-3, 3e-3):
    m = model_for(500, wp=w, wr=0.8 * w)
    J = jacobian(m, FREQ, names, h=2e-3, scheme="richardson")
    sb = {N: {s: fisher(weighted(J[:, [names.index(n) for n in SETS_CAT[N]]], sig_phi_deg=s))["sigma"]
              for s in (0.01, 0.1)} for N in ("N0", "N1", "N2", "N3", "N4")}
    res["large_beam"][f"wp={w * 1e6:g}um"] = {"sigma": sb,
                                             "rel_diff_to_1d": {N: {s: sb[N][s] / s1d[N][s] - 1 for s in sb[N]}
                                                                for N in sb}}
    print("beam", w, {N: round(sb[N][0.1] / s1d[N][0.1] - 1, 6) for N in sb}, flush=True)

# (2) design sweep with the independent model
dens = 80 / np.log10(2e8 / 1e4)
for L in (300, 500):
    m = model_for(L)
    for dAu in (80e-9, 40e-9, 20e-9):
        for fmax in (2e8, 1e9, 5e9, 2e10):
            nf = int(round(dens * np.log10(fmax / 1e4))) + 1
            f = np.geomspace(1e4, fmax, nf)
            ind = Indep(base_from_model(m, dAu=dAu), f)
            J = ind.jacobian(SETS_CAT["N3"])
            row = {"n_freq": nf, "max_abs_dphi_dlntau_deg": float(np.degrees(np.max(np.abs(J[:, 0].imag)))),
                   "f_at_max": float(f[np.argmax(np.abs(J[:, 0].imag))])}
            for N in ("N0", "N1", "N2", "N3"):
                idx = [SETS_CAT["N3"].index(n) for n in SETS_CAT[N]]
                for s in (0.01, 0.1):
                    r = fisher(weighted(J[:, idx], sig_phi_deg=s))
                    row[f"{N}|{s}"] = r["sigma"]
                    row[f"{N}|{s}|cond"] = r["cond"]
            res["design"][f"L{L}|dAu{dAu * 1e9:g}nm|fmax{fmax:g}"] = row
            print("design", L, dAu, fmax, nf, {k: (round(v, 4) if isinstance(v, float) else v)
                                               for k, v in row.items() if "|" in k and "cond" not in k}, flush=True)
# twin spot checks of the sweep (80 nm Au, 1 GHz and 5 GHz; 20 nm Au, 1 GHz) at 500 nm
chk = {}
for dAu, fmax in ((80e-9, 1e9), (80e-9, 5e9), (20e-9, 1e9)):
    nf = int(round(dens * np.log10(fmax / 1e4))) + 1
    f = np.geomspace(1e4, fmax, nf)
    m = model_for(500)
    mt = model_for(500)
    th = {"ln_dAu": np.log(dAu / 80e-9)}
    # twin Jacobian at the modified Au thickness: perturb around theta = th
    from common import Model  # noqa: F401
    cols = []
    for n in SETS_CAT["N3"]:
        if n == "g":
            cols.append(np.ones(nf, complex))
        elif n == "phi0":
            cols.append(np.full(nf, 1j))
        elif n == "td_ns":
            cols.append(-1j * 2 * np.pi * f * 1e-9)
        else:
            def cd(hh, n=n):
                tp, tm = dict(th), dict(th)
                tp[n] = tp.get(n, 0.0) + hh
                tm[n] = tm.get(n, 0.0) - hh
                return np.log(mt.thermal(f, mt.phys(tp)) / mt.thermal(f, mt.phys(tm))) / (2 * hh)
            cols.append((4 * cd(1e-3) - cd(2e-3)) / 3)
    Jt = np.stack(cols, axis=1)
    Ji = Indep(base_from_model(m, dAu=dAu), f).jacobian(SETS_CAT["N3"])
    d = {}
    for N in ("N0", "N1", "N2", "N3"):
        idx = [SETS_CAT["N3"].index(n) for n in SETS_CAT[N]]
        a = fisher(weighted(Jt[:, idx], sig_phi_deg=0.1))["sigma"]
        b = fisher(weighted(Ji[:, idx], sig_phi_deg=0.1))["sigma"]
        d[N] = {"twin": a, "indep": b, "rel_diff": abs(a / b - 1)}
    chk[f"dAu{dAu * 1e9:g}|fmax{fmax:g}"] = d
    print("twin check", dAu, fmax, {N: v["rel_diff"] for N, v in d.items()}, flush=True)
res["design_twin_check"] = chk

# (3) phase-only and amplitude-only information (the gain column is zero in phase-only data and the offset/delay
#     columns are zero in amplitude-only data, so they are removed from the respective sets)
po = {}
for L in (30, 100, 300, 500):
    m = model_for(L)
    nm5 = SETS_CAT["N5"]
    J = Indep(base_from_model(m), FREQ).jacobian(nm5)
    for N in ("N0", "N1", "N2", "N3", "N4", "N5"):
        nph = [n for n in SETS_CAT[N] if n != "g"]
        nam = [n for n in SETS_CAT[N] if n not in ("phi0", "td_ns")]
        for s in (0.01, 0.1):
            Jw = weighted(J[:, [nm5.index(n) for n in nph]], sig_phi_deg=s)
            rp = fisher(Jw[FREQ.size:])
            Jw2 = weighted(J[:, [nm5.index(n) for n in nam]], sig_phi_deg=s)
            ra = fisher(Jw2[:FREQ.size])
            rb = fisher(weighted(J[:, [nm5.index(n) for n in SETS_CAT[N]]], sig_phi_deg=s))
            po[f"L{L}|{N}|{s}"] = {"both": rb["sigma"], "phase_only": rp["sigma"], "amp_only": ra["sigma"],
                                   "phase_only_cond": rp["cond"], "amp_only_cond": ra["cond"]}
    print("phase/amp", L, {k: (round(v["both"], 3), round(v["phase_only"], 3) if v["phase_only"] else None,
                              round(v["amp_only"], 3) if v["amp_only"] else None)
                           for k, v in po.items() if k.startswith(f"L{L}|") and k.endswith("|0.1")}, flush=True)
res["phase_amp_split"] = po

# (4) robustness against the size of tau_R itself (the regulator is heuristic): sigma(ln tau_R) for tau_R scaled by
#     0.5 ... 35 at 300 and 500 nm (Cattaneo, beams, full band); linear regime predicts sigma proportional to 1/tau_R
ts = {}
for L in (300, 500):
    m = model_for(L)
    for fac in (0.5, 1.0, 2.0, 5.0, 10.0, 35.0):
        J = Indep(base_from_model(m, tauR=m.tauR * fac, taur_ip=m.tauR * fac), FREQ).jacobian(SETS_CAT["N5"])
        row = {}
        for N in ("N0", "N1", "N2", "N3", "N4", "N5"):
            idx = [SETS_CAT["N5"].index(n) for n in SETS_CAT[N]]
            for s in (0.01, 0.1):
                row[f"{N}|{s}"] = fisher(weighted(J[:, idx], sig_phi_deg=s))["sigma"]
        ts[f"L{L}|x{fac:g}"] = row
        print("tau scale", L, fac, {k: round(v, 4) for k, v in row.items() if k.endswith("|0.1") or k == "N0|0.01"},
              flush=True)
res["tau_scaling"] = ts
res["seconds"] = time.perf_counter() - T0
(OUT / "a6_design.json").write_text(json.dumps(jsonable(res), indent=1))
print("done", res["seconds"])
