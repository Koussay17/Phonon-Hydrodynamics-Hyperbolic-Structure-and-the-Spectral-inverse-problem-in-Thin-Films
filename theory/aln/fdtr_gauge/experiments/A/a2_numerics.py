"""A2: numerical verification of the twin-based Fisher analysis.

1. Hankel quadrature: fixed trapezoid step 0.16 ... 0.01 and x_low -22 / -30 against the adaptive rule
   (epsrel 1e-12); errors on H, on the tau_R column and on sigma(ln tau_R) for N0..N5.
2. Finite-difference step of the Jacobian: central h = 1e-1 ... 1e-5 and Richardson (h = 2e-3) against exact
   autodiff of the independent model (Cattaneo) or against Richardson (GK); errors on sigma(ln tau_R).
3. Independent 1D quadrupole (indep.py): twin kernel(k=0) vs closed-form chain (Cattaneo and 1D GK lam_eff),
   10 kHz - 10 GHz; large-beam convergence of the twin Gaussian average to the 1D limit (observed order);
   sigma(ln tau_R) in the 1D limit from both implementations.
4. Independent finite-beam model (Cattaneo): H, Jacobian and sigma(ln tau_R) at every thickness.
5. Frequency grid: sigma(ln tau_R) * sqrt(N) for N = 20 ... 320 log-spaced points (information density).
6. In-plane relaxation time: the twin uses one tau_R for q_z and q_r; effect of tau_r,in-plane = tau_x(L).
Usage: python a2_numerics.py   -> out/a2_numerics.json
"""
from __future__ import annotations

import json
import time

import numpy as np

from common import (FREQ, SETS_CAT, SETS_GK, SIG_PHI_DEG, OUT, Model, jacobian, weighted, fisher, model_for,
                    load_base, jsonable)
from indep import Indep, base_from_model, gl_rule

T0 = time.perf_counter()
base = load_base()
res = {}


def sigmas(J, names, sets, sig=(0.01, 0.1)):
    out = {}
    for N, nm in sets.items():
        idx = [names.index(n) for n in nm]
        out[N] = {str(s): fisher(weighted(J[:, idx], sig_phi_deg=s))["sigma"] for s in sig}
    return out


def maxrel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def sig_err(sa, sb):
    return {N: {s: abs(sa[N][s] / sb[N][s] - 1) for s in sa[N]} for N in sa}


# ---------------------------------------------------------------- 1. quadrature
q = {}
for L, law, tf in ((500, "cattaneo", 0.0), (30, "cattaneo", 0.0), (500, "gk", 2.0), (30, "gk", 0.5)):
    key = f"{law}_L{L}" + (f"_taul{tf}" if law == "gk" else "")
    sets = SETS_GK if law == "gk" else SETS_CAT
    names = sets["N5"]
    mref = model_for(L, law=law, taul_factor=tf, quad="adaptive", epsrel=1e-12)
    lref = mref.lnH(FREQ)
    entry = {"lnH_err_vs_adaptive": {}, "tau_col_err_vs_step0.01": {}, "sigma_rel_err_vs_step0.01": {}}
    Js = {}
    for st in (0.16, 0.08, 0.04, 0.02, 0.01):
        m = model_for(L, law=law, taul_factor=tf, step=st)
        entry["lnH_err_vs_adaptive"][str(st)] = float(np.max(np.abs(m.lnH(FREQ) - lref)))
        Js[st] = jacobian(m, FREQ, names, h=2e-3, scheme="richardson")
    m30 = model_for(L, law=law, taul_factor=tf, step=0.04, x_low=-30.0)
    entry["lnH_err_xlow30_vs_adaptive"] = float(np.max(np.abs(m30.lnH(FREQ) - lref)))
    s_ref = sigmas(Js[0.01], names, sets)
    for st in (0.16, 0.08, 0.04, 0.02):
        entry["tau_col_err_vs_step0.01"][str(st)] = maxrel(Js[st][:, 0], Js[0.01][:, 0])
        entry["sigma_rel_err_vs_step0.01"][str(st)] = sig_err(sigmas(Js[st], names, sets), s_ref)
    # FD Jacobian on the adaptive rule (nodes move between evaluations; h = 1e-2 to stay above its noise)
    if L == 500:
        Ja = jacobian(mref, FREQ, names, h=1e-2, scheme="central")
        entry["sigma_adaptive_fd_h1e-2_rel_err_vs_step0.01"] = sig_err(sigmas(Ja, names, sets), s_ref)
    entry["sigma_step0.01"] = s_ref
    q[key] = entry
    print("quad", key, entry["lnH_err_vs_adaptive"], f"{time.perf_counter() - T0:.0f}s", flush=True)
res["quadrature"] = q

# ---------------------------------------------------------------- 2. FD step
fd = {}
for L, law, tf in ((500, "cattaneo", 0.0), (30, "cattaneo", 0.0), (500, "gk", 2.0), (100, "gk", 0.5)):
    key = f"{law}_L{L}" + (f"_taul{tf}" if law == "gk" else "")
    sets = SETS_GK if law == "gk" else SETS_CAT
    names = sets["N5"]
    m = model_for(L, law=law, taul_factor=tf)
    if law == "cattaneo":
        Jex = Indep(base_from_model(m), FREQ).jacobian(names)
        ref_name = "autodiff (independent model)"
    else:
        Jex = jacobian(m, FREQ, names, h=1e-3, scheme="richardson")
        ref_name = "twin Richardson h=1e-3"
    s_ex = sigmas(Jex, names, sets)
    e = {"reference": ref_name, "col_err": {}, "sigma_rel_err": {}}
    for h in (1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-5):
        J = jacobian(m, FREQ, names, h=h, scheme="central")
        e["col_err"][f"central_{h:g}"] = {n: maxrel(J[:, j], Jex[:, j]) for j, n in enumerate(names)}
        e["sigma_rel_err"][f"central_{h:g}"] = sig_err(sigmas(J, names, sets), s_ex)
    for h in (2e-2, 2e-3):
        J = jacobian(m, FREQ, names, h=h, scheme="richardson")
        e["col_err"][f"richardson_{h:g}"] = {n: maxrel(J[:, j], Jex[:, j]) for j, n in enumerate(names)}
        e["sigma_rel_err"][f"richardson_{h:g}"] = sig_err(sigmas(J, names, sets), s_ex)
    fd[key] = e
    print("fd", key, {k: v["ln_tauR"] for k, v in e["col_err"].items()}, f"{time.perf_counter() - T0:.0f}s", flush=True)
res["fd_step"] = fd

# ---------------------------------------------------------------- 3. 1D limit and independent quadrupole
oned = {}
fw = np.geomspace(1e4, 1e10, 61)
for L in (30, 100, 300, 500):
    for law, tf in (("cattaneo", 0.0), ("gk", 0.5), ("gk", 2.0)):
        key = f"{law}_L{L}" + (f"_taul{tf}" if law == "gk" else "")
        m1 = model_for(L, law=law, taul_factor=tf, geometry="1d")
        ind = Indep(base_from_model(m1), fw, geometry="1d", law=law)
        d = {"lnH_max_abs_diff_10k_10G": float(np.max(np.abs(m1.lnH(fw) - ind.lnH())))}
        sets = SETS_GK if law == "gk" else SETS_CAT
        names = [n for n in sets["N4"]]
        Jt = jacobian(m1, FREQ, names, h=2e-3, scheme="richardson")
        nm_ind = [("taul_lin" if n == "ln_taul" else n) for n in names]
        Ji = Indep(base_from_model(m1), FREQ, geometry="1d", law=law).jacobian(nm_ind)
        if law == "gk":   # d/d ln tau_l = tau_l * d/d tau_l = (tau_l/tau_R) * d/d taul_lin
            j = names.index("ln_taul")
            Ji[:, j] *= tf
        s_t = sigmas(Jt, names, {k: v for k, v in sets.items() if k != "N5"})
        s_i = sigmas(Ji, names, {k: v for k, v in sets.items() if k != "N5"})
        d["sigma_twin_fd"] = s_t
        d["sigma_indep_autodiff"] = s_i
        d["sigma_rel_diff"] = sig_err(s_t, s_i)
        oned[key] = d
        print("1d", key, d["lnH_max_abs_diff_10k_10G"], f"{time.perf_counter() - T0:.0f}s", flush=True)
# large-beam convergence of the twin Gaussian average to the 1D limit (Cattaneo and GK, 500 nm)
conv = {}
fb = np.array([1e5, 1e6, 1e7, 1e8, 2e8])
for law, tf in (("cattaneo", 0.0), ("gk", 2.0)):
    errs = []
    radii = (30e-6, 100e-6, 300e-6, 1e-3, 3e-3)
    for w in radii:
        mb = model_for(500, law=law, taul_factor=tf, wp=w, wr=0.8 * w)
        m1 = model_for(500, law=law, taul_factor=tf, wp=w, wr=0.8 * w, geometry="1d")
        errs.append(float(np.max(np.abs(mb.lnH(fb) - m1.lnH(fb)))))
    errs = np.array(errs)
    order = np.log(errs[:-1] / errs[1:]) / np.log(np.array(radii[1:]) / np.array(radii[:-1]))
    conv[f"{law}_L500"] = {"radii_wp": radii, "max_abs_dlnH": errs.tolist(), "observed_order": order.tolist()}
    print("large-beam", law, errs, order, flush=True)
oned["large_beam_convergence"] = conv
res["one_d"] = oned

# ---------------------------------------------------------------- 4. independent finite-beam model
ib = {}
for L in (30, 100, 300, 500):
    for lm in ("rta", "iso"):
        m = model_for(L, lamr_mode=lm)
        names = SETS_CAT["N5"]
        ind = Indep(base_from_model(m), FREQ)
        Jt = jacobian(m, FREQ, names, h=2e-3, scheme="richardson")
        Ji = ind.jacobian(names)
        s_t, s_i = sigmas(Jt, names, SETS_CAT, SIG_PHI_DEG), sigmas(Ji, names, SETS_CAT, SIG_PHI_DEG)
        ib[f"L{L}_lamr_{lm}"] = {"lnH_max_abs_diff": float(np.max(np.abs(m.lnH(FREQ) - ind.lnH()))),
                                 "col_err": {n: maxrel(Jt[:, j], Ji[:, j]) for j, n in enumerate(names)},
                                 "sigma_twin": s_t, "sigma_indep": s_i, "sigma_rel_diff": sig_err(s_t, s_i)}
        print("indep-beam", L, lm, ib[f"L{L}_lamr_{lm}"]["lnH_max_abs_diff"],
              max(max(v.values()) for v in ib[f"L{L}_lamr_{lm}"]["sigma_rel_diff"].values()), flush=True)
res["indep_beam"] = ib

# ---------------------------------------------------------------- 5. frequency grid
fg = {}
for L in (100, 500):
    m = model_for(L)
    names = SETS_CAT["N5"]
    rows = {}
    for N in (20, 40, 80, 160, 320):
        f = np.geomspace(1e4, 2e8, N)
        J = Indep(base_from_model(m), f).jacobian(names)
        s = sigmas(J, names, SETS_CAT, (0.01, 0.1))
        rows[str(N)] = {k: {sg: v * np.sqrt(N) for sg, v in d.items()} for k, d in s.items()}
    # twin spot check at N = 160
    f = np.geomspace(1e4, 2e8, 160)
    Jt = jacobian(m, f, names, h=2e-3, scheme="richardson")
    rows["160_twin"] = {k: {sg: v * np.sqrt(160) for sg, v in d.items()}
                        for k, d in sigmas(Jt, names, SETS_CAT, (0.01, 0.1)).items()}
    fg[f"L{L}"] = rows
    print("grid", L, {N: rows[N]["N2"]["0.1"] for N in rows}, flush=True)
res["frequency_grid_sigma_sqrtN"] = fg

# ---------------------------------------------------------------- 6. in-plane relaxation time
ip = {}
b1 = json.loads((OUT / "a1_base_values.json").read_text())
for L in (30, 100, 300, 500):
    m = model_for(L)
    taux = b1["rta"]["m646435"]["tau_x"][str(L)] * 1e-12
    names = SETS_CAT["N5"]
    ind0 = Indep(base_from_model(m), FREQ)
    ind1 = Indep(base_from_model(m, taur_ip=taux), FREQ)
    J0, J1 = ind0.jacobian(names), ind1.jacobian(names)
    # z-only derivative: d/d ln tau_R minus the in-plane part
    Jz = J1.copy()
    Jz[:, 0] = J1[:, 0] - ind1.jacobian(["ln_taur_inplane"])[:, 0]
    s0, s1, sz = (sigmas(J, names, SETS_CAT, (0.01, 0.1)) for J in (J0, J1, Jz))
    ip[f"L{L}"] = {"tau_x_ps": taux * 1e12, "max_abs_dlnH": float(np.max(np.abs(ind1.lnH() - ind0.lnH()))),
                   "sigma_scalar_tau": s0, "sigma_tau_ip_eq_taux_dlnTauBoth": s1, "sigma_tau_ip_eq_taux_dlnTauZ": sz}
    print("in-plane", L, ip[f"L{L}"]["max_abs_dlnH"], s0["N2"], sz["N2"], flush=True)
res["inplane_tau"] = ip

res["seconds"] = time.perf_counter() - T0
(OUT / "a2_numerics.json").write_text(json.dumps(jsonable(res), indent=1))
print("done", res["seconds"])
