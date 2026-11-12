"""R-A5: objection A-2 -- the FDTR-visible combination delta = tau_R - tau_l and its static-data range.

1. Availability of the normal/resistive split: the stored film results (results_c/film3/*.json, film/*.json) contain
   tau, kappa and tauRTA only; the event lists carry a per-event normal/umklapp flag (ev_umk, BZ-representative test
   of aln_events.py, 'reporting only'), and the witnesses' amplitudes are stored, so the split is computable without
   new optimisation (R-A1).
2. Per witness (REF + 12 three-mesh optima), 7x7x5 and 9x9x6, L in {30, 100, 300, 500, 1000} nm:
   (a) kinetic closure: tau_l = exact k^2 coefficient of the regularised film kernel (ra_kinetic.py; reduces to
       4 tau_c / 5 for the gray isotropic Callaway operator, R-A0), delta_kin = tau_mem - tau_l;
   (b) gray closure of note 14 with a gray normal time from the N events: tau_l = (4/5)(1/tau_N + 1/tau_R)^-1,
       tau_R = tau_mem, tau_N in {G1 (kappa-weighted), G2 (capacity-weighted), RQ (Rayleigh quotient on Q V x)}.
   Ranges over the witnesses (inner bounds, robust over the two meshes as in R-A3), relative and in units of tau_R,0.
3. FDTR side, 1D large-beam limit, local interface-temperature convention (indep.impedance_1d_gk; Assumption M of B):
   sigma(delta) [ps] at the kinetic GK base point against sigma(tau_R) [ps] of the Cattaneo film (A); the margin
   m = sigma / (half-width), N1-N3, 0.01/0.03/0.1 deg, R1/R2 at the note-16 values and at the favourable corner.
   Linearity of H in delta up to the kinetic delta_0, and the Fourier-film test (A7 analogue): noiseless profile
   Delta chi^2 between the kinetic GK film and a Fourier film with the nuisances refitted (Gauss-Newton).
4. Gray-closure bound from (w_R, w_N, x = tau_N / tau_R): |Delta delta| / tau_R <= [1 - 0.8 x^2/(1+x)^2] w_R
   + [0.8 x/(1+x)^2] w_N.
Output: out/ra5_delta.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ["XLA_FLAGS"] = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE.parent / "A"))
from indep import impedance_1d_gk  # noqa: E402
import common  # noqa: E402

WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
CLASSES = {"orb0.1": WIT[:6], "all": WIT}
LS = (30, 100, 300, 500, 1000)
key = lambda L: f"{float(L):.6g}"  # noqa: E731
ra = {m: json.loads((OUT / f"ra1_film_{m}.json").read_text()) for m in ("m775", "m996")}
rc = json.loads((OUT / "ra1c_ref_rta.json").read_text())
r2 = json.loads((OUT / "ra2_rta_film_c.json").read_text())
a1 = json.loads((HERE.parent / "A" / "out" / "a1_base_values.json").read_text())["base"]
res = {"availability": {
    "stored_film_results_keys": ["kappa_x", "tau_x", "tauRTA_x", "kappa_z", "tau_z", "tauRTA_z"],
    "event_list_flag": "ev_umk (0 normal, 1 umklapp) per event; BZ-representative test, reporting only",
    "m775": {k: ra["m775"][k] for k in ("n_events", "n_normal", "n_event_orbits", "n_mixed_NU_orbits",
                                         "REF_rate_share_normal", "REF_rate_share_mixed")},
    "m996": {k: ra["m996"][k] for k in ("n_events", "n_normal", "n_event_orbits", "n_mixed_NU_orbits",
                                         "REF_rate_share_normal", "REF_rate_share_mixed")},
    "gamma_g": {"m775": ra["m775"]["gamma_g"], "m996": ra["m996"]["gamma_g"], "rta64": r2["gamma_g_RTA64"],
                "gray_isotropic": 0.8}}}


def gray_delta(rec, which):
    tN = rec[f"tauN_{which}"]
    tc = 1.0 / (1.0 / tN + 1.0 / rec["tau_z"])
    return rec["tau_z"] - 0.8 * tc


def q(m, W, L, name):
    rec = ra[m][W][key(L)]
    if name == "tau":
        return rec["tau_z"]
    if name == "taul":
        return rec["taul_kin"]
    if name == "delta":
        return rec["delta_kin"]
    if name.startswith("tauN_"):
        return rec[name]
    if name.startswith("dgray_"):
        return gray_delta(rec, name.split("_")[1])
    raise KeyError(name)


QUANT = ("tau", "taul", "delta", "tauN_G1", "tauN_G2", "tauN_RQ", "dgray_G1", "dgray_G2", "dgray_RQ")
tab = {}
for L in LS:
    row = {"REF": {}}
    for m in ra:
        rec = ra[m]["REF"][key(L)]
        row["REF"][m] = {nm: q(m, "REF", L, nm) for nm in QUANT}
        row["REF"][m]["tauc_RQ"] = rec["tauc_RQ"]
        row["REF"][m]["tauc_eff"] = rec["taul_kin"] / ra[m]["gamma_g"]
        row["REF"][m]["fracN"] = rec["fracN_kappa_weighted"]
        row["REF"][m]["taul_RTA_diag"] = rc[m][key(L)]["taul_RTA"]
        # normal/umklapp classification ambiguity (orbits mixing N and U labels): bounds on the gray normal time
        for G in ("G1", "G2"):
            lo, hi = rec.get(f"tauN_{G}_orbitN_any"), rec.get(f"tauN_{G}_orbitN_only")
            if lo is not None and hi is not None:
                tR = rec["tau_z"]
                dg = [tR - 0.8 / (1 / t + 1 / tR) for t in (lo, hi)]
                row["REF"][m][f"tauN_{G}_bounds"] = [lo, hi]
                row["REF"][m][f"dgray_{G}_bounds"] = dg
    row["REF"]["rta64_taul"] = r2["film"][key(L)]["taul_RTA"]
    exc = {}
    for W in WIT:
        exc[W] = {}
        for nm in QUANT:
            x7 = q("m775", W, L, nm) / q("m775", "REF", L, nm) - 1
            x9 = q("m996", W, L, nm) / q("m996", "REF", L, nm) - 1
            exc[W][nm] = {"m775": x7, "m996": x9, "robust": (x7 if abs(x7) < abs(x9) else x9) if x7 * x9 > 0 else 0.0}
    row["excursions"] = exc
    cls = {}
    for cl, ws in CLASSES.items():
        cls[cl] = {}
        for nm in QUANT:
            up = max([exc[W][nm]["robust"] for W in ws] + [0.0])
            lo = max([-exc[W][nm]["robust"] for W in ws] + [0.0])
            cls[cl][nm] = {"upper": up, "lower": lo, "w_half": 0.5 * (up + lo)}
        # absolute half-widths in units of tau_R,0 (9x9x6 REF ratios delta0/tau0)
        R9 = row["REF"]["m996"]
        for nm in ("delta", "dgray_G1", "dgray_G2", "dgray_RQ"):
            cls[cl][nm]["abs_over_tauR0_m996"] = cls[cl][nm]["w_half"] * abs(R9[nm]) / R9["tau"]
        cls[cl]["rho_kin_over_wtau"] = cls[cl]["delta"]["abs_over_tauR0_m996"] / cls[cl]["tau"]["w_half"]
        cls[cl]["rho_grayG1_over_wtau"] = cls[cl]["dgray_G1"]["abs_over_tauR0_m996"] / cls[cl]["tau"]["w_half"]
        cls[cl]["rho_grayG2_over_wtau"] = cls[cl]["dgray_G2"]["abs_over_tauR0_m996"] / cls[cl]["tau"]["w_half"]
        cls[cl]["rho_grayRQ_over_wtau"] = cls[cl]["dgray_RQ"]["abs_over_tauR0_m996"] / cls[cl]["tau"]["w_half"]
    row["classes"] = cls
    # REF tau_l estimates for the FDTR base point
    t9, t7 = row["REF"]["m996"]["taul"], row["REF"]["m775"]["taul"]
    corr = r2["film"][key(L)]["taul_RTA"] / rc["m996"][key(L)]["taul_RTA"]
    row["taul0_estimates_ps"] = {"m775": t7, "m996": t9, "m996_meshcorr_rta64": t9 * corr, "corr_factor": corr}
    tab[key(L)] = row
res["table"] = tab

# ------------------------------------------------------------------ FDTR 1D Fisher (kinetic GK base point)
FREQ = common.FREQ
F = jnp.asarray(FREQ)


def lnH1d(th, bv):
    """th = (dtauR [ps], dtaul [ps], ln lamz, ln C, ln R1, ln R2, g, phi0, td_ns); bv = (tauR, taul, lamz, C, R1, R2, d)."""
    e = jnp.exp
    P = {"tauR": bv[0] + th[0] * 1e-12, "taul": bv[1] + th[1] * 1e-12, "lamz": bv[2] * e(th[2]),
         "C": bv[3] * e(th[3]), "R1": bv[4] * e(th[4]), "R2": bv[5] * e(th[5]), "dAlN": bv[6],
         "lamS": common.SAP_LAM, "CS": common.SAP_C, "lamAu": common.AU_LAM, "CAu": common.AU_C, "dAu": common.AU_D}
    p = 2j * jnp.pi * F
    z = jnp.log(impedance_1d_gk(p, P)) + th[6] + 1j * (th[7] - 2 * jnp.pi * F * th[8] * 1e-9)
    return jnp.concatenate([z.real, z.imag])


jac1d = jax.jit(jax.jacfwd(lnH1d))
f1d = jax.jit(lnH1d)
SETS = {"N1": [6, 7, 8], "N2": [6, 7, 8, 4, 5], "N3": [6, 7, 8, 4, 5, 2, 3]}
SIGPHI = (0.01, 0.03, 0.1)


def whiten(J, sp):
    nf = FREQ.size
    return np.concatenate([J[:nf] / common.SIG_A, J[nf:] / np.deg2rad(sp)])


def proj_resid(v, Jn):
    sc = np.linalg.norm(Jn, axis=0)
    Q, _ = np.linalg.qr(Jn / sc)
    r = v - Q @ (Q.T @ v)
    return r - Q @ (Q.T @ r)


def sigma_dir(J, dvec, nuis, sp):
    """sigma along the parameter direction dvec (in the first two coordinates), nuisances profiled."""
    Jw = whiten(J, sp)
    v = Jw[:, :2] @ np.asarray(dvec, float)
    return 1.0 / np.linalg.norm(proj_resid(v, Jw[:, nuis]))


fd = {}
for L in (30, 100, 300, 500):
    b = a1[str(L)]
    tau0, lamz = b["tau_R_ps"] * 1e-12, b["lam_z"]
    est = tab[key(L)]["taul0_estimates_ps"]
    row = {"tauR0_ps": tau0 * 1e12, "taul0_ps": est}
    for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
        for ttag, taul0 in (("cattaneo", 0.0), ("kin_m996", est["m996"] * 1e-12),
                            ("kin_meshcorr", est["m996_meshcorr_rta64"] * 1e-12)):
            bv = jnp.asarray([tau0, taul0, lamz, common.C_ALN, R1, R2, L * 1e-9])
            J = np.asarray(jac1d(jnp.zeros(9), bv))
            ent = {}
            for N, nuis in SETS.items():
                for sp in SIGPHI:
                    sR = sigma_dir(J, (1, 0), nuis, sp)          # tau_R at fixed tau_l  (= +delta)
                    sL = sigma_dir(J, (0, -1), nuis, sp)         # -tau_l at fixed tau_R (= +delta)
                    ent[f"{N}|{sp}"] = {"sigma_delta_ps_via_tauR": sR, "sigma_delta_ps_via_taul": sL}
            # linearity in delta from tau_l = 0 to taul0 (exact change vs first-order prediction)
            if taul0 > 0:
                z0 = np.asarray(f1d(jnp.zeros(9), bv))
                z1 = np.asarray(f1d(jnp.asarray([0.0, -taul0 * 1e12, 0, 0, 0, 0, 0, 0, 0]), bv))
                lin = J[:, 1] * (-taul0 * 1e12)
                ent["linearity_rel_err"] = float(np.max(np.abs((z1 - z0) - lin)) / np.max(np.abs(lin)))
            row[f"{Rtag}|{ttag}"] = ent
    fd[key(L)] = row
res["fdtr_1d"] = fd


# ------------------------------------------------------------------ margins
def margins(L, cl="all"):
    t = tab[key(L)]
    c = t["classes"][cl]
    out = {}
    for Rtag in ("note16", "favourable"):
        for N in SETS:
            for sp in SIGPHI:
                sC = fd[key(L)][f"{Rtag}|cattaneo"][f"{N}|{sp}"]["sigma_delta_ps_via_tauR"]
                sK = fd[key(L)][f"{Rtag}|kin_meshcorr"][f"{N}|{sp}"]["sigma_delta_ps_via_tauR"]
                tau0 = fd[key(L)]["tauR0_ps"]
                m_tau = sC / (c["tau"]["w_half"] * tau0)
                m_kin = sK / (c["delta"]["abs_over_tauR0_m996"] * tau0)
                m_g1 = sC / (c["dgray_G1"]["abs_over_tauR0_m996"] * tau0)
                out[f"{Rtag}|{N}|{sp}"] = {"m_tau_cattaneo_1d": m_tau, "m_delta_kinetic_1d": m_kin,
                                           "m_delta_grayG1_1d": m_g1, "kinetic_over_tau": m_kin / m_tau}
    return out


res["margins_1d"] = {cl: {key(L): margins(L, cl) for L in (30, 100, 300, 500)} for cl in CLASSES}


# ------------------------------------------------------------------ Fourier film against the kinetic GK film (A7 analogue)
def gn_profile(bv_true, bv_alt, nuis, sp, iters=30):
    """min over nuisances of the whitened misfit between data = true film and the alternative film (noiseless)."""
    y = np.asarray(f1d(jnp.zeros(9), bv_true))
    th = np.zeros(9)
    for _ in range(iters):
        z = np.asarray(f1d(jnp.asarray(th), bv_alt))
        J = np.asarray(jac1d(jnp.asarray(th), bv_alt))
        rw = whiten(y - z, sp)
        Jw = whiten(J, sp)[:, nuis]
        step, *_ = np.linalg.lstsq(Jw, rw, rcond=None)
        th[nuis] += step
        if np.max(np.abs(step)) < 1e-12:
            break
    z = np.asarray(f1d(jnp.asarray(th), bv_alt))
    return float(np.sum(whiten(y - z, sp) ** 2))


a7 = {}
for L in (100, 300, 500):
    b = a1[str(L)]
    tau0, lamz = b["tau_R_ps"] * 1e-12, b["lam_z"]
    tl = tab[key(L)]["taul0_estimates_ps"]["m996_meshcorr_rta64"] * 1e-12
    for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8)}.items():
        bv_k = jnp.asarray([tau0, tl, lamz, common.C_ALN, R1, R2, L * 1e-9])
        bv_c = jnp.asarray([tau0, 0.0, lamz, common.C_ALN, R1, R2, L * 1e-9])
        bv_f = jnp.asarray([0.0, 0.0, lamz, common.C_ALN, R1, R2, L * 1e-9])
        for N, nuis in SETS.items():
            for sp in (0.01, 0.1):
                a7[f"{L}|{Rtag}|{N}|{sp}"] = {"dchi2_kineticGK_vs_Fourier": gn_profile(bv_k, bv_f, nuis, sp),
                                              "dchi2_cattaneo_vs_Fourier": gn_profile(bv_c, bv_f, nuis, sp)}
res["fourier_test_1d"] = a7

# ------------------------------------------------------------------ gray-closure bound
gb = {}
for L in (30, 100, 300, 500):
    t = tab[key(L)]
    for which in ("G1", "G2", "RQ"):
        x = t["REF"]["m996"][f"tauN_{which}"] / t["REF"]["m996"]["tau"]
        cR = 1 - 0.8 * x * x / (1 + x) ** 2
        cN = 0.8 * x / (1 + x) ** 2
        wR = t["classes"]["all"]["tau"]["w_half"]
        wN = t["classes"]["all"][f"tauN_{which}"]["w_half"]
        gb[f"{L}|{which}"] = {"x_tauN_over_tauR": x, "coef_R": cR, "coef_N": cN, "w_R": wR, "w_N_witnesses": wN,
                              "bound_abs_over_tauR0": cR * wR + cN * wN, "bound_over_wR": (cR * wR + cN * wN) / wR}
res["gray_bound"] = gb
(OUT / "ra5_delta.json").write_text(json.dumps(res, indent=1))

# ------------------------------------------------------------------ print
for L in LS:
    t = tab[key(L)]
    R9 = t["REF"]["m996"]
    print(f"L={L}: tau {R9['tau']:.3f} taul_kin {R9['taul']:.2f} (m775 {t['REF']['m775']['taul']:.2f}, "
          f"meshcorr {t['taul0_estimates_ps']['m996_meshcorr_rta64']:.2f}) delta {R9['delta']:.2f} "
          f"tauN G1/G2/RQ {R9['tauN_G1']:.2f}/{R9['tauN_G2']:.2f}/{R9['tauN_RQ']:.2f} fracN {R9['fracN']:.3f}")
    for cl in CLASSES:
        c = t["classes"][cl]
        print(f"   {cl}: w_tau {c['tau']['w_half']:.4f} w_taul {c['taul']['w_half']:.4f} w_delta {c['delta']['w_half']:.4f}"
              f" |Ddelta|/tau0 {c['delta']['abs_over_tauR0_m996']:.4f} (x{c['rho_kin_over_wtau']:.2f} of w_tau)"
              f" gray G1 x{c['rho_grayG1_over_wtau']:.3f} G2 x{c['rho_grayG2_over_wtau']:.3f} w_tauN(G1) {c['tauN_G1']['w_half']:.4f}")
for L in (30, 100, 300, 500):
    for k, v in res["margins_1d"]["all"][key(L)].items():
        if k.startswith("note16") and ("N2" in k or "N3" in k):
            print(L, k, {kk: round(vv, 3) for kk, vv in v.items()})
print(json.dumps(a7, indent=0))
print(json.dumps(gb, indent=0))
