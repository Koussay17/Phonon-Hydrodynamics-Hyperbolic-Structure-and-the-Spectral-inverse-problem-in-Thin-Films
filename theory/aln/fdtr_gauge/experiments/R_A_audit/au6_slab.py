"""AU6: kinetic slab with diffuse black walls (no Matthiessen regulator) against the regularised-kernel GK film.

Model (RTA with energy conservation, one local temperature, bulk 3ph lifetimes of the converged phono3py mesh,
frequency domain e^{pt}): for each mode class i (capacity c_i per direction, lifetime tau_i, |v_i| = |v_z|),
    +-|v_i| dT_i^+-/dz = -(p + 1/tau_i) T_i^+- + T_loc / tau_i,   0 < z < d,
    T_loc = [sum_i (c_i/tau_i)(T_i^+ + T_i^-)] / W(p),  W = 2 sum_i c_i/tau_i + 2 sum_loc c_0 p/(1 + p tau_0)
(classes with v_z = 0 are local: T_0 = T_loc/(1 + p tau_0)); walls black and diffuse: T_i^+(0) = T_A, T_i^-(d) = T_B.
Exact modal solution: S = T^+ + T^- obeys S'' = M S, M = diag(g^2/v^2) - diag(2g/(tau v^2)) 1 beta^T, g = p + 1/tau,
beta_j = (c_j/tau_j)/W; S = Phi [a e^{-lam z} + b e^{-lam (d-z)}], lam = sqrt(eig M) (Re > 0), D = -(|v|/g) S',
2N boundary conditions solved in even/odd form. Two-port (q(0), q(d)) = Y (T_A, T_B); checks: energy balance
q(0) - q(d) = p U, reciprocity Y12 = Y21, symmetry Y11 = -Y22, and the large-d limit G d -> lam_bulk.
Mode classes: 64x64x35 modes binned in (log tau, log |v_z|) (moving) and in log tau (v_z = 0), representative tau =
capacity mean, |v| preserving sum c v^2 tau; the regularised film parameters (lam_M, tau_R,M, tau_l,M; regulator
|v_z|/d, V1 completion) are computed on the same classes, so slab and GK film describe the same material.
FDTR (1D, au5 stack: Au 80 nm / R1 / film / R2 / sapphire), port temperatures = wall emission temperatures (R1, R2
are wall-emission interface resistances) and, as a convention check, the energy temperature at the walls.
Fits of the slab data by the GK film (au5.lnH, local convention): (a) tau_l fixed at tau_l,M (the claim's model);
(b) tau_l free; (c) Cattaneo (tau_l = 0); (d) Fourier; nuisance sets N1/N2/N3 (N2/N3 absorb the static boundary
resistance of the slab). Transfer test: tau_3ph x (1 +- 0.1) in both slab and GK film: fitted delta change against
the change of delta_M.
usage: python -B au6_slab.py [nbins]      Output: out/au6_slab[_nb].json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au2_rta as a2  # noqa: E402
import au5_fisher1d as a5  # noqa: E402

KB = 1.380649e-23
VCELL = 41.81335463721578e-30          # m^3 (events_m775 'volume', Angstrom^3)
H_OVER_KB = 47.99243073366221


def classes(md, nb, eta=0.0, vfloor=1.0):
    """Binned mode classes. md from au2.load_modes (cap in arbitrary units, tau3 [s], vz [m/s])."""
    Nq = float(np.prod(md["mesh"]))
    cap = md["cap"] * KB / (Nq * VCELL)          # J/(m^3 K), total over both directions
    tau = md["tau3"] * (1 + eta)
    v = md["vz"]
    mov = v > vfloor
    lt, lv = np.log10(tau), np.log10(np.maximum(v, vfloor))
    et = np.linspace(lt.min() - 1e-9, lt.max() + 1e-9, nb + 1)
    ev = np.linspace(lv[mov].min() - 1e-9, lv[mov].max() + 1e-9, nb + 1)
    it = np.clip(np.searchsorted(et, lt) - 1, 0, nb - 1)
    iv = np.clip(np.searchsorted(ev, lv) - 1, 0, nb - 1)
    k = it * nb + iv
    C = np.bincount(k[mov], weights=cap[mov], minlength=nb * nb)
    Ct = np.bincount(k[mov], weights=(cap * tau)[mov], minlength=nb * nb)
    Cv2t = np.bincount(k[mov], weights=(cap * v * v * tau)[mov], minlength=nb * nb)
    ok = C > 0
    c_m = C[ok] / 2                              # per direction
    t_m = Ct[ok] / C[ok]
    v_m = np.sqrt(Cv2t[ok] / Ct[ok])
    # local classes (v_z = 0)
    loc = ~mov
    if loc.any():
        k0 = it[loc]
        C0 = np.bincount(k0, weights=cap[loc], minlength=nb)
        Ct0 = np.bincount(k0, weights=(cap * tau)[loc], minlength=nb)
        ok0 = C0 > 0
        c_0, t_0 = C0[ok0], Ct0[ok0] / C0[ok0]           # total capacity of the local class
    else:
        c_0, t_0 = np.zeros(0), np.zeros(0)
    return {"c": c_m, "tau": t_m, "v": v_m, "c0": c_0, "tau0": t_0, "C_total": float(cap.sum()),
            "lam_bulk_exact": float((cap * v * v * tau).sum()), "n_classes": int(ok.sum()), "n_local": int(c_0.size)}


def matthiessen(cl, d):
    """Regularised film parameters on the classes: lam_M, tau_R,M, tau_l,M (V1 completion), C."""
    c2 = np.concatenate([2 * cl["c"], cl["c0"]])          # total capacity per class
    v = np.concatenate([cl["v"], np.zeros(cl["c0"].size)])
    t = np.concatenate([cl["tau"], cl["tau0"]])
    te = 1.0 / (1.0 / t + v / d)
    S = c2.sum()
    lam = float((c2 * v * v * te).sum())
    s = lam / S
    tauR = float((c2 * v * v * te ** 2).sum() / lam)
    taul = float((c2 * te * (v * v * te - s) ** 2).sum() / (S * s * s))
    return {"lam": lam, "tauR": tauR, "taul": taul, "C": float(S), "delta": tauR - taul}


def slab_twoport(cl, d, p):
    """Y such that (q(0), q(d)) = Y (T_A, T_B); also energy-balance and reciprocity residuals.
    Lengths scaled by d (z = d zeta)."""
    c, tau, v = cl["c"], cl["tau"], cl["v"]
    c0, t0 = cl["c0"], cl["tau0"]
    g = p + 1.0 / tau
    # c: capacity per direction of a moving class; c0: total capacity of a local (v_z = 0) class
    W = 2 * np.sum(c / tau) + np.sum(c0 * p / (1 + p * t0))
    beta = (c / tau) / W
    alpha = (g / v) ** 2 * d * d
    u = (2 * g / (tau * v * v)) * d * d
    Mm = np.diag(alpha) - np.outer(u, beta)
    mu, Phi = sla.eig(Mm, check_finite=False)
    # polish the small eigenvalues on the secular equation in cancellation-free form:
    # f(mu) = F0 - mu sum(beta u / (alpha (alpha - mu))),  F0 = 1 - sum(beta u/alpha) = p [sum 2c/(1+p tau) + sum c0/(1+p tau0)] / W
    F0 = p * (np.sum(2 * c / (1 + p * tau)) + np.sum(c0 / (1 + p * t0))) / W
    bu = beta * u
    Ctot = 2 * c.sum() + c0.sum()
    mu_d0 = p * d * d * Ctot / np.sum(2 * c * v * v * tau)       # first-order diffusive eigenvalue (scaled)
    m_ = int(np.argmin(np.abs(mu - mu_d0)))
    if abs(mu[m_]) < 1.0:
        x_ = mu[m_]
        for _ in range(30):
            f_ = F0 - x_ * np.sum(bu / (alpha * (alpha - x_)))
            fp_ = -np.sum(bu / (alpha - x_) ** 2)
            dx = f_ / fp_
            x_ = x_ - dx
            if abs(dx) < 1e-15 * max(abs(x_), 1e-300):
                break
        if abs(x_ - mu[m_]) < 1e-3 * abs(mu[m_]) + 1e-14:
            mu[m_] = x_
            ph = u / (alpha - x_)                               # exact eigenvector of diag(alpha) - u beta^T
            Phi[:, m_] = ph / np.linalg.norm(ph)
    lam = np.sqrt(mu.astype(complex))
    lam = np.where(lam.real < 0, -lam, lam)
    E = np.exp(-lam)                                   # e^{-lam d} in scaled units
    Vg = (v / g) / d                                   # |v|/g in units of d
    PhiL = Phi * lam[None, :]
    Pm = 0.5 * (Phi + Vg[:, None] * PhiL)
    Mn = 0.5 * (Phi - Vg[:, None] * PhiL)
    N = c.size
    one = np.ones(N, complex)
    A_plus = Pm + Mn * E[None, :]
    A_minus = Pm - Mn * E[None, :]
    luP = sla.lu_factor(A_plus, check_finite=False)
    luM = sla.lu_factor(A_minus, check_finite=False)
    wq = c * v * v / g                                 # q = - sum wq S'(z)  (S' in physical units = S'_scaled / d)

    def solve(TA, TB):
        sp_ = sla.lu_solve(luP, (TA + TB) * one, check_finite=False)
        sm_ = sla.lu_solve(luM, (TA - TB) * one, check_finite=False)
        a = 0.5 * (sp_ + sm_)
        b = 0.5 * (sp_ - sm_)
        dS0 = Phi @ (-lam * a + lam * E * b) / d
        dSd = Phi @ (-lam * E * a + lam * b) / d
        q0 = -wq @ dS0
        qd = -wq @ dSd
        intS = Phi @ ((a + b) * (1 - E) / lam) * d
        U = c @ intS + np.sum(c0 / (1 + p * t0)) * (beta @ intS)
        S0 = Phi @ (a + E * b)
        Sd = Phi @ (E * a + b)
        Tl0, Tld = beta @ S0, beta @ Sd
        Ctot = 2 * c.sum() + c0.sum()
        th0 = (c @ S0 + np.sum(c0 / (1 + p * t0)) * Tl0) / Ctot
        thd = (c @ Sd + np.sum(c0 / (1 + p * t0)) * Tld) / Ctot
        return q0, qd, U, th0, thd
    q0A, qdA, UA, th0A, thdA = solve(1.0, 0.0)
    q0B, qdB, UB, th0B, thdB = solve(0.0, 1.0)
    Y = np.array([[q0A, q0B], [qdA, qdB]])
    Th = np.array([[th0A, th0B], [thdA, thdB]])        # energy temperatures at the walls
    ebal = max(abs((q0A - qdA) - p * UA), abs((q0B - qdB) - p * UB)) / max(abs(q0A), 1e-300)
    recip = abs(Y[0, 1] + Y[1, 0]) / abs(Y[0, 1])     # with q(d) the flux leaving at the bottom: Y21 = -Y12
    symm = abs(Y[0, 0] + Y[1, 1]) / abs(Y[0, 0])
    return Y, Th, {"energy_balance_rel": float(ebal), "reciprocity_rel": float(recip), "symmetry_rel": float(symm)}


def Z_stack_slab(Ys, Ths, p, R1, R2, port="wall"):
    """Surface impedance of Au / R1 / slab / R2 / sapphire. port = 'wall' (R1, R2 on the wall emission temperatures)
    or 'energy' (R1, R2 on the energy temperatures at the walls)."""
    Zsub = 1.0 / np.sqrt(a5.SAP_LAM * a5.SAP_C * p)
    out = []
    for Y, Th, pp in zip(Ys, Ths, p):
        Zs = 1.0 / np.sqrt(a5.SAP_LAM * a5.SAP_C * pp)
        if port == "wall":
            # bottom: T_B = (Zs + R2) q(d); q(d) = Y10 T_A + Y11 T_B
            Zb = Zs + R2
            qd_per_TA = Y[1, 0] / (1 - Y[1, 1] * Zb)
            q0_per_TA = Y[0, 0] + Y[0, 1] * Zb * qd_per_TA
            Zt = 1.0 / q0_per_TA
        else:
            # energy temperatures: theta_d = Zs q(d) + ... with R2 between theta_d and the substrate:
            # theta_d - T_sub = R2 q(d), T_sub = Zs q(d)  ->  theta_d = (Zs + R2) q(d); unknowns T_A, T_B
            Zb = Zs + R2
            # theta = Th @ (T_A, T_B), q = Y @ (T_A, T_B); condition theta_d = Zb q_d, then Z = theta_0 / q_0
            row = Th[1] - Zb * Y[1]
            # T_B = -row[0]/row[1] T_A
            r = -row[0] / row[1]
            q0 = Y[0, 0] + Y[0, 1] * r
            th0 = Th[0, 0] + Th[0, 1] * r
            Zt = th0 / q0
        Zt = Zt + R1
        out.append(a5.layer(Zt, a5.AU_LAM, a5.AU_C, a5.AU_D, pp))
    return np.array(out)


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    Ls = tuple(int(x) for x in sys.argv[2].split(",")) if len(sys.argv) > 2 else (300, 500)
    etas = (0.0,) if (len(sys.argv) > 3 and sys.argv[3] == "eta0") else (0.0, -0.1, 0.1)
    tag = f"{nb}" + ("" if len(sys.argv) <= 2 else "_" + "_".join(str(x) for x in Ls)) + ("_eta0" if etas == (0.0,) else "")
    t0 = time.perf_counter()
    md = a2.load_modes("m646435")
    p = 2j * np.pi * a5.FREQ
    res = {"nbins": nb}
    # large-d check of the slab conductance against the bulk conductivity (static limit, small p)
    cl = classes(md, nb)
    res["classes"] = {k: cl[k] for k in ("n_classes", "n_local", "C_total", "lam_bulk_exact")}
    lam_b_cls = float(np.sum(2 * cl["c"] * cl["v"] ** 2 * cl["tau"]))
    res["classes"]["lam_bulk_classes"] = lam_b_cls
    chk = {}
    for dd in (1e-6, 1e-5, 1e-4):
        Y, Th, diag = slab_twoport(cl, dd, 2j * np.pi * 1.0)
        G = -Y[0, 1].real                          # q(0) per unit T_B (negative) -> conductance
        chk[f"{dd:g}"] = {"G_d_over_lam_bulk": float(G * dd / lam_b_cls), **diag}
    res["large_d_check"] = chk
    print("classes", res["classes"], "large-d", {k: round(v["G_d_over_lam_bulk"], 4) for k, v in chk.items()},
          f"{time.perf_counter() - t0:.0f}s", flush=True)
    for L in Ls:
        d = L * 1e-9
        ent = {}
        variants = {}
        for eta in etas:
            cle = classes(md, nb, eta=eta)
            Mt = matthiessen(cle, d)
            Ys, Ths, diags = [], [], []
            for pp in p:
                Y, Th, dg = slab_twoport(cle, d, pp)
                Ys.append(Y), Ths.append(Th), diags.append(dg)
            Y0, _, _ = slab_twoport(cle, d, 2j * np.pi * 1e2)
            G0 = float(Y0[0, 0].real)
            variants[eta] = {"M": Mt, "Ys": Ys, "Ths": Ths, "G0": G0,
                             "diag_max": {k: max(x[k] for x in diags) for k in diags[0]}}
        v0 = variants[0.0]
        Mt = v0["M"]
        ent["matthiessen"] = {k: (vv * 1e12 if k in ("tauR", "taul", "delta") else vv) for k, vv in Mt.items()}
        ent["slab_static"] = {"G_slab": v0["G0"], "R_slab": 1 / v0["G0"], "R_matthiessen": d / Mt["lam"],
                              "R_slab_minus_R_M": 1 / v0["G0"] - d / Mt["lam"], "lam_slab_eff": v0["G0"] * d}
        ent["slab_checks_max"] = v0["diag_max"]
        cl0 = classes(md, nb)
        om = 2 * np.pi * np.array([1e3, 2e3])
        Yl = [slab_twoport(cl0, d, 1j * w_)[0] for w_ in om]
        a1 = [Y_[0, 0].imag / w_ for Y_, w_ in zip(Yl, om)]
        b1 = [Y_[0, 1].imag / w_ for Y_, w_ in zip(Yl, om)]
        a1x, b1x = (4 * a1[0] - a1[1]) / 3, (4 * b1[0] - b1[1]) / 3     # Richardson in omega^2 (error O(omega^2))
        Gs = v0["G0"]
        ent["two_port_first_order"] = {"a1": a1x, "b1": b1x, "a1_plus_b1_over_Cd2": (a1x + b1x) / (Mt["C"] * d / 2),
                                       "delta_2port_wall_ps": (Mt["C"] * d / 6 - (a1x - b1x)) / (2 * Gs) * 1e12,
                                       "delta_lumped_wall_resistance_only_ps":
                                       (Mt["C"] * d / 12) * (1 / Gs) * (1 - ((d / Mt["lam"]) * Gs) ** 2) * 1e12}
        print(L, "Matthiessen", {k: round(vv, 4) for k, vv in ent["matthiessen"].items()}, "static",
              {k: f"{vv:.4e}" for k, vv in ent["slab_static"].items()}, "checks", ent["slab_checks_max"], flush=True)
        for port in ("wall", "energy"):
            for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
                Z = Z_stack_slab(v0["Ys"], v0["Ths"], p, R1, R2, port=port)
                y = np.log(Z)
                base = {"tauR": Mt["tauR"], "taul": Mt["taul"], "lamz": Mt["lam"], "C": Mt["C"], "R1": R1, "R2": R2, "d": d}
                row = {}
                base_s = dict(base)
                base_s["lamz"] = v0["G0"] * d                      # the slab's own static conductance x d
                SL = {"N2": (a5.SETS["N2"], base), "N2s": (a5.SETS["N2"], base_s),
                      "N2lam": (a5.SETS["N2"] + [2], base), "N3": (a5.SETS["N3"], base)}
                for N, (nuis, base) in SL.items():
                    for sp in (0.01, 0.1):
                        chiA, thA = a5.gn_refit(y, base, nuis, sp)
                        chiB, thB = a5.gn_refit(y, base, nuis, sp, free=(1,))
                        bC = dict(base)
                        bC["taul"] = 0.0
                        chiC, thC = a5.gn_refit(y, bC, nuis, sp, free=(0,))
                        bF = dict(base)
                        bF["tauR"], bF["taul"] = 0.0, 0.0
                        chiF, thF = a5.gn_refit(y, bF, nuis, sp)
                        dfit = (Mt["tauR"] - Mt["taul"] - thB[1] * 1e-12) * 1e12
                        tr = {}
                        for eta in [x for x in etas if x != 0.0]:
                            ve = variants[eta]
                            ye = np.log(Z_stack_slab(ve["Ys"], ve["Ths"], p, R1, R2, port=port))
                            chie, the = a5.gn_refit(ye, base, nuis, sp, free=(1,))
                            dfit_e = (Mt["tauR"] - Mt["taul"] - the[1] * 1e-12) * 1e12
                            dM = (ve["M"]["delta"] - Mt["delta"]) * 1e12
                            # Delta chi^2 between the two slabs (eta vs 0) with nuisances refitted, film model = GK
                            # with tau_l free (an experimenter's model): distinguishability of the two slabs
                            tr[str(eta)] = {"delta_M_change_ps": dM, "fitted_delta_change_ps": dfit_e - dfit,
                                            "transfer": (dfit_e - dfit) / dM if dM != 0 else None}
                        row[f"{N}|{sp}"] = {"dchi2_gk_taulM": chiA, "dchi2_gk_taul_free": chiB, "delta_fit_ps": dfit,
                                            "delta_M_ps": Mt["delta"] * 1e12, "dchi2_cattaneo_tauR_free": chiC,
                                            "tauR_fit_cattaneo_ps": (Mt["tauR"] + thC[0] * 1e-12) * 1e12,
                                            "dchi2_fourier": chiF, "lnR1_fit": thB[4], "lnR2_fit": thB[5],
                                            "ln_lamz_fit": thB[2], "transfer": tr}
                ent[f"{port}|{Rtag}"] = row
                print(L, port, Rtag, {k: {"chi(taulM)": round(vv["dchi2_gk_taulM"], 2), "chi(free)": round(vv["dchi2_gk_taul_free"], 3),
                                          "d_fit": round(vv["delta_fit_ps"], 1), "d_M": round(vv["delta_M_ps"], 1),
                                          "chiF": round(vv["dchi2_fourier"], 2),
                                          "transfer": [None if vv["transfer"][x]["transfer"] is None else round(vv["transfer"][x]["transfer"], 3) for x in vv["transfer"]]}
                                      for k, vv in row.items() if k.split("|")[1] == "0.01"},
                      f"{time.perf_counter() - t0:.0f}s", flush=True)
        res[str(L)] = ent
        (OUT / f"au6_slab_{tag}.json").write_text(json.dumps(res, indent=1, default=float))
    res["seconds"] = time.perf_counter() - t0
    (OUT / f"au6_slab_{tag}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
