"""AU2: independent recomputation of the RTA kinetic nonlocal time tau_l(L) on two converged phono3py meshes, with
three routes, the even-sector completion variants, and the k-truncation test on the exact diffusive root.

Mode data: fc3-memory E27 RTA chunks (read-only) for 48x48x26 and 64x64x35; capacities at 300 K, tau = 1/(4 pi gamma),
|v_z| from gv (THz*A -> m/s); each irreducible mode is split into two directed copies (+|v_z|, -|v_z|) with half the
capacity (only even functions of v_z enter). Film regulator 1/tau_b = |v_z| / L (c = 1).

Symmetrised variables e_i = sqrt(c_i), b_i = v_i e_i. Operators (all symmetric PSD with C e = 0):
  V1 'combined'  : D - (D e)(D e)^T / (e^T D e),  D = D_3 + D_b       (campaign RTA route, one local temperature)
  V2 'separate'  : D - sum_m (D_m e)(D_m e)^T / (e^T D_m e), m = 3ph, boundary (each process conserves energy)
  V3 'QDQ-wall'  : D_3 - (D_3 e)(D_3 e)^T / (e^T D_3 e) + Q D_b Q    (the event-route completion of the regulator)
  V4 'QDQ-all'   : Q D Q                                                (comparison only)
Routes: (i) closed form tau_R = |x|^2/K, tau_l = |e|^2 (u.C^-1 u)/K^2, u = Q V x, with C^-1 on e-perp through
Woodbury (C + gamma ehat ehat^T); for V1 also the explicit moment formula <tau (v^2 tau - s)^2>/s^2;
(ii) exact driven kernel kappa(k, p) = |e|^2 (b.A^-1 e)/(e.A^-1 e)/(-i k), A = p + i k V + C (Woodbury, O(n)),
first derivatives by Richardson; (iii) diffusive root s = k^2 of p|e|^2 + s kappa = 0 (Newton, continuation in
omega), lam_eff = -p |e|^2 / s, compared at FDTR frequencies with the first-order form 1 - delta p, the GK
resummation (1 + tau_l p)/(1 + tau_R p), and the first-order-in-k^2 form exact in p, [K(p) + p|e|^2 M(p)/K(p)]/K(0)
(the quantity of R-A1d). Also the k^4 coefficient of the static kernel against K L^4 (GK resummation: ratio 1).
Output: out/au2_rta.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
RTA = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/rta")
H_OVER_KB = 47.99243073366221        # h/k_B in K/THz
LS = (30, 100, 300, 500, 1000)
FREQ_TEST = np.array([1e6, 1e7, 2e7, 5e7, 1e8, 1.4e8, 2e8, 5e8, 1e9])


def load_modes(mesh):
    D0 = RTA / mesh
    meta = json.loads((D0 / "meta.json").read_text())
    fs = sorted(f for f in D0.glob("chunk_*.npz") if ".tmp" not in f.name)
    assert len(fs) == meta["n_chunks"]
    w, f, gam, gv = [], [], [], []
    for fp in fs:
        z = np.load(fp)
        w.append(z["w"]), f.append(z["f"]), gam.append(z["gamma"]), gv.append(z["gv"])
    w = np.concatenate(w).astype(float)
    f = np.concatenate(f)
    gam = np.concatenate(gam)
    gv = np.concatenate(gv)
    assert abs(w.sum() - np.prod(meta["mesh"])) < 1e-6
    ok = (gam > 0) & (f > 1e-3)
    x = H_OVER_KB * f[ok] / 300.0
    cap = np.broadcast_to(w[:, None], f.shape)[ok] * x * x * np.exp(x) / np.expm1(x) ** 2
    tau3 = 1.0 / (4 * np.pi * gam[ok]) * 1e-12            # s
    vz = np.abs(gv[..., 2][ok]) * 100.0                  # m/s
    return {"mesh": meta["mesh"], "cap": cap, "tau3": tau3, "vz": vz, "n_ir_modes": int(ok.sum())}


def directed(md):
    c = np.concatenate([md["cap"], md["cap"]]) / 2
    v = np.concatenate([md["vz"], -md["vz"]])
    t3 = np.concatenate([md["tau3"], md["tau3"]])
    return c, v, t3


class LowRankOp:
    """C = diag(d) + U S U^T (symmetric, U n x r, S r x r)."""

    def __init__(self, d, U, S):
        self.d, self.U, self.S = d, U, S

    def solve_shift(self, shift, rhs):
        """(diag(d + shift) + U S U^T)^-1 rhs (shift may be a complex vector), Woodbury."""
        Di = 1.0 / (self.d + shift)
        if self.U.shape[1] == 0:
            return Di[:, None] * rhs if rhs.ndim == 2 else Di * rhs
        R = rhs if rhs.ndim == 2 else rhs[:, None]
        DU = Di[:, None] * self.U
        cap = np.linalg.inv(self.S) + self.U.T @ DU
        Y = Di[:, None] * R - DU @ np.linalg.solve(cap, DU.T @ R)
        return Y if rhs.ndim == 2 else Y[:, 0]

    def matvec(self, x):
        return self.d * x + self.U @ (self.S @ (self.U.T @ x))


def build_variant(name, c, v, t3, L):
    e = np.sqrt(c)
    ehat = e / np.linalg.norm(e)
    d3 = 1.0 / t3
    db = np.abs(v) / (L * 1e-9)
    D = d3 + db
    if name == "V1":
        g = D * e
        U = g[:, None]
        S = np.array([[-1.0 / (e @ g)]])
    elif name == "V2":
        g3, gb = d3 * e, db * e
        U = np.stack([g3, gb], axis=1)
        S = np.diag([-1.0 / (e @ g3), -1.0 / (e @ gb)])
    elif name == "V3":
        g3 = d3 * e
        hb = db * ehat
        # Q D_b Q = D_b - ehat hb^T - hb ehat^T + (ehat.hb) ehat ehat^T ; plus -(g3 g3^T)/(e.g3)
        U = np.stack([g3, ehat, hb], axis=1)
        S = np.array([[-1.0 / (e @ g3), 0, 0], [0, ehat @ hb, -1.0], [0, -1.0, 0]])
    elif name == "V4":
        h = D * ehat
        U = np.stack([ehat, h], axis=1)
        S = np.array([[ehat @ h, -1.0], [-1.0, 0.0]])
    else:
        raise KeyError(name)
    op = LowRankOp(D, U, S)
    op.ie = 1 if name == "V3" else (0 if name == "V4" else None)     # column of ehat in U, if present
    return op, e, ehat


def regularised(op, ehat, gam):
    """C + gam ehat ehat^T as a LowRankOp (adds gam to the ehat entry of S if ehat is already a column of U)."""
    if op.ie is not None:
        S = op.S.copy()
        S[op.ie, op.ie] += gam
        return LowRankOp(op.d, op.U, S)
    r = op.S.shape[0]
    S = np.zeros((r + 1, r + 1))
    S[:r, :r] = op.S
    S[r, r] = gam
    return LowRankOp(op.d, np.concatenate([op.U, ehat[:, None]], axis=1), S)


def closed_form(op, e, ehat, v):
    """tau_R, tau_l via C + gamma ehat ehat^T (gamma-independence checked)."""
    out = {}
    for gam in (np.median(op.d), 10 * np.median(op.d)):
        reg = regularised(op, ehat, gam)
        b = v * e
        x = reg.solve_shift(0.0, b)
        K = b @ x
        tauR = x @ x / K
        u = v * x
        u = u - ehat * (ehat @ u)
        y = reg.solve_shift(0.0, u)
        M = u @ y
        out[gam] = (K, tauR, (e @ e) * M / K ** 2, M, abs(ehat @ x) / np.linalg.norm(x))
    (K, tauR, taul, M, ex), (K2, tauR2, taul2, M2, _) = out.values()
    return {"K": K, "tauR": tauR, "taul": taul, "M": M, "gamma_indep_rel": abs(taul2 / taul - 1),
            "ehat_dot_x_rel": ex}


def kernel(op, e, v, k, p):
    A_shift = p + 1j * k * v
    Ae = op.solve_shift(A_shift, e.astype(complex))
    num = (v * e) @ Ae
    den = e @ Ae
    return (e @ e) * num / den / (-1j * k)


def kernel_coefs(op, e, v, K_cf, tauR_cf, Lam_max):
    """Taylor coefficients of kappa in s = k^2 by a Cauchy integral on |s| = r (trapezoidal rule = FFT, spectrally
    accurate inside the convergence radius 1/Lam_max^2), at p = 0 and at small +-p for d/dp (Richardson in p)."""
    r = (0.5 / Lam_max) ** 2
    nth = 64
    th = 2 * np.pi * np.arange(nth) / nth
    s_pts = r * np.exp(1j * th)
    k_pts = np.sqrt(s_pts)

    def coefs(p):
        vals = np.array([kernel(op, e, v, kk, p) for kk in k_pts])
        c = np.fft.fft(vals) / nth
        return c[:5] / r ** np.arange(5)
    c = coefs(0.0)
    hp = 1e-3 / tauR_cf
    d1 = (coefs(hp)[0] - coefs(-hp)[0]) / (2 * hp)
    d2 = (coefs(hp / 2)[0] - coefs(-hp / 2)[0]) / hp
    dp = (4 * d2 - d1) / 3
    return {"K": c[0].real, "M": -c[1].real, "c4": c[2].real, "c6": -c[3].real, "imag_rel": float(abs(c[0].imag) / abs(c[0])),
            "tauR": -dp.real / c[0].real, "r_s": r}


def diffusive_root(op, e, v, K0, freqs):
    e2 = e @ e
    out = []
    s = None
    for f in freqs:
        p = 2j * np.pi * f
        if s is None:
            s = -p * e2 / K0
        else:
            s = s * (f / fprev)                       # continuation: s ~ p
        for it in range(100):
            k = np.sqrt(s)
            if k.imag > 0:
                k = -k
            F = None
            for _try in range(8):
                try:
                    F = p * e2 + s * kernel(op, e, v, k, p)
                    if np.isfinite(F):
                        break
                except np.linalg.LinAlgError:
                    pass
                s = s * (1 + 1e-7 * (1 + 1j))
                k = np.sqrt(s)
                if k.imag > 0:
                    k = -k
            ds = 1e-6 * s
            k2 = np.sqrt(s + ds)
            if k2.imag > 0:
                k2 = -k2
            F2 = p * e2 + (s + ds) * kernel(op, e, v, k2, p)
            step = F * ds / (F2 - F)
            s = s - step
            if abs(step) < 1e-13 * abs(s):
                break
        fprev = f
        out.append({"f": float(f), "lam_eff_rel": complex(-p * e2 / s / K0), "newton_iters": it + 1,
                    "k_times_1nm": complex(np.sqrt(s) * 1e-9)})
    return out


def first_order_in_k2_exact_in_p(op, e, ehat, v, K0, freqs):
    """[K(p) + p |e|^2 M(p)/K(p)]/K(0): K(p) = b.(p+C)^-1 b, M(p) = u(p).(p+C)^-1 u(p), u(p) = Q V (p+C)^-1 b."""
    reg = regularised(op, ehat, np.median(op.d))
    e2 = e @ e
    b = v * e
    out = []
    for f in freqs:
        p = 2j * np.pi * f
        # (p + C) on e-perp: the regulariser gamma ehat ehat^T becomes (p + gamma) on ehat; b and u are in e-perp
        x = reg.solve_shift(p, b.astype(complex))
        K = b @ x
        u = v * x
        u = u - ehat * (ehat @ u)
        y = reg.solve_shift(p, u)
        Mp = u @ y
        out.append(complex((K + p * e2 * Mp / K) / K0))
    return out


def main():
    res = {}
    for mesh in ("m484826", "m646435"):
        md = load_modes(mesh)
        c, v, t3 = directed(md)
        row = {"mesh": md["mesh"], "n_ir_modes": md["n_ir_modes"]}
        for L in LS:
            ent = {}
            # explicit V1 moment formula (independent of the operator algebra)
            tl = 1.0 / (1.0 / md["tau3"] + md["vz"] / (L * 1e-9))
            S = md["cap"].sum()
            s_ = (md["cap"] * md["vz"] ** 2 * tl).sum() / S
            taul_mom = (md["cap"] * tl * (md["vz"] ** 2 * tl - s_) ** 2).sum() / (S * s_ ** 2)
            tmem_mom = (md["cap"] * md["vz"] ** 2 * tl ** 2).sum() / (md["cap"] * md["vz"] ** 2 * tl).sum()
            mean_tau = (md["cap"] * tl).sum() / S
            ent["moment_V1"] = {"tau_mem_ps": tmem_mom * 1e12, "taul_ps": taul_mom * 1e12, "mean_tau_ps": mean_tau * 1e12,
                                "direct_term_ps": ((md["cap"] * md["vz"] ** 4 * tl ** 3).sum() / (S * s_ ** 2)) * 1e12}
            for var in ("V1", "V2", "V3", "V4"):
                op, e, ehat = build_variant(var, c, v, t3, L)
                cf = closed_form(op, e, ehat, v)
                ent[var] = {"tauR_ps": cf["tauR"] * 1e12, "taul_ps": cf["taul"] * 1e12,
                            "delta_ps": (cf["tauR"] - cf["taul"]) * 1e12, "gamma_indep_rel": cf["gamma_indep_rel"],
                            "ehat_dot_x_rel": cf["ehat_dot_x_rel"]}
                if var in ("V1", "V3") and L in (300, 500):
                    Lam_max = float(np.max(np.abs(v) / op.d))
                    kc = kernel_coefs(op, e, v, cf["K"], cf["tauR"], Lam_max)
                    L2 = cf["M"] / cf["K"]
                    ent[var]["kernel_route"] = {"rel_K": abs(kc["K"] / cf["K"] - 1), "rel_M": abs(kc["M"] / cf["M"] - 1),
                                                "rel_tauR": abs(kc["tauR"] / cf["tauR"] - 1),
                                                "c4_over_K_L4": kc["c4"] / (cf["K"] * L2 ** 2),
                                                "c6_over_K_L6": kc["c6"] / (cf["K"] * L2 ** 3),
                                                "L_GK_nm": float(np.sqrt(L2) * 1e9), "Lambda_max_nm": Lam_max * 1e9,
                                                "imag_rel": kc["imag_rel"]}
                    droot = diffusive_root(op, e, v, cf["K"], FREQ_TEST)
                    fo = first_order_in_k2_exact_in_p(op, e, ehat, v, cf["K"], FREQ_TEST)
                    cmp_ = []
                    for r_, fo_ in zip(droot, fo):
                        pp = 2j * np.pi * r_["f"]
                        exact = r_["lam_eff_rel"]
                        gk = (1 + cf["taul"] * pp) / (1 + cf["tauR"] * pp)
                        lin = 1 - (cf["tauR"] - cf["taul"]) * pp
                        deff = -(exact - 1) / pp * 1e12
                        cmp_.append({"f_Hz": r_["f"], "exact": [exact.real, exact.imag], "gk": [gk.real, gk.imag],
                                     "first_order": [lin.real, lin.imag], "k2_exact_in_p": [fo_.real, fo_.imag],
                                     "delta_eff_exact_ps": [deff.real, deff.imag],
                                     "rel_dev_exact_vs_gk": abs(exact / gk - 1),
                                     "rel_dev_exact_vs_k2_exact_in_p": abs(exact / fo_ - 1),
                                     "dev_of_nonFourier_part_vs_gk": abs((exact - 1) / (gk - 1) - 1),
                                     "abs_k_L_GK": abs(np.sqrt(-pp / (cf["K"] / (e @ e)))) * np.sqrt(L2),
                                     "newton_iters": r_["newton_iters"]})
                    ent[var]["dispersion_test"] = cmp_
            row[str(L)] = ent
            print(mesh, L, "V1 moment taul", round(ent["moment_V1"]["taul_ps"], 2), "tau_mem", round(ent["moment_V1"]["tau_mem_ps"], 3),
                  {var: (round(ent[var]["tauR_ps"], 3), round(ent[var]["taul_ps"], 2)) for var in ("V1", "V2", "V3", "V4")},
                  flush=True)
            if "kernel_route" in ent["V1"]:
                print("    kernel route V1:", {k_: (f"{v_:.2e}" if isinstance(v_, float) else v_) for k_, v_ in ent["V1"]["kernel_route"].items()})
                for r_ in ent["V1"]["dispersion_test"]:
                    print(f"    f={r_['f_Hz']:.1e} |k|L_GK={r_['abs_k_L_GK']:.3f} exact={r_['exact'][0]:.4f}{r_['exact'][1]:+.4f}i "
                          f"gk={r_['gk'][0]:.4f}{r_['gk'][1]:+.4f}i  dev(nonFourier part) {r_['dev_of_nonFourier_part_vs_gk']:.3f} "
                          f"dev vs k2-exact-in-p {r_['rel_dev_exact_vs_k2_exact_in_p']:.3f} delta_eff {r_['delta_eff_exact_ps'][0]:.1f}{r_['delta_eff_exact_ps'][1]:+.1f}i",
                          flush=True)
        res[mesh] = row
    (OUT / "au2_rta.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
