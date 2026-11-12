"""AU5: independent 1D (large-beam) FDTR Fisher check of sigma(delta) against the static-data range of delta.

Own implementation (NumPy, no JAX, no campaign module): impedance recursion Au(80 nm) / R1 / AlN film / R2 /
sapphire (semi-infinite), film lam_eff(p) = lam_z (1 + tau_l p)/(1 + tau_R p) (1D GK film of the twin, local
interface temperature), observable ln H = ln Z + g + i (phi0 - 2 pi f t_d), 80 log-spaced frequencies 10 kHz-200 MHz,
noise 1 % on ln|H| and sigma_phi on the phase. Thermal Jacobian columns by central differences with Richardson
(h, h/2), checked against h/4; instrument columns exact. Nuisance sets: N1 = {g, phi0, t_d}, N2 = N1 + {ln R1,
ln R2}, N3 = N2 + {ln lam_z, ln C}. sigma(delta) = 1/||(I - P_N) J_delta|| with J_delta = d lnH/d tau_R at fixed
tau_l (and, as a check, -d lnH/d tau_l at fixed tau_R).
Ranges of delta from au4 (robust 7x7x5/9x9x6 rule; also 9x9x6-only and the geometric held-out estimate).
Direct test: noiseless data of the REF film against each witness film (robust excursions of tau_R and tau_l applied
to the base point), nuisances refitted by damped Gauss-Newton (Levenberg-Marquardt); m_direct = 1/sqrt(max Delta chi^2).
The k-truncation test (exact diffusive-root kernel against the first-order GK film) is in au5b_ktrunc.py.
Output: out/au5_fisher1d.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
A1 = json.loads((HERE.parent / "A" / "out" / "a1_base_values.json").read_text())["base"]
AU4 = json.loads((OUT / "au4_delta_range.json").read_text())

AU_LAM, AU_C, AU_D = 150.0, 2.49e6, 80e-9
SAP_LAM, SAP_C = 35.0, 3.03e6
C_ALN = 2.41e6
FREQ = np.geomspace(1e4, 2e8, 80)
SIG_A = 0.01
PNAMES = ("tauR_ps", "taul_ps", "ln_lamz", "ln_C", "ln_R1", "ln_R2", "g", "phi0", "td_ns")
SETS = {"N1": [6, 7, 8], "N2": [6, 7, 8, 4, 5], "N3": [6, 7, 8, 4, 5, 2, 3]}
SETS_EXTRA = {"N2lam": [6, 7, 8, 4, 5, 2]}      # lam_z free, C fixed (not a campaign set)
SIGPHI = (0.01, 0.03, 0.1)


def layer(Z_below, lam, C, d, p):
    g = np.sqrt(C * p / lam)
    Y = lam * g
    t = np.tanh(g * d)
    return (Z_below + t / Y) / (1 + Y * Z_below * t)


def lnH(th, base, ratio=None):
    """th: increments in PNAMES order; base: dict with tauR, taul [s], lamz, C, R1, R2, d. ratio: optional complex
    lam_eff/lam_z over FREQ replacing the GK form (exact-kernel test)."""
    p = 2j * np.pi * FREQ
    tauR = base["tauR"] + th[0] * 1e-12
    taul = base["taul"] + th[1] * 1e-12
    lamz = base["lamz"] * np.exp(th[2])
    C = base["C"] * np.exp(th[3])
    R1 = base["R1"] * np.exp(th[4])
    R2 = base["R2"] * np.exp(th[5])
    r = (1 + taul * p) / (1 + tauR * p) if ratio is None else ratio
    Z = 1.0 / np.sqrt(SAP_LAM * SAP_C * p) + R2
    Z = layer(Z, lamz * r, C, base["d"], p) + R1
    Z = layer(Z, AU_LAM, AU_C, AU_D, p)
    return np.log(Z) + th[6] + 1j * (th[7] - 2 * np.pi * FREQ * th[8] * 1e-9)


def jac(base, ratio=None, h=(1e-3, 1e-3, 1e-5, 1e-5, 1e-5, 1e-5)):
    cols = []
    th0 = np.zeros(9)
    for i in range(6):
        def cd(hh):
            tp, tm = th0.copy(), th0.copy()
            tp[i] += hh
            tm[i] -= hh
            return (lnH(tp, base, ratio) - lnH(tm, base, ratio)) / (2 * hh)
        hi = h[i] * (max(base["tauR"] * 1e12, 1.0) if i == 0 else (max(base["taul"] * 1e12, 1.0) if i == 1 else 1.0))
        d1, d2 = cd(hi), cd(hi / 2)
        cols.append((4 * d2 - d1) / 3)
    cols.append(np.ones(FREQ.size, complex))
    cols.append(np.full(FREQ.size, 1j))
    cols.append(-1j * 2 * np.pi * FREQ * 1e-9)
    return np.stack(cols, axis=1)


def whiten(Jc, sp):
    return np.concatenate([Jc.real / SIG_A, Jc.imag / np.deg2rad(sp)], axis=0)


def sigma_dir(J, vec, nuis, sp):
    Jw = whiten(J, sp)
    v = Jw[:, :2] @ np.asarray(vec, float)
    N = Jw[:, nuis]
    N = N / np.linalg.norm(N, axis=0)
    Q, _ = np.linalg.qr(N)
    r = v - Q @ (Q.T @ v)
    r = r - Q @ (Q.T @ r)
    return 1.0 / np.linalg.norm(r)


def gn_refit(y, base_alt, nuis, sp, ratio=None, free=(), iters=40):
    """min over nuisances (and optional free target coordinates) of the whitened misfit to data y (complex lnH)."""
    th = np.zeros(9)
    idx = list(nuis) + list(free)

    def cost(t):
        z = lnH(t, base_alt, ratio)
        r = whiten((y - z)[:, None], sp)[:, 0]
        return float(r @ r) if np.all(np.isfinite(r)) else np.inf, r
    c0, rw = cost(th)
    mu = 1e-6
    for _ in range(iters):
        J = jac_at(th, base_alt, ratio)
        Jw = whiten(J, sp)[:, idx]
        sc = np.linalg.norm(Jw, axis=0)
        sc[sc == 0] = 1.0
        A = Jw / sc
        g_ = A.T @ rw
        H_ = A.T @ A
        accepted = False
        for _ls in range(30):
            step = np.linalg.solve(H_ + mu * np.diag(np.diag(H_) + 1e-12), g_) / sc
            # bound the log-parameter steps (positions 2..5) to one e-fold per iteration
            for jj, ii in enumerate(idx):
                if 2 <= ii <= 5:
                    step[jj] = np.clip(step[jj], -1.0, 1.0)
            tn = th.copy()
            tn[idx] += step
            cn, rn = cost(tn)
            if cn <= c0:
                th, c0, rw = tn, cn, rn
                mu = max(mu / 10, 1e-12)
                accepted = True
                break
            mu *= 10
        if not accepted or np.max(np.abs(step * sc)) < 1e-11:
            break
    return c0, th


def jac_at(th, base, ratio=None):
    """Jacobian at a shifted point (rebuild the base so that jac() differentiates about th)."""
    b2 = dict(base)
    b2["tauR"] = base["tauR"] + th[0] * 1e-12
    b2["taul"] = base["taul"] + th[1] * 1e-12
    b2["lamz"] = base["lamz"] * np.exp(th[2])
    b2["C"] = base["C"] * np.exp(th[3])
    b2["R1"] = base["R1"] * np.exp(th[4])
    b2["R2"] = base["R2"] * np.exp(th[5])
    J = jac(b2, ratio)
    return J


def main():
    res = {"fdtr": {}, "direct": {}}
    # Jacobian accuracy check (h vs h/2 Richardson against h/4)
    L = 500
    b = A1[str(L)]
    est = AU4["table"][f"{float(L):.6g}"]["all"]["base_point"]
    base = {"tauR": b["tau_R_ps"] * 1e-12, "taul": est["taul0_ps"] * 1e-12, "lamz": b["lam_z"], "C": C_ALN,
            "R1": 1e-8, "R2": 2e-8, "d": L * 1e-9}
    J1 = jac(base)
    J2 = jac(base, h=tuple(x / 4 for x in (1e-3, 1e-3, 1e-5, 1e-5, 1e-5, 1e-5)))
    res["jacobian_step_check_max_rel"] = float(np.max(np.abs(J1[:, :6] - J2[:, :6])) / np.max(np.abs(J1[:, :6])))
    print("Jacobian step check", res["jacobian_step_check_max_rel"], flush=True)
    for L in (100, 300, 500):
        b = A1[str(L)]
        t4 = AU4["table"][f"{float(L):.6g}"]["all"]
        est = t4["base_point"]
        for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
            for ttag, taul0 in (("cattaneo", 0.0), ("kinetic", est["taul0_ps"] * 1e-12)):
                base = {"tauR": b["tau_R_ps"] * 1e-12, "taul": taul0, "lamz": b["lam_z"], "C": C_ALN, "R1": R1, "R2": R2,
                        "d": L * 1e-9}
                J = jac(base)
                ent = {}
                for N, nuis in {**SETS, **SETS_EXTRA}.items():
                    for sp in SIGPHI:
                        sR = sigma_dir(J, (1, 0), nuis, sp)
                        sL = sigma_dir(J, (0, -1), nuis, sp)
                        ent[f"{N}|{sp}"] = {"sigma_delta_ps_via_tauR": sR, "sigma_delta_ps_via_taul": sL}
                res["fdtr"][f"{L}|{Rtag}|{ttag}"] = ent
        # margins against the half-ranges (campaign normalisation, and the alternative readings of au4)
        for Rtag in ("note16", "favourable"):
            m = {}
            for N in {**SETS, **SETS_EXTRA}:
                for sp in SIGPHI:
                    sK = res["fdtr"][f"{L}|{Rtag}|kinetic"][f"{N}|{sp}"]["sigma_delta_ps_via_tauR"]
                    sC = res["fdtr"][f"{L}|{Rtag}|cattaneo"][f"{N}|{sp}"]["sigma_delta_ps_via_tauR"]
                    tau0 = b["tau_R_ps"]
                    hr_c = est["halfrange_delta_ps_campaign_rule"]
                    hr_b = est["halfrange_delta_ps_base_delta0"]
                    hr_9 = t4["m996"]["halfrange_delta_over_tau0"] * tau0
                    hr_x = t4["x12_geometric"]["halfrange_delta_over_tau0"] * tau0
                    hr_tau = t4["robust"]["w_tau"] * tau0
                    m[f"{N}|{sp}"] = {"m_kin_campaign_range": sK / hr_c, "m_kin_base_delta0": sK / hr_b,
                                      "m_kin_996only": sK / hr_9, "m_kin_geometric_x12": sK / hr_x,
                                      "m_cattaneo_tau_range": sC / hr_tau}
            res["fdtr"][f"{L}|{Rtag}|margins"] = m
        mm = res["fdtr"][f"{L}|note16|margins"]
        print(L, "note16 margins N2 0.01/0.1:", {k: round(v, 3) for k, v in mm["N2|0.01"].items()},
              {k: round(v, 3) for k, v in mm["N2|0.1"].items()}, flush=True)
        print(L, "note16 margins N3 0.01/0.1:", {k: round(v, 3) for k, v in mm["N3|0.01"].items()},
              {k: round(v, 3) for k, v in mm["N3|0.1"].items()}, flush=True)
        # direct nonlinear test with the extreme witnesses (robust excursions applied to the base point)
        exc = AU4["table"][f"{float(L):.6g}"]["excursions"]
        for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
            base = {"tauR": b["tau_R_ps"] * 1e-12, "taul": est["taul0_ps"] * 1e-12, "lamz": b["lam_z"], "C": C_ALN,
                    "R1": R1, "R2": R2, "d": L * 1e-9}
            y = lnH(np.zeros(9), base)
            per = {}
            for W, ex in exc.items():
                bw = dict(base)
                bw["tauR"] = base["tauR"] * (1 + ex["tau"]["robust"])
                bw["taul"] = base["taul"] * (1 + ex["taul"]["robust"])
                per[W] = {}
                for N, nuis in SETS.items():
                    for sp in (0.01, 0.1):
                        chi, _ = gn_refit(y, bw, nuis, sp)
                        per[W][f"{N}|{sp}"] = chi
            summ = {}
            for N in SETS:
                for sp in (0.01, 0.1):
                    k = f"{N}|{sp}"
                    Wm = max(per, key=lambda x: per[x][k])
                    summ[k] = {"max_dchi2": per[Wm][k], "at": Wm, "m_direct": 1 / np.sqrt(per[Wm][k])}
            res["direct"][f"{L}|{Rtag}"] = {"per_witness": per, "summary": summ}
            print(L, Rtag, "direct:", {k: (round(v["m_direct"], 3), v["at"]) for k, v in summ.items()}, flush=True)
        (OUT / "au5_fisher1d.json").write_text(json.dumps(res, indent=1))

    (OUT / "au5_fisher1d.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
