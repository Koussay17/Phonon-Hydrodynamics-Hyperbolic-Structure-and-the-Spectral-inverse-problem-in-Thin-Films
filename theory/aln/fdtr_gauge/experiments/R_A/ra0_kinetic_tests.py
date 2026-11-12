"""R-A0: manufactured tests of the kinetic nonlocal time (ra_kinetic.py) before it is applied to the AlN operator.

T1  gray isotropic Callaway operator: tau_R(kinetic) = tau_R, tau_l = (4/5) tau_c for tau_N/tau_R in [0.01, 100];
    the unprojected (historical) variant |e|^2 tau_c |V b|^2 / |b|^4 = (9/5) tau_c.
T2  exact static kernel kappa(k, 0) by complex solves on e-perp: its k^2 coefficient (Richardson in k) equals
    -K L^2 from the closed formula, for the gray operator and for a random non-gray symmetric operator with the same
    conservation laws (energy exactly, momentum by N only), parity-symmetric mode set.
T3  invariance: tau_l does not depend on the regularisation constant gamma of the padded Cholesky route used for the
    event operator (C + gamma ehat ehat^T), compared with the Householder reduction to e-perp.
Output: out/ra0_kinetic_tests.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from ra_kinetic import kinetic_times_dense, kernel_k2_numeric, gray_isotropic_operator  # noqa: E402

OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261112)
res = {}

# ---------------------------------------------------------------- T1
t1 = []
for x in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
    Ch, e, vz = gray_isotropic_operator(n_mu=8, n_phi=16, v=1.0, tauR=1.0, tauN=x)
    r = kinetic_times_dense(Ch, e, vz)
    tauc = 1.0 / (1.0 / x + 1.0)
    b = vz * e
    hist = float(e @ e) * tauc * float((vz * b) @ (vz * b)) / float(b @ b) ** 2
    t1.append({"tauN_over_tauR": x, "tauR_kin": r["tauR"], "taul_over_tauc": r["taul"] / tauc,
               "gamma_g": r["gamma_g"], "hist_over_tauc": hist / tauc,
               "delta_over_tauR": (r["tauR"] - r["taul"]) / 1.0,
               "delta_closure_formula": 1.0 - 0.8 * tauc})
res["T1_gray_isotropic"] = t1
err1 = max(max(abs(d["taul_over_tauc"] - 0.8), abs(d["tauR_kin"] - 1.0), abs(d["hist_over_tauc"] - 1.8),
               abs(d["gamma_g"] - 0.8)) for d in t1)
res["T1_max_abs_error"] = err1
print("T1 max |error|", err1)


# ---------------------------------------------------------------- T2
def k2_coefficient(Ch, e, vz, k0):
    """-(d kappa / d k^2) at k = 0 by Richardson on the even part, from complex solves at k0, k0/2."""
    def c2(k):
        kap = kernel_k2_numeric(Ch, e, vz, [k, -k])
        K0 = kernel_k2_numeric(Ch, e, vz, [0.0])[0]
        even = 0.5 * (kap[0] + kap[1])
        return (K0 - even).real / k ** 2, K0.real
    a1, K0 = c2(k0)
    a2, _ = c2(k0 / 2)
    return (4 * a2 - a1) / 3, K0


def random_nongray(n_half=40, seed=1):
    """Parity-symmetric mode set (pairs +v, -v), random capacities and velocities; operator
    Chat = Q_e S_R Q_e + Q_em S_N Q_em with S random SPD, Q_e = 1 - P_e, Q_em = 1 - P_e - P_m (z-momentum mode
    m ~ e * k_z with random 'wavevectors'); both parts commute with the parity permutation."""
    r = np.random.default_rng(seed)
    w = r.uniform(0.2, 1.0, n_half)
    v = r.normal(size=n_half)
    kz = v * r.uniform(0.5, 1.5, n_half)            # momentum not proportional to the current (non-gray)
    e = np.sqrt(np.concatenate([w, w]))
    vz = np.concatenate([v, -v])
    mz = e * np.concatenate([kz, -kz])
    n = 2 * n_half
    Pp = np.zeros((n, n))
    Pp[np.arange(n_half), np.arange(n_half) + n_half] = 1
    Pp[np.arange(n_half) + n_half, np.arange(n_half)] = 1

    def sym_spd(scale):
        A = r.normal(size=(n, n))
        S = A @ A.T / n + np.diag(r.uniform(0.5, 2.0, n))
        S = 0.5 * (S + Pp @ S @ Pp)                   # parity symmetric
        return scale * S
    ehat = e / np.linalg.norm(e)
    mhat = mz / np.linalg.norm(mz)
    Qe = np.eye(n) - np.outer(ehat, ehat)
    # mz is odd, ehat even, so they are orthogonal
    Qem = Qe - np.outer(mhat, mhat)
    Ch = Qe @ sym_spd(1.0) @ Qe + Qem @ sym_spd(3.0) @ Qem
    return 0.5 * (Ch + Ch.T), e, vz


t2 = {}
for name, (Ch, e, vz) in {"gray_x1": gray_isotropic_operator(6, 12, 1.0, 1.0, 1.0),
                          "gray_x0.1": gray_isotropic_operator(6, 12, 1.0, 1.0, 0.1),
                          "nongray_1": random_nongray(40, 1), "nongray_2": random_nongray(60, 7)}.items():
    r = kinetic_times_dense(Ch, e, vz)
    c2, K0 = k2_coefficient(Ch, e, vz, 1e-3)
    t2[name] = {"K_formula": r["K"], "K_numeric": K0, "KL2_formula": r["K"] * r["L2"], "k2_coef_numeric": c2,
                "rel_err_L2": abs(c2 / (r["K"] * r["L2"]) - 1), "taul": r["taul"], "tauR": r["tauR"],
                "gamma_g": r["gamma_g"]}
    print(name, t2[name]["rel_err_L2"], r["taul"] / r["tauR"], r["gamma_g"])
res["T2_kernel_expansion"] = t2

# ---------------------------------------------------------------- T3 (gamma independence, padded-Cholesky route)
Ch, e, vz = random_nongray(50, 3)
ehat = e / np.linalg.norm(e)
ref = kinetic_times_dense(Ch, e, vz)
t3 = []
for gam in (1e-3, 1.0, np.trace(Ch) / len(e), 1e3):
    A = Ch + gam * np.outer(ehat, ehat)
    cf = sla.cho_factor(A, lower=True)
    b = vz * e
    x = sla.cho_solve(cf, b)
    K = b @ x
    u = vz * x
    u = u - ehat * (ehat @ u)
    y = sla.cho_solve(cf, u)
    taul = (e @ e) * (u @ y) / K ** 2
    t3.append({"gamma": float(gam), "taul_rel_diff": abs(taul / ref["taul"] - 1),
               "tauR_rel_diff": abs((x @ x) / K / ref["tauR"] - 1), "ehat_dot_x": float(abs(ehat @ x)),
               "ehat_dot_y": float(abs(ehat @ y))})
res["T3_regularisation_independence"] = t3
print("T3", max(d["taul_rel_diff"] for d in t3))
(OUT / "ra0_kinetic_tests.json").write_text(json.dumps(res, indent=1))
