"""Audit, part 2: the twin itself (src/twin_fdtr.py, imported read-only) through an own wrapper (common.py of the
branch is not used).

(a) Hankel quadrature of gaussian_response: trapezoid in x = ln(sqrt(eta) k) at steps 0.32 ... 0.01, x_low = -14, -22,
    -30, and the upper end x_high = 2 (fixed inside gaussian_response) against x_high = 3 (rebuilt here from the
    twin's own trapezoid_nodes and solve_kernel); the twin's default per-frequency step (step=None) as well.
(b) Production Jacobian of the branch (central differences of ln H in the log parameters, Richardson h = 2e-3,
    rule frozen at step 0.04, x_low -22) for all 13 thermal/beam parameters at 30/100/300/500 nm; also plain
    central differences at h = 1e-1 ... 1e-5 for the tau_R column (500 nm, 30 nm) to show the FD error curve.
Outputs out/aud2_twin_checks.json (complex numbers as [re, im]).
Usage: python aud2_twin_checks.py
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np  # noqa: E402

sys.path.insert(0, r"C:/Users/Koussay/these/src")
import twin_fdtr as tw  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
INP = json.loads((OUT / "aud_inputs.json").read_text())
F = np.array(INP["freq80"])
S = INP["stack"]
PN = ["ln_tauR", "ln_lamz", "ln_lamr", "ln_C", "ln_R1", "ln_R2", "ln_lamAu", "ln_CAu", "ln_dAu", "ln_dAlN",
      "ln_lamS", "ln_CS", "ln_wscale"]
T0 = time.perf_counter()


def phys(Lkey, theta=None):
    b = INP["base"][Lkey]
    L = float(Lkey.split("_")[0]) * 1e-9
    th = dict(theta or {})
    e = lambda n: float(np.exp(th.get(n, 0.0)))  # noqa: E731
    return {"tauR": b["tauR"] * e("ln_tauR"), "lamz": b["lamz"] * e("ln_lamz"), "lamr": b["lamr"] * e("ln_lamr"),
            "C": S["C"] * e("ln_C"), "R1": S["R1"] * e("ln_R1"), "R2": S["R2"] * e("ln_R2"),
            "lamAu": S["lamAu"] * e("ln_lamAu"), "CAu": S["CAu"] * e("ln_CAu"), "dAu": S["dAu"] * e("ln_dAu"),
            "dAlN": L * e("ln_dAlN"), "lamS": S["lamS"] * e("ln_lamS"), "CS": S["CS"] * e("ln_CS"),
            "wp": S["wp"] * e("ln_wscale"), "wr": S["wr"] * e("ln_wscale")}


def stack_from(P):
    aln = tw.Material(P["lamz"], P["C"], lam_r=P["lamr"], tau_r=P["tauR"])
    return tw.au_aln_sapphire(aln, aln_thickness=P["dAlN"], au_thickness=P["dAu"],
                              au=tw.Material(P["lamAu"], P["CAu"]), sapphire=tw.Material(P["lamS"], P["CS"]),
                              r_au_aln=P["R1"], r_aln_sapphire=P["R2"])


def H(P, step=0.04, x_low=-22.0, f=F):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tw.gaussian_response(f, stack_from(P), P["wp"], P["wr"], step=step, x_low=x_low, passivity="ignore")


def H_xhigh(P, step, x_low, x_high, f=F):
    st = stack_from(P)
    kn, w = tw.trapezoid_nodes(f, st, P["wp"], P["wr"], step, x_low, x_high)
    out = np.empty(f.size, complex)
    for i in range(0, f.size, 8):
        sol = tw.solve_kernel(f[i:i + 8][:, None], kn[None, :], st)
        out[i:i + 8] = sol.impedance @ w
    return out


def cz(a):
    return [[float(z.real), float(z.imag)] for z in np.asarray(a)]


def jac_fd(Lkey, names, h=2e-3, scheme="richardson", step=0.04, x_low=-22.0):
    cols = []
    for n in names:
        def cd(hh):
            hp = H(phys(Lkey, {n: hh}), step, x_low)
            hm = H(phys(Lkey, {n: -hh}), step, x_low)
            return np.log(hp / hm) / (2 * hh)
        if scheme == "richardson":
            d1, d2 = cd(h), cd(h / 2)
            cols.append((4 * d2 - d1) / 3)
        else:
            cols.append(cd(h))
    return np.stack(cols, axis=1)


res = {"freq": F.tolist(), "quadrature": {}, "jacobian": {}, "fd_curve": {}}
for Lk in ("30", "100", "300", "500"):
    P0 = phys(Lk)
    q = {}
    for st in (0.32, 0.16, 0.08, 0.04, 0.02, 0.01):
        q[f"step{st}_xlow-22"] = cz(np.log(H(P0, st, -22.0)))
    for xl in (-14.0, -30.0):
        q[f"step0.04_xlow{xl:g}"] = cz(np.log(H(P0, 0.04, xl)))
    q["step0.04_xlow-22_xhigh3"] = cz(np.log(H_xhigh(P0, 0.04, -22.0, 3.0)))
    q["step0.04_xlow-22_xhigh1.5"] = cz(np.log(H_xhigh(P0, 0.04, -22.0, 1.5)))
    q["default_step"] = cz(np.log(H(P0, None, None)))
    res["quadrature"][Lk] = q
    t1 = time.perf_counter()
    J = jac_fd(Lk, PN)
    res["jacobian"][Lk] = {"names": PN, "J": [cz(J[:, j]) for j in range(len(PN))]}
    print(Lk, "quadrature and production Jacobian done", f"{time.perf_counter() - T0:.0f}s (J {time.perf_counter() - t1:.0f}s)",
          flush=True)
    (OUT / "aud2_twin_checks.json").write_text(json.dumps(res))
for Lk in ("500", "30"):
    curve = {}
    for h in (1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-5):
        curve[f"central_{h:g}"] = cz(jac_fd(Lk, ["ln_tauR"], h=h, scheme="central")[:, 0])
    curve["richardson_2e-3"] = cz(jac_fd(Lk, ["ln_tauR"], h=2e-3)[:, 0])
    res["fd_curve"][Lk] = curve
    print(Lk, "fd curve done", f"{time.perf_counter() - T0:.0f}s", flush=True)
res["seconds"] = time.perf_counter() - T0
(OUT / "aud2_twin_checks.json").write_text(json.dumps(res))
print("done", res["seconds"])
