"""A7: is the film memory time itself visible?  Noiseless profile test of the Fourier film (tau_R = 0) and of
large changes of tau_R against the Cattaneo base point, every thickness, nuisance sets N0..N3 refitted (twin).

Data = twin at the base point (Cattaneo, beams 10/8 um, 80 frequencies 10 kHz - 200 MHz). Alternatives:
tau_R = 0 (Fourier film), tau_R x 2, tau_R x 10. For each, Delta chi2_min = min over the set's nuisances of the
whitened misfit (instrument parameters by exact linear projection, thermal ones by Levenberg-Marquardt).
Linear prediction for a relative change w of tau_R when H is linear in tau_R (omega tau << 1): (w / sigma_F)^2.
Usage: python a7_fourier_limit.py -> out/a7_fourier_limit.json
"""
from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import least_squares

from common import FREQ, OUT, SETS_CAT, model_for, jacobian, weighted, fisher, jsonable

T0 = time.perf_counter()
INSTR = ("g", "phi0", "td_ns")


def qbasis(names, f, sa, sp):
    cols = []
    for n in names:
        if n == "g":
            cols.append(np.concatenate([np.ones(f.size) / sa, np.zeros(f.size)]))
        elif n == "phi0":
            cols.append(np.concatenate([np.zeros(f.size), np.ones(f.size) / sp]))
        elif n == "td_ns":
            cols.append(np.concatenate([np.zeros(f.size), -2 * np.pi * f * 1e-9 / sp]))
    if not cols:
        return None
    Q, _ = np.linalg.qr(np.stack(cols, axis=1))
    return Q


res = {}
for L in (30, 100, 300, 500):
    m = model_for(L)
    lnHd = m.lnH(FREQ)
    J = jacobian(m, FREQ, SETS_CAT["N3"], h=2e-3, scheme="richardson")
    for N in ("N0", "N1", "N2", "N3"):
        names = SETS_CAT[N]
        thn = [n for n in names if n not in INSTR and n != "ln_tauR"]
        for s in (0.01, 0.1):
            sa, sp = 0.01, np.deg2rad(s)
            Q = qbasis([n for n in names if n in INSTR], FREQ, sa, sp)
            sF = fisher(weighted(J[:, [SETS_CAT["N3"].index(n) for n in names]], sig_phi_deg=s))["sigma"]
            row = {"sigma_F": sF}
            for lab, fac in (("fourier", 0.0), ("x2", 2.0), ("x10", 10.0)):
                def resid(x, fac=fac):
                    P = m.phys(dict(zip(thn, x)))
                    P["tauR"] = m.tauR * fac
                    h = m.thermal(FREQ, P)
                    d = lnHd - (np.log(np.abs(h)) + 1j * np.unwrap(np.angle(h)))
                    d = np.log(np.abs(np.exp(d))) + 1j * np.angle(np.exp(d))
                    r = np.concatenate([d.real / sa, d.imag / sp])
                    return r if Q is None else r - Q @ (Q.T @ r)
                if thn:
                    r = least_squares(resid, np.zeros(len(thn)), method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                                      diff_step=1e-4)
                    c2, x = float(2 * r.cost), r.x.tolist()
                else:
                    rr = resid(np.zeros(0))
                    c2, x = float(rr @ rr), []
                w = fac - 1.0
                row[lab] = {"dchi2": c2, "linear_pred_(w/sigmaF)^2": (w / sF) ** 2, "nuisance": dict(zip(thn, x))}
            res[f"L{L}|{N}|{s}"] = row
            print(L, N, s, "sigma_F", round(sF, 4), {k: (round(v["dchi2"], 5), round(v["linear_pred_(w/sigmaF)^2"], 5))
                                                    for k, v in row.items() if k != "sigma_F"},
                  f"{time.perf_counter() - T0:.0f}s", flush=True)
            (OUT / "a7_fourier_limit.json").write_text(json.dumps(jsonable(res), indent=1))
print("done", time.perf_counter() - T0)
