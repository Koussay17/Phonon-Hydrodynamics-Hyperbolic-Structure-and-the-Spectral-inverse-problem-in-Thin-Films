"""AU3b: exact diffusive root (all orders in k) of the regularised conserving event operator on 7x7x5 for REF and the
witnesses that set the ends of the delta range, and the FDTR distinguishability with the exact kernel.

Operator from au3_event775 (independent assembly; E-QDQ completion, as in the campaign). Kernel on e-perp:
  kappa(k, p) = b.(p + C_hat + g ehat ehat^T + i k Q V Q)^-1 b      (dense complex LU per evaluation),
root s = k^2 of p|e|^2 + s kappa = 0 by secant iteration in s from the first-order-in-k^2 guess, at 12 frequencies
1-200 MHz; lam_eff/lam0 = -p|e|^2/(s K0). Below 1 MHz the first-order form is used (|k| L_GK < 0.05 there).
FDTR (au5 stack, note-16 and favourable R): data of REF vs each witness, both with the exact kernel, nuisances refitted
(N1/N2/N3, 0.01 and 0.1 deg), against the same test with the first-order GK films of the static times (lam_z, C at
the A1 base for every tensor, as in the campaign). Ratio of the Delta chi^2 values = how much of the first-order
distinguishability survives the k^2 truncation.
Output: out/au3b_exact_root775.json
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
from scipy.interpolate import CubicSpline  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au3_event775 as a3  # noqa: E402
import au5_fisher1d as a5  # noqa: E402

FOURPI = 4 * np.pi
NAMES = ["REF", "cal3w_m3L100max", "cal3w_m3L100min", "orb0.1_m3L30max", "cal1nz_m3L100min"]
FR = np.array([1e6, 3e6, 1e7, 2e7, 3.5e7, 5e7, 6.5e7, 8e7, 1e8, 1.4e8, 1.7e8, 2e8])


def operator(C, e, vz, L):
    ehat = e / np.linalg.norm(e)
    d = np.abs(vz) / (FOURPI * L * 10.0)
    h = d * ehat
    Ch = C + np.diag(d) - np.outer(ehat, h) - np.outer(h, ehat) + (ehat @ h) * np.outer(ehat, ehat)
    gam = np.trace(Ch) / e.size
    return Ch + gam * np.outer(ehat, ehat), ehat


def static_times(Cr, e, ehat, vz):
    b = vz * e
    cf = sla.cho_factor(Cr, lower=True, check_finite=False)
    x = sla.cho_solve(cf, b, check_finite=False)
    K = b @ x
    u = vz * x
    u = u - ehat * (ehat @ u)
    M = u @ sla.cho_solve(cf, u, check_finite=False)
    return K, (x @ x) / K / FOURPI, (e @ e) * M / (FOURPI * K * K), M


def kernel(Cr, ehat, QVQ, b, k, pc):
    A = Cr + 1j * k * QVQ
    A[np.diag_indices_from(A)] += pc
    x = sla.solve(A, b.astype(complex), assume_a="sym", check_finite=False)   # complex symmetric (LDL^T)
    return b @ x


def root(Cr, e, ehat, QVQ, vz, K0, M0, f):
    """s = k^2 in operator units (velocity THz*A, rate THz/(4 pi)): p_c |e|^2 + s kappa(sqrt s, p_c) = 0."""
    pc = 1j * 2 * np.pi * f * 1e-12 / FOURPI
    e2 = e @ e
    b = vz * e
    # first-order-in-k^2 guess: kappa ~ K0 - s M0  ->  p e2 + s K0 - s^2 M0 = 0 (small root)
    s = -pc * e2 / K0
    s = s * (1 + pc * e2 * M0 / K0 ** 2)
    def F(s_):
        k_ = np.sqrt(s_)
        return pc * e2 + s_ * kernel(Cr, ehat, QVQ, b, k_, pc)
    s1, F1 = s, F(s)
    s2 = s * (1 + 1e-3)
    F2 = F(s2)
    for it in range(30):
        s3 = s2 - F2 * (s2 - s1) / (F2 - F1)
        s1, F1 = s2, F2
        s2, F2 = s3, F(s3)
        if abs(s2 - s1) < 1e-12 * abs(s2):
            break
    return complex(-pc * e2 / s2 / K0), it + 3


def main():
    t0 = time.perf_counter()
    G = a3.geometry()
    ref = np.load(a3.PPD / "pp_r01_REF.npz")
    conv, g0w = float(ref["conv"]), ref["g0w"]
    e, vz = G["e"], G["vz"]
    ehat = e / np.linalg.norm(e)
    Qm = np.eye(e.size) - np.outer(ehat, ehat)
    QVQ = Qm @ (vz[:, None] * Qm)
    res = {"freqs_Hz": FR.tolist(), "names": NAMES}
    for nm in NAMES:
        g = a3.rates(G, np.load(a3.PPD / f"pp_r01_{nm}.npz")["pp"], conv, g0w)
        C = a3.assemble(G, g)
        res[nm] = {}
        for L in (300.0, 500.0):
            Cr, _ = operator(C, e, vz, L)
            K0, tauR, taul, M0 = static_times(Cr, e, ehat, vz)
            rows = []
            for f in FR:
                lam, its = root(Cr, e, ehat, QVQ, vz, K0, M0, f)
                rows.append({"f": float(f), "lam_eff_rel": [lam.real, lam.imag], "iters": its})
            res[nm][f"{L:g}"] = {"tauR_ps": tauR, "taul_ps": taul, "rows": rows}
            r200 = rows[-1]["lam_eff_rel"]
            pp = 2j * np.pi * 2e8
            gk = (1 + taul * 1e-12 * pp) / (1 + tauR * 1e-12 * pp)
            print(nm, L, f"tauR {tauR:.3f} taul {taul:.2f} | 200 MHz exact {r200[0]:.4f}{r200[1]:+.4f}i gk {gk.real:.4f}{gk.imag:+.4f}i",
                  f"({time.perf_counter() - t0:.0f}s)", flush=True)
        (OUT / "au3b_exact_root775.json").write_text(json.dumps(res, indent=1))
    # ---------------------------------------------------------------- FDTR distinguishability
    A1 = a5.A1
    fd = {}
    for L in (300, 500):
        def ratio_exact(nm):
            rows = res[nm][f"{L:g}"]["rows"]
            f = np.array([r["f"] for r in rows])
            lam = np.array([complex(*r["lam_eff_rel"]) for r in rows])
            p = 2j * np.pi * f
            deff = -(lam - 1) / p
            tR, tl = res[nm][f"{L:g}"]["tauR_ps"] * 1e-12, res[nm][f"{L:g}"]["taul_ps"] * 1e-12
            # anchor below 1 MHz on the static first-order value
            fx = np.concatenate([[1e3, 1e4, 1e5], f])
            dx = np.concatenate([[tR - tl] * 3, deff])
            sre, sim = CubicSpline(np.log(fx), dx.real), CubicSpline(np.log(fx), dx.imag)
            dI = sre(np.log(a5.FREQ)) + 1j * sim(np.log(a5.FREQ))
            return 1 - 2j * np.pi * a5.FREQ * dI
        for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
            b = A1[str(L)]
            r0 = res["REF"][f"{L:g}"]
            base = {"tauR": r0["tauR_ps"] * 1e-12, "taul": r0["taul_ps"] * 1e-12, "lamz": b["lam_z"], "C": a5.C_ALN,
                    "R1": R1, "R2": R2, "d": L * 1e-9}
            rr = ratio_exact("REF")
            y_ref_ex = a5.lnH(np.zeros(9), base, rr)
            y_ref_gk = a5.lnH(np.zeros(9), base)
            ent = {"truncation_dchi2": {}}
            for N, nuis in a5.SETS.items():
                for sp in (0.01, 0.1):
                    c_, _ = a5.gn_refit(y_ref_ex, base, nuis, sp)
                    ent["truncation_dchi2"][f"{N}|{sp}"] = c_
            for nm in NAMES[1:]:
                rw = res[nm][f"{L:g}"]
                bw = dict(base)
                bw["tauR"], bw["taul"] = rw["tauR_ps"] * 1e-12, rw["taul_ps"] * 1e-12
                y_w_ex = a5.lnH(np.zeros(9), bw, ratio_exact(nm))
                y_w_gk = a5.lnH(np.zeros(9), bw)
                per = {}
                for N, nuis in a5.SETS.items():
                    for sp in (0.01, 0.1):
                        cex, _ = a5.gn_refit(y_w_ex, base, nuis, sp, ratio=rr)
                        cgk, _ = a5.gn_refit(y_w_gk, base, nuis, sp)
                        per[f"{N}|{sp}"] = {"dchi2_exact": cex, "dchi2_firstorder": cgk,
                                            "ratio": cex / cgk if cgk > 0 else None,
                                            "m_exact": 1 / np.sqrt(cex) if cex > 0 else None,
                                            "m_firstorder": 1 / np.sqrt(cgk) if cgk > 0 else None}
                ent[nm] = per
            fd[f"{L}|{Rtag}"] = ent
            print(L, Rtag, {nm: {k: (round(v["m_exact"], 3), round(v["m_firstorder"], 3)) for k, v in ent[nm].items()
                                 if k.startswith("N2") or k.startswith("N3")} for nm in NAMES[1:]}, flush=True)
    res["fdtr"] = fd
    res["seconds"] = time.perf_counter() - t0
    (OUT / "au3b_exact_root775.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
