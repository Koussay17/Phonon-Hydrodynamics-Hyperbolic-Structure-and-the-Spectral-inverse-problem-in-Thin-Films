"""R-A1d: frequency dependence of the first-order film response of the regularised event operator (objection A-2).

The static times (tau_R = tau_mem, tau_l = k^2 coefficient, R-A1) describe the film to first order in p only if the
modes that carry tau_l respond quasi-statically in the FDTR band. Test with the exact kernel at finite p = i omega:
  K(p) = b.(p + C')^-1 b,  M(p) = u(p).(p + C')^-1 u(p),  u(p) = Q V (p + C')^-1 b      (complex symmetric solves),
  C' = C + Q D Q + gamma ehat ehat^T (energy-conserving regulator, as R-A1),
  effective conductivity of the dispersion relation p |e|^2 + k^2 kappa(k, p) = 0 to first order in k^2:
  lam_eff(p) / lam(0) = [K(p) + p |e|^2 M(p) / K(p)] / K(0)     (code units; p in operator units p_phys / 4 pi),
  delta_eff(omega) = -(lam_eff/lam0 - 1) / (i omega)   (complex; -> tau_R - tau_l as omega -> 0),
compared with Guyer-Krumhansl (1 + tau_l p)/(1 + tau_R p) with the static kinetic times and Cattaneo 1/(1 + tau_R p).
Check (7x7x5, REF, 500 nm): the k^2 coefficient M(0) against the static kernel computed directly at small k by
complex solves on e-perp (Richardson in k).
Memory: padded (n+1)^2 buffer as R-A1; one real and one complex (n+1)^2 array at a time.
usage: python -B ra1d_freq_kernel.py m775|m996 L1,L2,... [name ...]
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

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
sys.path.insert(0, str(FC3 / "scripts_c"))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ra1_film_event_ext import _C_padded, WITNESSES  # noqa: E402

OUT = HERE / "out"
mesh = sys.argv[1]
LS = tuple(float(x) for x in sys.argv[2].split(","))
MESH = {"m775": ("results_c/m775x/events_m775_s0.1.npz", "results_c/m775x/film2_pp"),
        "m996": ("results_c/m996/events_m996_s0.1.npz", "results_c/m996/film2_pp")}
NAMES = sys.argv[3:] or (["REF"] + WITNESSES)
FREQS = np.array([1e4, 1e5, 1e6, 3e6, 1e7, 2e7, 4e7, 6.5e7, 1e8, 1.4e8, 2e8, 5e8, 1e9, 3e9])

t0 = time.perf_counter()
M = FC3Class(str(FC3 / MESH[mesh][0]), FC3 / MESH[mesh][1] / "pp_r01_REF.npz", allow_missing=True)
for attr in ("wpair_flat", "wpair_val", "rp", "Wg_dense"):
    if hasattr(M.ops, attr):
        setattr(M.ops, attr, None)
for attr in ("A", "W"):
    if hasattr(M.geom, attr):
        setattr(M.geom, attr, None)
geom = M.geom
n = geom.n
vz, b, e, ehat = geom.vel[:, 2], geom.b[:, 2], geom.e, geom.ehat
e2 = float(e @ e)
bp = np.append(b, 0.0)
outp = OUT / f"ra1d_freq_kernel_{mesh}.json"
res = json.loads(outp.read_text()) if outp.exists() else {"mesh": mesh, "freqs_Hz": FREQS.tolist()}
print("setup", round(time.perf_counter() - t0), "s", flush=True)


def build(g, L):
    """C + Q D Q + gamma ehat ehat^T on the padded buffer (dummy slot decoupled), by row blocks."""
    r = M.ops.diag_full(g)
    gam = float(r.sum()) / n
    d = np.abs(vz) / (FOURPI * L * 10.0)
    dt = d * ehat
    c0 = float(ehat @ dt) + gam
    CL = _C_padded(M, g)
    di = np.arange(n)
    CL[di, di] += d
    for i0 in range(0, n, 1024):
        i1 = min(i0 + 1024, n)
        CL[i0:i1, :n] += (-ehat[i0:i1, None] * dt[None, :] - dt[i0:i1, None] * ehat[None, :]
                          + c0 * ehat[i0:i1, None] * ehat[None, :])
    return CL


for nm in NAMES:
    g = M.rates(np.load(FC3 / MESH[mesh][1] / f"pp_r01_{nm}.npz")["pp"])
    res.setdefault(nm, {})
    for L in LS:
        t1 = time.perf_counter()
        CL = build(g, L)
        cf = sla.cho_factor(CL, lower=True, overwrite_a=True, check_finite=False)
        x0 = sla.cho_solve(cf, bp, check_finite=False)[:n]
        K0, N0 = float(b @ x0), float(x0 @ x0)
        u0 = vz * x0
        u0 = u0 - ehat * float(ehat @ u0)
        M0 = float(u0 @ sla.cho_solve(cf, np.append(u0, 0.0), check_finite=False)[:n])
        del cf, CL
        CL = build(g, L)                                   # rebuilt (the factorisation overwrote it); kept real
        tauR, taul = N0 / K0 / FOURPI, e2 * M0 / (FOURPI * K0 * K0)
        k2chk = None
        if nm == "REF" and n < 4000 and abs(L - 500.0) < 1e-9:
            Cr_full = CL[:n, :n].copy()
            v = ehat.copy()
            v[0] += np.sign(v[0])
            v /= np.linalg.norm(v)
            Bh = np.eye(n) - 2.0 * np.outer(v, v)
            Bh = Bh[:, 1:]
            Cr = Bh.T @ Cr_full @ Bh
            Vr = Bh.T @ (vz[:, None] * Bh)
            br = Bh.T @ b

            def kap(k):
                return complex(br @ np.linalg.solve(Cr + 1j * k * Vr, br.astype(complex)))
            k0v = 1e-4 * float(np.min(np.diag(Cr_full))) / float(np.max(np.abs(vz)))
            Kk0 = kap(0.0).real
            c1 = (Kk0 - 0.5 * (kap(k0v) + kap(-k0v)).real) / k0v ** 2
            c2 = (Kk0 - 0.5 * (kap(k0v / 2) + kap(-k0v / 2)).real) / (k0v / 2) ** 2
            coef = (4 * c2 - c1) / 3
            k2chk = {"K_direct": Kk0, "K_formula": K0, "M_direct": coef, "M_formula": M0,
                     "rel_diff_K": abs(Kk0 / K0 - 1), "rel_diff_M": abs(coef / M0 - 1), "k_used": k0v}
            del Bh, Cr, Vr, Cr_full
        rows = []
        idx = np.arange(n)
        for jf, f in enumerate(FREQS):
            pc = 1j * 2 * np.pi * f * 1e-12 / FOURPI
            Af = CL.astype(complex)
            Af[idx, idx] += pc
            lu = sla.lu_factor(Af, overwrite_a=True, check_finite=False)
            x = sla.lu_solve(lu, bp.astype(complex), check_finite=False)[:n]
            K = complex(b @ x)
            u = vz * x
            u = u - ehat * complex(ehat @ u)
            y = sla.lu_solve(lu, np.append(u, 0.0), check_finite=False)[:n]
            del lu, Af
            Mp = complex(u @ y)
            lam = (K + pc * e2 * Mp / K) / K0
            pphys = 1j * 2 * np.pi * f * 1e-12
            gk = (1 + taul * pphys) / (1 + tauR * pphys)
            cat = 1 / (1 + tauR * pphys)
            deff = -(lam - 1) / pphys
            rows.append({"f": float(f), "lam_eff_rel": [lam.real, lam.imag], "gk_static": [gk.real, gk.imag],
                         "cattaneo": [cat.real, cat.imag], "delta_eff_ps": [deff.real, deff.imag],
                         "K_rel": [(K / K0).real, (K / K0).imag]})
        del CL
        res[nm][f"{L:.6g}"] = {"tauR_static": tauR, "taul_static": taul, "delta_static": tauR - taul, "rows": rows,
                               "k2_direct_check": k2chk}
        print(nm, L, f"tauR {tauR:.3f} taul {taul:.2f} ({time.perf_counter() - t1:.0f}s)", " ".join(
            f"{r['f']:.0e}:{r['delta_eff_ps'][0]:.1f}{r['delta_eff_ps'][1]:+.1f}i" for r in rows[::2]),
            ("k2chk %.2e/%.2e" % (k2chk["rel_diff_K"], k2chk["rel_diff_M"])) if k2chk else "", flush=True)
        tmp = outp.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(res, indent=1))
        tmp.replace(outp)
print("total", round(time.perf_counter() - t0), "s")
