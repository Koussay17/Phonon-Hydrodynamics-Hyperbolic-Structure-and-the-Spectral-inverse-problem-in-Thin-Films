"""R-A5b: direct test of the observable for objection A-2 (cf. referee decisive test 7), 1D large-beam limit.

Each tensor W (REF and the three-mesh optima) defines the dynamic part of the film conductivity
r_W(omega) = lam_eff,W(omega) / lam_eff,W(0) from the exact first-order-in-k^2 kernel of its regularised event
operator (R-A1d, 7x7x5 and, where computed, 9x9x6), with no Guyer-Krumhansl or Cattaneo fit. The film of the twin
stack (Au 80 nm / R1 / AlN L / R2 / sapphire, A4 values, lam_z and C fixed at the A1 base) is evaluated with
lam_eff(omega) = lam_z r_W(omega) (local interface-temperature convention, as the twin's GK film in 1D), and
  Delta chi^2_W = || (I - P_nuis) w [ln H(r_W) - ln H(r_REF)] ||^2   (whitened, nuisances N1/N2/N3 profiled to first
  order; noise 1 % on ln|H|, sigma_phi on arg H, the 80 A4 frequencies).
Compared with the same test for (a) Cattaneo r = 1/(1 + i omega tau_R,W) with tau_R,W = tau_mem of W (the quantity A
compares), (b) static Guyer-Krumhansl r = (1 + i omega tau_l)/(1 + i omega tau_R) with the static kinetic times, and
(c) the same with the entropic interface-temperature convention (theta = (1 + tau_l p) T at both AlN interfaces).
Margin analogue of A: m = 1 / sqrt(max_W Delta chi^2_W)  (for Cattaneo this is sigma(ln tau_R)/w at the extreme W).
r_W(omega) is interpolated from the R-A1d frequency grid to the 80 frequencies by cubic splines of delta_eff in ln f.
Output: out/ra5b_observable.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE.parent / "A"))
import common  # noqa: E402

FREQ = common.FREQ
a1 = json.loads((HERE.parent / "A" / "out" / "a1_base_values.json").read_text())["base"]
SIGPHI = (0.01, 0.03, 0.1)


def zfilm(f, r, P):
    """1D impedance Au/R1/AlN(lam_z r(omega))/R2/sapphire; r complex array over f.
    Local interface-temperature convention unless P['g'] is given: entropic convention at k = 0, the interfaces see
    theta = g T with g = 1 + tau_l p (twin GKBoundary 'entropic'), i.e. load Z_b/g below and g Z_film above."""
    p = 2j * np.pi * f
    lam = P["lamz"] * r
    s = np.sqrt(P["C"] * p / lam)
    Y = lam * s
    th = np.tanh(s * P["d"])
    g = P.get("g", 1.0)
    Zb = (1 / np.sqrt(P["lamS"] * P["CS"] * p) + P["R2"]) / g
    Z = g * (Zb + th / Y) / (1 + Y * Zb * th) + P["R1"]
    sa = np.sqrt(P["CAu"] * p / P["lamAu"])
    Ya = P["lamAu"] * sa
    ta = np.tanh(sa * P["dAu"])
    return (Z + ta / Ya) / (1 + Ya * Z * ta)


def lnH(r, P, inst=(0.0, 0.0, 0.0)):
    z = np.log(zfilm(FREQ, r, P))
    return z + inst[0] + 1j * (inst[1] - 2 * np.pi * FREQ * inst[2] * 1e-9)


def nuis_jac(r, P):
    """Columns g, phi0, td_ns (exact), ln R1, ln R2, ln lam_z, ln C (central differences, h = 1e-5)."""
    cols = [np.ones(FREQ.size, complex), np.full(FREQ.size, 1j), -1j * 2 * np.pi * FREQ * 1e-9]
    for k in ("R1", "R2", "lamz", "C"):
        h = 1e-5
        Pp, Pm = dict(P), dict(P)
        Pp[k] *= np.exp(h)
        Pm[k] *= np.exp(-h)
        cols.append((lnH(r, Pp) - lnH(r, Pm)) / (2 * h))
    return np.stack(cols, axis=1)


SETS = {"N1": [0, 1, 2], "N2": [0, 1, 2, 3, 4], "N3": [0, 1, 2, 3, 4, 5, 6]}


def dchi2(dz, Jn, sp):
    w = np.concatenate([np.full(FREQ.size, 1 / common.SIG_A), np.full(FREQ.size, 1 / np.deg2rad(sp))])
    v = w * np.concatenate([dz.real, dz.imag])
    J = w[:, None] * np.vstack([Jn.real, Jn.imag])
    sc = np.linalg.norm(J, axis=0)
    Q, _ = np.linalg.qr(J / sc)
    res = v - Q @ (Q.T @ v)
    res = res - Q @ (Q.T @ res)
    return float(res @ res)


def r_kinetic(rows):
    f = np.array([x["f"] for x in rows])
    de = np.array([complex(*x["delta_eff_ps"]) for x in rows])
    sre, sim = CubicSpline(np.log(f), de.real), CubicSpline(np.log(f), de.imag)
    dI = sre(np.log(FREQ)) + 1j * sim(np.log(FREQ))
    return 1 - 2j * np.pi * FREQ * 1e-12 * dI


out = {}
# check of the closed-form 1D film (local and entropic conventions) against the twin's k = 0 kernel (GK AlN layer)
chk = {}
for conv in ("local", "entropic"):
    for L, tR, tl in ((500, 25e-12, 180e-12), (100, 11e-12, 30e-12)):
        m = common.Model(L_nm=float(L), tauR=tR, lamz=a1[str(L)]["lam_z"], lamr=a1[str(L)]["lam_z"], law="gk",
                         taul=tl, geometry="1d", extra={"temperature": conv})
        h = m.thermal(FREQ, m.phys())
        P = {"lamz": a1[str(L)]["lam_z"], "C": common.C_ALN, "d": L * 1e-9, "lamS": common.SAP_LAM,
             "CS": common.SAP_C, "R1": common.R1, "R2": common.R2, "lamAu": common.AU_LAM, "CAu": common.AU_C,
             "dAu": common.AU_D}
        pp = 2j * np.pi * FREQ
        if conv == "entropic":
            P["g"] = 1 + pp * tl
        zz = zfilm(FREQ, (1 + pp * tl) / (1 + pp * tR), P)
        ratio = h / zz
        chk[f"{conv}|{L}"] = float(np.max(np.abs(ratio / ratio[0] - 1)))
out["check_closed_form_vs_twin_k0"] = chk
print("closed-form 1D film vs twin kernel(k=0), max |ratio/ratio0 - 1|:", chk, flush=True)
for mesh in ("m775", "m996"):
    fp = OUT / f"ra1d_freq_kernel_{mesh}.json"
    if not fp.exists():
        continue
    rd = json.loads(fp.read_text())
    names = [k for k in rd if k not in ("mesh", "freqs_Hz")]
    for Lk in rd["REF"]:
        L = int(float(Lk))
        P = {"lamz": a1[str(L)]["lam_z"], "C": common.C_ALN, "d": L * 1e-9, "lamS": common.SAP_LAM,
             "CS": common.SAP_C, "R1": common.R1, "R2": common.R2, "lamAu": common.AU_LAM, "CAu": common.AU_C,
             "dAu": common.AU_D}
        for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
            P["R1"], P["R2"] = R1, R2
            ref = rd["REF"][Lk]
            gref = 1 + 2j * np.pi * FREQ * 1e-12 * ref["taul_static"]
            r_ref = {"kinetic": r_kinetic(ref["rows"]),
                     "gk_static_entropic": (1 + 2j * np.pi * FREQ * 1e-12 * ref["taul_static"]) /
                                           (1 + 2j * np.pi * FREQ * 1e-12 * ref["tauR_static"]),
                     "gk_static": (1 + 2j * np.pi * FREQ * 1e-12 * ref["taul_static"]) /
                                  (1 + 2j * np.pi * FREQ * 1e-12 * ref["tauR_static"]),
                     "cattaneo": 1 / (1 + 2j * np.pi * FREQ * 1e-12 * ref["tauR_static"])}
            blk = {}
            for model, r0 in r_ref.items():
                P.pop("g", None)
                if model.endswith("entropic"):
                    P["g"] = gref
                z0 = lnH(r0, P)
                Jn = nuis_jac(r0, P)
                per = {}
                for W in names:
                    if W == "REF" or Lk not in rd[W]:
                        continue
                    rw = rd[W][Lk]
                    if model == "kinetic":
                        r1 = r_kinetic(rw["rows"])
                    elif model.startswith("gk_static"):
                        r1 = (1 + 2j * np.pi * FREQ * 1e-12 * rw["taul_static"]) / \
                             (1 + 2j * np.pi * FREQ * 1e-12 * rw["tauR_static"])
                    else:
                        r1 = 1 / (1 + 2j * np.pi * FREQ * 1e-12 * rw["tauR_static"])
                    if model.endswith("entropic"):
                        P1 = dict(P)
                        P1["g"] = 1 + 2j * np.pi * FREQ * 1e-12 * rw["taul_static"]
                        dz = lnH(r1, P1) - z0
                    else:
                        dz = lnH(r1, P) - z0
                    per[W] = {f"{N}|{sp}": dchi2(dz, Jn[:, cols], sp) for N, cols in SETS.items() for sp in SIGPHI}
                    per[W]["max_abs_dphi_deg"] = float(np.degrees(np.max(np.abs(dz.imag))))
                summ = {}
                for N in SETS:
                    for sp in SIGPHI:
                        k = f"{N}|{sp}"
                        Wmax = max(per, key=lambda x: per[x][k])
                        summ[k] = {"max_dchi2": per[Wmax][k], "at": Wmax, "margin": 1 / np.sqrt(per[Wmax][k])}
                blk[model] = {"per_witness": per, "summary": summ}
            P.pop("g", None)
            out[f"{mesh}|{L}|{Rtag}"] = blk
            print(mesh, L, Rtag, {m: {k: round(v["margin"], 3) for k, v in blk[m]["summary"].items()
                                      if k.endswith("|0.01") or k.endswith("|0.1")} for m in blk}, flush=True)
(OUT / "ra5b_observable.json").write_text(json.dumps(out, indent=1))
