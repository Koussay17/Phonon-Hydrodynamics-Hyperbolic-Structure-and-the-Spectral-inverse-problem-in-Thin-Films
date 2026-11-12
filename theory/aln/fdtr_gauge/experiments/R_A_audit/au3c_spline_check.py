"""AU3c: is the N3 amplitude-channel distinguishability of AU3b physical or an interpolation artefact?
For REF and cal1nz_m3L100min (500 nm, favourable stack): (i) smoothness of the spline of delta_eff(omega) between the
12 computed frequencies (leave-one-out prediction error at interior nodes); (ii) split of the N3 Delta chi^2 into the
amplitude and phase channels; (iii) the same test with the spline replaced by a monotone (PCHIP) interpolant.
Output: out/au3c_spline_check.json"""
import json, os, sys
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import au5_fisher1d as a5
d = json.loads((HERE / "out" / "au3b_exact_root775.json").read_text())
L = 500
res = {}
def deff_nodes(nm):
    r = d[nm][str(L)]
    f = np.array([x["f"] for x in r["rows"]])
    lam = np.array([complex(*x["lam_eff_rel"]) for x in r["rows"]])
    de = -(lam - 1) / (2j * np.pi * f)
    tR, tl = r["tauR_ps"] * 1e-12, r["taul_ps"] * 1e-12
    fx = np.concatenate([[1e3, 1e4, 1e5], f])
    dx = np.concatenate([[tR - tl] * 3, de])
    return fx, dx, tR, tl
def ratio(fx, dx, kind):
    I = CubicSpline if kind == "cubic" else PchipInterpolator
    sre, sim = I(np.log(fx), dx.real), I(np.log(fx), dx.imag)
    dI = sre(np.log(a5.FREQ)) + 1j * sim(np.log(a5.FREQ))
    return 1 - 2j * np.pi * a5.FREQ * dI
for nm in ("REF", "cal1nz_m3L100min"):
    fx, dx, _, _ = deff_nodes(nm)
    loo = []
    for i in range(4, len(fx) - 1):
        m = np.ones(len(fx), bool); m[i] = False
        s = CubicSpline(np.log(fx[m]), dx[m].real)(np.log(fx[i])) + 1j * CubicSpline(np.log(fx[m]), dx[m].imag)(np.log(fx[i]))
        loo.append(float(abs(s - dx[i]) / abs(dx[i])))
    res[f"{nm}_leave_one_out_rel_err"] = loo
b = a5.A1[str(L)]
r0 = d["REF"][str(L)]
for Rtag, (R1, R2) in {"favourable": (2e-9, 2e-9), "note16": (1e-8, 2e-8)}.items():
    base = {"tauR": r0["tauR_ps"] * 1e-12, "taul": r0["taul_ps"] * 1e-12, "lamz": b["lam_z"], "C": a5.C_ALN, "R1": R1, "R2": R2, "d": L * 1e-9}
    for kind in ("cubic", "pchip"):
        fx0, dx0, _, _ = deff_nodes("REF")
        rr = ratio(fx0, dx0, kind)
        fx1, dx1, tR1, tl1 = deff_nodes("cal1nz_m3L100min")
        rw = ratio(fx1, dx1, kind)
        bw = dict(base); bw["tauR"], bw["taul"] = tR1, tl1
        y = a5.lnH(np.zeros(9), bw, rw)
        for sp in (0.01, 0.1):
            chi, th = a5.gn_refit(y, base, a5.SETS["N3"], sp, ratio=rr)
            z = a5.lnH(th, base, rr)
            amp = float(np.sum(((y - z).real / a5.SIG_A) ** 2))
            ph = float(np.sum(((y - z).imag / np.deg2rad(sp)) ** 2))
            res[f"{Rtag}|{kind}|{sp}"] = {"dchi2": chi, "amplitude_part": amp, "phase_part": ph,
                                          "max_abs_dlnamp": float(np.max(np.abs((y - z).real))),
                                          "max_abs_dphase_deg": float(np.degrees(np.max(np.abs((y - z).imag))))}
print(json.dumps(res, indent=1))
(HERE / "out" / "au3c_spline_check.json").write_text(json.dumps(res, indent=1))
