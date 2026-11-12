"""Audit, part 3: arbitrary-precision reference for H and d ln H / d ln tau_R at selected grid frequencies.

mpmath, 50 significant digits. Same physics as aud_model.jl but a third implementation: unscaled cosh/sinh
quadrupoles (no overflow in mpmath), tanh-sinh quadrature of H = 1/(4 pi eta) int_0^inf exp(-s) Z(s/eta) ds split
at s = 1e-3, 1e-2, 0.1, 1, 10; the tau_R derivative is taken under the integral by a central difference with
h = 1e-20 in ln tau_R (truncation ~1e-40, roundoff ~1e-30 at 50 digits).
Grid indices 0, 20, 40, 55, 65, 70, 75, 79 of the 80-point grid; thicknesses 30 and 500 nm.
Usage: python aud3_mp_check.py -> out/aud3_mp_check.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
INP = json.loads((OUT / "aud_inputs.json").read_text())
mp.mp.dps = 50
S = {k: mp.mpf(repr(v)) for k, v in INP["stack"].items()}
IDX = (0, 20, 40, 55, 65, 70, 75, 79)


def Z(k2, p, P):
    gs = mp.sqrt((P["CS"] * p + P["lamS"] * k2) / P["lamS"])
    Ys = P["lamS"] * gs
    g = 1 + P["tauR"] * p
    gn = mp.sqrt((P["C"] * p * g + P["lamr"] * k2) / P["lamz"])
    Yn = P["lamz"] * gn / g
    ga = mp.sqrt((P["CAu"] * p + P["lamAu"] * k2) / P["lamAu"])
    Ya = P["lamAu"] * ga

    def M(gam, Y, d):
        c, s = mp.cosh(gam * d), mp.sinh(gam * d)
        return mp.matrix([[c, s / Y], [Y * s, c]])

    def R(r):
        return mp.matrix([[1, r], [0, 1]])

    T = M(ga, Ya, P["dAu"]) * R(P["R1"]) * M(gn, Yn, P["dAlN"]) * R(P["R2"])
    return (T[0, 0] + T[0, 1] * Ys) / (T[1, 0] + T[1, 1] * Ys)


def run(Lkey):
    b = INP["base"][Lkey]
    P = dict(S)
    P.update({"tauR": mp.mpf(repr(b["tauR"])), "lamz": mp.mpf(repr(b["lamz"])), "lamr": mp.mpf(repr(b["lamr"])),
              "dAlN": mp.mpf(Lkey.split("_")[0]) * mp.mpf("1e-9")})
    eta = (P["wp"] ** 2 + P["wr"] ** 2) / 8
    h = mp.mpf("1e-20")
    Pp, Pm = dict(P), dict(P)
    Pp["tauR"] = P["tauR"] * mp.exp(h)
    Pm["tauR"] = P["tauR"] * mp.exp(-h)
    out = {}
    for i in IDX:
        f = mp.mpf(repr(INP["freq80"][i]))
        p = mp.mpc(0, 2 * mp.pi * f)
        pts = [0, mp.mpf("1e-3"), mp.mpf("1e-2"), mp.mpf("0.1"), 1, 10, mp.inf]
        t = time.perf_counter()
        Hv = mp.quad(lambda s: mp.exp(-s) * Z(s / eta, p, P), pts) / (4 * mp.pi * eta)
        dH = mp.quad(lambda s: mp.exp(-s) * (Z(s / eta, p, Pp) - Z(s / eta, p, Pm)) / (2 * h), pts) / (4 * mp.pi * eta)
        # error estimate of the quadrature: same integral at 10 more digits
        with mp.workdps(60):
            Hv2 = mp.quad(lambda s: mp.exp(-s) * Z(s / eta, p, P), pts) / (4 * mp.pi * eta)
        lnH = mp.log(Hv)
        J = dH / Hv
        out[str(i)] = {"f": float(f), "lnH": [mp.nstr(lnH.real, 30), mp.nstr(lnH.imag, 30)],
                       "Jtau": [mp.nstr(J.real, 30), mp.nstr(J.imag, 30)],
                       "quad_self_check_rel": float(abs(Hv2 / Hv - 1)), "seconds": time.perf_counter() - t}
        print(Lkey, i, float(f), mp.nstr(J, 20), "self-check", out[str(i)]["quad_self_check_rel"],
              f"{out[str(i)]['seconds']:.1f}s", flush=True)
    return out


res = {"dps": 50, "indices": IDX}
for Lk in ("500", "30"):
    res[Lk] = run(Lk)
    (OUT / "aud3_mp_check.json").write_text(json.dumps(res, indent=1))
print("done")
