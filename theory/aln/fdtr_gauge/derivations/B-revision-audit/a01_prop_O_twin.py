"""Audit of Proposition O (proofs/B-revision.md, Sec. 1.2): which temperature the twin reads, and where the
interface convention is applied.  Independent closed forms (reflection form of the theta-quadrupole), not the
Mobius-matrix code of b2_01.

Checks, at k = 0, against src/twin_fdtr.py (imported read-only):
 (1) bare GK film / Fourier substrate: twin T(0)  ==  (1/B) (1 + r E)/(1 - r E),  r = (B L - N)/(B L + N),
     N = 1 + tau_n p, tau_n = c_n C   (theta-jump at the bottom face only; nothing at the heated face);
 (2) twin theta_e(0) == (1 + tau_l p) * (1)  for every interface convention;
 (3) Au / R1 / GK film / R2 / sapphire:  twin == Mob(T_Au; R1 + (N/B)(1 + r E)/(1 - r E))  with N_t = N_b = N;
 (4) the alternative readings are excluded: (a) theta read at the heated face of the bare film,
     (b) convention applied at the bottom face only in the stack, (c) at the top face only;
 (5) S-family decoupling at k = 0: T(0) does not depend on the slip coefficient C_s (0, 1, inf);
 (6) defaults: observable='T' in kernel / gaussian_response / gaussian_response_adaptive / plane_response, and
     synthetic_measurement passes no observable.
Double precision (the twin is double precision).  Output: a01_prop_O_twin.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
import json
import inspect
import warnings
import numpy as np

sys.path.insert(0, r"C:\Users\Koussay\these\src")
import twin_fdtr as tw  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}


def film_reflection(lam, C, tR, tl, e, p, load, N_top, N_bot):
    """N_top * T(0)/q(0) of a film whose bottom sees theta_b = load * q with theta_b = N_bot T (reflection form)."""
    g = (1 + tR * p) / (1 + tl * p)
    a = lam / C
    A = e * np.sqrt(p * g / a)                  # xi sqrt(p g), xi = e / sqrt(a)
    B = np.sqrt(lam * C) * np.sqrt(p / g)       # b sqrt(p/g)
    u = B * load / N_bot
    r = (u - 1) / (u + 1)
    E = np.exp(-2 * A)
    return N_top / B * (1 + r * E) / (1 - r * E)


def fourier_T(lam, C, d, p):
    A = d * np.sqrt(p * C / lam)
    B = np.sqrt(lam * C) * np.sqrt(p)
    return ((np.cosh(A), np.sinh(A) / B), (B * np.sinh(A), np.cosh(A)))


def mob(M, x):
    (a, b_), (c, d) = M
    return (a * x + b_) / (c * x + d)


C_f, lam_f, e_f = 2.41e6, 188.6, 500e-9
sap = tw.Material(35.0, 3.03e6)
au = tw.Material(150.0, 2.49e6)
f = np.concatenate([np.logspace(4, np.log10(2e8), 25), np.logspace(9, 11, 7)])
p = 2j * np.pi * f
S = np.sqrt(sap.lam_z * sap.rho_c) * np.sqrt(p)

rows = []
worst = {"bare_T": 0.0, "bare_theta_e": 0.0, "stack": 0.0}
alt_min = {"bare_theta_read": np.inf, "stack_bottom_only": np.inf, "stack_top_only": np.inf}
slip_dev = 0.0
for alpha in (1.0 / 3.0, 2.0):
    for (tR, tl) in ((25e-12, 10e-12), (25e-12, 25e-12), (0.0, 20e-12), (30e-12, 0.0)):
        if tl == 0.0:
            ell = 0.0
            film = tw.Material(lam_f, C_f, tau_r=tR)          # Cattaneo: c_n must be 0 in every convention
        else:
            ell = np.sqrt(tl * lam_f / ((1 + alpha) * C_f))
            film = tw.Material(lam_f, C_f, tau_r=tR, ell=ell, alpha=alpha)
        w = ell ** 2 / lam_f
        bcs = {"local": tw.GKBoundary(temperature="local"),
               "entropic": tw.GKBoundary(temperature="entropic"),
               "kinetic+0.4": tw.GKBoundary(temperature="kinetic", c_n=0.4 * w if w else 1e-17),
               "kinetic-0.7": tw.GKBoundary(temperature="kinetic", c_n=-0.7 * w if w else -1e-17)}
        for name, bc in bcs.items():
            c_n, _ = bc.coefficients(film)
            tau_n = c_n * C_f
            N = 1 + tau_n * p
            for R in (0.0, 3e-8):
                st = tw.Stack((tw.Layer(film, e_f, R, bc),), sap)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    zT = tw.kernel(f, 0.0, st)
                    zth = tw.kernel(f, 0.0, st, observable="theta_e")
                cT = film_reflection(lam_f, C_f, tR, film.tau_l, e_f, p, R + 1 / S, 1.0, N)
                dT = float(np.max(np.abs(zT / cT - 1)))
                dth = float(np.max(np.abs(zth / ((1 + film.tau_l * p) * cT) - 1)))
                worst["bare_T"] = max(worst["bare_T"], dT)
                worst["bare_theta_e"] = max(worst["bare_theta_e"], dth)
                if tau_n != 0:
                    alt = float(np.max(np.abs(zT / (N * cT) - 1)))
                    alt_min["bare_theta_read"] = min(alt_min["bare_theta_read"], alt)
                rows.append({"alpha": alpha, "tR": tR, "tl": tl, "conv": name, "R": R, "tau_n/tau_l":
                             (tau_n / film.tau_l if film.tau_l else None), "bare_T_dev": dT, "bare_theta_e_dev": dth})
            # stack
            R1, R2 = 1e-8, 2e-8
            st = tw.Stack((tw.Layer(au, 80e-9, R1), tw.Layer(film, e_f, R2, bc)), sap)
            ztw = tw.kernel(f, 0.0, st)
            TA = fourier_T(au.lam_z, au.rho_c, 80e-9, p)
            both = mob(TA, R1 + film_reflection(lam_f, C_f, tR, film.tau_l, e_f, p, R2 + 1 / S, N, N))
            bot = mob(TA, R1 + film_reflection(lam_f, C_f, tR, film.tau_l, e_f, p, R2 + 1 / S, 1.0, N))
            top = mob(TA, R1 + film_reflection(lam_f, C_f, tR, film.tau_l, e_f, p, R2 + 1 / S, N, 1.0))
            ds = float(np.max(np.abs(ztw / both - 1)))
            worst["stack"] = max(worst["stack"], ds)
            if tau_n != 0:
                alt_min["stack_bottom_only"] = min(alt_min["stack_bottom_only"], float(np.max(np.abs(ztw / bot - 1))))
                alt_min["stack_top_only"] = min(alt_min["stack_top_only"], float(np.max(np.abs(ztw / top - 1))))
            # (5) slip coefficient must not matter at k = 0
            if ell > 0:
                for slip in (0.0, 1.0, np.inf):
                    bc2 = tw.GKBoundary(slip=slip, temperature=bc.temperature, c_n=bc.c_n, c_t=bc.c_t)
                    st2 = tw.Stack((tw.Layer(au, 80e-9, R1), tw.Layer(film, e_f, R2, bc2)), sap)
                    slip_dev = max(slip_dev, float(np.max(np.abs(tw.kernel(f, 0.0, st2) / ztw - 1))))
            rows.append({"alpha": alpha, "tR": tR, "tl": tl, "conv": name, "stack_dev": ds,
                         "tau_n": tau_n})

OUT["max_deviation_twin_vs_closed_form"] = worst
OUT["min_deviation_of_alternative_readings_when_c_n_nonzero"] = alt_min
OUT["max_change_of_stack_T0_with_slip_coefficient_at_k0"] = slip_dev
OUT["cases"] = rows

# (6) defaults
sig = {nm: inspect.signature(getattr(tw, nm)).parameters["observable"].default
       for nm in ("kernel", "gaussian_response", "gaussian_response_adaptive", "plane_response")}
src = inspect.getsource(tw.synthetic_measurement)
OUT["observable_defaults"] = sig
OUT["synthetic_measurement_mentions_observable"] = "observable" in src
OUT["cattaneo_coefficients_all_conventions"] = [
    tw.GKBoundary(temperature=t, c_n=(1e-17 if t == "kinetic" else 0.0)).coefficients(tw.Material(lam_f, C_f, tau_r=25e-12))
    for t in ("local", "entropic", "kinetic")]

print(json.dumps({k: OUT[k] for k in OUT if k != "cases"}, indent=1, default=str))
with open(os.path.join(HERE, "a01_prop_O_twin.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
