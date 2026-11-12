"""B-2 / B-4, step 1: which temperature the twin observes, and the k = 0 closed forms under every interface convention.

Claims tested against src/twin_fdtr.py (read-only import):
 (a) The observable of `kernel`, `gaussian_response`, `gaussian_response_adaptive` and `plane_response` defaults to
     T(0) (energy-density temperature at the heated face); `synthetic_measurement` calls `gaussian_response` without an
     `observable` argument, i.e. always with T(0). The alternative output `theta_e` is the entropic temperature
     (1 + tau_l p) T(0) at k = 0 when the top layer is GK.
 (b) Bare GK film on a Fourier substrate, k = 0, convention c_n (local 0, entropic (1+alpha) ell^2/lam, kinetic given):
        Z_T = Mobius(M; L/N),   N = 1 + tau_n p,   tau_n = c_n C,   L = R + 1/S,
     i.e. the film sees the load L/(1 + tau_n p); theta_e(0) = (1 + tau_l p) Z_T.
 (c) Known Fourier transducer on top (Au/R1/GK film/R2/sapphire), k = 0:
        Z_tot = Mobius(T_Au; R1 + N * Mobius(M; L2/N)),   L2 = R2 + 1/S,
     i.e. the film enters through the theta-quadrupole N M N^-1 (theta-jumps at both faces, same c_n).
Double precision (the twin is double precision). Output: b2_01_twin_observable.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "2"
import sys
import json
import inspect
import warnings
import numpy as np

sys.path.insert(0, r"C:\Users\Koussay\these\src")
import twin_fdtr as tw  # noqa: E402

OUT = {}

# ---------------------------------------------------------------- (a) defaults of the observable
sig = {}
for name in ("kernel", "gaussian_response", "gaussian_response_adaptive", "plane_response"):
    sig[name] = inspect.signature(getattr(tw, name)).parameters["observable"].default
src = inspect.getsource(tw.synthetic_measurement)
call = [ln.strip() for ln in src.splitlines() if "gaussian_response(" in ln]
OUT["a_observable_defaults"] = {"defaults": sig,
                                "synthetic_measurement_call": call,
                                "word_observable_in_synthetic_measurement_source": ("observable" in src)}
src_k = inspect.getsource(tw.solve_kernel)
OUT["a_impedance_is_T_component"] = ("T0 = T0 + amps[0][..., i] * comp[0]" in src_k)


# ---------------------------------------------------------------- closed forms (k = 0)
def mobius(M, load):
    (a, b_), (c, d) = M
    return (a * load + b_) / (c * load + d)


def quad(m, thick, p):
    """1D lambda_eff quadrupole of a Fourier/Cattaneo/GK layer at k = 0: [[cosh, sinh/B], [B sinh, cosh]]."""
    g = (1 + m.tau_r * p) / (1 + m.tau_l * p) if m.law == "gk" else (1 + m.tau_r * p)
    lam_eff = m.lam_z / g
    sig_ = np.sqrt(m.rho_c * p / lam_eff)
    A = sig_ * thick
    B = lam_eff * sig_
    return ((np.cosh(A), np.sinh(A) / B), (B * np.sinh(A), np.cosh(A)))


def tau_n_of(m, bc):
    c_n, _ = bc.coefficients(m)
    return c_n * m.rho_c


def bare_closed(f, film, e, sub, R, bc):
    p = 2j * np.pi * f
    S = np.sqrt(sub.lam_z * sub.rho_c) * np.sqrt(p)
    L = R + 1 / S
    N = 1 + tau_n_of(film, bc) * p
    ZT = mobius(quad(film, e, p), L / N)
    theta_e = (1 + film.tau_l * p) * ZT
    return ZT, theta_e


def stack_closed(f, au, d_au, R1, film, e, R2, sub, bc):
    p = 2j * np.pi * f
    S = np.sqrt(sub.lam_z * sub.rho_c) * np.sqrt(p)
    L2 = R2 + 1 / S
    N = 1 + tau_n_of(film, bc) * p
    W = R1 + N * mobius(quad(film, e, p), L2 / N)
    return mobius(quad(au, d_au, p), W)


# ---------------------------------------------------------------- (b), (c): comparisons
C_aln, lam, e = 2.41e6, 188.6, 500e-9
sap = tw.Material(35.0, 3.03e6)
au = tw.Material(150.0, 2.49e6)
f = np.concatenate([np.logspace(4, np.log10(2e8), 40), np.logspace(9, 11, 9)])
res_b, res_c, sep = {}, {}, {}
for alpha in (1.0 / 3.0, 2.0):
    for (tR, tl) in ((25e-12, 10e-12), (25e-12, 25e-12), (0.0, 20e-12)):
        ell = np.sqrt(tl * lam / ((1 + alpha) * C_aln))
        film = tw.Material(lam, C_aln, tau_r=tR, ell=ell, alpha=alpha)
        w = ell ** 2 / lam
        bcs = {"local": tw.GKBoundary(temperature="local"),
               "entropic": tw.GKBoundary(temperature="entropic"),
               "kinetic_k=+0.3": tw.GKBoundary(temperature="kinetic", c_n=0.3 * w),
               "kinetic_k=-0.5": tw.GKBoundary(temperature="kinetic", c_n=-0.5 * w),
               "kinetic_k=+2.0": tw.GKBoundary(temperature="kinetic", c_n=2.0 * w),
               "beardo_Si311K": tw.GKBoundary.beardo(2.71 / 1068e6, -21e-9, -31e-9, -16e-9)}
        for name, bc in bcs.items():
            key = f"alpha={alpha:.3g},tR={tR:.0e},tl={tl:.0e},{name}"
            for R in (0.0, 2e-8):
                st = tw.Stack((tw.Layer(film, e, R, bc),), sap)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    zT = tw.kernel(f, 0.0, st)
                    zth = tw.kernel(f, 0.0, st, observable="theta_e")
                cT, cth = bare_closed(f, film, e, sap, R, bc)
                res_b[key + f",R={R:.0e}"] = {
                    "tau_n/tau_l": float(tau_n_of(film, bc) / film.tau_l),
                    "T_vs_ZT": float(np.max(np.abs(zT / cT - 1))),
                    "theta_e_vs_(1+tl p)ZT": float(np.max(np.abs(zth / cth - 1))),
                    # the alternative reading 'observable = theta at the front with the interface coefficient'
                    "T_vs_(1+tau_n p)ZT": float(np.max(np.abs(zT / ((1 + tau_n_of(film, bc) * 2j * np.pi * f) * cT) - 1))),
                }
            st = tw.Stack((tw.Layer(au, 80e-9, 1e-8), tw.Layer(film, e, 2e-8, bc)), sap)
            ztw = tw.kernel(f, 0.0, st)
            zc = stack_closed(f, au, 80e-9, 1e-8, film, e, 2e-8, sap, bc)
            # alternative: interface convention applied only at the bottom face (load L2/N, T at the film top)
            p = 2j * np.pi * f
            S = np.sqrt(sap.lam_z * sap.rho_c) * np.sqrt(p)
            N = 1 + tau_n_of(film, bc) * p
            z_alt = mobius(quad(au, 80e-9, p), 1e-8 + mobius(quad(film, e, p), (2e-8 + 1 / S) / N))
            res_c[key] = {"twin_vs_theta_quadrupole": float(np.max(np.abs(ztw / zc - 1))),
                          "twin_vs_bottom_only": float(np.max(np.abs(ztw / z_alt - 1)))}
OUT["b_bare_film_k0"] = res_b
OUT["c_transducer_stack_k0"] = res_c
OUT["summary"] = {
    "max_dev_bare_T": max(v["T_vs_ZT"] for v in res_b.values()),
    "max_dev_bare_theta_e": max(v["theta_e_vs_(1+tl p)ZT"] for v in res_b.values()),
    "min_dev_bare_alternative_nonlocal": min(v["T_vs_(1+tau_n p)ZT"] for k, v in res_b.items()
                                             if v["tau_n/tau_l"] != 0),
    "max_dev_stack": max(v["twin_vs_theta_quadrupole"] for v in res_c.values()),
    "min_dev_stack_alternative_nonlocal": min(v["twin_vs_bottom_only"] for k, v in res_c.items()
                                              if "local" not in k.split(",")[-1]),
}
print(json.dumps(OUT["a_observable_defaults"], indent=1))
print("impedance built from the T component:", OUT["a_impedance_is_T_component"])
print(json.dumps(OUT["summary"], indent=1))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_01_twin_observable.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
