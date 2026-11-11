"""Workstream B: two-frequency local identifiability.

(The low-frequency leading-order check of Theorem N2 is in 08_low_frequency_two_point.py.)
(B) Scan: for random theta (log-uniform) and on the boundary strata (b = b_s, tau_l = 0, tau_R = 0),
    the largest inverse condition number of the 4x4 two-frequency Jacobian over a 31x31 grid of
    frequency pairs. A value far above double rounding (1e-16) means N = 2 gives rank 4 at that theta.
(C) Conditioning near the resonance: smallest singular value of the function-level Jacobian versus
    the relative distance delta = (tau_l - tau_R)/tau_R.
Run:  python 06_two_frequency_scan.py > 06_two_frequency_scan.out
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import mpmath as mp
from model import jac_analytic, real_jacobian

mp.mp.dps = 80
bs = mp.mpf(1)

print("\n=== (B) two-frequency scan (double precision) ===")
BS = 1.0


def J_np(w, xi, b, tR, tl):
    p = 1j*w
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*np.sqrt(p*g); B = b*np.sqrt(p/g); S = BS*np.sqrt(p)
    L = 1/S
    r = (B*L - 1)/(B*L + 1); E = np.exp(-2*A); den = 1 - r**2*E**2
    ul = -4*A*r*E/den
    vl = -1 + 2*E*(2*B*L/(B*L + 1)**2)/den
    hR, hl = p/(1 + tR*p), p/(1 + tl*p)
    W = ul - vl
    sc = [xi, b, tR if tR > 0 else 1.0, tl if tl > 0 else 1.0]
    J = np.array([ul/xi, vl/b, W*hR/2, -W*hl/2])
    return J*np.array(sc)


grid = np.logspace(-2, 2, 31)


def best_invcond(th):
    Js = np.array([J_np(w, *th) for w in grid])          # (31, 4)
    R = np.stack([Js.real, Js.imag], axis=1)              # (31, 2, 4)
    best = 0.0
    for i in range(len(grid)):
        for j in range(i + 1, len(grid)):
            M = np.vstack([R[i], R[j]])
            s = np.linalg.svd(M, compute_uv=False)
            best = max(best, s[-1]/s[0])
    return best


rng = np.random.default_rng(7)
vals = []
for k in range(150):
    th = [10**rng.uniform(-0.5, 0.5), 10**rng.uniform(-0.7, 0.7), 10**rng.uniform(-1.5, 1), 10**rng.uniform(-1.5, 1)]
    vals.append((best_invcond(th), th))
vals.sort(key=lambda t: t[0])
print("  random interior points: min over theta of (max over pairs of sigma_min/sigma_max) =", f"{vals[0][0]:.3e}", "at", np.round(vals[0][1], 4))
print("  quartiles:", [f"{vals[int(q*(len(vals)-1))][0]:.2e}" for q in (0, 0.25, 0.5, 0.75, 1)])
for name, gen in [("b = b_s", lambda: [10**rng.uniform(-0.5, 0.5), 1.0, 10**rng.uniform(-1.5, 1), 10**rng.uniform(-1.5, 1)]),
                  ("tau_l = 0", lambda: [10**rng.uniform(-0.5, 0.5), 10**rng.uniform(-0.7, 0.7), 10**rng.uniform(-1.5, 1), 0.0]),
                  ("tau_R = 0", lambda: [10**rng.uniform(-0.5, 0.5), 10**rng.uniform(-0.7, 0.7), 0.0, 10**rng.uniform(-1.5, 1)]),
                  ("b = b_s, tau_l = 0", lambda: [10**rng.uniform(-0.5, 0.5), 1.0, 10**rng.uniform(-1.5, 1), 0.0]),
                  ("b = b_s, tau_R = 0", lambda: [10**rng.uniform(-0.5, 0.5), 1.0, 0.0, 10**rng.uniform(-1.5, 1)])]:
    v = sorted([(best_invcond(th), th) for th in [gen() for _ in range(60)]], key=lambda t: t[0])
    print(f"  stratum {name}: min over 60 points = {v[0][0]:.3e} at {np.round(v[0][1], 4)}; median = {v[30][0]:.2e}")
print("  resonance (control): ", f"{best_invcond([1.0, 2.0, 0.5, 0.5]):.3e}")

print("\n=== (C) conditioning near the resonance (function level, 61 frequencies) ===")
mp.mp.dps = 40
omegas = [mp.mpf(10)**(mp.mpf(k)/10) for k in range(-30, 31)]
for dlt in ['1e-1', '1e-2', '1e-3', '1e-4', '1e-5']:
    d = mp.mpf(dlt)
    th = [mp.mpf(1), mp.mpf(2), mp.mpf('0.5'), mp.mpf('0.5')*(1 + d)]
    Jm = real_jacobian(omegas, th, bs)
    U, Sv, V = mp.svd_r(Jm)
    print(f"  delta = {dlt}: sigma_min/sigma_max = {mp.nstr(Sv[3]/Sv[0], 6)}; ratio to delta = {mp.nstr(Sv[3]/Sv[0]/d, 6)}")
