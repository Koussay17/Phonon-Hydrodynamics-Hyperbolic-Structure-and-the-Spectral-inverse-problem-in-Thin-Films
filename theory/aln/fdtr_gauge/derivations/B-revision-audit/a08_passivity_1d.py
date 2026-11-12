"""Probe of the remark in proofs/B-revision.md Sec. 7, B-2(iv): 'the entropic one is the one with a sign-definite 1D
dissipation identity'.

At k = 0 the film obeys q = -(lam/g) T' and C p T + q' = 0, so
    conj(T(0)) q(0) = C p int|T|^2 + conj(T(e)) q(e) + (conj(g)/lam) int|q|^2,
and with the local jump T(e) = L q(e):  Re Z_T = Re(L)|q(e)|^2 + (Re g/lam) int|q|^2 > 0  (Re g > 0 on i R).
With the entropic variable theta = (1 + tl p) T and the jump theta(e) = L q(e) (eps = 1) the analogous identity bounds
Re[(1 + tl p) Z_T].  For other eps both identities carry an interface term without sign.
Check (double precision, nondimensional, b_s = 1): sample random (xi, b, tR, tl, R, omega) and report
min Re Z_T / |Z_T| (read-out T(0)) and min Re[(1 + tl p) Z_T] / |.| for eps in {0, 1, 0.3, -0.5, 2}.
Output: a08_passivity_1d.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261112)
n = 400000
xi = np.exp(rng.uniform(np.log(0.05), np.log(5), n))
b = np.exp(rng.uniform(np.log(0.1), np.log(10), n))
tR = np.exp(rng.uniform(np.log(1e-3), np.log(1e2), n)) * (rng.uniform(size=n) > 0.1)
tl = np.exp(rng.uniform(np.log(1e-3), np.log(1e2), n)) * (rng.uniform(size=n) > 0.1)
R = np.exp(rng.uniform(np.log(1e-3), np.log(10), n)) * (rng.uniform(size=n) > 0.3)
w = np.exp(rng.uniform(np.log(1e-3), np.log(1e3), n))
p = 1j * w
g = (1 + tR * p) / (1 + tl * p)
A = xi * np.sqrt(p * g)
B = b * np.sqrt(p / g)
L = R + 1 / np.sqrt(p)
E = np.exp(-2 * A)
OUT = {}
for eps in (0.0, 1.0, 0.3, -0.5, 2.0):
    N = 1 + eps * tl * p
    u = B * L / N
    r = (u - 1) / (u + 1)
    Z = (1 + r * E) / (B * (1 - r * E))
    th = (1 + tl * p) * Z
    ok = np.isfinite(Z) & np.isfinite(th)
    mT = (Z.real / np.abs(Z))[ok]
    mth = (th.real / np.abs(th))[ok]
    OUT[f"eps={eps}"] = {"min Re Z_T/|Z_T|": float(mT.min()), "fraction Re Z_T < 0": float((mT < 0).mean()),
                         "min Re theta_e/|theta_e|": float(mth.min()), "fraction Re theta_e < 0": float((mth < 0).mean())}
print(json.dumps(OUT, indent=1))
with open(os.path.join(HERE, "a08_passivity_1d.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)

# ---------------------------------------------------------------- stack: Fourier transducer + R1 + film (same eps at both faces)
xiA = np.exp(rng.uniform(np.log(0.02), np.log(2), n))
bA = np.exp(rng.uniform(np.log(0.2), np.log(5), n))
R1 = np.exp(rng.uniform(np.log(1e-3), np.log(10), n)) * (rng.uniform(size=n) > 0.3)
sq = np.sqrt(p)
EA = np.exp(-2 * xiA * sq)
thA = (1 - EA) / (1 + EA)
BA = bA * sq
OUT_S = {}
for eps in (0.0, 1.0, 0.3, -0.5, 2.0):
    N = 1 + eps * tl * p
    u = B * L / N
    r = (u - 1) / (u + 1)
    Phi = N * (1 + r * E) / (B * (1 - r * E))
    W = R1 + Phi
    Zt = (W + thA / BA) / (1 + BA * thA * W)
    ok = np.isfinite(Zt)
    m = (Zt.real / np.abs(Zt))[ok]
    OUT_S[f"eps={eps}"] = {"min Re Z_tot/|Z_tot|": float(m.min()), "fraction Re Z_tot < 0": float((m < 0).mean())}
print(json.dumps(OUT_S, indent=1))
OUT["stack"] = OUT_S
with open(os.path.join(HERE, "a08_passivity_1d.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
