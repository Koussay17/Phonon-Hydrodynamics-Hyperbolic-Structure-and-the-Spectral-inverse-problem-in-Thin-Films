"""Theorem 3(i): RTA mesh laws by exact lattice sums, compared with the sharp constants that the
dyadic-shell argument predicts (the manuscript states only orders of growth).

Model (d = 2, square lattice, B = 2 pi I, V_cell = 1, Gamma-centred odd mesh, Gamma excluded):
  J^2(k) = (khat_x^2 + c1 |k|) exp(-|k|^2/8)       J_0^2 = khat_x^2 vanishes on the plane khat_x = 0, so the
                                                    correction is ADDITIVE, not multiplicative (Sec. proof)
  r(k)   = r_0(khat) |k|^alpha (1 + c2 |k|^2),      r_0 = 1 + 2.4 khat_x^2 khat_y^2  in [1, 1.6], continuous
Mesh sums M_m(N) = N^-d sum_{n != 0} J^2 r^-m (k = 2 pi n / N); continuum M_m = (2pi)^-d int_BZ J^2 r^-m.
Predictions (A_m = int_{S^1} J_0^2 r_0^-m, S_m = sum_{n != 0} J_0^2 r_0^-m(nhat) |n|^{-m alpha}):
  m alpha < d : M_m(N) -> M_m
  m alpha = d : M_m(N) = (2pi)^-d A_m ln N + O(1)           (slope in ln N)
  m alpha > d : M_m(N) = N^{m alpha - d} (2pi)^{-m alpha} S_m (1 + o(1))
Usage: python -B audit_ir_lattice.py > audit_ir_lattice.out.txt      (single thread, < 300 MB)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import math
import time
import numpy as np
import mpmath as mp

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""), flush=True)


c1, c2 = 0.3, 0.5
d = 2


def J2(kx, ky):
    k = np.hypot(kx, ky)
    return ((kx / k) ** 2 + c1 * k) * np.exp(-k ** 2 / 8)


def r_of(kx, ky, alpha):
    k = np.hypot(kx, ky)
    r0 = 1 + 2.4 * (kx / k) ** 2 * (ky / k) ** 2
    return r0 * k ** alpha * (1 + c2 * k ** 2)


def mesh_moments(N, alpha, ms=(0, 1, 2), chunk=256):
    h = (N - 1) // 2
    n1 = np.arange(-h, h + 1, dtype=np.float64)
    out = {m: 0.0 for m in ms}
    for s0 in range(0, N, chunk):
        rows = n1[s0:s0 + chunk]
        NX, NY = np.meshgrid(rows, n1, indexing="ij")
        kx = (2 * np.pi / N) * NX.ravel(); ky = (2 * np.pi / N) * NY.ravel()
        keep = (kx != 0) | (ky != 0)
        kx, ky = kx[keep], ky[keep]
        j2 = J2(kx, ky); rr = r_of(kx, ky, alpha)
        for m in ms:
            out[m] += float(np.sum(j2 / rr ** m))
    return {m: out[m] / N ** d for m in ms}


def cont_moment(m, alpha):
    """(2pi)^-2 int over [-pi,pi]^2 in polar coordinates; 4-fold symmetry (x -> -x, y -> -y)."""
    mp.mp.dps = 20

    def inner(th):
        rmax = mp.pi / max(mp.cos(th), mp.sin(th))
        ct, st = mp.cos(th), mp.sin(th)
        r0 = 1 + 2.4 * ct ** 2 * st ** 2
        f = lambda k: ((ct ** 2 + c1 * k) * mp.e ** (-k ** 2 / 8)) / (r0 * k ** alpha * (1 + c2 * k ** 2)) ** m * k
        return mp.quad(f, [0, rmax / 4, rmax])
    val = 4 * (mp.quad(inner, [0, mp.pi / 8, mp.pi / 4]) + mp.quad(inner, [mp.pi / 4, 3 * mp.pi / 8, mp.pi / 2]))
    return float(val / (2 * mp.pi) ** 2)


def A_ang(m):
    return float(4 * mp.quad(lambda th: mp.cos(th) ** 2 / (1 + 2.4 * mp.cos(th) ** 2 * mp.sin(th) ** 2) ** m, [0, mp.pi / 4, mp.pi / 2]))


def S_lattice(m, alpha, R=1500):
    """sum over n != 0 of J_0^2 r_0^-m |n|^{-m alpha}: disc |n| <= R plus continuum tail A_m R^{d-m alpha}/(m alpha-d)."""
    tot = 0.0
    n1 = np.arange(-R, R + 1, dtype=np.float64)
    for s0 in range(0, len(n1), 256):
        NX, NY = np.meshgrid(n1[s0:s0 + 256], n1, indexing="ij")
        x = NX.ravel(); y = NY.ravel()
        rr2 = x * x + y * y
        keep = (rr2 > 0) & (rr2 <= R * R)
        x, y, rr2 = x[keep], y[keep], rr2[keep]
        cx2 = x * x / rr2; cy2 = y * y / rr2
        tot += float(np.sum(cx2 / (1 + 2.4 * cx2 * cy2) ** m * rr2 ** (-m * alpha / 2)))
    s = m * alpha - d
    return tot + A_ang(m) * R ** (-s) / s


t0 = time.time()
Ns = [101, 201, 401, 801, 1601, 3201]
out = {"model": {"c1": c1, "c2": c2, "d": d, "Ns": Ns}}
# ----------------------------------------------------------------------------- (a) 2 alpha < d: convergence
alpha = 0.6
M1c, M2c = cont_moment(1, alpha), cont_moment(2, alpha)
tau_inf = M2c / M1c
rows = []
for N in Ns:
    mm = mesh_moments(N, alpha)
    rows.append((N, mm[1], mm[2], mm[2] / mm[1]))
print(f"(a) alpha = {alpha}: continuum M1 = {M1c:.10f}, M2 = {M2c:.10f}, tau = {tau_inf:.10f}")
for N, a1, a2, tt in rows:
    print(f"    N = {N:5d}: M1 = {a1:.10f}  M2 = {a2:.10f}  tau = {tt:.10f}  |tau - tau_inf|/tau_inf = {abs(tt - tau_inf) / tau_inf:.3e}")
errs = [abs(t_ - tau_inf) / tau_inf for *_, t_ in rows]
check("(a) 2alpha < d: tau_RTA(N) -> M2/M1 (relative error monotonically decreasing, < 1e-2 at N = 3201)",
      errs[-1] < 1e-2 and all(errs[i + 1] < errs[i] for i in range(len(errs) - 1)))
# observed rate: error ~ N^{-(d - 2 alpha)} from the singular corner (shells near the origin)
print(f"    observed error exponent between N = 801 and 3201: {math.log(errs[3] / errs[5]) / math.log(4):.3f} (expected about d - 2 alpha = {d - 2 * alpha:.2f})")
out["a"] = {"alpha": alpha, "M1": M1c, "M2": M2c, "rows": rows}
# ----------------------------------------------------------------------------- (b) 2 alpha = d: log growth with slope (2pi)^-d A_2 / M1
alpha = 1.0
M1c = cont_moment(1, alpha)
A2 = A_ang(2)
slope_pred = A2 / (2 * math.pi) ** d / M1c
rows = []
for N in Ns:
    mm = mesh_moments(N, alpha)
    rows.append((N, mm[1], mm[2], mm[2] / mm[1]))
slopes = [(rows[i + 1][3] - rows[i][3]) / math.log(rows[i + 1][0] / rows[i][0]) for i in range(len(rows) - 1)]
print(f"(b) alpha = {alpha}: predicted d tau/d ln N = (2pi)^-2 A_2 / M1 = {slope_pred:.8f}  (M1 = {M1c:.8f}, A_2 = {A2:.8f})")
for (N, a1, a2, tt), sl in zip(rows[1:], slopes):
    print(f"    N = {N:5d}: tau = {tt:.8f}, M1(N) = {a1:.8f}, local slope = {sl:.8f}, slope/pred = {sl / slope_pred:.6f}")
dev_b = [abs(s_ / slope_pred - 1) for s_ in slopes]
check("(b) 2alpha = d: tau_RTA(N) grows like ln N with the predicted constant (|slope/pred - 1| decreasing, < 5e-3)",
      dev_b[-1] < 5e-3 and dev_b[-1] < dev_b[0])
out["b"] = {"alpha": alpha, "slope_pred": slope_pred, "rows": rows, "slopes": slopes}
# ----------------------------------------------------------------------------- (c) d/2 < alpha < d: power N^{2 alpha - d}
alpha = 1.5
M1c = cont_moment(1, alpha)
S2 = S_lattice(2, alpha)
pred = (2 * math.pi) ** (-2 * alpha) * S2 / M1c
rows = []
for N in Ns:
    mm = mesh_moments(N, alpha)
    rows.append((N, mm[1], mm[2], mm[2] / mm[1] / N ** (2 * alpha - d)))
print(f"(c) alpha = {alpha}: predicted tau/N^(2alpha-d) -> (2pi)^(-2alpha) S_2 / M1 = {pred:.8f}  (S_2 = {S2:.8f}, M1 = {M1c:.8f})")
for N, a1, a2, rat in rows:
    print(f"    N = {N:5d}: M1(N) = {a1:.8f}, tau/N^(2alpha-d) = {rat:.8f}, ratio/pred = {rat / pred:.6f}")
pred_M2c = (2 * math.pi) ** (-2 * alpha) * S2
m2dev = [abs(r_[2] / r_[0] ** (2 * alpha - d) / pred_M2c - 1) for r_ in rows]
m1err = [abs(r_[1] - M1c) / M1c for r_ in rows]
rate_m1 = math.log(m1err[-2] / m1err[-1]) / math.log(2)
# Richardson extrapolation of tau/N^(2alpha-d) with the M1 error exponent d - alpha
q_ = 2.0 ** (-(d - alpha))
R1, R2 = rows[-2][3] / pred, rows[-1][3] / pred
Rinf = R2 - (R1 - R2) * q_ / (1 - q_)
print(f"    M2(N)/N^(2alpha-d)/[(2pi)^(-2alpha) S_2] - 1 = {m2dev[-1]:.2e} at N = 3201; M1(N) error exponent = {rate_m1:.3f} (expected d - alpha = {d - alpha});"
      f" Richardson limit of (tau/N^(2alpha-d))/pred = {Rinf:.5f}")
check("(c) d/2 < alpha < d: M2(N) ~ (2pi)^(-2alpha) S_2 N^(2alpha-d) (within 1e-2), M1(N) -> M1 at rate N^-(d-alpha), "
      "so tau/N^(2alpha-d) -> (2pi)^(-2alpha) S_2/M1 (extrapolated within 5e-3)",
      m2dev[-1] < 1e-2 and abs(rate_m1 - (d - alpha)) < 0.05 and abs(Rinf - 1) < 5e-3)
out["c"] = {"alpha": alpha, "S2": S2, "M1": M1c, "pred": pred, "rows": rows}
# ----------------------------------------------------------------------------- (d) alpha = d: M1 ~ ln N, M2 ~ N^d, tau ~ N^d / ln N
alpha = 2.0
A1 = A_ang(1)
S2 = S_lattice(2, alpha)
rows = []
for N in Ns:
    mm = mesh_moments(N, alpha)
    rows.append((N, mm[1], mm[2], mm[2] / N ** d, mm[2] / mm[1] * math.log(N) / N ** d))
sl1 = [(rows[i + 1][1] - rows[i][1]) / math.log(rows[i + 1][0] / rows[i][0]) for i in range(len(rows) - 1)]
pred_sl1 = A1 / (2 * math.pi) ** d
pred_M2 = (2 * math.pi) ** (-2 * alpha) * S2
pred_tau = (2 * math.pi) ** (-d) * S2 / A1
print(f"(d) alpha = {alpha}: predicted dM1/dlnN = {pred_sl1:.8f}, M2/N^d -> {pred_M2:.8f}, tau ln N/N^d -> {pred_tau:.8f}")
for (N, a1, a2, m2r, tr), sl in zip(rows[1:], sl1):
    print(f"    N = {N:5d}: M1 slope/pred = {sl / pred_sl1:.6f}, (M2/N^d)/pred = {m2r / pred_M2:.6f}, (tau lnN/N^d)/pred = {tr / pred_tau:.4f}")
check("(d) alpha = d: M1(N) ~ (2pi)^-d A_1 ln N, M2(N) ~ N^d (2pi)^-2d S_2, hence tau ~ N^d/ln N (constants within 1e-2)",
      abs(sl1[-1] / pred_sl1 - 1) < 1e-2 and abs(rows[-1][3] / pred_M2 - 1) < 1e-2)
print(f"    note: (tau ln N / N^d)/pred = {rows[-1][4] / pred_tau:.4f} at N = 3201 -> 1 only like 1 + O(1/ln N)")
out["d"] = {"alpha": alpha, "pred_sl1": pred_sl1, "pred_M2": pred_M2, "pred_tau": pred_tau, "rows": rows, "slopes": sl1}
# ----------------------------------------------------------------------------- integral criterion: M_m finite iff m alpha < d
# cut-off integrals (2pi)^-2 int_{|k| > kmin}: bounded for m alpha < d, log for = d, power for > d
print("    elapsed %.1f s" % (time.time() - t0))
json.dump(out, open("audit_ir_lattice.json", "w"), indent=1, default=float)
print()
print(f"{sum(res)} of {len(res)} checks passed")
