"""Independent checks of the mesh laws of Theorem 3(i) (Supplemental Sec. S4) and of Table S2.

Part A recomputes the entries of Table S2 for the paper's test model by routes that differ from
checks/audit/audit_ir_lattice.py:
  * A_1 and A_2 in closed form (A_m = (1/2) int_0^{2pi} (1 + a sin^2 2t)^-m dt, a = 0.6);
  * continuum moments by QUADPACK with the algebraic endpoint weight k^{1 - m alpha} (scipy quad, weight='alg');
  * S_2 by a smoothed cutoff and Poisson summation (super-algebraic convergence in the cutoff radius), instead of
    a sharp cutoff with a continuum tail;
  * mesh sums folded onto one quadrant (the model is even in k_x and in k_y).
Part B tests the mesh laws on a hexagonal lattice with angular functions that are discontinuous on a closed
null set of directions which contains lattice directions (the situation that (H5) is meant to admit), with the
values on the discontinuity directions fixed by exact integer tests. It also shows that, for m alpha > d, the
constant S_m depends on the values of Phi_m on those lattice directions (a null set).
Usage (from this directory): python -B audit_s4_mesh_laws.py > audit_s4_mesh_laws.out.txt
Single thread; < 400 MB.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import json
import math
import time
import numpy as np
from scipy import integrate

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""), flush=True)


t0 = time.time()
TP = 2 * math.pi
Ns = [101, 201, 401, 801, 1601, 3201]
out = {}

# =============================================================================== Part A: the paper's test model
c1, c2, aa = 0.3, 0.5, 2.4
A1_cf = math.pi / math.sqrt(1 + 0.6)                      # (1/2) 2pi/sqrt(1+a)
A2_cf = 0.5 * math.pi * (2 + 0.6) / (1 + 0.6) ** 1.5      # (1/2) pi (2+a)/(1+a)^{3/2}
A1_q = integrate.quad(lambda t: math.cos(t) ** 2 / (1 + aa * math.cos(t) ** 2 * math.sin(t) ** 2), 0, TP, limit=200)[0]
A2_q = integrate.quad(lambda t: math.cos(t) ** 2 / (1 + aa * math.cos(t) ** 2 * math.sin(t) ** 2) ** 2, 0, TP, limit=200)[0]
print(f"A: A_1 = {A1_cf:.12f} (closed form) vs {A1_q:.12f} (quad); A_2 = {A2_cf:.12f} vs {A2_q:.12f}")
check("A: closed forms A_1 = pi/sqrt(1.6), A_2 = 1.3 pi/1.6^1.5 agree with quadrature",
      abs(A1_cf - A1_q) < 1e-11 and abs(A2_cf - A2_q) < 1e-11)


def cont_moment_A(m, alpha):
    """(2pi)^-2 int_{[-pi,pi]^2} J^2 r^-m d^2k = 4 (2pi)^-2 int_0^{pi/2} dtheta int_0^{K(theta)} ... k dk."""
    if m == 0:
        alpha = 0.0

    def radial(th):
        ct, st = math.cos(th), math.sin(th)
        K = math.pi / max(ct, st)
        r0 = 1 + aa * ct * ct * st * st
        f = lambda k: (ct * ct + c1 * k) * math.exp(-k * k / 8) / (r0 * (1 + c2 * k * k)) ** m
        val = integrate.quad(f, 0, K, weight="alg", wvar=(1 - m * alpha, 0.0), epsabs=1e-14, epsrel=1e-12, limit=200)[0]
        return val

    tot = integrate.quad(radial, 0, math.pi / 4, epsabs=1e-14, epsrel=1e-12, limit=200)[0] \
        + integrate.quad(radial, math.pi / 4, math.pi / 2, epsabs=1e-14, epsrel=1e-12, limit=200)[0]
    return 4 * tot / TP ** 2


def smooth_step(t):
    """C-infinity step: 0 for t <= 1/2, 1 for t >= 1."""
    x = np.clip((t - 0.5) / 0.5, 0.0, 1.0)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        a = np.where(x > 0, np.exp(-1 / np.where(x > 0, x, 1)), 0.0)
        b = np.where(x < 1, np.exp(-1 / np.where(x < 1, 1 - x, 1)), 0.0)
    return a / (a + b)


def S_smooth(Phi_xy, A_Phi, s, R):
    """sum_{n != 0} Phi(nhat) |n|^-s by a smoothed cutoff: sum of Phi |n|^-s (1 - g(|n|/R)) over |n| < R plus
    A_Phi R^{2-s} int_0^inf t^{1-s} g(t) dt (Poisson summation; the remainder is super-algebraically small in R)."""
    L = int(R) + 1
    n1 = np.arange(-L, L + 1, dtype=np.float64)
    NX, NY = np.meshgrid(n1, n1, indexing="ij")
    x, y = NX.ravel(), NY.ravel()
    r = np.hypot(x, y)
    keep = (r > 0) & (r < R)
    x, y, r = x[keep], y[keep], r[keep]
    part = float(np.sum(Phi_xy(x, y) * r ** (-s) * (1 - smooth_step(r / R))))
    gint = integrate.quad(lambda t: t ** (1 - s) * float(smooth_step(np.array([t]))[0]), 0.5, 1.0, epsabs=1e-15, limit=200)[0] \
        + 1.0 / (s - 2)
    return part + A_Phi * R ** (2 - s) * gint


Phi2_xy = lambda x, y: (x * x / (x * x + y * y)) / (1 + aa * x * x * y * y / (x * x + y * y) ** 2) ** 2


def mesh_moments_quadrant(N, alpha, chunk=200):
    """M_m(N) = N^-2 sum over the odd Gamma-centred mesh (Gamma excluded), folded onto n_x, n_y >= 0."""
    h = (N - 1) // 2
    n = np.arange(0, h + 1, dtype=np.float64)
    wt1 = np.where(n == 0, 1.0, 2.0)
    out_m = {0: 0.0, 1: 0.0, 2: 0.0}
    for s0 in range(0, h + 1, chunk):
        nx = n[s0:s0 + chunk]
        NX, NY = np.meshgrid(nx, n, indexing="ij")
        WX, WY = np.meshgrid(wt1[s0:s0 + chunk], wt1, indexing="ij")
        kx = TP * NX.ravel() / N
        ky = TP * NY.ravel() / N
        w = (WX * WY).ravel()
        keep = (kx > 0) | (ky > 0)
        kx, ky, w = kx[keep], ky[keep], w[keep]
        k2 = kx * kx + ky * ky
        k = np.sqrt(k2)
        cx2 = kx * kx / k2
        cy2 = ky * ky / k2
        J2 = (cx2 + c1 * k) * np.exp(-k2 / 8)
        r = (1 + aa * cx2 * cy2) * k ** alpha * (1 + c2 * k2)
        for m in (0, 1, 2):
            out_m[m] += float(np.sum(w * J2 / r ** m))
    return {m: out_m[m] / N ** 2 for m in out_m}


# S_2 for alpha = 1.5 and 2 (paper normalisation: S_2 = sum Phi_2 |Bn|^{-2 alpha} = (2pi)^{-2alpha} sum Phi_2 |n|^{-2alpha})
S2n = {}
for alpha in (1.5, 2.0):
    s = 2 * alpha
    vals = [S_smooth(Phi2_xy, A2_cf, s, R) for R in (40.0, 80.0, 160.0)]
    S2n[alpha] = vals[-1]
    print(f"A: alpha = {alpha}: sum_n Phi_2 |n|^-{s:g} by smoothed cutoff R = 40, 80, 160: "
          + ", ".join(f"{v:.13f}" for v in vals))
    check(f"A: smoothed-cutoff lattice sum for alpha = {alpha} is converged in R (change < 1e-10 between R = 80 and 160)",
          abs(vals[-1] - vals[-2]) < 1e-10)
print(f"    for comparison, the sharp-cutoff value of the earlier script (alpha = 1.5): 3.48512991")
check("A: S_2(alpha = 1.5) agrees with the sharp-cutoff value 3.48512991 used for Table S2 (to 1e-7)",
      abs(S2n[1.5] - 3.48512991) < 1e-7, f"{S2n[1.5]:.10f}")
S2p = {a_: (TP ** (-2 * a_)) * S2n[a_] for a_ in S2n}   # paper normalisation, V_cell = 1

# continuum moments
M = {}
for alpha in (0.6, 1.0, 1.5):
    for m in (1, 2):
        if m * alpha < 2:
            M[(m, alpha)] = cont_moment_A(m, alpha)
print("A: continuum moments (QUADPACK, algebraic weight): " + ", ".join(f"M_{m}(alpha={a_}) = {v:.10f}" for (m, a_), v in M.items()))
tau_inf = M[(2, 0.6)] / M[(1, 0.6)]
slope1 = (A2_cf / TP ** 2) / M[(1, 1.0)]
print(f"A: M_2/M_1 (alpha = 0.6) = {tau_inf:.10f} (Table S2: 0.6291354); (2pi)^-2 A_2/M_1 (alpha = 1) = {slope1:.10f} (Table S2: 0.3885388)")
check("A: Table S2 predictions 0.6291354 and 0.3885388 reproduced by an independent quadrature",
      abs(tau_inf - 0.6291354) < 6e-8 and abs(slope1 - 0.3885388) < 6e-8)

# mesh sums
mesh = {}
for alpha in (0.6, 1.0, 1.5, 2.0):
    mesh[alpha] = [mesh_moments_quadrant(N, alpha) for N in Ns]
# alpha = 0.6
errs = [abs(mm[2] / mm[1] - tau_inf) / tau_inf for mm in mesh[0.6]]
expo = math.log(errs[3] / errs[5]) / math.log(4)
print("A: alpha = 0.6: |tau(N) - M2/M1|/(M2/M1) = " + ", ".join(f"{e:.3e}" for e in errs) + f"; exponent (801 -> 3201) = {expo:.3f}")
check("A: alpha = 0.6 row of Table S2: relative error 2.6e-3 at N1 = 3201, decreasing like N1^-0.79",
      abs(errs[-1] - 2.6e-3) < 0.05e-3 and abs(expo - 0.79) < 0.006)
# alpha = 1
taus = [mm[2] / mm[1] for mm in mesh[1.0]]
ls = [(taus[i + 1] - taus[i]) / math.log(Ns[i + 1] / Ns[i]) / slope1 for i in range(len(Ns) - 1)]
print("A: alpha = 1: local slope / prediction = " + ", ".join(f"{x:.6f}" for x in ls))
check("A: alpha = 1 row: local slope 0.99657 x prediction at N1 = 3201, increasing", abs(ls[-1] - 0.99657) < 6e-6
      and all(ls[i + 1] > ls[i] for i in range(len(ls) - 1)))
# alpha = 1.5
M1c = M[(1, 1.5)]
m2dev = [mm[2] / N ** 1 / S2p[1.5] - 1 for mm, N in zip(mesh[1.5], Ns)]
m1err = [abs(mm[1] - M1c) / M1c for mm in mesh[1.5]]
r_m1 = math.log(m1err[-2] / m1err[-1]) / math.log(2)
rat = [mm[2] / mm[1] / N / (S2p[1.5] / M1c) for mm, N in zip(mesh[1.5], Ns)]
q = 2.0 ** (-0.5)
Rinf = rat[-1] - (rat[-2] - rat[-1]) * q / (1 - q)
print("A: alpha = 1.5: M2(N)/(N S2) - 1 = " + ", ".join(f"{x:+.3e}" for x in m2dev)
      + f"; M1 error exponent = {r_m1:.4f}; tau/[N S2/M1] = " + ", ".join(f"{x:.5f}" for x in rat)
      + f"; Richardson (exponent 1/2) = {Rinf:.5f}")
# a two-term fit 1 + a N^-1/2 + b N^-1 log N on the last three points, as a cross-check of the extrapolation
Xf = np.array([[1, N ** -0.5, math.log(N) / N] for N in Ns[-3:]])
coef = np.linalg.solve(Xf, np.array(rat[-3:]))
print(f"    three-term fit 1 + a N^-1/2 + b N^-1 log N on N = 801, 1601, 3201: limit = {coef[0]:.5f}")
check("A: alpha = 1.5 row: M2 deviation 2.5e-3, M1 exponent 0.500, tau ratio 1.026 at 3201, Richardson value 0.99785",
      abs(m2dev[-1] - 2.5e-3) < 0.06e-3 and abs(r_m1 - 0.5) < 0.0006 and abs(rat[-1] - 1.026) < 6e-4
      and abs(Rinf - 0.99785) < 6e-6)
# alpha = 2
sl1 = [(mesh[2.0][i + 1][1] - mesh[2.0][i][1]) / math.log(Ns[i + 1] / Ns[i]) / (A1_cf / TP ** 2) for i in range(len(Ns) - 1)]
m2r = [mm[2] / N ** 2 / S2p[2.0] for mm, N in zip(mesh[2.0], Ns)]
tl = [mm[2] / mm[1] * math.log(N) / N ** 2 / (TP ** 2 * S2p[2.0] / A1_cf) for mm, N in zip(mesh[2.0], Ns)]
print("A: alpha = 2: M1 slope/pred = " + ", ".join(f"{x:.6f}" for x in sl1) + "; (M2/N^2)/S2 = "
      + ", ".join(f"{x:.6f}" for x in m2r) + "; tau log N/N^2 / limit = " + ", ".join(f"{x:.4f}" for x in tl))
check("A: alpha = 2 row: ratios 1.00087 and 1.00158, tau log N1/N1^d = 1.057 x limit at 3201",
      abs(sl1[-1] - 1.00087) < 6e-6 and abs(m2r[-1] - 1.00158) < 6e-6 and abs(tl[-1] - 1.057) < 6e-4)
out["A"] = {"A1": A1_cf, "A2": A2_cf, "S2_n": S2n, "M": {f"{k}": v for k, v in M.items()},
            "mesh": {str(a_): v for a_, v in mesh.items()}}

# =============================================================================== Part B: hexagonal lattice, jumps
print("B: hexagonal lattice, angular functions with jumps on lattice directions")
B = TP * np.array([[1.0, 0.0], [-1 / math.sqrt(3), 2 / math.sqrt(3)]])   # columns b1, b2
detB = abs(np.linalg.det(B))
Vc = TP ** 2 / detB
Binv = np.linalg.inv(B)


def sector_flags(n1, n2):
    """Exact integer tests on n (direction of Bn): F = theta mod pi in [0, pi/3); U = k_y > 0."""
    t = 2 * n2 - n1          # k_y = (2 pi/sqrt 3)(2 n2 - n1), k_x = 2 pi n1
    in1 = (t >= 0) & (n2 < 2 * n1)                # 0 <= theta < pi/3
    in2 = (t <= 0) & (n2 > 2 * n1)                # pi <= theta < 4 pi/3
    U = t > 0
    return (in1 | in2), U


def model_B(kx, ky, F, U, alpha, alt=False):
    k2 = kx * kx + ky * ky
    k = np.sqrt(k2)
    cx2 = kx * kx / k2
    PhiJ = cx2 + 0.5 * F
    r0 = 1 + 0.4 * U
    if alt:   # the other one-sided values on the jump directions theta = 0 and pi (k_y = 0): r0 = 1.4 there
        r0 = np.where(ky == 0, 1.4, r0)
    J2 = (PhiJ + c1 * k) * np.exp(-k2 / 8)
    r = r0 * k ** alpha * (1 + c2 * k2)
    return J2, r, PhiJ, r0


def mesh_B(N, alpha, chunk=256):
    h = (N - 1) // 2
    n = np.arange(-h, h + 1, dtype=np.int64)
    acc = {0: 0.0, 1: 0.0, 2: 0.0}
    for s0 in range(0, len(n), chunk):
        N1, N2 = np.meshgrid(n[s0:s0 + chunk], n, indexing="ij")
        n1, n2 = N1.ravel(), N2.ravel()
        keep = (n1 != 0) | (n2 != 0)
        n1, n2 = n1[keep], n2[keep]
        kx = (B[0, 0] * n1 + B[0, 1] * n2) / N
        ky = (B[1, 0] * n1 + B[1, 1] * n2) / N
        t = 2 * n2 - n1
        ky = np.where(t == 0, 0.0, ky)
        F, U = sector_flags(n1, n2)
        J2, r, _, _ = model_B(kx, ky, F, U, alpha)
        for m in (0, 1, 2):
            acc[m] += float(np.sum(J2 / r ** m))
    return {m: acc[m] / (N ** 2 * Vc) for m in acc}


def S_B(m, alpha, R, alt=False):
    """S_m = sum_{n != 0} Phi_m(Bn/|Bn|) |Bn|^{-m alpha}, |Bn| < R, plus the continuum tail A_m R^{2-s}/((s-2)|det B|)."""
    s = m * alpha
    L = int(R * np.abs(Binv).sum(1).max()) + 2
    n = np.arange(-L, L + 1, dtype=np.int64)
    tot = 0.0
    for s0 in range(0, len(n), 256):
        N1, N2 = np.meshgrid(n[s0:s0 + 256], n, indexing="ij")
        n1, n2 = N1.ravel(), N2.ravel()
        kx = (B[0, 0] * n1 + B[0, 1] * n2).astype(float)
        ky = (B[1, 0] * n1 + B[1, 1] * n2).astype(float)
        t = 2 * n2 - n1
        ky = np.where(t == 0, 0.0, ky)
        kk = np.hypot(kx, ky)
        keep = (kk > 0) & (kk < R)
        F, U = sector_flags(n1[keep], n2[keep])
        _, _, PhiJ, r0 = model_B(kx[keep], ky[keep], F, U, alpha, alt)
        tot += float(np.sum(PhiJ / r0 ** m * kk[keep] ** (-s)))
    return tot + A_B(m) * R ** (2 - s) / ((s - 2) * detB)


def PhiB_theta(th, m):
    thm = math.fmod(th, math.pi)
    F = 1.0 if (0 <= thm < math.pi / 3) else 0.0
    U = 1.0 if math.sin(th) > 0 and 0 < th < math.pi else 0.0
    return (math.cos(th) ** 2 + 0.5 * F) / (1 + 0.4 * U) ** m


def A_B(m):
    pts = [math.pi / 3, math.pi, 4 * math.pi / 3]
    return integrate.quad(lambda t: PhiB_theta(t, m), 0, TP, points=pts, limit=400, epsabs=1e-13)[0]


def A_B_closed(m):
    # int cos^2 over the upper (U=1) and lower half: pi/2 each; jump 0.5 on [0, pi/3) (U=1) and [pi, 4pi/3) (U=0)
    return (math.pi / 2) / 1.4 ** m + (math.pi / 2) + 0.5 * (math.pi / 3) / 1.4 ** m + 0.5 * (math.pi / 3)


def K_of(th):
    w = np.array([math.cos(th), math.sin(th)])
    u = Binv @ w
    return min(0.5 / abs(u[0]) if u[0] != 0 else np.inf, 0.5 / abs(u[1]) if u[1] != 0 else np.inf)


def cont_B(m, alpha):
    corners = []
    for sx in (-0.5, 0.5):
        for sy in (-0.5, 0.5):
            v = B @ np.array([sx, sy])
            corners.append(math.atan2(v[1], v[0]) % TP)
    brk = sorted(set([0.0, math.pi / 3, math.pi, 4 * math.pi / 3, TP] + corners))

    def radial(th):
        ct = math.cos(th)
        thm = math.fmod(th, math.pi)
        F = 1.0 if (0 <= thm < math.pi / 3) else 0.0
        U = 1.0 if 0 < th < math.pi else 0.0
        r0 = 1 + 0.4 * U
        f = lambda k: (ct * ct + 0.5 * F + c1 * k) * math.exp(-k * k / 8) / (r0 * (1 + c2 * k * k)) ** m
        return integrate.quad(f, 0, K_of(th), weight="alg", wvar=(1 - m * alpha, 0.0), epsabs=1e-14, epsrel=1e-12,
                              limit=200)[0]

    tot = 0.0
    for a_, b_ in zip(brk[:-1], brk[1:]):
        tot += integrate.quad(radial, a_, b_, epsabs=1e-13, epsrel=1e-11, limit=200)[0]
    return tot / TP ** 2


okA = all(abs(A_B(m) - A_B_closed(m)) < 1e-10 for m in (1, 2))
print(f"    V_cell = {Vc:.10f}; A_1 = {A_B_closed(1):.12f}, A_2 = {A_B_closed(2):.12f} (closed form = quadrature: {okA})")
meshB = {alpha: [mesh_B(N, alpha) for N in Ns] for alpha in (0.6, 1.0, 1.5)}
# (i) m alpha < d: convergence to the continuum
Mb1, Mb2 = cont_B(1, 0.6), cont_B(2, 0.6)
e1 = [abs(mm[1] - Mb1) / Mb1 for mm in meshB[0.6]]
e2 = [abs(mm[2] - Mb2) / Mb2 for mm in meshB[0.6]]
print("    alpha = 0.6: |M1(N) - M1|/M1 = " + ", ".join(f"{x:.2e}" for x in e1) + "; |M2(N) - M2|/M2 = "
      + ", ".join(f"{x:.2e}" for x in e2))
check("B: m alpha < d with jumps on lattice directions: M_1(N1) -> M_1 and M_2(N1) -> M_2 (alpha = 0.6, errors decreasing, "
      "< 1e-2 at 3201)", okA and e1[-1] < 1e-2 and e2[-1] < 1e-2 and e2[-1] < e2[2] and e1[-1] < e1[2])
# (ii) m alpha = d: slope (2pi)^-2 A_2 for alpha = 1, m = 2
slB = [(meshB[1.0][i + 1][2] - meshB[1.0][i][2]) / math.log(Ns[i + 1] / Ns[i]) / (A_B_closed(2) / TP ** 2)
       for i in range(len(Ns) - 1)]
print("    alpha = 1: slope of M2(N) in log N / ((2pi)^-2 A_2) = " + ", ".join(f"{x:.5f}" for x in slB))
check("B: m alpha = d with jumps: slope of M_2(N1) in log N1 -> (2pi)^-2 A_2 (within 5e-3 at 3201)", abs(slB[-1] - 1) < 5e-3)
# (iii) m alpha > d: M_2(N)/N -> S_2/V_cell, with the values on the jump directions as in the mesh data
S2B = S_B(2, 1.5, 1500.0)
S2B_half = S_B(2, 1.5, 750.0)
S2B_alt = S_B(2, 1.5, 1500.0, alt=True)
rB = [mm[2] / N / (S2B / Vc) for mm, N in zip(meshB[1.5], Ns)]
print(f"    alpha = 1.5: S_2 = {S2B:.10f} (R = 1500; R = 750: {S2B_half:.10f}); with the other one-sided values of r_0 on "
      f"theta = 0, pi: S_2 = {S2B_alt:.10f} (relative change {S2B_alt / S2B - 1:+.4f})")
print("    alpha = 1.5: M2(N) V_cell/(N S_2) = " + ", ".join(f"{x:.5f}" for x in rB))
dev_rB = [abs(x - 1) for x in rB]
check("B: m alpha > d with jumps: M_2(N1)/N1^{2alpha-d} -> S_2/V_cell, where S_2 uses the values of Phi_2 on the lattice "
      "directions inside the jump set (deviation < 3e-3 at 3201 and decreasing over N1 = 801, 1601, 3201, after a "
      "pre-asymptotic crossing of 1)",
      dev_rB[-1] < 3e-3 and dev_rB[-1] < dev_rB[-2] < dev_rB[-3] and abs(S2B - S2B_half) < 1e-3 * S2B)
check("B: the constant S_2 changes by a finite amount (> 5%) if Phi_2 is changed only on the null set of jump directions (theta = 0, pi)",
      abs(S2B_alt / S2B - 1) > 0.05)
out["B"] = {"Vc": Vc, "A1": A_B_closed(1), "A2": A_B_closed(2), "M1_06": Mb1, "M2_06": Mb2, "S2": S2B,
            "S2_alt": S2B_alt, "mesh": {str(a_): v for a_, v in meshB.items()}}
json.dump(out, open("audit_s4_mesh_laws.json", "w"), indent=1, default=float)
print(f"    elapsed {time.time() - t0:.1f} s")
print()
print(f"{sum(res)} of {len(res)} checks passed")
