"""Independent checks of Lemma S1 (Supplemental Sec. S4, lattice sums), away from the paper's test case.

The paper's numerical test uses the square lattice B = 2 pi I and a continuous angular function. Here:
  * skew and non-orthogonal lattices (hexagonal, a random skew 2D lattice, the bcc reciprocal lattice of fcc in 3D);
  * angular functions Phi that are discontinuous on a closed null set of directions (jumps), and one that vanishes
    outside a small arc (the positivity argument of part (c));
  * adversarial error terms eps_N(n) = +theta(|Bn|/N) with slowly decaying theta.
Checked:
  1. counting bound #{0 < |Bn| < R} <= (3R/rho0)^d for R >= rho0, and the shell constant C1 = (6/rho0)^d of (S5);
  2. the shell Riemann sums a_l / 2^{l(d-s)} -> c* = |det B|^-1 A int_1^2 rho^{d-1-s} d rho, (S6);
  3. part (a): the uniform small-ball bound with C = C1/(1 - 2^{s-d}), and the limit of N^{s-d} S_N;
  4. part (b): S_N = |det B|^-1 A log N + o(log N), with an error term whose sum grows like log log N;
  5. part (c): S_N -> S_Phi > 0, including Phi supported on a small arc; error terms -> 0.
Usage (from this directory): python -B audit_s4_lattice_lemma.py > audit_s4_lattice_lemma.out.txt
Single thread; < 400 MB.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import math
import time
import numpy as np
from scipy import integrate

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""), flush=True)


t0 = time.time()
TWO_PI = 2 * math.pi
lattices = {
    "square 2pi I": np.array([[TWO_PI, 0.0], [0.0, TWO_PI]]),
    "hexagonal": TWO_PI * np.array([[1.0, 0.0], [-1 / math.sqrt(3), 2 / math.sqrt(3)]]),
    "skew 2D": np.array([[1.0, 0.37], [0.0, 0.83]]),
    "bcc (fcc reciprocal)": TWO_PI * np.array([[-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]]),
}


def lattice_points(B, Rmax, chunk=2_000_000):
    """Yield arrays (m x d) of points Bn with 0 < |Bn| < Rmax (n integer), in chunks."""
    d = B.shape[0]
    Binv = np.linalg.inv(B)
    # |n_i| <= Rmax * ||row i of B^-1||
    lim = [int(math.ceil(Rmax * np.linalg.norm(Binv[i]))) + 1 for i in range(d)]
    rng = [np.arange(-L, L + 1, dtype=np.float64) for L in lim]
    if d == 2:
        for a0 in range(0, len(rng[0]), max(1, chunk // len(rng[1]))):
            A0 = rng[0][a0:a0 + max(1, chunk // len(rng[1]))]
            NX, NY = np.meshgrid(A0, rng[1], indexing="ij")
            n = np.stack([NX.ravel(), NY.ravel()], axis=1)
            k = n @ B.T
            r = np.sqrt((k * k).sum(1))
            keep = (r > 0) & (r < Rmax)
            yield k[keep], r[keep]
    else:
        for a in rng[0]:
            NY, NZ = np.meshgrid(rng[1], rng[2], indexing="ij")
            n = np.stack([np.full(NY.size, a), NY.ravel(), NZ.ravel()], axis=1)
            k = n @ B.T
            r = np.sqrt((k * k).sum(1))
            keep = (r > 0) & (r < Rmax)
            yield k[keep], r[keep]


def rho0_of(B):
    best = np.inf
    for k, r in lattice_points(B, 4 * max(np.linalg.norm(B, axis=0))):
        if r.size:
            best = min(best, r.min())
    return best


# ---------------------------------------------------------------- 1. counting and shell constants
print("=== 1. counting bound and shell constant C1 = (6/rho0)^d")
ok1, worst_count, worst_shell = True, 0.0, 0.0
for name, B in lattices.items():
    d = B.shape[0]
    rho0 = rho0_of(B)
    Rtop = (40 if d == 2 else 14) * rho0
    radii = []
    for k, r in lattice_points(B, Rtop):
        radii.append(r)
    radii = np.sort(np.concatenate(radii))
    Rs = rho0 * np.geomspace(1.0, Rtop / rho0 * 0.999, 60)
    counts = np.searchsorted(radii, Rs, side="left")  # #{0 < |Bn| < R}
    ratio_c = (counts / (3 * Rs / rho0) ** d).max()
    worst_count = max(worst_count, ratio_c)
    C1 = (6 / rho0) ** d
    shell_ratios = []
    for s in (0.0, 0.5 * d, float(d), 1.5 * d):
        lmin = int(math.floor(math.log2(rho0))) - 1
        lmax = int(math.floor(math.log2(Rtop))) - 1
        for l in range(lmin, lmax + 1):
            sel = radii[(radii >= 2.0 ** l) & (radii < 2.0 ** (l + 1))]
            if sel.size == 0:
                if 2.0 ** (l + 1) > rho0 * 1.0000001 and 2.0 ** l <= rho0:
                    pass
                continue
            shell_ratios.append(float((sel ** (-s)).sum() / (C1 * 2.0 ** (l * (d - s)))))
            # emptiness rule: nonempty shells have 2^{l+1} > rho0
            ok1 = ok1 and 2.0 ** (l + 1) > rho0
    mr = max(shell_ratios)
    worst_shell = max(worst_shell, mr)
    ok1 = ok1 and ratio_c <= 1.0 and mr <= 1.0
    print(f"    {name:22s}: rho0 = {rho0:.6f}; max count/(3R/rho0)^d = {ratio_c:.4f}; "
          f"max shell sum/(C1 2^(l(d-s))) over s in (0, d/2, d, 3d/2) = {mr:.4f}")
check("counting bound (3R/rho0)^d and shell bound C1 2^{l(d-s)}, C1 = (6/rho0)^d, hold on square, hexagonal, skew 2D and bcc 3D "
      "lattices (all ratios <= 1)", ok1, f"worst count ratio {worst_count:.3f}, worst shell ratio {worst_shell:.3f}")

# ---------------------------------------------------------------- angular functions
A_CACHE = {}


def Phi_jump(k):
    """d = 2: discontinuous on the closed null set of directions {theta = 0, pi/3, pi, 4pi/3} (four points)."""
    th = np.mod(np.arctan2(k[:, 1], k[:, 0]), math.pi)  # period pi
    c2 = (k[:, 0] ** 2) / (k * k).sum(1)
    return c2 + 0.5 * (th < math.pi / 3) + 0.25


def Phi_jump_theta(th):
    th = math.fmod(th, math.pi)
    if th < 0:
        th += math.pi
    return math.cos(th) ** 2 + 0.5 * (th < math.pi / 3) + 0.25


def Phi_arc(k):
    """Phi = 1 on the arc 0.10 < theta < 0.12, else 0 (continuous outside a closed null set)."""
    th = np.arctan2(k[:, 1], k[:, 0])
    return ((th > 0.10) & (th < 0.12)).astype(float)


A_jump = integrate.quad(Phi_jump_theta, 0, 2 * math.pi, points=[math.pi / 3, math.pi, 4 * math.pi / 3], limit=200)[0]
A_jump_exact = math.pi + 0.5 * (2 * math.pi / 3) + 0.25 * 2 * math.pi  # int cos^2 + jump part + const
A_arc = 0.02
print(f"    A(Phi_jump) = {A_jump:.12f} (closed form {A_jump_exact:.12f}); A(Phi_arc) = {A_arc}")


def phi_sum(B, Phi, s, lo, hi, eps_fun=None, N=None):
    """sum over lo <= |Bn| < hi of Phi(dir) |Bn|^-s (+ eps_N(n) |Bn|^-s if eps_fun given)."""
    tot, tot_e = 0.0, 0.0
    for k, r in lattice_points(B, hi):
        sel = r >= lo
        k, r = k[sel], r[sel]
        w = r ** (-s)
        tot += float(np.sum(Phi(k) * w))
        if eps_fun is not None:
            tot_e += float(np.sum(eps_fun(r / N) * w))
    return tot, tot_e


# ---------------------------------------------------------------- 2. shell Riemann sums
print("=== 2. shell Riemann sums a_l / 2^{l(d-s)} -> c* (hexagonal lattice, Phi with jumps)")
B = lattices["hexagonal"]
detB = abs(np.linalg.det(B))
ok2 = True
for s in (1.0, 2.0, 3.0):
    cstar = A_jump / detB * (math.log(2) if s == 2 else (2 ** (2 - s) - 1) / (2 - s))
    rows = []
    for l in range(4, 12):
        a_l, _ = phi_sum(B, Phi_jump, s, 2.0 ** l, 2.0 ** (l + 1))
        rows.append((l, a_l / 2.0 ** (l * (2 - s)) / cstar - 1))
    errs = [abs(e) for _, e in rows]
    print(f"    s = {s}: c* = {cstar:.10f}; relative deviation by l = 4..11: " + ", ".join(f"{e:+.2e}" for _, e in rows))
    ok2 = ok2 and errs[-1] < 2e-3 and errs[-1] < errs[0] / 8
check("(S6): a_l/2^{l(d-s)} -> c* for s = 1, 2, 3 on the hexagonal lattice with an angular function discontinuous on "
      "four directions (deviation < 2e-3 at l = 11 and decreasing)", ok2)

# ---------------------------------------------------------------- 3. part (a)
print("=== 3. part (a), s = 1 < d = 2 (hexagonal and skew lattices)")
ok3 = True
for name in ("hexagonal", "skew 2D"):
    B = lattices[name]
    detB = abs(np.linalg.det(B))
    rho0 = rho0_of(B)
    s, d, k0 = 1.0, 2, rho0 / 4
    C = (6 / rho0) ** d / (1 - 2 ** (s - d))
    # uniform bound N^{s-d} sum_{0<|Bn|<delta N} |Bn|^{-s} <= C delta^{d-s}
    worst = 0.0
    for N in (50, 200, 800):
        for delta in rho0 * np.array([0.02, 0.1, 0.5, 2.0]):
            tot, _ = phi_sum(B, lambda k: np.ones(len(k)), s, 0.0, delta * N)
            worst = max(worst, N ** (s - d) * tot / (C * delta ** (d - s)))
    target = A_jump * k0 ** (d - s) / ((d - s) * detB)
    th_eps = lambda t: 0.2 * np.sqrt(np.minimum(t / k0, 1.0))
    pred_err = 0.2 * 2 * math.pi * integrate.quad(lambda q: math.sqrt(q / k0) * q ** (d - 1 - s), 0, k0)[0] / detB
    bound_err = 0.2 * 2 * math.pi * k0 ** (d - s) / ((d - s) * detB)
    rows = []
    for N in (250, 1000, 4000):
        tot, tot_e = phi_sum(B, Phi_jump, s, 0.0, k0 * N, th_eps, N)
        rows.append((N, N ** (s - d) * tot / target - 1, N ** (s - d) * tot_e))
    print(f"    {name}: max ratio in uniform bound = {worst:.4f} (C = {C:.4f}); N^(s-d) sum phi / target - 1 = "
          + ", ".join(f"{e:+.2e}" for _, e, _ in rows)
          + f"; error-term sum -> {rows[-1][2]:.6f} (integral prediction {pred_err:.6f}, bound (S4) {bound_err:.6f})")
    ok3 = ok3 and worst <= 1 and abs(rows[-1][1]) < 2e-3 and abs(rows[-1][1]) < abs(rows[0][1]) \
        and abs(rows[-1][2] - pred_err) < 2e-3 * pred_err + 1e-3 and rows[-1][2] <= bound_err
check("part (a): uniform small-ball bound with C = C1/(1-2^{s-d}); N^{s-d} sum phi -> A k0^{d-s}/((d-s)|det B|); "
      "error terms within the bound of (S4)", ok3)

# ---------------------------------------------------------------- 4. part (b)
print("=== 4. part (b), s = d = 2 (hexagonal lattice; Phi with jumps; error theta(t) = 0.5/(1+|log t|))")
B = lattices["hexagonal"]
detB = abs(np.linalg.det(B))
k0 = 1.0
slope_pred = A_jump / detB
th_b = lambda t: 0.5 / (1 + np.abs(np.log(t)))
Ns = [2 ** j for j in range(6, 15)]
vals = []
for N in Ns:
    tot, tot_e = phi_sum(B, Phi_jump, 2.0, 0.0, k0 * N, th_b, N)
    vals.append((N, tot, tot_e))
sl = [(vals[i + 1][1] - vals[i][1]) / math.log(vals[i + 1][0] / vals[i][0]) / slope_pred for i in range(len(vals) - 1)]
offs = [v[1] - slope_pred * math.log(v[0]) for v in vals]
print("    local slope / (|det B|^-1 A) by doubling: " + ", ".join(f"{x:.5f}" for x in sl))
print("    S_N - |det B|^-1 A log N: " + ", ".join(f"{x:.5f}" for x in offs))
print("    error-term sum / log N: " + ", ".join(f"{v[2] / math.log(v[0]):.5f}" for v in vals)
      + "   (sum itself: " + ", ".join(f"{v[2]:.4f}" for v in vals) + ")")
er = [v[2] / math.log(v[0]) for v in vals]
check("part (b): slope of S_N in log N -> |det B|^-1 A (within 5e-3 at the last doubling); "
      "error-term sum grows but its ratio to log N decreases",
      abs(sl[-1] - 1) < 5e-3 and all(er[i + 1] < er[i] for i in range(2, len(er) - 1)))

# ---------------------------------------------------------------- 5. part (c)
print("=== 5. part (c), s = 4 > d = 2 (square lattice; Phi supported on the arc 0.10 < theta < 0.12)")
B = lattices["square 2pi I"]
th_c = lambda t: np.minimum(t, 1.0) ** 0.3
vals = []
for N in (100, 400, 1600, 6400):
    tot, tot_e = phi_sum(B, Phi_arc, 4.0, 0.0, 1.0 * N, th_c, N)
    vals.append((N, tot, tot_e))
# first lattice direction in the arc: n = (9, 1), angle atan(1/9) = 0.1107
first = (2 * math.pi * math.hypot(9, 1)) ** (-4)
print("    S_N: " + ", ".join(f"{v[1]:.12e}" for v in vals) + f"; term n = (9,1): {first:.6e}")
print("    error-term sums: " + ", ".join(f"{v[2]:.4e}" for v in vals))
check("part (c): S_N converges to S_Phi > 0 for Phi supported on a small arc (positivity through a lattice direction "
      "in the arc); error terms -> 0",
      vals[-1][1] >= first > 0 and abs(vals[-1][1] - vals[-2][1]) < 1e-3 * vals[-1][1]
      and all(vals[i + 1][2] < vals[i][2] for i in range(len(vals) - 1)))
# 3D: part (c) and part (b) on bcc
print("=== 5b. d = 3, bcc reciprocal lattice: part (b) slope and part (c) limit with a jump function")
B = lattices["bcc (fcc reciprocal)"]
detB = abs(np.linalg.det(B))


def Phi3(k):
    c = k[:, 2] / np.sqrt((k * k).sum(1))
    return c * c + 0.5 * (k[:, 0] > 0)  # jump on the plane k_x = 0 of directions (closed, null)


A3 = 4 * math.pi / 3 + 0.5 * 2 * math.pi
vals3 = []
for N in (40, 80, 160, 320, 640):
    tot, _ = phi_sum(B, Phi3, 3.0, 0.0, 1.0 * N)
    vals3.append((N, tot))
sl3 = [(vals3[i + 1][1] - vals3[i][1]) / math.log(2) / (A3 / detB) for i in range(len(vals3) - 1)]
vals3c = []
for N in (40, 80, 160, 320):
    tot, _ = phi_sum(B, Phi3, 6.0, 0.0, 1.0 * N)
    vals3c.append(tot)
print("    s = d = 3: local slope / (|det B|^-1 A): " + ", ".join(f"{x:.5f}" for x in sl3))
print("    s = 6: S_N = " + ", ".join(f"{x:.10e}" for x in vals3c))
check("d = 3 (bcc, jump on a plane of directions): log slope -> |det B|^-1 A (within 2e-2); s > d sums converge",
      abs(sl3[-1] - 1) < 2e-2 and abs(vals3c[-1] - vals3c[-2]) < 1e-3 * vals3c[-1])

print(f"    elapsed {time.time() - t0:.1f} s")
print()
print(f"{sum(res)} of {len(res)} checks passed")
