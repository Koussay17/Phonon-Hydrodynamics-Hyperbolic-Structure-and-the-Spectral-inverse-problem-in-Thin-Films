"""Independent check of Proposition 1 (grey Callaway operator) on full cubic meshes with Bose weights.

Mesh: simple-cubic reciprocal lattice, k = (2 pi/N) n, n in {-(N-1)/2..(N-1)/2}^3, N odd (invariant
under the full cubic group O_h), Gamma excluded; one branch eps = hbar v |k| (three identical
branches only triple every sum). Units hbar = v = k_B = 1, V_cell = 1, T fixed.
Matrix-free: C y = (y - Pi_N y)/tau_N + (y - Pi_e y)/tau_R with Pi_N = U (U^T U)^{-1} U^T, U = [e,Px,Py,Pz].

Checks: j_x = v^2 P_x; C j = j/tau_R; e, P_i mutually orthogonal; |j_x|^2 = v^2 |e|^2/3
(lambda = rho c v^2 tau_R/3); spectrum of C on H = {1/tau_R (x3), 1/tau_c}; diag deviation
n * max(1/tau_N + 1/tau_R - r_nu) bounded under refinement; lifetimes -> tau_c; support of j_x.

Usage: python -B audit_callaway.py > audit_callaway.out.txt
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np

res = []


def check(name, ok, extra=""):
    res.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("   [" + extra + "]") if extra else ""))


tauN, tauR = 0.3, 2.0
tauc = 1.0 / (1.0 / tauN + 1.0 / tauR)
T = 1.2  # k_B T in units of hbar v (2 pi / a); max eps on the mesh is about pi*sqrt(3)
rows = []
for N in (5, 9, 13, 17, 25, 33):
    r1 = np.arange(-(N - 1) // 2, (N - 1) // 2 + 1)
    n1, n2, n3 = np.meshgrid(r1, r1, r1, indexing="ij")
    K = (2 * np.pi / N) * np.stack([n1.ravel(), n2.ravel(), n3.ravel()], axis=1).astype(float)
    kk = np.linalg.norm(K, axis=1)
    keep = kk > 0
    K, kk = K[keep], kk[keep]
    n = len(kk)
    Vol = float(N ** 3)  # V = N_k V_cell
    eps = kk  # hbar v |k|
    f0 = 1.0 / np.expm1(eps / T)
    wB = np.sqrt(f0 * (1 + f0))
    vel = K / kk[:, None]  # v * k_hat with v = 1
    e = wB * eps / np.sqrt(Vol)
    P = (wB[:, None] * K) / np.sqrt(Vol)  # hbar k
    jx = wB * eps * vel[:, 0] / np.sqrt(Vol)
    U = np.column_stack([e, P])
    G = U.T @ U
    Ginv = np.linalg.inv(G)

    def Cop(y):
        piN = U @ (Ginv @ (U.T @ y))
        pie = e * (e @ y) / (e @ e)
        return (y - piN) / tauN + (y - pie) / tauR

    off = np.max(np.abs(G - np.diag(np.diag(G)))) / np.max(np.abs(G))
    jx_err = np.max(np.abs(jx - P[:, 0])) / np.max(np.abs(jx))
    eig_err = np.max(np.abs(Cop(jx) - jx / tauR)) / np.max(np.abs(jx / tauR))
    ce = np.max(np.abs(Cop(e))) / np.max(np.abs(e))
    third = (jx @ jx) / (e @ e)  # should be v^2/3
    # diagonal of C: (1 - (Pi_N)_nn)/tauN + (1 - (Pi_e)_nn)/tauR
    dPiN = np.einsum("ij,jk,ik->i", U, Ginv, U)
    dPie = e ** 2 / (e @ e)
    r = (1 - dPiN) / tauN + (1 - dPie) / tauR
    dev = (1 / tauN + 1 / tauR) - r
    # spectrum on H: test a random vector orthogonal to span{e,P}: eigenvalue 1/tau_c
    rng = np.random.default_rng(N)
    z = rng.standard_normal(n)
    z -= U @ (Ginv @ (U.T @ z))
    eig_c = np.max(np.abs(Cop(z) - z / tauc)) / np.max(np.abs(z / tauc))
    # support of the physical current j_x: modes with k_x = 0 have v_x = 0
    zeros = int(np.sum(np.abs(jx) < 1e-14 * np.max(np.abs(jx))))
    rows.append((N, n, off, jx_err, eig_err, ce, third, n * dev.max(), dev.min() >= 0, eig_c, np.max(np.abs(1 / r - tauc)) / tauc, zeros))

print(" N     n    offdiag(U^TU)  |jx-v^2Px|  |Cj-j/tR|  |Ce|     |jx|^2/|e|^2   n*max dev  dev>=0  |Cz-z/tc|  max|1/r-tc|/tc  #zeros(jx)")
for row in rows:
    print("%3d %6d   %9.2e    %9.2e   %9.2e  %8.1e   %.12f   %8.4f   %s   %9.2e   %9.3e   %d" % row)
check("e, P_x, P_y, P_z mutually orthogonal on O_h-symmetric meshes", all(r_[2] < 1e-12 for r_ in rows))
check("j_x = v^2 P_x (linear isotropic dispersion)", all(r_[3] < 1e-13 for r_ in rows))
check("C j = j / tau_R and C e = 0", all(r_[4] < 1e-12 and r_[5] < 1e-12 for r_ in rows))
check("|j_x|^2 = v^2 |e|^2 / 3 (gives lambda = rho c v^2 tau_R / 3)", all(abs(r_[6] - 1 / 3) < 1e-12 for r_ in rows))
check("on span{e,P}^perp, C = 1/tau_c (so C|_H has rates 1/tau_R (x3) and 1/tau_c; ker C = span{e})", all(r_[9] < 1e-12 for r_ in rows))
check("diag C = 1/tau_N + 1/tau_R - dev with 0 <= dev and n*max(dev) bounded under refinement",
      all(r_[8] for r_ in rows) and max(r_[7] for r_ in rows) < 2 * min(r_[7] for r_ in rows[2:]))
check("lifetimes 1/r_nu -> tau_c (relative deviation decreasing to 0)", rows[-1][10] < rows[0][10] and rows[-1][10] < 1e-2)
check("physical current j_x vanishes on the plane k_x = 0 (it does NOT have full support)", all(r_[11] > 0 for r_ in rows),
      "zero entries: " + ", ".join(str(r_[11]) for r_ in rows))

print()
print(f"{sum(res)} of {len(res)} checks passed")
