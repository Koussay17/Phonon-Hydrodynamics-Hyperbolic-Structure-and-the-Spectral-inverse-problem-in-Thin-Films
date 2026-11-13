"""Independent audit core for the full-operator kinetic slab (campaign 20261111-fdtr-gauge-memory, branch R-S audit).

Written without reference to rs_core.py solvers. Only the raw upstream data are shared with the audited code:
event list, stored fc3 amplitudes (pp), symmetry maps of fc3-memory (read-only).

1. Operator (own assembly).  Modes (q, j) with f > cutoff; x = hbar w / kT; n = 1/(e^x - 1);
   capacity c_i = k_B x^2 n(n+1) / (N_q V)  [J m^-3 K^-1];  e_i = sqrt(c_i);  v_i = v_z [m/s].
   Event vectors (campaign convention): -lp/s_P at P, ld/s_A at A (2 ld/s_A if A = B), ld/s_B at B, with
   lp = sqrt((w_A + w_B)/w_P), ld = 1/lp, s = sqrt(n(n+1)).  Rate g = W pp, W = conv g0w / (4 sh_P sh_A sh_B)
   (A = B: / (8 sh_P sh_A^2)), sh = sinh(x/2).  Orbit tying: g averaged over the orbit of each event under the 24
   stored grid maps; orbit label = minimum event index over the 24 images (single pass, the maps form a group).
   C = 4 pi 1e12 sum g a a^T  [1/s], dense symmetric; C e = 0.

2. Slab solver (no parity reduction, no Schur complement, no eigendecomposition, no invariant imbedding).
   V f' + (p + C) f = 0 on 0 < z < d; black diffuse walls f_i(0) = T_A e_i (v_i > 0), f_i(d) = T_B e_i (v_i < 0).
   Splitting C = C_R - K', C_R = D - (De)(De)^T/W (energy-conserving RTA on the diagonal D = diag C, W = e^T D e),
   K' = C_R - C (e^T K' = 0, K' e = 0).
   Inner problem (RTA slab with source s): v_i f_i' + (p + D_i) f_i = D_i e_i theta(z) + s_i(z),
   theta = (De)^T f / W; discretised on a z mesh by exact characteristic integration of each mode with the total
   source represented as piecewise linear between nodes (linear-characteristic scheme, O(h^2)); local modes
   (v = 0): f_i = (D_i e_i theta + s_i)/(p + D_i) pointwise. The scalar field theta on the nodes is obtained from a
   dense (N+1) x (N+1) system (I - A) theta = r.
   Outer problem: f = R0[K' f] + F_bc, solved by GMRES on the node values of f (all modes, both directions).
   Two-port: q(0), q(d) (q = sum v e f) and energy temperatures e.f/C_tot at the walls for (T_A, T_B) = (1,0), (0,1).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
from scipy.sparse.linalg import LinearOperator, gmres  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c")
EV = {"m553": FC3 / "m553/events_m553_s0.1.npz", "m775": FC3 / "m775x/events_m775_s0.1.npz",
      "m996": FC3 / "m996/events_m996_s0.1.npz"}
PP = {("m553", "filmw"): FC3 / "m553/filmw_pp", ("m775", "film2"): FC3 / "m775x/film2_pp",
      ("m775", "filmw"): FC3 / "m775x/filmw_pp", ("m996", "film2"): FC3 / "m996/film2_pp"}
KB = 1.380649e-23


# ------------------------------------------------------------------------------------------------ operator
def load_modes(mesh):
    Z = np.load(EV[mesh])
    fr = Z["freqs"]
    Nq, nb = fr.shape
    valid = (fr > float(Z["cutoff"])).ravel()
    new = -np.ones(Nq * nb, dtype=np.int64)
    new[valid] = np.arange(valid.sum())
    w = fr.ravel()[valid]
    kT_THz = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    x = w / kT_THz
    n = 1.0 / np.expm1(x)
    s = np.sqrt(n * (n + 1))
    cap = KB * x * x * n * (n + 1) / (Nq * float(Z["volume"]) * 1e-30)
    vz = Z["gv"].reshape(Nq * nb, 3)[valid][:, 2] * 100.0
    pf, af, bf = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    P, A, B = new[pf], new[af], new[bf]
    if not ((P >= 0).all() and (A >= 0).all() and (B >= 0).all()):
        raise RuntimeError("event refers to a mode below cutoff")
    rep = A == B
    lp = np.sqrt((w[A] + w[B]) / w[P])
    ld = 1.0 / lp
    aP = -lp / s[P]
    aA = np.where(rep, 2.0 * ld, ld) / s[A]
    aB = np.where(rep, 0.0, ld / s[B])
    sh = np.sinh(x / 2)
    # orbits: minimum event index over the images under the 24 grid maps
    NF = Nq * nb
    key = (pf * NF + np.minimum(af, bf)) * NF + np.maximum(af, bf)
    order = np.argsort(key, kind="stable")
    skey = key[order]
    m = key.size
    lab = np.arange(m)
    missing = 0
    for R in Z["rot_maps"]:
        ip = R[pf // nb] * nb + pf % nb
        ia = R[af // nb] * nb + af % nb
        ib = R[bf // nb] * nb + bf % nb
        k2 = (ip * NF + np.minimum(ia, ib)) * NF + np.maximum(ia, ib)
        pos = np.clip(np.searchsorted(skey, k2), 0, m - 1)
        ok = skey[pos] == k2
        missing += int((~ok).sum())
        img = np.where(ok, order[pos], lab)
        lab = np.minimum(lab, img)
        del ip, ia, ib, k2, pos, ok, img
    # one pass suffices if the maps form a group; verify idempotence on one more pass
    _, orb = np.unique(lab, return_inverse=True)
    return {"mesh": mesh, "w": w, "x": x, "s": s, "sh": sh, "cap": cap, "e": np.sqrt(cap), "vz": vz, "P": P, "A": A,
            "B": B, "rep": rep, "aP": aP, "aA": aA, "aB": aB, "orb": orb, "n_orb": int(orb.max() + 1),
            "missing_images": missing, "n": int(w.size), "Nq": Nq, "C_tot": float(cap.sum()), "addr": Z["addr"],
            "meshdims": Z["mesh"], "rot_maps": Z["rot_maps"], "nb": nb, "valid": valid, "new": new,
            "unit_to_WmK": float(Z["unit_to_WmK"])}


def operator(G, group, name):
    Pz = np.load(PP[(G["mesh"], group)] / f"pp_r01_{name}.npz")
    pp, conv, g0w = Pz["pp"], float(Pz["conv"]), Pz["g0w"]
    P, A, B, rep, sh = G["P"], G["A"], G["B"], G["rep"], G["sh"]
    W = np.where(rep, conv * g0w / (8 * sh[P] * sh[A] * sh[A]),
                 conv * g0w / (4 * sh[P] * sh[A] * sh[np.where(rep, A, B)]))
    g = W * pp
    gs = np.bincount(G["orb"], weights=g, minlength=G["n_orb"]) / np.bincount(G["orb"], minlength=G["n_orb"])
    g = gs[G["orb"]]
    n = G["n"]
    idx = (P, A, B)
    val = (G["aP"], G["aA"], G["aB"])
    flat = np.zeros(n * n)
    for i in range(3):
        for j in range(3):
            np.add.at(flat, idx[i] * n + idx[j], g * val[i] * val[j]) if n > 4000 else                 flat.__iadd__(np.bincount(idx[i] * n + idx[j], weights=g * val[i] * val[j], minlength=n * n))
    C = flat.reshape(n, n)
    # symmetrise in place by blocks (C = (C + C^T)/2)
    bs = 512
    for i0 in range(0, n, bs):
        for j0 in range(i0, n, bs):
            a = C[i0:i0 + bs, j0:j0 + bs]
            b = C[j0:j0 + bs, i0:i0 + bs]
            m = 0.5 * (a + b.T)
            C[i0:i0 + bs, j0:j0 + bs] = m
            C[j0:j0 + bs, i0:i0 + bs] = m.T
    C *= 4.0 * np.pi * 1e12
    return C


def kinetic_bulk(C, e, v, L=None, chunk=512):
    """K (W/mK), tau_R, tau_l, delta (s) of a conserving operator (bordered solves; own implementation).
    L given: operator C + Q D_L Q with D_L = |v|/L (energy-conserving regulator, campaign's QDQ completion), built in
    place in the bordered matrix (memory-lean)."""
    n = e.size
    M = np.zeros((n + 1, n + 1))
    M[:n, :n] = C
    if L is not None:
        dL = np.abs(v) / L
        eh = e / np.linalg.norm(e)
        h = dL * eh
        M[np.arange(n), np.arange(n)] += dL
        s_ = float(eh @ h)
        for i0 in range(0, n, chunk):
            sl = slice(i0, min(i0 + chunk, n))
            M[sl, :n] += -np.outer(eh[sl], h) - np.outer(h[sl], eh) + s_ * np.outer(eh[sl], eh)
    M[:n, n] = e
    M[n, :n] = e
    lu = sla.lu_factor(M, check_finite=False, overwrite_a=True)
    b = v * e
    xx = sla.lu_solve(lu, np.append(b, 0.0), check_finite=False)[:n]
    K = float(b @ xx)
    eh = e / np.linalg.norm(e)
    u = v * xx
    u -= eh * (eh @ u)
    y = sla.lu_solve(lu, np.append(u, 0.0), check_finite=False)[:n]
    tR = float(xx @ xx) / K
    tl = float(e @ e) * float(u @ y) / K ** 2
    return {"K": K, "tauR": tR, "taul": tl, "delta": tR - tl}


# ------------------------------------------------------------------------------------------------ z mesh
def make_mesh(d, hmin=0.02e-9, gamma=0.15, hmax=12e-9, level=0):
    """Symmetric graded mesh: h(s) = smooth-min(hmin + gamma s, hmax), s = distance to the nearer wall; then each
    cell bisected 'level' times."""
    half = [0.0]
    while half[-1] < d / 2:
        s_ = half[-1]
        h = 1.0 / np.sqrt((hmin + gamma * s_) ** -2 + hmax ** -2)
        half.append(s_ + h)
    half = np.array(half)
    half *= (d / 2) / half[-1]
    z = np.concatenate([half, d - half[-2::-1]])
    for _ in range(level):
        z = np.sort(np.concatenate([z, 0.5 * (z[1:] + z[:-1])]))
    return z


# ------------------------------------------------------------------------------------------------ slab solver
def _phi12(x):
    """phi1 = (1 - e^-x)/x, phi2 = (x - 1 + e^-x)/x^2, cancellation-free."""
    small = np.abs(x) < 1e-3
    xs = np.where(small, 1.0, x)
    em = np.expm1(-xs)
    p1 = -em / xs
    p2 = (xs + em) / (xs * xs)
    xx = x[small]
    p1[small] = 1 - xx / 2 + xx ** 2 / 6 - xx ** 3 / 24 + xx ** 4 / 120
    p2[small] = 0.5 - xx / 6 + xx ** 2 / 24 - xx ** 3 / 120 + xx ** 4 / 720
    return p1, p2


def _psi(x, jmax=3):
    """psi_j(x) = int_0^1 exp(-x (1 - t)) t^j dt, j = 0..jmax (series for |x| < 1, stable recursion otherwise)."""
    x = np.asarray(x, complex)
    out = np.zeros((jmax + 1,) + x.shape, complex)
    small = np.abs(x) < 1.0
    xs = x[small]
    from math import factorial
    for j in range(jmax + 1):
        acc = np.zeros(xs.shape, complex)
        term = np.ones(xs.shape, complex)
        for n_ in range(32):
            acc = acc + term / factorial(n_ + j + 1)
            term = term * (-xs)
        out[j][small] = factorial(j) * acc
    xl = x[~small]
    p0 = -np.expm1(-xl) / xl
    out[0][~small] = p0
    prev = p0
    for j in range(1, jmax + 1):
        prev = (1.0 - j * prev) / xl
        out[j][~small] = prev
    return out


def _stencils(z, order):
    """Per cell k: stencil node indices (order 2: k, k+1; order 4: 4 nodes around the cell, shifted at the ends)."""
    N = z.size - 1
    if order == 2:
        return np.stack([np.arange(N), np.arange(1, N + 1)], 1)
    st = np.stack([np.arange(N) - 1 + i for i in range(4)], 1)
    st[0] = [0, 1, 2, 3]
    st[-1] = [N - 3, N - 2, N - 1, N]
    return st


def _cell_weights(z, st, direction, x, hv):
    """w[i, k, m] = (h_k/|v_i|) int_0^1 exp(-x_ik (1 - t)) l_m(t) dt, with t the normalised distance from the
    upstream node of cell k and l_m the Lagrange basis on the stencil nodes."""
    N = z.size - 1
    h = np.diff(z)
    nm = st.shape[1]
    up = z[:-1] if direction == "pos" else z[1:]
    sgn = 1.0 if direction == "pos" else -1.0
    tnodes = sgn * (z[st] - up[:, None]) / h[:, None]                 # (N, nm), signed along the propagation
    # l_m(t) = sum_j c[k, j, m] t^j with c = inv(V), V[k, m, j] = t_m^j
    c = np.linalg.inv(tnodes[:, :, None] ** np.arange(nm)[None, None, :])     # c[k, j, m]
    ps = _psi(x, nm - 1)                                               # (nm, n_modes, N)
    w = np.einsum("jik,kjm->ikm", ps, c) * hv[:, :, None]
    return w


def build_Kp(C, e, inplace=False, chunk=512):
    """K' = C_R - C = (D - C) - (De)(De)^T/W (frequency independent). inplace=True overwrites C (memory-lean)."""
    n = e.size
    D = np.diag(C).copy()
    De = D * e
    W = float(e @ De)
    Kp = C if inplace else C.copy()
    Kp *= -1.0
    Kp[np.diag_indices(n)] += D
    for i0 in range(0, n, chunk):
        Kp[i0:i0 + chunk] -= np.outer(De[i0:i0 + chunk], De) / W
    return Kp, D


class Slab:
    def __init__(self, C, v, e, p, z, vtol=1e-12, order=2, Kp=None, D=None):
        self.v, self.e, self.p, self.z = v, e, p, z
        n = e.size
        self.n = n
        self.N = z.size - 1
        self.order = order
        h = np.diff(z)
        if Kp is None:
            Kp, D = build_Kp(C, e)
        self.D = D
        self.W = float(e @ (D * e))
        self.De = D * e
        self.ipos = np.where(v > vtol)[0]
        self.ineg = np.where(v < -vtol)[0]
        self.iloc = np.where(np.abs(v) <= vtol)[0]
        self.Kp = Kp
        self.st = _stencils(z, order)
        self.dir = {}
        for tag, ii in (("pos", self.ipos), ("neg", self.ineg)):
            av = np.abs(v[ii])
            sig = (p + D[ii]) / av
            xk = sig[:, None] * h[None, :]
            hv = h[None, :] / av[:, None]
            E = np.exp(-xk)
            if order == 2:
                p1, p2 = _phi12(xk)
                w = np.stack([hv * (p1 - p2), hv * p2], 2)            # upstream node first
                stl = self.st if tag == "pos" else self.st[:, ::-1]
                self.dir[tag] = {"E": E, "w": w, "st": stl}
            else:
                w = _cell_weights(z, self.st, tag, xk, hv)
                self.dir[tag] = {"E": E, "w": w, "st": self.st}
        self.loc_fac = 1.0 / (p + D[self.iloc])
        # quadrature weights for the integral of a nodal field over the film (same interpolation, x = 0)
        qw = np.zeros(self.N + 1)
        wq = _cell_weights(z, self.st, "pos", np.zeros((1, self.N), complex), h[None, :]).real[0]
        for k in range(self.N):
            np.add.at(qw, self.st[k], wq[k])
        self.qw = qw
        self._build_A()

    # sweep of a given total source sig (n x (N+1)) with wall values (bc for pos at z=0, neg at z=d)
    def sweep(self, sig, TA=0.0, TB=0.0):
        N = self.N
        f = np.zeros((self.n, N + 1), complex)
        e = self.e
        for tag, ii in (("pos", self.ipos), ("neg", self.ineg)):
            dd = self.dir[tag]
            E, w, st = dd["E"], dd["w"], dd["st"]
            s_ = sig[ii]
            fd = np.zeros((ii.size, N + 1), complex)
            if tag == "pos":
                fd[:, 0] = TA * e[ii]
                for k in range(N):
                    fd[:, k + 1] = E[:, k] * fd[:, k] + np.einsum("im,im->i", w[:, k], s_[:, st[k]])
            else:
                fd[:, N] = TB * e[ii]
                for k in range(N - 1, -1, -1):
                    fd[:, k] = E[:, k] * fd[:, k + 1] + np.einsum("im,im->i", w[:, k], s_[:, st[k]])
            f[ii] = fd
        f[self.iloc] = self.loc_fac[:, None] * sig[self.iloc]
        return f

    def _build_A(self):
        """A[k, j] = (1/W) sum_i D_i e_i [sweep of D_i e_i hat_j]_i(z_k) (zero wall values)."""
        N = self.N
        A = np.zeros((N + 1, N + 1), complex)
        De = self.De
        for tag, ii in (("pos", self.ipos), ("neg", self.ineg)):
            dd = self.dir[tag]
            E, w, st = dd["E"], dd["w"], dd["st"]
            wt = De[ii] ** 2 / self.W
            F = np.zeros((ii.size, N + 1), complex)          # columns: source node j
            rng = range(N) if tag == "pos" else range(N - 1, -1, -1)
            for k in rng:
                F *= E[:, k][:, None]
                for m in range(st.shape[1]):
                    F[:, st[k, m]] += w[:, k, m]
                A[k + 1 if tag == "pos" else k] += wt @ F
        A[np.diag_indices(N + 1)] += np.sum(De[self.iloc] ** 2 * self.loc_fac) / self.W
        self.A = A
        self.luA = sla.lu_factor(np.eye(N + 1) - A, check_finite=False)

    def rta_solve(self, s, TA=0.0, TB=0.0):
        f0 = self.sweep(s, TA, TB)
        r = (self.De @ f0) / self.W
        th = sla.lu_solve(self.luA, r, check_finite=False)
        return f0 + self.sweep(np.outer(self.De, th)), th

    def _Kp_apply(self, f):
        fr = np.ascontiguousarray(f).view(np.float64)        # (n, 2(N+1)) interleaved
        return (self.Kp @ fr).view(np.complex128)

    def solve(self, TA, TB, rtol=1e-12, restart=60, maxiter=40, x0=None, rta_only=False):
        n, N1 = self.n, self.N + 1
        Fbc, _ = self.rta_solve(np.zeros((n, N1), complex), TA, TB)
        if rta_only:
            return Fbc, {"iters": 0}
        def mv(xv):
            f = xv.reshape(n, N1)
            g, _ = self.rta_solve(self._Kp_apply(f))
            return (f - g).ravel()
        Lop = LinearOperator((n * N1, n * N1), matvec=mv, dtype=complex)
        hist = []
        b = Fbc.ravel()
        sol, info = gmres(Lop, b, x0=None if x0 is None else x0.ravel(), rtol=rtol, atol=0.0, restart=restart,
                          maxiter=maxiter, callback=lambda r: hist.append(float(r)), callback_type="pr_norm")
        f = sol.reshape(n, N1)
        res = np.linalg.norm(mv(sol) - b) / np.linalg.norm(b)
        return f, {"info": int(info), "iters": len(hist), "true_relres": float(res),
                   "hist_last": hist[-5:] if hist else []}

    def outputs(self, f):
        v, e = self.v, self.e
        q = (v * e) @ f                                       # flux at every node
        U = e @ f
        Ct = float(e @ e)
        # energy balance with the trapezoidal integral of U (diagnostic of the z discretisation)
        intU = self.qw @ U
        return {"q0": q[0], "qd": q[-1], "th0": U[0] / Ct, "thd": U[-1] / Ct,
                "ebal": abs(q[0] - q[-1] - self.p * intU) / abs(q[0])}


def prolong(f_coarse, z_coarse, z_fine):
    """Initial guess on a bisected mesh: copy at the old nodes, linear interpolation at the midpoints."""
    out = np.empty((f_coarse.shape[0], z_fine.size), complex)
    for k in range(f_coarse.shape[0]):
        out[k] = np.interp(z_fine, z_coarse, f_coarse[k].real) + 1j * np.interp(z_fine, z_coarse, f_coarse[k].imag)
    return out


def twoport(C, v, e, p, z, rtol=1e-12, restart=60, maxiter=40, rta_only=False, x0s=None, order=2, Kp=None, D=None):
    t0 = time.perf_counter()
    S = Slab(C, v, e, p, z, order=order, Kp=Kp, D=D)
    tb = time.perf_counter() - t0
    out, infos, fs = [], [], []
    for k, (TA, TB) in enumerate(((1.0, 0.0), (0.0, 1.0))):
        f, inf = S.solve(TA, TB, rtol=rtol, restart=restart, maxiter=maxiter, rta_only=rta_only,
                         x0=None if x0s is None else x0s[k])
        out.append(S.outputs(f))
        infos.append(inf)
        fs.append(f)
    Y = np.array([[out[0]["q0"], out[1]["q0"]], [out[0]["qd"], out[1]["qd"]]])
    Th = np.array([[out[0]["th0"], out[1]["th0"]], [out[0]["thd"], out[1]["thd"]]])
    diag = {"ebal": max(o["ebal"] for o in out), "gmres": infos, "t_build": tb, "t_total": time.perf_counter() - t0,
            "N": S.N, "reciprocity": float(abs(Y[0, 1] + Y[1, 0]) / abs(Y[0, 1])),
            "symmetry": float(abs(Y[0, 0] + Y[1, 1]) / abs(Y[0, 0]))}
    return Y, Th, diag, fs
