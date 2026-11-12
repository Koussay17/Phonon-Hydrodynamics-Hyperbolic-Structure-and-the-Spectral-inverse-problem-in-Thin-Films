"""Twin-search tools on the independent model (cxmodel).

Equations: F(th') = [Re; Im](Z(omega_j; th')/Z(omega_j; th*) - 1), j = 1..N  (relative residual).
Unknowns in log variables y = log th' (positivity) unless stated.
"""
import numpy as np
import mpmath as mp
import cxmodel as cx


def resid_jac(spec, th, omega, Zt):
    """Relative residual and Jacobian w.r.t. th (not log)."""
    z, J = cx.ZJ_np(spec, th, omega)
    F = z / Zt - 1
    G = J / Zt[None, :]
    return np.concatenate([F.real, F.imag]), np.concatenate([G.real.T, G.imag.T], axis=0)


def newton(spec, th0, omega, Zt, tol=1e-14, maxit=60, logvars=True):
    th = np.array(th0, dtype=float)
    for it in range(maxit):
        F, J = resid_jac(spec, th, omega, Zt)
        nf = np.linalg.norm(F)
        if nf < tol:
            return th, nf, it, True
        if logvars:
            Jl = J * th[None, :]
            try:
                dy = np.linalg.lstsq(Jl, -F, rcond=None)[0]
            except np.linalg.LinAlgError:
                return th, nf, it, False
            step = 1.0
            while step > 1e-4:
                tn = th * np.exp(step * dy)
                Fn, _ = resid_jac(spec, tn, omega, Zt)
                if np.linalg.norm(Fn) < nf:
                    break
                step /= 2
            th = tn
        else:
            dx = np.linalg.lstsq(J, -F, rcond=None)[0]
            th = th + dx
    F, J = resid_jac(spec, th, omega, Zt)
    return th, np.linalg.norm(F), maxit, np.linalg.norm(F) < tol


# ---------------------------------------------------------------------------------------------
# batched Levenberg-Marquardt in log variables (many starts at once)
# ---------------------------------------------------------------------------------------------
def batch_lm(spec, Y0, omega, Zt, iters=200, lam0=1e-2, lo=-12.0, hi=6.0):
    """Y0: (m, n) start points in log variables. Returns Y, residual norms."""
    Y = np.array(Y0, dtype=float)
    m, n = Y.shape
    lam = np.full(m, lam0)

    def ev(Y):
        TH = np.exp(Y)
        Z, J = cx.ZJ_np_batch(spec, TH, omega)       # Z (m, nf), J (m, n, nf)
        F = Z / Zt[None, :] - 1
        G = J / Zt[None, None, :] * TH[:, :, None]   # d/dy
        Fr = np.concatenate([F.real, F.imag], axis=1)            # (m, 2nf)
        Gr = np.concatenate([G.real, G.imag], axis=2)            # (m, n, 2nf)
        return Fr, Gr

    with np.errstate(all="ignore"):
        F, G = ev(Y)
        cost = np.sum(F * F, axis=1)
        for it in range(iters):
            JTJ = np.einsum("mik,mjk->mij", G, G)
            JTF = np.einsum("mik,mk->mi", G, F)
            dia = np.einsum("mii->mi", JTJ)
            A = JTJ + lam[:, None, None] * (np.eye(n)[None] * (dia[:, :, None] + 1e-12))
            try:
                dY = -np.linalg.solve(A, JTF[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                dY = np.zeros_like(Y)
            dY = np.clip(dY, -2.0, 2.0)
            Yn = np.clip(Y + dY, lo, hi)
            Fn, Gn = ev(Yn)
            cn = np.sum(Fn * Fn, axis=1)
            bad = ~np.isfinite(cn)
            cn[bad] = np.inf
            acc = cn < cost
            Y[acc] = Yn[acc]
            F[acc] = Fn[acc]
            G[acc] = Gn[acc]
            cost[acc] = cn[acc]
            lam[acc] = np.maximum(lam[acc] / 3, 1e-12)
            lam[~acc] = np.minimum(lam[~acc] * 4, 1e8)
    return Y, np.sqrt(cost)


# ---------------------------------------------------------------------------------------------
# mpmath refinement and residual
# ---------------------------------------------------------------------------------------------
def refine_mp(spec, th0, omega, th_target, dps=60, maxit=40):
    mp.mp.dps = dps
    th = [mp.mpf(float(x)) for x in th0]
    Zt = [cx.Z_mp(spec, th_target, mp.mpf(w)) for w in omega]
    for it in range(maxit):
        F, J = [], []
        for w, zt in zip(omega, Zt):
            z, d = cx.ZJ_mp(spec, th, mp.mpf(w))
            F.append((z / zt - 1))
            J.append([x / zt for x in d])
        Fr = [f.real for f in F] + [f.imag for f in F]
        Jr = [[x.real for x in row] for row in J] + [[x.imag for x in row] for row in J]
        nF = mp.sqrt(mp.fsum(x * x for x in Fr))
        if nF < mp.mpf(10) ** (-dps + 8):
            break
        Jm = mp.matrix(Jr)
        if Jm.rows == Jm.cols:
            dx = mp.lu_solve(Jm, mp.matrix([-x for x in Fr]))
        else:
            dx = mp.qr_solve(Jm, mp.matrix([-x for x in Fr]))[0]
        th = [th[k] + dx[k] for k in range(len(th))]
    # final
    F = []
    for w, zt in zip(omega, Zt):
        z = cx.Z_mp(spec, th, mp.mpf(w))
        F.append(z / zt - 1)
    nF = mp.sqrt(mp.fsum(abs(f) ** 2 for f in F))
    return th, nF


def jac_cond_mp(spec, th, omega, Zt):
    J = []
    for w, zt in zip(omega, Zt):
        z, d = cx.ZJ_mp(spec, th, mp.mpf(w))
        J.append([x / zt for x in d])
    Jr = [[x.real for x in row] for row in J] + [[x.imag for x in row] for row in J]
    S = mp.svd_r(mp.matrix(Jr), compute_uv=False)
    sv = sorted([S[i] for i in range(len(S))], reverse=True)
    return sv
