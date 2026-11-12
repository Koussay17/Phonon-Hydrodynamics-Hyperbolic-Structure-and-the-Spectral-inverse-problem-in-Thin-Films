"""Pseudo-arclength continuation of twin roots F(th', lam) = 0 (6 real equations, th' in R^6, lam scalar).

The homotopy parameter lam can move the convention (eps0, eps), the target, or the frequencies.
Signs of th' are not constrained, so a root can be followed through tau < 0 and back; physical
admissibility (all entries >= 0, xi, b > 0) is recorded along the path.
"""
import numpy as np
import cxmodel as cx


class Homotopy:
    """F(th', lam) = [Re; Im](Z(omega(lam); th', spec(lam)) / Z(omega(lam); target(lam), spec(lam)) - 1)."""

    def __init__(self, spec_of, target_of, omega_of):
        self.spec_of, self.target_of, self.omega_of = spec_of, target_of, omega_of

    def F(self, th, lam):
        sp = self.spec_of(lam)
        om = np.asarray(self.omega_of(lam), dtype=float)
        Zt = cx.Z_np(sp, self.target_of(lam), om)
        z, J = cx.ZJ_np(sp, th, om)
        R = z / Zt - 1
        G = J / Zt[None, :]
        return np.concatenate([R.real, R.imag]), np.concatenate([G.real.T, G.imag.T], axis=0)

    def FJ(self, x, hl=1e-6):
        th, lam = x[:-1], x[-1]
        F, Jt = self.F(th, lam)
        Fp, _ = self.F(th, lam + hl)
        Fm, _ = self.F(th, lam - hl)
        Jl = (Fp - Fm) / (2 * hl)
        return F, np.concatenate([Jt, Jl[:, None]], axis=1)


def tangent(J, prev=None):
    # null vector of the 6 x 7 Jacobian
    u, s, vt = np.linalg.svd(J)
    t = vt[-1]
    if prev is not None and np.dot(t, prev) < 0:
        t = -t
    return t / np.linalg.norm(t)


def correct(H, xp, t, tol=1e-12, maxit=12):
    x = xp.copy()
    for it in range(maxit):
        F, J = H.FJ(x)
        g = np.concatenate([F, [np.dot(t, x - xp)]])
        A = np.vstack([J, t[None, :]])
        try:
            dx = np.linalg.solve(A, -g)
        except np.linalg.LinAlgError:
            return x, False
        x = x + dx
        if np.linalg.norm(dx) < tol * (1 + np.linalg.norm(x)) and np.linalg.norm(F) < 1e-10:
            F, _ = H.FJ(x)
            return x, np.linalg.norm(F) < 1e-10
    F, _ = H.FJ(x)
    return x, np.linalg.norm(F) < 1e-10


def follow(H, th0, lam0, lam_end, h0=0.01, hmax=0.05, hmin=1e-7, maxsteps=20000, target_tol=1e-6):
    """Follow the branch through (th0, lam0) in the direction of increasing (lam_end > lam0) lam.

    Returns dict with the list of (lam, th) samples, events, and the end state.
    """
    x = np.concatenate([np.asarray(th0, float), [lam0]])
    F, J = H.FJ(x)
    t = tangent(J)
    sgn = np.sign(lam_end - lam0)
    if t[-1] * sgn < 0:
        t = -t
    h = h0
    samples = [(float(lam0), x[:-1].tolist())]
    events = []
    lam_dir = np.sign(t[-1])
    for step in range(maxsteps):
        xp = x + h * t
        xn, ok = correct(H, xp, t)
        if not ok or np.linalg.norm(xn - x) > 3 * h:
            h /= 2
            if h < hmin:
                events.append(("stalled", float(x[-1])))
                break
            continue
        Fn, Jn = H.FJ(xn)
        tn = tangent(Jn, t)
        # fold in lam
        if np.sign(tn[-1]) != lam_dir and abs(tn[-1]) > 0:
            events.append(("fold", float(xn[-1])))
            lam_dir = np.sign(tn[-1])
        # crossing of the trivial branch (twin passes through the target)
        tgt = np.asarray(H.target_of(xn[-1]), float)
        dist = np.max(np.abs(xn[:-1] - tgt))
        if dist < target_tol:
            events.append(("touches target", float(xn[-1])))
        x, t = xn, tn
        samples.append((float(x[-1]), x[:-1].tolist()))
        h = min(h * 1.4, hmax)
        lam = x[-1]
        if (lam - lam_end) * sgn >= 0:
            break
        if np.max(np.abs(x[:-1])) > 1e4:
            events.append(("escapes to infinity", float(lam)))
            break
        if lam_dir * sgn < 0 and len([e for e in events if e[0] == "fold"]) > 6:
            events.append(("too many folds", float(lam)))
            break
    return dict(samples=samples, events=events, end=x.tolist())


def admissible(th, spec):
    th = np.asarray(th)
    return bool(np.all(th >= 0) and th[0] > 0 and th[1] > 0)


def land(H, path, lam_star, tol=1e-12):
    """Newton at lam = lam_star from every sample bracketing lam_star (all crossings)."""
    import cxtwin as ct
    sols = []
    S = path["samples"]
    for (l1, t1), (l2, t2) in zip(S[:-1], S[1:]):
        if (l1 - lam_star) * (l2 - lam_star) <= 0 and l1 != l2:
            w = (lam_star - l1) / (l2 - l1)
            th0 = (1 - w) * np.array(t1) + w * np.array(t2)
            sp = H.spec_of(lam_star)
            om = np.asarray(H.omega_of(lam_star), float)
            Zt = cx.Z_np(sp, H.target_of(lam_star), om)
            th, nf, it, ok = ct.newton(sp, th0, om, Zt, tol=tol, logvars=False)
            if ok:
                sols.append(th.tolist())
    return sols
