"""Augmented-Lagrangian L-BFGS on theta = log(gamma) for extremizing log N on the feasible slice.

f(theta) = -sign*log N + sum_i lam_i c_i + (mu/2) sum_i c_i^2,
c = (W gamma / r - 1, (K_k - K0_k)/scale_k).  Outer multiplier updates; the final iterate is
projected exactly onto F by MultiSlice.restore (Newton), so the reported value is attained by a
feasible witness. Works with optim_multi.MultiSlice.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize


def alb_run(S, gam0, sign, outer=8, inner=300, mu0=1e2, verbose=False):
    theta = np.log(np.maximum(gam0, 1e-300))
    nr = S.nr
    kp = S.kpairs
    lam = np.zeros(nr + len(kp))
    mu = mu0
    best = None

    def fg(th):
        gam = np.exp(np.clip(th, -700, 700))
        try:
            R = S.evaluate(gam, grad=True)
        except Exception:
            return 1e30, np.zeros_like(th)
        c_d = S.ops.Wv(gam) / S.rs - 1.0
        c_k = S.kres(R)
        c = np.concatenate([c_d, c_k])
        val = -sign * math.log(R["N"]) + lam @ c + 0.5 * mu * (c @ c)
        w = lam + mu * c
        # gradient wrt gamma
        gg = -sign * R["dN"] / R["N"]
        gg = gg + S.ops.WTv(w[:nr] / S.rs)
        for j, k in enumerate(kp):
            gg = gg + w[nr + j] * R["dK"][k] / S.kscale[k]
        return val, gg * gam

    for o in range(outer):
        res = minimize(fg, theta, jac=True, method="L-BFGS-B", options={"maxiter": inner, "maxcor": 30})
        theta = res.x
        gam = np.exp(theta)
        R = S.evaluate(gam, grad=False)
        c = np.concatenate([S.ops.Wv(gam) / S.rs - 1.0, S.kres(R)])
        cn = float(np.max(np.abs(c)))
        if verbose:
            print(f"outer {o}: tau={R['tau']:.6g} viol={cn:.2e} mu={mu:.1e} it={res.nit}", flush=True)
        # try exact restoration from the current point
        g2, ok, R2 = S.restore(gam)
        if ok:
            if best is None or sign * (R2["tau"] - best["tau"]) > 0:
                best = {"g": g2, "tau": R2["tau"], "N": R2["N"], "K": R2["K"], "feas": S.feas(g2, R2), "outer": o}
        lam = lam + mu * c
        if cn > 1e-4:
            mu *= 10.0
    return best
