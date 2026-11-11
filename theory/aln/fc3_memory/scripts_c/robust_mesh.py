"""Block B (T5, T6): robust-objective optimisation over several samplings of the operator.

One fc3 (scaled coordinates u of a per-orbit class, orbit_scaling.py) drives the exact linearised models of every mesh
in --meshes (5x5x3, 7x7x5, 9x9x6; 9x9x6 stored in float32). Constraints on every constraint mesh: all mode-orbit
lifetimes within +-dr and kappa_xx, kappa_zz within +-dK (block-averaged convention), per-orbit bound ||(W u)_o|| <= eps.
Objective (sense max or min) on the objective meshes: the smooth minimum (max) or smooth maximum (min) of the
log tau ratios, -1/beta log sum exp(-beta s_M), s_M = sgn log(tau_M / tau_M,0), i.e. the change that holds on every
objective mesh simultaneously. Augmented Lagrangian (PHR) with L-BFGS-B inner solves, as optimize_fc3.run; only iterates
feasible for the nominal tolerances are reported (inner bounds). Optional random start inside the box (multistart).

usage: python -B robust_mesh.py --tag orb|cal1nz|cal3w|... --eps E --col 0|2 --sense max|min
       --meshes m553,m775,m996 --objective m553,m775,m996 --dr 0.1 --dK 0.01 --out PREFIX
       [--start zero|rand:SEED|file:U.npy] [--maxiter 300] [--maxouter 15] [--beta 40] [--film L_nm] [--rtol 1e-4]
Warm start: with --start zero, an existing PREFIX_start.npy is used as the initial point (recorded as start_used).
Stopping: after outer >= 3, stop when two consecutive iterates are feasible for the nominal tolerances and the robust
log ratio changed by less than rtol (or, as before, violation < 1e-7 and ratio change < 1e-6 after outer >= 6).
--film L: the objective is the film-regularised memory time tau_mem(L) of C_L = C + diag(|v_z| / (4 pi L)) (E30,
film_event.py); constraints stay the bulk static data. D does not depend on the fc3, so dK_L/dg and dN_L/dg take the
bulk form with C_L in place of the regularised C.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "4")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from fc3class import FC3Class, FOURPI  # noqa: E402
from optimize_fc3 import FC3Model  # noqa: E402

RES = HERE.parent / "results_c"
MESH = {"m553": RES / "m553", "m775": RES / "m775x", "m996": RES / "m996"}
EVENTS = {"m553": "events_m553_s0.1.npz", "m775": "events_m775_s0.1.npz", "m996": "events_m996_s0.1.npz"}


class MeshModel(FC3Model):
    """Linearised model of one mesh in scaled coordinates (float32 storage allowed for large meshes)."""

    def __init__(self, mesh, tag, film=None):
        rdir = MESH[mesh]
        self.M = M = FC3Class(str(rdir / EVENTS[mesh]), rdir / "pp_r01_REF.npz")
        lin = "linearization.npz" if (rdir / "linearization.npz").exists() else "lin_meta.npz"
        Z = np.load(rdir / lin)
        pre = "scaled_" if tag == "orb" else f"scaled_{tag}_"
        self.Gr = np.load(rdir / f"{pre}Gr.npy", mmap_mode=None)
        self.Gi = np.load(rdir / f"{pre}Gi.npy", mmap_mode=None)
        self.n_var = self.Gr.shape[1]
        self.Phi0 = Z["Phi0"]
        self.row_orb = Z["row_orb"]
        self.allowed = Z["allowed"]
        self.uc = float(Z["uc"])
        self.w_row = Z["W_orb"][self.row_orb] * self.uc / Z["blk"][self.row_orb]
        self.g_fixed = np.where(self.allowed, 0.0, Z["g_ref_orb"])
        self.n_orb = M.n_orb
        self.ops = M.ops
        self.reps = M.tie["reps"]
        self.Wr = M.geom.W[self.reps].tocsr()
        self.e = M.geom.ehat
        self.b = M.geom.b
        self.dfilm = None if film is None else np.abs(M.geom.vel[:, 2]) / (FOURPI * film * 10.0)
        base = self.evaluate(np.zeros(self.n_var))
        self.r0, self.K0, self.base = base["r"], {0: base["K"][0], 2: base["K"][2]}, base

    def evaluate(self, y, cols=(0, 2), grad=False):
        out = super().evaluate(y, cols, grad)
        if self.dfilm is None:
            return out
        C = self.ops.C(out["g_orb"][self.M.orb])
        C[np.diag_indices_from(C)] += self.dfilm
        cf = sla.cho_factor(C, lower=True, check_finite=False)
        for k in ("Kf", "Nf", "tauf", "dKf", "dNf"):
            out[k] = {}
        for c in cols:
            b = self.b[:, c]
            x = sla.cho_solve(cf, b, check_finite=False)
            K, N = float(b @ x), float(x @ x)
            out["Kf"][c], out["Nf"][c], out["tauf"][c] = K, N, N / K
            if grad:
                w = sla.cho_solve(cf, x, check_finite=False)
                p, q = self.ops.ATx(x), self.ops.ATx(w)
                out["dKf"][c] = np.bincount(self.M.orb, weights=-(p * p), minlength=self.n_orb)
                out["dNf"][c] = np.bincount(self.M.orb, weights=-2 * p * q, minlength=self.n_orb)
        return out

    def orbit_rates(self, y):
        yy = y.astype(self.Gr.dtype, copy=False)
        Phi = self.Phi0 + (self.Gr @ yy).astype(float) + 1j * (self.Gi @ yy).astype(float)
        g = self.g_fixed + np.bincount(self.row_orb, weights=self.w_row * (Phi.real ** 2 + Phi.imag ** 2),
                                       minlength=self.n_orb)
        return g, Phi

    def chain(self, Dorb, Phi):
        v = self.w_row * Dorb[self.row_orb] * np.conj(Phi)
        vr = v.real.astype(self.Gr.dtype)
        vi = v.imag.astype(self.Gr.dtype)
        return 2.0 * ((vr @ self.Gr).astype(float) - (vi @ self.Gi).astype(float))


def run(models, objm, consm, Wmat, col_block, eps, col, sense, dr, dK, beta, u0, maxouter, maxiter, logfile, film=False,
        rtol=1e-4):
    T, NN, KK, DN, DK = ("tauf", "Nf", "Kf", "dNf", "dKf") if film else ("tau", "N", "K", "dN", "dK")
    n = len(u0)
    nb = int(col_block.max()) + 1
    u = u0.copy()
    dr_i, dK_i, oe_i = dr * (1 - 2e-3), dK * (1 - 2e-2), eps * (1 - 1e-3)
    lam = {m: {"rp": np.zeros(len(models[m].r0)), "rm": np.zeros(len(models[m].r0)), "Kp": np.zeros(2),
               "Km": np.zeros(2)} for m in consm}
    lam_o = np.zeros(nb)
    mu, sgn = 50.0, (1.0 if sense == "max" else -1.0)
    hist, best, t0 = [], None, time.perf_counter()

    def rel2_of(u):
        Wu = Wmat @ u
        return np.bincount(col_block, weights=Wu * Wu, minlength=nb), Wu

    def cons(m, ev):
        md = models[m]
        s = ev["r"] / md.r0 - 1.0
        kk = np.array([ev["K"][0] / md.K0[0] - 1.0, ev["K"][2] / md.K0[2] - 1.0])
        return {"rp": s - dr_i, "rm": -s - dr_i, "Kp": kk - dK_i, "Km": -kk - dK_i}, s, kk

    def robust(svals):
        # smooth min of s_M (s already signed); returns value and weights
        a = -beta * np.array(svals)
        amax = a.max()
        w = np.exp(a - amax)
        val = -(amax + np.log(w.sum())) / beta
        return val, w / w.sum()

    def fun(u):
        f, gu = 0.0, np.zeros(n)
        evs = {m: models[m].evaluate(u, grad=True) for m in models}
        svals = [sgn * np.log(evs[m][T][col] / models[m].base[T][col]) for m in objm]
        val, wts = robust(svals)
        f -= val
        for m in models:
            md, ev = models[m], evs[m]
            Dorb = np.zeros(md.n_orb)
            if m in objm:
                w = wts[objm.index(m)]
                Dorb -= w * sgn * (ev[DN][col] / ev[NN][col] - ev[DK][col] / ev[KK][col])
            if m in consm:
                h, _, _ = cons(m, ev)
                L = lam[m]
                wp, wm = np.maximum(0, L["rp"] + mu * h["rp"]), np.maximum(0, L["rm"] + mu * h["rm"])
                kp, km = np.maximum(0, L["Kp"] + mu * h["Kp"]), np.maximum(0, L["Km"] + mu * h["Km"])
                f += (np.sum(wp ** 2 - L["rp"] ** 2) + np.sum(wm ** 2 - L["rm"] ** 2)
                      + np.sum(kp ** 2 - L["Kp"] ** 2) + np.sum(km ** 2 - L["Km"] ** 2)) / (2 * mu)
                Dorb += md.dr_dorb((wp - wm) / md.r0)
                for i, c in enumerate((0, 2)):
                    Dorb += (kp[i] - km[i]) * ev["dK"][c] / md.K0[c]
            gu += md.chain(Dorb, ev["Phi"])
        rel2, Wu = rel2_of(u)
        ho = rel2 / oe_i ** 2 - 1.0
        wo = np.maximum(0.0, lam_o + mu * ho)
        f += np.sum(wo ** 2 - lam_o ** 2) / (2 * mu)
        gu += 2.0 * (Wmat.T @ ((wo / oe_i ** 2)[col_block] * Wu))
        return f, gu

    for it in range(maxouter):
        res = minimize(fun, u, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-13, "gtol": 1e-9})
        u = res.x
        rec = {"outer": it, "nit": int(res.nit), "mu": mu, "mesh": {}}
        viol, feas, svals = 0.0, True, []
        for m in models:
            ev = models[m].evaluate(u)
            h, s, kk = cons(m, ev)
            tr = [ev["tau"][c] / models[m].base["tau"][c] for c in (0, 2)]
            rec["mesh"][m] = {"tau_ratio": tr, "kappa_ratio": [1 + kk[0], 1 + kk[1]], "max_dr": float(np.max(np.abs(s)))}
            if film:
                rec["mesh"][m]["film_tau_ratio"] = [ev[T][c] / models[m].base[T][c] for c in (0, 2)]
            if m in objm:
                svals.append(np.log(ev[T][col] / models[m].base[T][col]))
            if m in consm:
                viol = max(viol, max(float(np.max(h[k])) for k in h))
                feas &= bool(np.max(np.abs(s)) <= dr and np.max(np.abs(kk)) <= dK)
                for k in h:
                    lam[m][k] = np.maximum(0.0, lam[m][k] + mu * h[k])
        rec["objective_logratio_minmax"] = [float(min(svals)), float(max(svals))]
        rec["robust_ratio"] = float(np.exp(min(svals) if sense == "max" else max(svals)))
        rel2 = rel2_of(u)[0]
        ho = rel2 / oe_i ** 2 - 1.0
        viol = max(viol, float(np.max(ho)))
        rec["orbit_rel_over_eps_max"] = float(np.sqrt(np.max(rel2)) / eps)
        feas &= rec["orbit_rel_over_eps_max"] <= 1.0
        lam_o[:] = np.maximum(0.0, lam_o + mu * ho)
        rec.update({"viol": viol, "feasible": bool(feas), "t": time.perf_counter() - t0})
        hist.append(rec)
        print(json.dumps(rec), flush=True)
        with open(logfile, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        score = np.log(rec["robust_ratio"]) * (1 if sense == "max" else -1)
        if feas and (best is None or score > best[0]):
            best = (score, u.copy(), rec)
        if len(hist) > 1 and viol > 0.25 * hist[-2]["viol"] and viol > 1e-7:
            mu = min(mu * 4, 1e8)
        if it >= 6 and viol < 1e-7 and abs(hist[-1]["robust_ratio"] - hist[-2]["robust_ratio"]) < 1e-6:
            break
        if (it >= 3 and feas and hist[-2]["feasible"]
                and abs(np.log(hist[-1]["robust_ratio"] / hist[-2]["robust_ratio"])) < rtol):
            break
    return best, hist


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--eps", type=float, required=True)
    ap.add_argument("--col", type=int, required=True)
    ap.add_argument("--sense", choices=["max", "min"], required=True)
    ap.add_argument("--meshes", required=True)
    ap.add_argument("--objective", required=True)
    ap.add_argument("--dr", type=float, default=0.1)
    ap.add_argument("--dK", type=float, default=0.01)
    ap.add_argument("--beta", type=float, default=40.0)
    ap.add_argument("--start", default="zero")
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--maxouter", type=int, default=15)
    ap.add_argument("--out", required=True)
    ap.add_argument("--film", type=float, default=None)
    ap.add_argument("--rtol", type=float, default=1e-4)
    a = ap.parse_args()
    consm = a.meshes.split(",")
    objm = a.objective.split(",")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    models = {m: MeshModel(m, a.tag, a.film) for m in sorted(set(consm) | set(objm))}
    pre = "scaled_" if a.tag == "orb" else f"scaled_{a.tag}_"
    Wmat = np.load(RES / "m553" / f"{pre}W.npy")
    col_block = np.load(RES / "m553" / "cluster_orbits.npz")["col_block"]
    n = models[objm[0]].n_var
    start_used = a.start
    if a.start == "zero" and Path(str(out) + "_start.npy").exists():
        start_used = "file:" + str(out) + "_start.npy"
    if start_used.startswith("file:"):
        u0 = np.load(start_used[5:]).astype(float)
        assert u0.shape == (n,), (u0.shape, n)
    elif a.start == "zero":
        u0 = np.zeros(n)
    else:
        rng = np.random.default_rng(int(a.start.split(":")[1]))
        u0 = rng.standard_normal(n)
        nb = int(col_block.max()) + 1
        r = np.sqrt(np.bincount(col_block, weights=(Wmat @ u0) ** 2, minlength=nb))
        u0 *= 0.5 * a.eps / r.max()            # inside the box, largest orbit at half its bound
    best, hist = run(models, objm, consm, Wmat, col_block, a.eps, a.col, a.sense, a.dr, a.dK, a.beta, u0,
                     a.maxouter, a.maxiter, str(out) + ".log", film=a.film is not None, rtol=a.rtol)
    summ = {"args": vars(a), "start_used": start_used, "base_tau_ps": {m: [models[m].base["tau"][0] / FOURPI, models[m].base["tau"][2] / FOURPI]
                                             for m in models},
            "best": best[2] if best else None, "history": hist}
    if a.film is not None:
        summ["base_film_tau_ps"] = {m: [models[m].base["tauf"][0] / FOURPI, models[m].base["tauf"][2] / FOURPI]
                                    for m in models}
    if best:
        np.save(str(out) + "_u.npy", best[1])
    Path(str(out) + ".json").write_text(json.dumps(summ, indent=1))
