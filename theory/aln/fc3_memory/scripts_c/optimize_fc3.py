"""Task 3: range of tau_mem over symmetric fc3 perturbations that keep the 300 K lifetimes and kappa
(C: mirror, consolidated version of the D: optimiser).

y = c/|c0|; Phi_m(y) = Phi0_m + |c0| G_m y (exact, linearize.py); g_O(y) = W_O uc mean_block |Phi_m|^2 on allowed
orbits, selection-rule orbits fixed. Base point y = 0 (block-averaged reference).
Problem: extremise tau_col(y) subject to |r_mu/r_mu(0) - 1| <= dr (mode-orbit representatives),
|K_aa/K_aa(0) - 1| <= dK (xx, zz), optional |y| <= rho (relative Frobenius size), optional per-cluster-orbit
relative bound ||d fc3_o|| <= eps ||fc3_0,o|| for each of the 456 symmetry-distinct cluster orbits.
Method: augmented Lagrangian (PHR) with L-BFGS-B inner solves; internal margins 0.998 dr, 0.98 dK, (1 - 1e-4) rho,
(1 - 1e-3) eps; only iterates feasible for the NOMINAL tolerances are reported (inner bounds on the range).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "1")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from fc3class import FC3Class, FOURPI  # noqa: E402


class FC3Model:
    def __init__(self, ev_path, rdir):
        self.M = M = FC3Class(ev_path, Path(rdir) / "pp_r01_REF.npz")
        Z = np.load(Path(rdir) / "linearization.npz")
        self.c0n = float(np.linalg.norm(Z["c0"]))
        Gs = Z["G"] * self.c0n
        self.Gr = np.ascontiguousarray(Gs.real)
        self.Gi = np.ascontiguousarray(Gs.imag)
        del Gs
        self.n_var = self.Gr.shape[1]
        self.Phi0 = Z["Phi0"]
        self.row_orb = Z["row_orb"]
        blk = Z["blk"]
        self.allowed = Z["allowed"]
        self.uc = float(Z["uc"])
        self.w_row = Z["W_orb"][self.row_orb] * self.uc / blk[self.row_orb]
        self.g_fixed = np.where(self.allowed, 0.0, Z["g_ref_orb"])
        self.n_orb = M.n_orb
        self.ops = M.ops
        self.reps = M.tie["reps"]
        self.Wr = M.geom.W[self.reps].tocsr()
        self.e = M.geom.ehat
        self.b = M.geom.b
        self.Bmat = None
        base = self.evaluate(np.zeros(self.n_var))
        self.r0 = base["r"]
        self.K0 = {0: base["K"][0], 2: base["K"][2]}
        self.base = base

    def load_orbit_structure(self, basisdir):
        CO = np.load(Path(basisdir) / "cluster_orbits.npz")
        self.col_block = CO["col_block"]
        self.n_blocks = int(self.col_block.max()) + 1
        self.block_refnorm2 = CO["block_refnorm"] ** 2
        self.Bmat = np.load(Path(basisdir) / "basis_B.npy")

    def orbit_rel2(self, y):
        By = self.Bmat @ y
        return self.c0n ** 2 * np.bincount(self.col_block, weights=By * By, minlength=self.n_blocks) / self.block_refnorm2, By

    def orbit_rates(self, y):
        Phi = self.Phi0 + (self.Gr @ y) + 1j * (self.Gi @ y)
        g = self.g_fixed + np.bincount(self.row_orb, weights=self.w_row * (Phi.real ** 2 + Phi.imag ** 2),
                                       minlength=self.n_orb)
        return g, Phi

    def evaluate(self, y, cols=(0, 2), grad=False):
        g_orb, Phi = self.orbit_rates(y)
        g = g_orb[self.M.orb]
        C = self.ops.C(g)
        sc = np.trace(C) / len(C)
        cf = sla.cho_factor(C + sc * np.outer(self.e, self.e), lower=True, check_finite=False)
        out = {"g_orb": g_orb, "Phi": Phi, "K": {}, "N": {}, "tau": {}, "dK": {}, "dN": {}}
        for c in cols:
            b = self.b[:, c]
            x = sla.cho_solve(cf, b, check_finite=False)
            K, N = float(b @ x), float(x @ x)
            out["K"][c], out["N"][c], out["tau"][c] = K, N, N / K
            if grad:
                u = sla.cho_solve(cf, x, check_finite=False)
                p, q = self.ops.ATx(x), self.ops.ATx(u)
                out["dK"][c] = np.bincount(self.M.orb, weights=-(p * p), minlength=self.n_orb)
                out["dN"][c] = np.bincount(self.M.orb, weights=-2 * p * q, minlength=self.n_orb)
        out["r"] = self.Wr @ g
        return out

    def chain(self, Dorb, Phi):
        v = self.w_row * Dorb[self.row_orb] * np.conj(Phi)
        return 2.0 * (v.real @ self.Gr - v.imag @ self.Gi)

    def dr_dorb(self, wmode):
        return np.bincount(self.M.orb, weights=self.Wr.T @ wmode, minlength=self.n_orb)


def run(model, col, sense, rho, dr, dK, orbit_eps=None, maxouter=15, maxiter=200, logfile=None, y_init=None):
    n = model.n_var
    y = np.zeros(n) if y_init is None else np.array(y_init, float)
    dr_i, dK_i = dr * (1 - 2e-3), dK * (1 - 2e-2)
    rho_i = rho * (1 - 1e-4) if np.isfinite(rho) else np.inf
    use_orb = orbit_eps is not None
    oe_i = orbit_eps * (1 - 1e-3) if use_orb else None
    nm = len(model.r0)
    lam = {"rp": np.zeros(nm), "rm": np.zeros(nm), "Kp": np.zeros(2), "Km": np.zeros(2), "rho": 0.0,
           "orb": np.zeros(model.n_blocks) if use_orb else None}
    mu = 50.0
    sgn = -1.0 if sense == "max" else 1.0
    other = 2 if col == 0 else 0
    hist, t0, best = [], time.perf_counter(), None

    def cons(ev, y, rel2):
        s = ev["r"] / model.r0 - 1.0
        kk = np.array([ev["K"][0] / model.K0[0] - 1.0, ev["K"][2] / model.K0[2] - 1.0])
        h = {"rp": s - dr_i, "rm": -s - dr_i, "Kp": kk - dK_i, "Km": -kk - dK_i,
             "rho": (y @ y - rho_i ** 2) if np.isfinite(rho_i) else -1.0}
        if use_orb:
            h["orb"] = rel2 / oe_i ** 2 - 1.0
        return h, s, kk

    def fun(y):
        ev = model.evaluate(y, grad=True)
        rel2, By = model.orbit_rel2(y) if use_orb else (None, None)
        f = sgn * np.log(ev["tau"][col])
        Dorb = sgn * (ev["dN"][col] / ev["N"][col] - ev["dK"][col] / ev["K"][col])
        h, s, kk = cons(ev, y, rel2)
        wp = np.maximum(0.0, lam["rp"] + mu * h["rp"])
        wm = np.maximum(0.0, lam["rm"] + mu * h["rm"])
        f += (np.sum(wp ** 2 - lam["rp"] ** 2) + np.sum(wm ** 2 - lam["rm"] ** 2)) / (2 * mu)
        Dorb = Dorb + model.dr_dorb((wp - wm) / model.r0)
        kp = np.maximum(0.0, lam["Kp"] + mu * h["Kp"])
        km = np.maximum(0.0, lam["Km"] + mu * h["Km"])
        f += (np.sum(kp ** 2 - lam["Kp"] ** 2) + np.sum(km ** 2 - lam["Km"] ** 2)) / (2 * mu)
        for i, c in enumerate((0, 2)):
            Dorb = Dorb + (kp[i] - km[i]) * ev["dK"][c] / model.K0[c]
        gy = model.chain(Dorb, ev["Phi"])
        if np.isfinite(rho_i):
            wr = max(0.0, lam["rho"] + mu * h["rho"])
            f += (wr ** 2 - lam["rho"] ** 2) / (2 * mu)
            gy += wr * 2 * y
        if use_orb:
            wo = np.maximum(0.0, lam["orb"] + mu * h["orb"])
            f += np.sum(wo ** 2 - lam["orb"] ** 2) / (2 * mu)
            gy += 2.0 * (model.Bmat.T @ ((wo * model.c0n ** 2 / (model.block_refnorm2 * oe_i ** 2))[model.col_block] * By))
        return f, gy

    for it in range(maxouter):
        res = minimize(fun, y, jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-13, "gtol": 1e-9})
        y = res.x
        ev = model.evaluate(y)
        rel2 = model.orbit_rel2(y)[0] if use_orb else None
        h, s, kk = cons(ev, y, rel2)
        viol = max([float(np.max(h[k])) for k in ("rp", "rm", "Kp", "Km")]
                   + [float(h["rho"]) if np.isfinite(rho_i) else -1.0, float(np.max(h["orb"])) if use_orb else -1.0, 0.0])
        orb_max = float(np.sqrt(np.max(rel2))) if use_orb else None
        # orbit_eps may be a scalar or a per-orbit array (calibrated class); identical test for a scalar
        orb_ratio = float(np.max(np.sqrt(rel2) / orbit_eps)) if use_orb else None
        feas = (np.max(np.abs(s)) <= dr) and (np.max(np.abs(kk)) <= dK) and \
            ((not np.isfinite(rho)) or np.linalg.norm(y) <= rho) and ((not use_orb) or orb_ratio <= 1.0)
        rec = {"outer": it, "tau": ev["tau"][col] / FOURPI, "tau_other": ev["tau"][other] / FOURPI, "viol": viol,
               "max_dr": float(np.max(np.abs(s))), "max_dK": float(np.max(np.abs(kk))), "frob": float(np.linalg.norm(y)),
               "orbit_rel_max": orb_max, "orbit_rel_over_eps_max": orb_ratio, "mu": mu, "nit": int(res.nit), "feasible": bool(feas), "t": time.perf_counter() - t0}
        hist.append(rec)
        if logfile:
            with open(logfile, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
        if feas and (best is None or sgn * np.log(ev["tau"][col]) < sgn * np.log(best[0])):
            best = (ev["tau"][col], y.copy(), rec)
        for k in ("rp", "rm", "Kp", "Km"):
            lam[k] = np.maximum(0.0, lam[k] + mu * h[k])
        if np.isfinite(rho_i):
            lam["rho"] = max(0.0, lam["rho"] + mu * h["rho"])
        if use_orb:
            lam["orb"] = np.maximum(0.0, lam["orb"] + mu * h["orb"])
        if len(hist) > 1 and viol > 0.25 * hist[-2]["viol"] and viol > 1e-7:
            mu = min(mu * 4, 1e8)
        if it >= 6 and viol < 1e-7 and abs(hist[-1]["tau"] - hist[-2]["tau"]) < 1e-6 * hist[-1]["tau"]:
            break
    return best, hist


class FC3ModelScaled(FC3Model):
    """Model in the well-conditioned coordinates u (orbit_scaling.py): Phi = Phi0 + G M^-1/2 u, per-orbit relative
    change ||(W u)_o||. The orbit terms of run() are reused with Bmat = W, block_refnorm2 = 1, c0n = 1.
    Unscaled coordinates: y = M^-1/2 u / |c0| (needed only to rebuild the fc3 tensor)."""

    def __init__(self, ev_path, rdir, tag=""):
        super().__init__(ev_path, rdir)
        rdir = Path(rdir)
        pre = f"scaled_{tag}_" if tag else "scaled_"     # tag: per-orbit weighted scaling (orbit_scaling.py --eps)
        self.Gr = np.load(rdir / f"{pre}Gr.npy")
        self.Gi = np.load(rdir / f"{pre}Gi.npy")
        self.n_var = self.Gr.shape[1]
        CO = np.load(rdir / "cluster_orbits.npz")
        self.col_block = CO["col_block"]
        self.n_blocks = int(self.col_block.max()) + 1
        self.block_refnorm2 = np.ones(self.n_blocks)
        self.Bmat = np.load(rdir / f"{pre}W.npy")
        self.c0n_true = self.c0n
        self.c0n = 1.0
