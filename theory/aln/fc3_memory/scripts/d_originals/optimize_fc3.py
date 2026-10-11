"""Task 3: range of tau_mem over symmetric fc3 perturbations that keep the 300 K lifetimes and kappa.

Model (exact for the fc3 class, linearize.py): y = c / |c0| (relative perturbation coordinates in the
7752-dim symmetric + ASR fc3 space), Phi_m(y) = Phi0_m + |c0| G_m y,
g_O(y) = W_O uc mean_{m in block(O)} |Phi_m|^2 on allowed orbits; selection-rule orbits keep reference rates.
Base point y = 0: block-averaged reference (kappa 214.30 / 201.90 W/(m K), tau 34.253 / 27.441 ps).

Problem:  extremise tau_col(y)  subject to
   |r_mu(y)/r_mu(0) - 1| <= dr   for every mode (mode-orbit representatives; rates are symmetric),
   |K_xx(y)/K_xx(0) - 1| <= dK,  |K_zz(y)/K_zz(0) - 1| <= dK,
   |y| <= rho   (relative Frobenius size of the fc3 change; rho = inf allowed).
Method: augmented Lagrangian (PHR, inequality form) with L-BFGS-B inner solves; analytic gradients via
dK/dg = -(a^T x)^2, dN/dg = -2 (a^T x)(a^T u), u = C^+ x, and dg_O/dy = 2 W_O uc |c0| Re(conj(Phi) G)/blk.
Every reported value is that of an explicit fc3 (an inner bound on the attainable range).
usage: python -B optimize_fc3.py <events.npz> <resultdir> <outdir> col sense rho dr dK seed [init_scale]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "2")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from fc3class import FC3Class, FOURPI  # noqa: E402


class FC3Model:
    def __init__(self, ev_path, rdir, subspace=None, embed_dir=None):
        self.M = M = FC3Class(ev_path, Path(rdir) / "pp_r01_REF.npz")
        Z = np.load(Path(rdir) / "linearization.npz")
        self.c0n = float(np.linalg.norm(Z["c0"]))
        self.subspace = subspace
        T = None
        if subspace is not None:
            # restriction to the symmetric + ASR fc3 supported within a cutoff (isometric embedding T)
            tdir = Path(rdir) if (Path(rdir) / f"embed_T_r{float(subspace):.2f}.npy").exists() else Path(embed_dir or rdir)
            T = np.load(tdir / f"embed_T_r{float(subspace):.2f}.npy")
            self.T = T
        if "G" in Z.files:
            Gs = Z["G"] * self.c0n
            if T is not None:
                Gs = Gs @ T
            self.Gr = np.ascontiguousarray(Gs.real)
            self.Gi = np.ascontiguousarray(Gs.imag)
            del Gs
        else:                                   # chunked parts (linearize_big.py)
            parts = int(Z["parts"])
            Gr, Gi = [], []
            for k in range(parts):
                a_ = np.load(Path(rdir) / f"lin_Gr_{k}.npy")
                b_ = np.load(Path(rdir) / f"lin_Gi_{k}.npy")
                if T is not None:
                    a_, b_ = a_ @ T, b_ @ T
                Gr.append(a_ * self.c0n); Gi.append(b_ * self.c0n)
            self.Gr = np.ascontiguousarray(np.vstack(Gr)); self.Gi = np.ascontiguousarray(np.vstack(Gi))
            del Gr, Gi
        self.n_var = self.Gr.shape[1]
        self.Bmat = None
        self.Phi0 = Z["Phi0"]
        self.row_orb = Z["row_orb"]
        blk = Z["blk"]
        self.allowed = Z["allowed"]
        self.uc = float(Z["uc"])
        self.w_row = Z["W_orb"][self.row_orb] * self.uc / blk[self.row_orb]
        self.g_fixed = np.where(self.allowed, 0.0, Z["g_ref_orb"])
        self.n_orb = M.n_orb
        self.ops = M.ops
        geom = M.geom
        self.reps = M.tie["reps"]
        self.Wr = geom.W[self.reps].tocsr()                         # diag C rows of mode-orbit reps
        self.e = geom.ehat
        self.b = geom.b
        # base point
        base = self.evaluate(np.zeros(self.n_var), cols=(0, 2))
        self.r0 = base["r"]
        self.K0 = {0: base["K"][0], 2: base["K"][2]}
        self.base = base

    def load_orbit_structure(self, basisdir):
        """Cluster-orbit blocks of the symmetric basis (cluster_orbits.npz) for per-orbit relative bounds."""
        if self.subspace is not None:
            raise ValueError("orbit-relative bounds implemented for the full space only")
        CO = np.load(Path(basisdir) / "cluster_orbits.npz")
        self.col_block = CO["col_block"]
        self.n_blocks = int(self.col_block.max()) + 1
        self.block_refnorm2 = CO["block_refnorm"] ** 2
        self.Bmat = np.load(Path(basisdir) / "basis_B.npy")

    def orbit_rel2(self, y):
        """Per-orbit relative squared change ||d eta_o||^2 / ||eta0_o||^2 and B y."""
        By = self.Bmat @ y
        rel2 = self.c0n ** 2 * np.bincount(self.col_block, weights=By * By, minlength=self.n_blocks) / self.block_refnorm2
        return rel2, By

    def orbit_rates(self, y):
        Phi = self.Phi0 + (self.Gr @ y) + 1j * (self.Gi @ y)
        g = self.g_fixed.copy()
        g += np.bincount(self.row_orb, weights=self.w_row * (Phi.real ** 2 + Phi.imag ** 2), minlength=self.n_orb)
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
            K = float(b @ x)
            N = float(x @ x)
            out["K"][c], out["N"][c], out["tau"][c] = K, N, N / K
            if grad:
                u = sla.cho_solve(cf, x, check_finite=False)
                p = self.ops.ATx(x)
                q = self.ops.ATx(u)
                out["dK"][c] = np.bincount(self.M.orb, weights=-(p * p), minlength=self.n_orb)
                out["dN"][c] = np.bincount(self.M.orb, weights=-2 * p * q, minlength=self.n_orb)
        out["r"] = self.Wr @ g
        return out

    def chain(self, Dorb, Phi):
        """d/dy of sum_O Dorb[O] g_O(y)."""
        v = self.w_row * Dorb[self.row_orb] * np.conj(Phi)
        return 2.0 * (v.real @ self.Gr - v.imag @ self.Gi)

    def dr_dorb(self, wmode):
        """Orbit derivative of sum_mu wmode[mu] r_mu."""
        ge = self.Wr.T @ wmode                                     # per event
        return np.bincount(self.M.orb, weights=ge, minlength=self.n_orb)


def run(model, col, sense, rho, dr, dK, seed, init_scale=0.0, maxouter=25, log=None, orbit_eps=None, y_init=None,
        logfile=None, maxiter=300):
    rng = np.random.default_rng(seed)
    n = model.n_var
    y = np.zeros(n) if y_init is None else np.array(y_init, dtype=float)
    if init_scale > 0 and y_init is None:
        u = rng.standard_normal(n)
        y = u / np.linalg.norm(u) * init_scale * (rho if np.isfinite(rho) else 1.0)
    dr_i, dK_i = dr * (1 - 2e-3), dK * (1 - 2e-2)
    rho_i = rho * (1 - 1e-4) if np.isfinite(rho) else np.inf
    nmode = len(model.r0)
    lam = {"rp": np.zeros(nmode), "rm": np.zeros(nmode), "Kp": np.zeros(2), "Km": np.zeros(2), "rho": 0.0}
    use_orb = orbit_eps is not None
    if use_orb:
        oe_i = orbit_eps * (1 - 1e-3)
        lam["orb"] = np.zeros(model.n_blocks)
    mu = 50.0
    sgn = -1.0 if sense == "max" else 1.0
    other = 2 if col == 0 else 0
    hist = []
    t0 = time.perf_counter()

    def constraints(ev, y, rel2=None):
        s = ev["r"] / model.r0 - 1.0
        kk = np.array([ev["K"][0] / model.K0[0] - 1.0, ev["K"][2] / model.K0[2] - 1.0])
        h = {"rp": s - dr_i, "rm": -s - dr_i, "Kp": kk - dK_i, "Km": -kk - dK_i,
             "rho": (y @ y - rho_i ** 2) if np.isfinite(rho_i) else -1.0}
        if use_orb:
            if rel2 is None:
                rel2, _ = model.orbit_rel2(y)
            h["orb"] = rel2 / oe_i ** 2 - 1.0
        return h, s, kk

    def fun(y):
        ev = model.evaluate(y, cols=(0, 2), grad=True)
        tau = ev["tau"][col]
        f = sgn * np.log(tau)
        Dorb = sgn * ((ev["dN"][col] / ev["N"][col]) - (ev["dK"][col] / ev["K"][col]))   # d log tau / dg
        if use_orb:
            rel2_, By = model.orbit_rel2(y)
        else:
            rel2_, By = None, None
        h, s, kk = constraints(ev, y, rel2_)
        gy = np.zeros_like(y)
        # AL terms: (1/2mu) [max(0, lam + mu h)^2 - lam^2]
        wp = np.maximum(0.0, lam["rp"] + mu * h["rp"])
        wm = np.maximum(0.0, lam["rm"] + mu * h["rm"])
        f += (np.sum(wp ** 2 - lam["rp"] ** 2) + np.sum(wm ** 2 - lam["rm"] ** 2)) / (2 * mu)
        wmode = (wp - wm) / model.r0                                # d/dr of the AL terms
        Dorb = Dorb + model.dr_dorb(wmode)
        kp = np.maximum(0.0, lam["Kp"] + mu * h["Kp"])
        km = np.maximum(0.0, lam["Km"] + mu * h["Km"])
        f += (np.sum(kp ** 2 - lam["Kp"] ** 2) + np.sum(km ** 2 - lam["Km"] ** 2)) / (2 * mu)
        for i, c in enumerate((0, 2)):
            Dorb = Dorb + (kp[i] - km[i]) * ev["dK"][c] / model.K0[c]
        gy += model.chain(Dorb, ev["Phi"])
        if np.isfinite(rho_i):
            wr = max(0.0, lam["rho"] + mu * h["rho"])
            f += (wr ** 2 - lam["rho"] ** 2) / (2 * mu)
            gy += wr * 2 * y
        if use_orb:
            wo = np.maximum(0.0, lam["orb"] + mu * h["orb"])
            f += np.sum(wo ** 2 - lam["orb"] ** 2) / (2 * mu)
            coef = wo * model.c0n ** 2 / (model.block_refnorm2 * oe_i ** 2)
            gy += 2.0 * (model.Bmat.T @ (coef[model.col_block] * By))
        return f, gy

    best_feasible = None
    for it in range(maxouter):
        res = minimize(fun, y, jac=True, method="L-BFGS-B", options={"maxiter": maxiter, "maxcor": 30, "ftol": 1e-13, "gtol": 1e-9})
        y = res.x
        ev = model.evaluate(y, cols=(0, 2))
        h, s, kk = constraints(ev, y)
        viol = max(float(np.max(h["rp"])), float(np.max(h["rm"])), float(np.max(h["Kp"])), float(np.max(h["Km"])),
                   float(h["rho"]) if np.isfinite(rho_i) else -1.0, float(np.max(h["orb"])) if use_orb else -1.0, 0.0)
        # feasibility w.r.t. the NOMINAL tolerances
        feas = (np.max(np.abs(s)) <= dr) and (np.max(np.abs(kk)) <= dK) and (
            (not np.isfinite(rho)) or (np.linalg.norm(y) <= rho))
        if use_orb:
            orb_max = float(np.sqrt(np.max(model.orbit_rel2(y)[0])))
            feas = feas and (orb_max <= orbit_eps)
        rec = {"outer": it, "tau": ev["tau"][col] / FOURPI, "tau_other": ev["tau"][other] / FOURPI, "viol": viol,
               "max_dr": float(np.max(np.abs(s))), "max_dK": float(np.max(np.abs(kk))), "rho": float(np.linalg.norm(y)),
               "mu": mu, "nit": int(res.nit), "feasible": bool(feas), "t": time.perf_counter() - t0}
        if use_orb:
            rec["orbit_rel_max"] = orb_max
        hist.append(rec)
        if log:
            print(json.dumps(rec), flush=True)
        if logfile:
            with open(logfile, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
        if feas:
            val = ev["tau"][col]
            if best_feasible is None or (sgn * np.log(val) < sgn * np.log(best_feasible[0])):
                best_feasible = (val, y.copy(), rec)
        # multiplier update
        lam["rp"] = np.maximum(0.0, lam["rp"] + mu * h["rp"])
        lam["rm"] = np.maximum(0.0, lam["rm"] + mu * h["rm"])
        lam["Kp"] = np.maximum(0.0, lam["Kp"] + mu * h["Kp"])
        lam["Km"] = np.maximum(0.0, lam["Km"] + mu * h["Km"])
        if np.isfinite(rho_i):
            lam["rho"] = max(0.0, lam["rho"] + mu * h["rho"])
        if use_orb:
            lam["orb"] = np.maximum(0.0, lam["orb"] + mu * h["orb"])
        if len(hist) > 1 and viol > 0.25 * hist[-2]["viol"] and viol > 1e-7:
            mu = min(mu * 4, 1e8)
        if it >= 6 and viol < 1e-7 and abs(hist[-1]["tau"] - hist[-2]["tau"]) < 1e-6 * hist[-1]["tau"]:
            break
    return best_feasible, hist


def main():
    ev_path, rdir, outdir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    col = int(sys.argv[4]); sense = sys.argv[5]
    rho = float(sys.argv[6]); dr = float(sys.argv[7]); dK = float(sys.argv[8]); seed = int(sys.argv[9])
    init_scale = float(sys.argv[10]) if len(sys.argv) > 10 else 0.0
    y_init = np.load(sys.argv[11]) if len(sys.argv) > 11 else None
    outdir.mkdir(parents=True, exist_ok=True)
    name = f"c{col}_{sense}_rho{rho:g}_dr{dr:g}_dK{dK:g}_s{seed}_i{init_scale:g}" + ("_warm" if y_init is not None else "")
    if (outdir / f"{name}.json").exists():
        return
    model = FC3Model(ev_path, rdir)
    best, hist = run(model, col, sense, rho, dr, dK, seed, init_scale, log=True, y_init=y_init)
    out = {"name": name, "col": col, "sense": sense, "rho": rho, "dr": dr, "dK": dK, "seed": seed,
           "init_scale": init_scale, "tau_base_ps": model.base["tau"][col] / FOURPI, "history": hist}
    if best is not None:
        val, y, rec = best
        ev = model.evaluate(y, cols=(0, 2))
        g = ev["g_orb"]
        gb = model.base["g_orb"]
        al = model.allowed
        ratio = g[al] / gb[al]
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": ev["tau"][0] / FOURPI, "tau_z_ps": ev["tau"][2] / FOURPI,
                       "kappa_x_ratio": ev["K"][0] / model.K0[0], "kappa_z_ratio": ev["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(ev["r"] / model.r0 - 1))), "rho": float(np.linalg.norm(y)),
                       "orbit_rate_ratio_min": float(ratio.min()), "orbit_rate_ratio_max": float(ratio.max()),
                       "n_orbits_suppressed_10x": int((ratio < 0.1).sum()), "n_orbits_suppressed_100x": int((ratio < 0.01).sum()),
                       "record": rec}
        np.save(outdir / f"{name}_y.npy", y)
    (outdir / f"{name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out.get("best"), indent=1))


if __name__ == "__main__":
    main()
