"""E20: range of tau_mem with the static data imposed on TWO samplings of the same smeared operator (5x5x3 and
7x7x5 event sets, sigma 0.1 THz, 4 sigma, 300 K), for one fc3 shared by both.

A witness for a statement about AlN must keep the data and the memory effect when the Brillouin-zone sampling is
refined (review, section 4.2: the 5x5x3 witnesses do not). Here the unknown u (scaled coordinates of a per-orbit
class, orbit_scaling.py) drives both linearised models, Phi_mesh = Phi0_mesh + G_mesh M^-1/2 u, and every mesh in
the constraint list must keep all mode-orbit lifetimes within dr and kappa_xx, kappa_zz within dK (block-averaged
convention, as optimize_fc3.py). The objective is log tau on one mesh (default: the finer one).
Same augmented-Lagrangian / L-BFGS-B scheme, margins and feasibility test as optimize_fc3.run.

usage: python -B joint_mesh.py <tag|orb> <eps> <col 0|2> <max|min> <objective mesh m553|m775> <constraint meshes,
       comma separated> <out prefix> [maxiter] [maxouter]
  tag orb: uniform per-cluster class, bound eps; tag cal1, cal1nz, ...: calibrated class (bound 1 in its coordinates)
"""
from __future__ import annotations

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
from scipy.optimize import minimize  # noqa: E402

from fc3class import FC3Class, FOURPI  # noqa: E402
from optimize_fc3 import FC3Model  # noqa: E402

RES = HERE.parent / "results"
MESH = {"m553": (RES / "m553" / "events_m553_s0.1.npz", RES / "m553"),
        "m775": (RES / "m775x" / "events_m775_s0.1.npz", RES / "m775x")}


class MeshModel(FC3Model):
    """Linearised model of one mesh in the scaled coordinates of a class (no unscaled G is loaded)."""

    def __init__(self, ev_path, rdir, tag, scaledir):
        self.M = M = FC3Class(ev_path, Path(rdir) / "pp_r01_REF.npz")
        lin = "linearization.npz" if (Path(rdir) / "linearization.npz").exists() else "lin_meta.npz"
        Z = np.load(Path(rdir) / lin)
        pre = "scaled_" if tag == "orb" else f"scaled_{tag}_"
        self.Gr = np.load(Path(rdir) / f"{pre}Gr.npy")
        self.Gi = np.load(Path(rdir) / f"{pre}Gi.npy")
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
        base = self.evaluate(np.zeros(self.n_var))
        self.r0, self.K0, self.base = base["r"], {0: base["K"][0], 2: base["K"][2]}, base


def run_joint(models, obj, cmeshes, Wmat, col_block, eps, col, sense, dr=0.1, dK=0.01, maxouter=15, maxiter=300, logfile=None):
    n = models[obj].n_var
    nb = int(col_block.max()) + 1
    u = np.zeros(n)
    dr_i, dK_i, oe_i = dr * (1 - 2e-3), dK * (1 - 2e-2), eps * (1 - 1e-3)
    names = list(models)
    lam = {m: {"rp": np.zeros(len(models[m].r0)), "rm": np.zeros(len(models[m].r0)), "Kp": np.zeros(2),
               "Km": np.zeros(2)} for m in names}
    lam_o = np.zeros(nb)
    mu, sgn = 50.0, (-1.0 if sense == "max" else 1.0)
    hist, best, t0 = [], None, time.perf_counter()

    def rel2_of(u):
        Wu = Wmat @ u
        return np.bincount(col_block, weights=Wu * Wu, minlength=nb), Wu

    def cons(m, ev):
        md = models[m]
        s = ev["r"] / md.r0 - 1.0
        kk = np.array([ev["K"][0] / md.K0[0] - 1.0, ev["K"][2] / md.K0[2] - 1.0])
        return {"rp": s - dr_i, "rm": -s - dr_i, "Kp": kk - dK_i, "Km": -kk - dK_i}, s, kk

    def fun(u):
        f, gu = 0.0, np.zeros(n)
        for m in names:
            md = models[m]
            ev = md.evaluate(u, grad=True)
            Dorb = np.zeros(md.n_orb)
            if m == obj:
                f += sgn * np.log(ev["tau"][col])
                Dorb += sgn * (ev["dN"][col] / ev["N"][col] - ev["dK"][col] / ev["K"][col])
            if m in cmeshes:
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
        viol, feas = 0.0, True
        for m in names:
            ev = models[m].evaluate(u)
            h, s, kk = cons(m, ev)
            rec["mesh"][m] = {"tau_x": ev["tau"][0] / FOURPI, "tau_z": ev["tau"][2] / FOURPI,
                              "tau_ratio": [ev["tau"][0] / models[m].base["tau"][0], ev["tau"][2] / models[m].base["tau"][2]],
                              "kappa_ratio": [1 + kk[0], 1 + kk[1]], "max_dr": float(np.max(np.abs(s)))}
            if m in cmeshes:
                viol = max(viol, max(float(np.max(h[k])) for k in h))
                feas &= bool(np.max(np.abs(s)) <= dr and np.max(np.abs(kk)) <= dK)
                for k in h:
                    lam[m][k] = np.maximum(0.0, lam[m][k] + mu * h[k])
            if m == obj:
                rec["tau"] = ev["tau"][col] / FOURPI
        rel2 = rel2_of(u)[0]
        ho = rel2 / oe_i ** 2 - 1.0
        viol = max(viol, float(np.max(ho)))
        rec["orbit_rel_over_eps_max"] = float(np.sqrt(np.max(rel2)) / eps)
        feas &= rec["orbit_rel_over_eps_max"] <= 1.0
        lam_o[:] = np.maximum(0.0, lam_o + mu * ho)
        rec.update({"viol": viol, "feasible": bool(feas), "t": time.perf_counter() - t0})
        hist.append(rec)
        print(json.dumps(rec), flush=True)
        if logfile:
            with open(logfile, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
        if feas and (best is None or sgn * np.log(rec["tau"]) < sgn * np.log(best[0])):
            best = (rec["tau"], u.copy(), rec)
        if len(hist) > 1 and viol > 0.25 * hist[-2]["viol"] and viol > 1e-7:
            mu = min(mu * 4, 1e8)
        if it >= 6 and viol < 1e-7 and abs(hist[-1]["tau"] - hist[-2]["tau"]) < 1e-6 * hist[-1]["tau"]:
            break
    return best, hist


if __name__ == "__main__":
    tag, eps, col, sense, obj = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
    cmeshes = sys.argv[6].split(",")
    out = Path(sys.argv[7])
    maxiter = int(sys.argv[8]) if len(sys.argv) > 8 else 300
    maxouter = int(sys.argv[9]) if len(sys.argv) > 9 else 15
    out.parent.mkdir(parents=True, exist_ok=True)
    pre = "scaled_" if tag == "orb" else f"scaled_{tag}_"
    used = sorted(set(cmeshes) | {obj})
    models = {m: MeshModel(*MESH[m], tag, RES / "m553") for m in used}
    Wmat = np.load(RES / "m553" / f"{pre}W.npy")
    col_block = np.load(RES / "m553" / "cluster_orbits.npz")["col_block"]
    best, hist = run_joint(models, obj, cmeshes, Wmat, col_block, eps, col, sense, maxiter=maxiter, maxouter=maxouter,
                           logfile=str(out) + ".log")
    summ = {"tag": tag, "eps": eps, "col": col, "sense": sense, "objective_mesh": obj, "constraint_meshes": cmeshes,
            "maxiter": maxiter, "maxouter": maxouter,
            "base_tau_ps": {m: [models[m].base["tau"][0] / FOURPI, models[m].base["tau"][2] / FOURPI] for m in used},
            "best": best[2] if best else None, "history": hist}
    if best:
        np.save(str(out) + "_u.npy", best[1])
    Path(str(out) + ".json").write_text(json.dumps(summ, indent=1))
