"""Parallel job runner for the review attacks (box priors, forbidden-event exclusion, multi-T, extras).

usage: python cx_jobs.py <jobs.json> <nproc>
Each job (dict): id, col, sign, lo, hi, forbid_zero (bool), temps (list, first = 300 K objective block),
kcols_T (columns of K fixed at every temperature), extras (list), start, maxit, time_limit, seed.
start: "ref" | "rand:<sigma>" | "witness:<name>:<t_mix>" | "file:<path.npy>"
Outputs: review/cx_runs/<id>.json and <id>_u.npy (orbit log-rates relative to the base rates).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = os.environ.get("CX_THREADS", "2")
sys.dont_write_bytecode = True
HERE = Path(r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

OUT = HERE / "cx_runs"
_SETUP = {}


def setup():
    if "M" in _SETUP:
        return _SETUP
    from cx_common import AlN
    M = AlN()
    o = M.orbits()
    gref = M.g_ref()
    _SETUP.update(M=M, o=o, gref=gref, forb=gref < 1e-12 * gref.max(), w0=gref / M.Bf)
    return _SETUP


def build(job):
    from cx_common import AlN
    from cx_opt import Block, Problem
    S = setup()
    o, gref, forb, w0 = S["o"], S["gref"], S["forb"], S["w0"]
    temps = job.get("temps", [300.0])
    kcols = tuple(job.get("kcols_T", [0, 2]))
    blocks = []
    for T in temps:
        MT = S["M"] if T == 300.0 else AlN(T=T)
        g_T = w0 * MT.Bf                      # reference rates at T (T-independent |Phi|^2 delta)
        base = g_T * (~forb) if job.get("forbid_zero") else g_T
        r_T = MT.diag(g_T)                    # data: reference lifetimes at T
        moT = MT.moments(g_T, cols=kcols)
        K0 = {c: moT[c]["K"] for c in kcols}
        fixL = (T == 300.0) or (T in [float(t) for t in job.get("lifetimes_T", temps)])
        blocks.append(Block(MT, base, r_T[o["mode_reps"]], K0, kcols=kcols, fix_lifetimes=fixL))
    P = Problem(blocks, o["ev_orb"], o["n_ev_orb"], o["mode_reps"], job["lo"], job["hi"],
                obj_col=job["col"], sign=job["sign"], extras=[tuple(e) for e in job.get("extras", [])])
    if P.extras:
        P.set_extra_targets(np.zeros(P.p))     # extras fixed at their reference values
    P.rcond = job.get("rcond")                 # truncated pseudo-inverse for nearly dependent constraints
    P.tol = job.get("tol", 1e-11)              # feasibility tolerance (relative residuals)
    return P


def start_u(job, P):
    S = setup()
    st = job["start"]
    rng = np.random.default_rng(job.get("seed", 0))
    if st == "ref":
        return np.zeros(P.p)
    if st.startswith("rand:"):
        sig = float(st.split(":")[1])
        return np.clip(rng.normal(0, sig, P.p), P.lo, P.hi)
    if st.startswith("witness:"):
        from cx_common import their_ev_orb, load_witness
        _, name, tmix = st.split(":")
        tmix = float(tmix)
        ev_orb_theirs, _ = their_ev_orb()
        g = load_witness(name, ev_orb_theirs)
        gref = S["gref"]
        gm = (1 - tmix) * g + tmix * gref
        o = S["o"]
        ratio = np.bincount(o["ev_orb"], weights=gm, minlength=o["n_ev_orb"]) / \
            np.bincount(o["ev_orb"], weights=gref, minlength=o["n_ev_orb"])
        return np.clip(np.log(np.maximum(ratio, 1e-300)), P.lo, P.hi)
    if st.startswith("file:"):
        return np.clip(np.load(st[5:]), P.lo, P.hi)
    if st.startswith("filemix:"):
        # rates (1-t) * g_file + t * g_base, then clipped into the box (u relative to the base rates)
        _, path, tmix = st.rsplit(":", 2) if st.count(":") > 2 else (None, st.split(":")[1], st.split(":")[2])
        path = st[len("filemix:"):st.rfind(":")]
        tmix = float(st[st.rfind(":") + 1:])
        uf = np.load(path)
        return np.clip(np.log((1 - tmix) * np.exp(np.clip(uf, -700, 700)) + tmix), P.lo, P.hi)
    raise ValueError(st)


def run_job(job):
    OUT.mkdir(exist_ok=True)
    fj = OUT / f"{job['id']}.json"
    if fj.exists():
        return str(fj)
    t0 = time.time()
    try:
        P = build(job)
        u0 = start_u(job, P)
        hom = None
        if job.get("homotopy"):
            u0, hom = P.homotopy(u0, time_limit=0.4 * job.get("time_limit", 3600))
            if u0 is None:
                raise RuntimeError("homotopy failed at the reference")
        res = P.run(u0, maxit=job.get("maxit", 400), verbose=False, time_limit=job.get("time_limit"),
                    max_du=job.get("max_du", 3.0))
        if res is not None:
            res["homotopy_fraction"] = hom
    except Exception as exc:  # noqa: BLE001
        rec = {"job": job, "error": repr(exc), "seconds": time.time() - t0}
        fj.write_text(json.dumps(rec, indent=1))
        return str(fj)
    rec = {"job": job, "seconds": time.time() - t0}
    if res is None:
        rec["result"] = None
    else:
        u = res["u"]
        np.save(OUT / f"{job['id']}_u.npy", u)
        I = P.evaluate(u, need_grad=False)
        # report tau for all K columns at 300 K and constraint residual
        S = setup()
        g300 = P.rates(u, 0)
        mo = P.blocks[0].M.moments(g300, cols=(0, 2))
        rec["result"] = {"tau_ps": res["tau_ps"], "c_max": res["c_max"], "iters": res["iters"], "nfev": res["nfev"],
                         "tau_x_ps": mo[0]["tau_ps"], "tau_z_ps": mo[2]["tau_ps"],
                         "frac_at_lo": float(np.mean(u <= P.lo + 1e-9)), "frac_at_hi": float(np.mean(u >= P.hi - 1e-9)),
                         "homotopy_fraction": res.get("homotopy_fraction"),
                         "hist_tail": [float(np.exp(h * job["sign"]) / (4 * np.pi)) for h in res["hist"][-5:]]}
    fj.write_text(json.dumps(rec, indent=1))
    return str(fj)


def main():
    jobs = json.loads(Path(sys.argv[1]).read_text())
    nproc = int(sys.argv[2])
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, p in enumerate(pool.imap_unordered(run_job, jobs)):
            print(f"[{time.time()-t0:8.0f}s] {i+1}/{len(jobs)} {p}", flush=True)


if __name__ == "__main__":
    main()
