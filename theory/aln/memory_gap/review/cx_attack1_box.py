"""Attack 1: how much of the AlN memory gap survives a factor-F prior on the event rates?

Prior P_F: g_ref/F <= g <= F g_ref per event orbit (6mm x time reversal symmetric rates), plus the
campaign's data: all 897 lifetimes (117 orbit rows) and K_xx, K_zz fixed.
(a) RIGOROUS bound (Loewner order): g >= g_ref/F  =>  C(g) >= C(g_ref)/F on H  =>
    tau(g) = |x|^2/(x^T C x) <= 1/lambda_min(C(g)|H) <= F / lambda_min(C(g_ref)|H).
    Only the lower side of the box is used; lifetimes and K are not needed.
(b) inner values: maximise / minimise tau_x, tau_z over P_F with the data fixed (gradient projection).
usage: python cx_attack1_box.py F [nstarts] [maxit]
Output: review/cx_attack1_box_F{F}.json
"""
import json
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI  # noqa: E402
from cx_opt import Block, Problem  # noqa: E402

F = float(sys.argv[1])
nstarts = int(sys.argv[2]) if len(sys.argv) > 2 else 3
maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 400
tasks = sys.argv[4].split(",") if len(sys.argv) > 4 else ["0max", "2max", "0min", "2min"]

M = AlN()
o = M.orbits()
gref = M.g_ref()
r = M.diag(gref)
mo = M.moments(gref, cols=(0, 2))
K0 = {c: mo[c]["K"] for c in (0, 2)}
lamH = M.moments_eig(gref, cols=(0, 2))["lam_min_H"]
bound_ps = F / lamH / FOURPI
out = {"F": F, "lam_min_H_ref": lamH, "rigorous_upper_bound_tau_ps": bound_ps,
       "tau_ref_ps": {c: mo[c]["tau_ps"] for c in (0, 2)},
       "tau_CS_ps": {c: K0[c] / float(M.b[:, c] @ M.b[:, c]) / FOURPI for c in (0, 2)}, "runs": {}}
print(json.dumps({k: v for k, v in out.items() if k != "runs"}), flush=True)

L = np.log(F)
rng = np.random.default_rng(int(1000 * F) + 7)
for task in tasks:
    col = int(task[0])
    sign = +1 if task.endswith("max") else -1
    best = None
    runs = []
    for s in range(nstarts):
        blk = Block(M, gref, r[o["mode_reps"]], K0, kcols=(0, 2))
        P = Problem([blk], o["ev_orb"], o["n_ev_orb"], o["mode_reps"], -L, L, obj_col=col, sign=sign)
        if s == 0:
            u0 = np.zeros(P.p)
        else:
            u0 = np.clip(rng.normal(0.0, 0.5 * L, P.p), -L, L)
        t0 = time.time()
        res = P.run(u0, maxit=maxit, verbose=False, time_limit=1500)
        if res is None:
            runs.append({"start": s, "result": None})
            print(task, s, "restore failed", flush=True)
            continue
        u = res["u"]
        frac_lo = float(np.mean(u <= -L + 1e-9))
        frac_hi = float(np.mean(u >= L - 1e-9))
        rec = {"start": s, "tau_ps": res["tau_ps"], "c_max": res["c_max"], "iters": res["iters"],
               "nfev": res["nfev"], "frac_orbits_at_lower": frac_lo, "frac_orbits_at_upper": frac_hi,
               "seconds": time.time() - t0}
        runs.append(rec)
        print(task, s, json.dumps(rec), flush=True)
        if best is None or sign * res["tau_ps"] > sign * best["tau_ps"]:
            best = dict(rec)
            np.save(CAMP / "review" / f"cx_attack1_box_F{int(F)}_{task}_u.npy", u)
    out["runs"][task] = {"best": best, "all": runs}
    (CAMP / "review" / f"cx_attack1_box_F{int(F)}.json").write_text(json.dumps(out, indent=1))
print("done")
