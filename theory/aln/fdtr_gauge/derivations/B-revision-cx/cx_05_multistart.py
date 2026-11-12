"""cx_05: batched multistart for three-frequency twins Z(omega_j; th') = Z(omega_j; th*), j = 1..3.

Log variables (positive parameters); starts uniform in a log-box of half-width W around the target
(plus a wide absolute box for the times); batched Levenberg-Marquardt; roots with relative residual
< 1e-11 and log-distance > 1e-5 from the target are refined at 60 digits (Newton) and kept if the
refined residual is < 1e-50 and the Jacobian is regular.

Usage: python cx_05_multistart.py <case-set> [nstarts]
"""
import json
import sys
import time
sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle
throttle.wait_slot(2)          # light wait before any numerical import
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
case_set = sys.argv[1] if len(sys.argv) > 1 else "branch"
NST = int(sys.argv[2]) if len(sys.argv) > 2 else 20000

th_aln10, (xiA_aln, bA_aln), t0 = cx.aln_stack_nd(tau_l=10e-12)
th_aln0, _, _ = cx.aln_stack_nd(tau_l=0.0)
th_aln25, _, _ = cx.aln_stack_nd(tau_l=25e-12)
f_MHz = np.array([1.0, 20.0, 150.0])
om_band = 2 * np.pi * f_MHz * 1e6 * t0
f_MHz2 = np.array([0.1, 5.0, 60.0])
om_band2 = 2 * np.pi * f_MHz2 * 1e6 * t0

TARGETS = {
    "generic": ([0.8, 2.0, 0.5, 0.1, 0.3, 0.2], (0.2, 3.0)),
    "b=bs": ([0.8, 1.0, 0.5, 0.1, 0.3, 0.2], (0.2, 3.0)),
    "AlN-like(branch)": ([1.0, 2.07, 0.03, 0.012, 0.05, 0.1], (0.2, 3.0)),
    "AlN note16 tl=10ps": (list(th_aln10), (xiA_aln, bA_aln)),
    "AlN note16 tl=0": (list(th_aln0), (xiA_aln, bA_aln)),
    "AlN note16 tl=25ps": (list(th_aln25), (xiA_aln, bA_aln)),
}
TRIPLES = {"T1": [0.5, 4.0, 1.7], "T2": [0.2, 1.1, 6.0], "band1": list(om_band), "band2": list(om_band2)}
CONVS = {"entropic": (1.0, 1.0), "local": (0.0, 0.0), "kin0.084": (0.084, 0.084), "kin-0.5": (-0.5, -0.5)}

if case_set == "branch":
    runs = [(t, tr, cv) for t in ["generic", "b=bs", "AlN-like(branch)"] for tr in ["T1", "T2"] for cv in ["entropic", "local"]]
elif case_set == "aln":
    runs = [(t, tr, cv) for t in ["AlN note16 tl=10ps", "AlN note16 tl=0", "AlN note16 tl=25ps"]
            for tr in ["band1", "band2", "T1"] for cv in ["entropic", "local"]]
elif case_set == "aln2":
    runs = [(t, tr, cv) for t in ["AlN note16 tl=0", "AlN note16 tl=25ps"] for tr in ["band1", "band2", "T1"]
            for cv in ["entropic", "local"]]
elif case_set == "kin":
    runs = [(t, tr, cv) for t in ["generic", "b=bs"] for tr in ["T1", "T2"] for cv in ["kin0.084", "kin-0.5"]]
else:
    raise SystemExit("unknown case set")

rng = np.random.default_rng(7)
results = {}
for (tn, trn, cvn) in runs:
    tgt, (xiA, bA) = TARGETS[tn]
    om = np.array(TRIPLES[trn], float)
    e0, e = CONVS[cvn]
    sp = cx.Spec("stack", "tied", e0, e, xiA=xiA, bA=bA)
    tgt = np.array(tgt, float)
    Zt = cx.Z_np(sp, tgt, om)
    t_start = time.time()
    lt = np.log(np.maximum(tgt, 1e-3))          # box centre (a zero target entry, e.g. tau_l = 0, is centred at 1e-3)
    roots = []
    nb = 2000
    for batch in range(NST // nb):
        Y0 = lt[None, :] + rng.uniform(-3.0, 3.0, (nb, 6))
        # half of the starts: times drawn from an absolute wide box
        k = nb // 2
        Y0[:k, 2] = rng.uniform(np.log(1e-4), np.log(10.0), k)
        Y0[:k, 3] = rng.uniform(np.log(1e-4), np.log(10.0), k)
        Y, res = ct.batch_lm(sp, Y0, om, Zt, iters=250)
        good = np.isfinite(res) & (res < 1e-11)
        for y in Y[good]:
            d = np.max(np.abs(np.exp(y) - tgt) / np.maximum(tgt, 1e-3))
            if d > 1e-5 and not any(np.max(np.abs(y - r)) < 1e-6 for r in roots):
                roots.append(y)
    # refine at 60 digits
    mp.mp.dps = 60
    tg_mp = [mp.mpf(repr(float(x))) for x in tgt]
    refined = []
    for y in roots:
        thm, nF = ct.refine_mp(sp, np.exp(y), list(om), tg_mp, dps=60)
        if nF > mp.mpf("1e-45") or min(thm) <= 0:
            continue
        Ztm = [cx.Z_mp(sp, tg_mp, mp.mpf(repr(float(w)))) for w in om]
        sv = ct.jac_cond_mp(sp, thm, list(om), Ztm)
        mis = abs(cx.Z_mp(sp, thm, mp.mpf(repr(float(om[0] * 1.37)))) / cx.Z_mp(sp, tg_mp, mp.mpf(repr(float(om[0] * 1.37)))) - 1)
        dist = max(abs(mp.log(thm[i] / tg_mp[i])) if tg_mp[i] > 0 else abs(thm[i]) for i in range(6))
        if dist < mp.mpf("1e-8"):
            continue
        if any(max(abs(mp.log(thm[i] / r["_mp"][i])) for i in range(6)) < 1e-20 for r in refined):
            continue
        refined.append(dict(_mp=thm, theta=[mp.nstr(x, 25) for x in thm], residual_60=float(nF),
                            smin_over_smax=float(sv[-1] / sv[0]), max_log_dist=float(dist),
                            mismatch_off_data=float(mis)))
    for r in refined:
        r.pop("_mp")
    key = f"{tn} | {trn} {list(np.round(om, 6))} | {cvn}"
    results[key] = dict(target=list(tgt), transducer=(xiA, bA), omega=list(om), conv=(e0, e), nstarts=NST,
                        n_double_roots=len(roots), twins=refined, seconds=time.time() - t_start)
    print(f"{key}: starts {NST}, distinct double roots {len(roots)}, certified-quality twins {len(refined)}, "
          f"{time.time() - t_start:.0f}s", flush=True)
    for r in refined:
        print("    twin", r["theta"][:6], "res %.1e smin/smax %.1e dist %.3f mis %.1e" %
              (r["residual_60"], r["smin_over_smax"], r["max_log_dist"], r["mismatch_off_data"]), flush=True)
    with open(DIR + f"cx_05_multistart_{case_set}.json", "w") as f:
        json.dump(results, f, indent=1)
