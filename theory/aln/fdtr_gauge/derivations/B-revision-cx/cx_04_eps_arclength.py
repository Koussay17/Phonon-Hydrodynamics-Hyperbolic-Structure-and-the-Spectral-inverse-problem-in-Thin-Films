"""cx_04: pseudo-arclength continuation in eps (eps0 = eps, twin convention) of the three local
three-frequency twins, target and frequencies fixed.  Lands on eps in {1, 0.084, 0.5, 2.5, -0.5, -1}.
Any landed root that is admissible (all entries >= 0) and distinct from the target is a twin of that
convention; it is refined at 60 digits.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct
import cxcont as cc

CASES = json.load(open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cases_tw3.json"))
EPS_LAND = [1.0, 0.084, 0.5, 2.5, -0.5, -1.0]
out = {}
for name, c in CASES.items():
    target = np.array(c["target"], float)
    om = np.array(c["om"], float)
    H = cc.Homotopy(lambda e: cx.Spec("stack", "tied", e, e, xiA=0.2, bA=3.0), lambda e: target, lambda e: om)
    th0 = np.array([float(x) for x in c["twin"]])
    rec = {}
    for lam_end in (3.0, -3.0):
        path = cc.follow(H, th0, 0.0, lam_end, h0=0.005, hmax=0.02, maxsteps=6000)
        lams = [s[0] for s in path["samples"]]
        adm = [cc.admissible(s[1], None) for s in path["samples"]]
        print(name, "to", lam_end, "events", path["events"][:12], "lam range [%.3f, %.3f]" % (min(lams), max(lams)),
              "end", np.round(path["end"], 5), flush=True)
        landed = {}
        for e in EPS_LAND:
            sols = cc.land(H, path, e)
            good = []
            for th in sols:
                d = float(np.max(np.abs(np.array(th) - target)))
                good.append(dict(theta=th, admissible=cc.admissible(th, None), dist_to_target=d))
            if good:
                landed["%g" % e] = good
        rec["path_to_%g" % lam_end] = dict(events=path["events"], lam_min=min(lams), lam_max=max(lams),
                                            n_samples=len(lams), frac_admissible=float(np.mean(adm)),
                                            samples=path["samples"][:: max(1, len(lams) // 60)], landed=landed)
        for e, lst in landed.items():
            for g in lst:
                tag = "TWIN" if (g["admissible"] and g["dist_to_target"] > 1e-6) else ("target" if g["dist_to_target"] <= 1e-6 else "inadmissible")
                print("   eps=%s" % e, tag, np.round(g["theta"], 6))
                if tag == "TWIN":
                    sp = cx.Spec("stack", "tied", float(e), float(e), xiA=0.2, bA=3.0)
                    mp.mp.dps = 60
                    tg = [mp.mpf(x) for x in c["target"]]
                    thm, nF = ct.refine_mp(sp, g["theta"], c["om"], tg, dps=60)
                    Zt = [cx.Z_mp(sp, tg, mp.mpf(w)) for w in c["om"]]
                    sv = ct.jac_cond_mp(sp, thm, c["om"], Zt)
                    mis = abs(cx.Z_mp(sp, thm, mp.mpf("0.9")) / cx.Z_mp(sp, tg, mp.mpf("0.9")) - 1)
                    g["theta_60"] = [mp.nstr(x, 40) for x in thm]
                    g["residual_60"] = float(nF)
                    g["smin_over_smax"] = float(sv[-1] / sv[0])
                    g["mismatch_omega_0p9"] = float(mis)
                    print("      refined res %.2e  smin/smax %.2e  mismatch(0.9) %.2e" % (float(nF), float(sv[-1] / sv[0]), float(mis)))
    out[name] = rec

with open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_04_eps_arclength.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
