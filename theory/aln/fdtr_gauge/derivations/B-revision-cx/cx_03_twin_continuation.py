"""cx_03: (1) reproduce the three certified local three-frequency twins (TW3) with the independent model;
(2) continue each twin in the face coefficient eps (same at both film faces, as in the twin) from 0 to 1
    (entropic) and beyond / below, for fixed target and frequencies.  A twin that survives the continuation
    is a three-frequency twin of the entropic (or kinetic) convention.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct

CASES = {
    "N3a generic": dict(target=[0.8, 2.0, 0.5, 0.1, 0.3, 0.2], om=[0.5, 4.0, 1.7],
                        twin=["0.9858503942665071106433822115673965998921", "1.819856744816197917019237601594725930427",
                              "0.3801740956228803828235321290266019457612", "0.05427629762193134489886040884168456454904",
                              "0.274726741561021621413692925038094213953", "0.2915564967990104512584506216640756828577"]),
    "N3b generic": dict(target=[0.8, 2.0, 0.5, 0.1, 0.3, 0.2], om=[0.2, 1.1, 6.0],
                        twin=["0.7196743272911544303131335376428813807069", "2.112433652816248966617088964429045679345",
                              "0.5739888670843872666640827662008432719659", "0.1242013061577606432303457746016427968476",
                              "0.312705172762064213989233307685753353126", "0.165359300192196034630728615239604413143"]),
    "N3b b=bs": dict(target=[0.8, 1.0, 0.5, 0.1, 0.3, 0.2], om=[0.2, 1.1, 6.0],
                     twin=["0.3875508326944282554972740782225254071126", "1.545673439828704871335581129752434850496",
                           "1.358582444181814981033815993382740900089", "0.2350530113522599065743906315783995799671",
                           "0.4277541214805285925045666966086423250086", "0.415310782921830296055312722130337125253"]),
}

out = {}
mp.mp.dps = 60
for name, c in CASES.items():
    sp = cx.Spec("stack", "tied", 0.0, 0.0, xiA=0.2, bA=3.0)
    tw = [mp.mpf(x) for x in c["twin"]]
    tg = [mp.mpf(x) for x in c["target"]]
    res = mp.sqrt(mp.fsum(abs(cx.Z_mp(sp, tw, mp.mpf(w)) / cx.Z_mp(sp, tg, mp.mpf(w)) - 1) ** 2 for w in c["om"]))
    th_ref, nF = ct.refine_mp(sp, tw, c["om"], tg, dps=60)
    mov = max(abs(th_ref[k] - tw[k]) for k in range(6))
    Zt = [cx.Z_mp(sp, tg, mp.mpf(w)) for w in c["om"]]
    sv = ct.jac_cond_mp(sp, th_ref, c["om"], Zt)
    mism = abs(cx.Z_mp(sp, th_ref, mp.mpf("0.9")) / cx.Z_mp(sp, tg, mp.mpf("0.9")) - 1)
    print(name, "residual of printed twin %.2e" % float(res), "after Newton %.2e" % float(nF),
          "moved %.1e" % float(mov), "smin/smax %.2e" % float(sv[-1] / sv[0]), "mismatch at 0.9: %.3e" % float(mism))
    rec = dict(printed_residual=float(res), refined_residual=float(nF), moved=float(mov),
               jac_smin_over_smax=float(sv[-1] / sv[0]), mismatch_omega_0p9=float(mism))

    # continuation in eps (eps0 = eps), double precision with Newton, then mp refinement at the end points
    paths = {}
    for direction, eps_end in (("up", 3.0), ("down", -3.0)):
        th = np.array([float(x) for x in th_ref])
        eps = 0.0
        h = 0.01 if direction == "up" else -0.01
        hist = []
        status = "reached end"
        while (eps < eps_end - 1e-12) if direction == "up" else (eps > eps_end + 1e-12):
            e_new = eps + h
            if (direction == "up" and e_new > eps_end) or (direction == "down" and e_new < eps_end):
                e_new = eps_end
            spn = cx.Spec("stack", "tied", e_new, e_new, xiA=0.2, bA=3.0)
            Ztn = cx.Z_np(spn, c["target"], np.array(c["om"]))
            thn, nf, it, ok = ct.newton(spn, th, np.array(c["om"]), Ztn, tol=1e-13)
            far = np.max(np.abs(np.log(thn / np.array(c["target"])))) if np.all(thn > 0) else 0
            if ok and np.all(thn > 0) and np.max(np.abs(np.log(thn / th))) < 0.2 and far > 1e-6:
                eps = e_new
                th = thn
                hist.append((eps, [float(x) for x in th]))
                h = h * 1.3 if abs(h) < 0.05 else h
            else:
                h = h / 2
                if abs(h) < 1e-6:
                    status = "stopped at eps=%.6f (fold, exit or merge with target)" % eps
                    break
        paths[direction] = dict(status=status, eps_last=eps, theta_last=[float(x) for x in th],
                                samples=hist[:: max(1, len(hist) // 25)])
        print("   continuation", direction, status, "eps_last=%.4f" % eps, "theta=", np.round(th, 6))
    rec["continuation"] = paths

    # refine and characterise the twin at the conventions of interest
    conv_pts = {}
    for eps_q in [1.0, 0.084, -0.5, 0.5, 2.0]:
        hist = paths["up"]["samples"] if eps_q >= 0 else paths["down"]["samples"]
        allh = paths["up"]["samples"] + paths["down"]["samples"]
        # nearest sample in eps
        cand = [hh for hh in allh if abs(hh[0] - eps_q) < 0.2]
        if not cand:
            continue
        e0, th0 = min(cand, key=lambda hh: abs(hh[0] - eps_q))
        spq = cx.Spec("stack", "tied", eps_q, eps_q, xiA=0.2, bA=3.0)
        Ztq = cx.Z_np(spq, c["target"], np.array(c["om"]))
        thq, nf, it, ok = ct.newton(spq, np.array(th0), np.array(c["om"]), Ztq, tol=1e-13)
        if not ok:
            continue
        thm, nFm = ct.refine_mp(spq, thq, c["om"], [mp.mpf(x) for x in c["target"]], dps=60)
        Ztm = [cx.Z_mp(spq, [mp.mpf(x) for x in c["target"]], mp.mpf(w)) for w in c["om"]]
        svq = ct.jac_cond_mp(spq, thm, c["om"], Ztm)
        misq = abs(cx.Z_mp(spq, thm, mp.mpf("0.9")) / cx.Z_mp(spq, [mp.mpf(x) for x in c["target"]], mp.mpf("0.9")) - 1)
        dist = max(abs(mp.log(thm[k] / c["target"][k])) for k in range(6))
        conv_pts["eps=%g" % eps_q] = dict(theta_twin=[mp.nstr(x, 30) for x in thm], residual_60=float(nFm),
                                          jac_smin_over_smax=float(svq[-1] / svq[0]),
                                          max_log_distance_to_target=float(dist), mismatch_omega_0p9=float(misq))
        print("   eps=%g twin" % eps_q, [mp.nstr(x, 12) for x in thm], "res %.1e" % float(nFm),
              "smin/smax %.2e" % float(svq[-1] / svq[0]), "dist %.3f" % float(dist), "mism0.9 %.2e" % float(misq))
    rec["twins_by_convention"] = conv_pts
    out[name] = rec

with open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_03_twin_continuation.json", "w") as f:
    json.dump(out, f, indent=1)
