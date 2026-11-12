"""cx_09: Krawczyk certification (python-flint balls, 400 bits, box radius 1e-40) of twin roots.

Input: a JSON list of cases {conv: [eps0, eps], transducer: [xiA, bA] (decimal strings), target: [...],
omega: [...] (decimal strings), twin: [...] (>= 30 digits)}.  Each twin is first refined at 110 digits
by Newton (mpmath, exact decimal data), then certified; a control with the centre moved by 1e-30 must fail.
Usage: python cx_09_certify.py <cases.json> <out.json>
"""
import json
import sys
sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle
throttle.wait_slot(2)          # light wait before any numerical import
import flint
flint.ctx.prec = 400
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct
import cxkraw as kr

cases = json.load(open(sys.argv[1]))
res = []
for c in cases:
    e0, e = c["conv"]
    xiA, bA = c["transducer"]
    om = c["omega"]
    mp.mp.dps = 110
    spm = cx.Spec("stack", "tied", mp.mpf(e0), mp.mpf(e), xiA=mp.mpf(xiA), bA=mp.mpf(bA))
    tg = [mp.mpf(x) for x in c["target"]]
    th, nF = ct.refine_mp(spm, [mp.mpf(x) for x in c["twin"]], om, tg, dps=110)
    mid = [mp.nstr(x, 100) for x in th]
    flint.ctx.prec = 400
    spf = cx.Spec("stack", "tied", flint.arb(e0), flint.arb(e), xiA=flint.arb(xiA), bA=flint.arb(bA))
    r = kr.krawczyk_target(spf, c["target"], om, mid, rad="1e-40", prec=400)
    mid2 = [mp.nstr(th[0] + mp.mpf("1e-30"), 100)] + mid[1:]
    ctrl = kr.krawczyk_target(spf, c["target"], om, mid2, rad="1e-40", prec=400)["certified"]
    dist = max(abs(mp.log(th[k] / tg[k])) for k in range(6))
    rec = dict(name=c.get("name", ""), conv=c["conv"], transducer=c["transducer"], target=c["target"], omega=om,
               twin_40=[mp.nstr(x, 40) for x in th], newton_residual_110=float(nF), certified=r["certified"],
               K_dev_over_radius=r["max_K_dev_over_radius"], control_fails=not ctrl, max_log_distance_to_target=float(dist),
               admissible=bool(min(th) > 0))
    res.append(rec)
    print(rec["name"], "certified" if r["certified"] else "NOT certified", "K-dev/r %.1e" % r["max_K_dev_over_radius"],
          "control fails" if not ctrl else "CONTROL PASSED (bad)", "dist %.3f" % float(dist), flush=True)
json.dump(res, open(sys.argv[2], "w"), indent=1)
