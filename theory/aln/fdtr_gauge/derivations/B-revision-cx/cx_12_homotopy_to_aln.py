"""cx_12: homotopy of three-frequency twins towards the Au/AlN/sapphire note-16 target.

Start: twin pairs (th0, th0') found by cx_06 (transducer (0.2, 3), omega = (0.5, 4, 1.7)).
Path lam in [0, 1]: log-linear interpolation of the target th*(lam), of the transducer (xiA, bA)(lam) and of the
three frequencies omega(lam), ending at the note-16 target (tau_l = 10 ps) with its Au transducer and three
frequencies in the FDTR band.  The twin th'(lam) is followed by pseudo-arclength continuation; reaching lam = 1
with an admissible, distinct root gives a three-frequency twin of the AlN target.
Usage: python cx_12_homotopy_to_aln.py <conv: entropic|local> <band: band1|band2> [npairs]
"""
import sys as _s
_s.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle as _t
_t.wait_slot(2)          # light wait before any numerical import
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct
import cxcont as cc

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
conv = sys.argv[1]
band = sys.argv[2]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 20
eps = 1.0 if conv == "entropic" else 0.0
pairs = json.load(open(DIR + f"cx_06_joint_{conv}_N3.json"))["pairs"]
th_aln, (xiA1, bA1), t0 = cx.aln_stack_nd(tau_l=10e-12)
fMHz = {"band1": [1.0, 20.0, 150.0], "band2": [0.1, 5.0, 60.0]}[band]
om1 = 2 * np.pi * np.array(fMHz) * 1e6 * t0
om0 = np.array([0.5, 4.0, 1.7])
# order the end frequencies to match the start ordering (by magnitude)
order0 = np.argsort(om0)
om1s = np.empty(3); om1s[order0] = np.sort(om1)
out = []
for ip, P in enumerate(pairs[:NP]):
    th0 = np.array([float(x) for x in P["theta"]])
    tw0 = np.array([float(x) for x in P["theta_twin"]])
    lt0, lt1 = np.log(th0), np.log(np.array(th_aln))
    H = cc.Homotopy(lambda l: cx.Spec("stack", "tied", eps, eps, xiA=float(np.exp((1 - l) * np.log(0.2) + l * np.log(xiA1))),
                                      bA=float(np.exp((1 - l) * np.log(3.0) + l * np.log(bA1)))),
                    lambda l: np.exp((1 - l) * lt0 + l * lt1),
                    lambda l: np.exp((1 - l) * np.log(om0) + l * np.log(om1s)))
    path = cc.follow(H, tw0, 0.0, 1.0, h0=0.002, hmax=0.01, maxsteps=4000)
    lams = [s[0] for s in path["samples"]]
    landed = cc.land(H, path, 1.0)
    rec = dict(pair=ip, lam_max=max(lams), lam_min=min(lams), events=path["events"][:10], landed=[])
    for th in landed:
        d = float(np.max(np.abs(np.log(np.abs(np.array(th)) / np.array(th_aln)))))
        adm = cc.admissible(th, None)
        rec["landed"].append(dict(theta=th, admissible=adm, max_log_dist=d))
    out.append(rec)
    print(f"pair {ip}: lam reached [{min(lams):.3f}, {max(lams):.3f}], events {path['events'][:4]}, landed at lam=1: "
          f"{[(np.round(l['theta'], 5).tolist(), l['admissible'], round(l['max_log_dist'], 3)) for l in rec['landed']]}",
          flush=True)
    json.dump(dict(conv=conv, band=band, omega_end=list(om1s), target_end=list(th_aln), transducer_end=(xiA1, bA1),
                   results=out), open(DIR + f"cx_12_homotopy_{conv}_{band}.json", "w"), indent=1)
