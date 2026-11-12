"""B-4, step 7: certification of the three-frequency twins found by b4_06 (Krawczyk test in Arb ball arithmetic).

For each twin theta' of a target theta* at frequencies (w1, w2, w3) listed in b4_06_twin_search.json:
 (1) refine theta' by Newton's method at 80 digits (mpmath, models_mp.response);
 (2) Krawczyk operator on the box X = theta~ + [-r, r]^6 (r = 1e-40) for F(theta) = [Re, Im](Z_tot(w_j; theta) - Z_tot(w_j; theta*)),
     K(X) = theta~ - Y F(theta~) + (I - Y J(X)) (X - theta~),  Y = mid(J(theta~))^-1, J(X) by forward-mode AD over balls;
     K(X) in the interior of X proves a unique root in X;
 (3) report the separation from theta* and the response mismatch at a fourth frequency (certified nonzero).
Output: b4_07_certify_twin.json
"""
import os
import sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp  # noqa: E402
from flint import arb as _arb, arb_mat, ctx  # noqa: E402
from ad_arb import stack_response, arb  # noqa: E402
from models_mp import response  # noqa: E402

ctx.prec = 320
mp.mp.dps = 80
NAMES = ["xi", "b", "tR", "tl", "R1", "R2"]
HERE = os.path.dirname(os.path.abspath(__file__))


def spec_of(conv):
    d = dict(transducer=True, xiA=mp.mpf("0.2"), bA=mp.mpf(3), obs="theta_q")
    if conv == "local":
        d["iface"] = "local"
    else:
        d["iface"] = "eps"
        d["eps"] = mp.mpf(1)
    return d


def refine(th0, th_star, ws, conv):
    sp_ = spec_of(conv)
    star = dict(zip(NAMES, [mp.mpf(str(v)) for v in th_star]))
    tgt = [response(star, mp.mpf(str(w)), sp_) for w in ws]

    def F(*x):
        par = dict(zip(NAMES, x))
        out = []
        for w, z0 in zip(ws, tgt):
            z = response(par, mp.mpf(str(w)), sp_)
            out += [mp.re(z - z0), mp.im(z - z0)]
        return out
    x = mp.findroot(F, [mp.mpf(repr(v)) for v in th0], tol=mp.mpf(10) ** (-70), maxsteps=200)
    return [x[i] for i in range(6)]


def F_and_J(vals, names, ws, conv, target):
    Fv, rows = [], []
    for w, z0 in zip(ws, target):
        Z = stack_response(vals, names, str(w), conv)      # decimal frequency, as in the refinement
        dz = Z.v - z0
        Fv += [dz.real, dz.imag]
        rows.append([x.real for x in Z.d])
        rows.append([x.imag for x in Z.d])
    return Fv, rows


def krawczyk(th_mid, th_star, ws, conv, r):
    star_vals = dict(zip(NAMES, [str(v) for v in th_star]))
    target = [stack_response(star_vals, [], str(w), conv).v for w in ws]
    mids = [_arb(mp.nstr(v, 75)) for v in th_mid]
    Fm, Jm = F_and_J(dict(zip(NAMES, mids)), NAMES, ws, conv, target)
    Jmid = arb_mat([[_arb(x.mid()) for x in row] for row in Jm])
    Y = Jmid.inv()
    Y = arb_mat([[_arb(Y[i, j].mid()) for j in range(6)] for i in range(6)])
    box = [_arb(m.mid(), r) for m in mids]
    _, JX = F_and_J(dict(zip(NAMES, box)), NAMES, ws, conv, target)
    JX = arb_mat(JX)
    I6 = arb_mat([[1 if i == j else 0 for j in range(6)] for i in range(6)])
    dX = arb_mat([[box[i] - mids[i]] for i in range(6)])
    Fcol = arb_mat([[f] for f in Fm])
    K = arb_mat([[m] for m in mids]) - Y * Fcol + (I6 - Y * JX) * dX
    inside = all(((K[i, 0] - mids[i]).abs_upper() < _arb(r)) for i in range(6))
    widths = [float((K[i, 0] - mids[i]).abs_upper()) for i in range(6)]
    return bool(inside), widths, [float(abs(_arb(f.mid()))) for f in Fm]


data = json.load(open(os.path.join(HERE, "b4_06_twin_search.json")))
OUT = {}
for key, rec in data.items():
    conv, tname, wname = key.split(":")
    twins = [r for r in rec["roots"] if not r["is_target"]]
    for k, tw in enumerate(twins):
        ws = rec["frequencies"]
        if len(ws) != 3:
            OUT[f"{key}#{k}"] = {"note": "overdetermined case; refinement only", "theta": tw["theta"]}
            continue
        th = refine(tw["theta"], rec["target"], ws, conv)
        ok, widths, res = krawczyk(th, rec["target"], ws, conv, 1e-40)
        # mismatch at a fourth frequency
        sp_ = spec_of(conv)
        w4 = mp.mpf("0.9")
        z1 = response(dict(zip(NAMES, th)), w4, sp_)
        z0 = response(dict(zip(NAMES, [mp.mpf(str(v)) for v in rec["target"]])), w4, sp_)
        OUT[f"{key}#{k}"] = {"target": rec["target"], "frequencies": ws, "twin_40_digits": [mp.nstr(v, 40) for v in th],
                            "krawczyk_box_radius": 1e-40, "certified_unique_root_in_box": ok,
                            "K_minus_mid_abs_upper": widths, "F_at_mid": res,
                            "relative_mismatch_at_omega_0.9": float(abs(z1 - z0) / abs(z0))}
        print(key, k, "certified:", ok, "twin:", [mp.nstr(v, 12) for v in th],
              "mismatch at 0.9: %.3e" % OUT[f"{key}#{k}"]["relative_mismatch_at_omega_0.9"], flush=True)
with open(os.path.join(HERE, "b4_07_certify_twin.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
