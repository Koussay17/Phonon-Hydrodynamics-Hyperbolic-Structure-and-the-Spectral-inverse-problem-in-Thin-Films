"""cx_07: adversarial minimisation of the r-th singular value (r = claimed generic rank of the face) of
the column-normalised Jacobian of log Z (25 frequencies in [1e-3, 1e3]), over every face of the
parameter domain, for many conventions (local, tied eps incl. negative, free c_n, mixed face
coefficients (eps0, eps)), stack and bare film.

A face fixes some coordinates (tau_R = 0, tau_l = 0, R2 = 0, b = b_s, tau_R = tau_l, ...); the other
coordinates are free in log variables in the well-conditioned box xi, b in [0.2, 5], times in [0.01, 5],
resistances in [0.01, 3] (outside it the Jacobian is asymptotically degenerate in double precision);
minima on the box boundary are reported separately.  Batched walkers, Adam on
log10(s_r / s_1) with central finite-difference gradients.  Objective: log10(s_r/s_{r-1}) (isolated drop). A minimum below 1e-8 away from the known
sub-strata would be a degenerate set missing from the tables.

Usage: python cx_07_adversarial_smin.py <group>
"""
import json
import sys
import time
import numpy as np

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
group = sys.argv[1] if len(sys.argv) > 1 else "stack_tied"
OM = np.logspace(-3, 3, 25)
rng = np.random.default_rng(31)


def srank(spec, TH, r):
    Z, J = cx.ZJ_np_batch(spec, TH, OM)
    G = J / Z[:, None, :]
    G = np.concatenate([G.real, G.imag], axis=2).transpose(0, 2, 1)       # (m, 2nf, n)
    nrm = np.linalg.norm(G, axis=1, keepdims=True)
    nrm = np.where(nrm > 0, nrm, 1.0)
    s = np.linalg.svd(G / nrm, compute_uv=False)
    if r >= 2:
        return s[:, r - 1] / s[:, r - 2], s      # isolated drop of the r-th singular value
    return s[:, r - 1] / s[:, 0], s


def make_face(n, rules):
    """rules: list of length n, entries 'f' (free), float constant, or ('=', j) tie to coordinate j."""
    free = [i for i, r in enumerate(rules) if r == "f"]

    def th_of(Y):
        m = Y.shape[0]
        TH = np.zeros((m, n))
        for k, i in enumerate(free):
            TH[:, i] = np.exp(Y[:, k])
        for i, r in enumerate(rules):
            if isinstance(r, (int, float)) and r != "f":
                TH[:, i] = r
        for i, r in enumerate(rules):
            if isinstance(r, tuple):
                TH[:, i] = TH[:, r[1]]
        return TH
    return free, th_of


def box(spec, free):
    lo, hi = [], []
    for i in free:
        nm = spec.names[i]
        if nm in ("xi", "b"):
            lo.append(np.log(0.2)); hi.append(np.log(5.0))
        elif nm in ("tR", "tl", "tn"):
            lo.append(np.log(0.01)); hi.append(np.log(5.0))
        else:
            lo.append(np.log(0.01)); hi.append(np.log(3.0))
    return np.array(lo), np.array(hi)


def search(spec, rules, r, M=64, iters=200, lr=0.06):
    free, th_of = make_face(spec.n, rules)
    k = len(free)
    LO, HI = box(spec, free)
    Y = LO + (HI - LO) * rng.uniform(0, 1, (M, k))
    mom = np.zeros_like(Y); vel = np.zeros_like(Y)
    h = 1e-4
    best = np.full(M, np.inf); bestY = Y.copy()
    with np.errstate(all="ignore"):
        for it in range(iters):
            P = [Y]
            for j in range(k):
                e = np.zeros(k); e[j] = h
                P += [Y + e, Y - e]
            YY = np.concatenate(P, axis=0)
            f, _ = srank(spec, th_of(YY), r)
            f = np.log10(np.maximum(f, 1e-300))
            f0 = f[:M]
            imp = f0 < best
            best[imp] = f0[imp]; bestY[imp] = Y[imp]
            g = np.stack([(f[M * (1 + 2 * j):M * (2 + 2 * j)] - f[M * (2 + 2 * j):M * (3 + 2 * j)]) / (2 * h)
                          for j in range(k)], axis=1)
            g = np.where(np.isfinite(g), g, 0.0)
            mom = 0.9 * mom + 0.1 * g
            vel = 0.999 * vel + 0.001 * g * g
            Y = np.clip(Y - lr * mom / (np.sqrt(vel / (1 - 0.999 ** (it + 1))) + 1e-8) / (1 - 0.9 ** (it + 1)), LO, HI)
    order = np.argsort(best)
    THb = th_of(bestY[order[:6]])
    onb = [bool(np.any(np.abs(bestY[i] - LO) < 1e-3) or np.any(np.abs(bestY[i] - HI) < 1e-3)) for i in order[:6]]
    return [(float(best[i]), th.tolist(), ob) for i, th, ob in zip(order[:6], THb, onb)]


# faces of the stack parameter domain (xi, b, tR, tl, R1, R2); free c_n inserts tn after tl
F = "f"
STACK_FACES = {
    "interior": [F, F, F, F, F, F],
    "R2=0": [F, F, F, F, F, 0.0],
    "tR=0": [F, F, 0.0, F, F, F],
    "tl=0": [F, F, F, 0.0, F, F],
    "Delta": [F, F, F, ("=", 2), F, F],
    "Delta&R2=0": [F, F, F, ("=", 2), F, 0.0],
    "b=bs": [F, 1.0, F, F, F, F],
    "b=bs&R2=0": [F, 1.0, F, F, F, 0.0],
    "tl=0&R2=0": [F, F, F, 0.0, F, 0.0],
    "tR=0&R2=0": [F, F, 0.0, F, F, 0.0],
    "Pi": [F, F, 0.0, 0.0, F, F],
    "Pi&R2=0": [F, F, 0.0, 0.0, F, 0.0],
    "D2": [F, 1.0, F, ("=", 2), F, F],
    "D2&R2=0": [F, 1.0, F, ("=", 2), F, 0.0],
    "R1=0&R2=0": [F, F, F, F, 0.0, 0.0],
}


def claimed_rank_stack(e0, e, face):
    """Claimed generic rank on the face (Tables 0.3 and the mixed-face paragraph of §4)."""
    if e0 == 0 and e == 0:          # local
        table = {"Delta": 5, "Delta&R2=0": 4, "Pi": 5, "Pi&R2=0": 4, "D2": 5, "D2&R2=0": 4}
        return table.get(face, 6)
    if e0 == e:                     # tied eps != 0
        return {"Pi&R2=0": 5}.get(face, 6)
    # mixed
    r = 6
    if face == "Pi&R2=0":
        r = 5
    if face in ("Delta&R2=0", "D2&R2=0") and e == 0:
        r = 5
    if e0 == 0 and e == 0 and face in ("Delta", "D2"):
        r = 5
    return r


results = {}
t0 = time.time()
if group in ("stack_tied", "stack_mixed"):
    if group == "stack_tied":
        convs = [(0.0, 0.0), (1.0, 1.0), (0.084, 0.084), (-0.5, -0.5), (2.5, 2.5), (0.5, 0.5), (-2.0, -2.0)]
    else:
        convs = [(1.0, 0.0), (1.0, -1.0), (0.3, 0.7), (0.0, 1.0), (0.5, 0.0), (2.0, -1.0), (-0.5, 1.0),
                 (1.0, 0.5), (0.5, 1.0), (0.0, -0.5), (1.5, 0.0)]
    for (e0, e) in convs:
        sp = cx.Spec("stack", "tied", e0, e, xiA=0.2, bA=3.0)
        for face, rules in STACK_FACES.items():
            r = claimed_rank_stack(e0, e, face)
            best = search(sp, rules, r)
            key = f"stack eps0={e0:g} eps={e:g} | {face} | claimed rank {r}"
            results[key] = best
            flag = "  <-- LOW" if best[0][0] < -8 else ""
            print(f"{key}: min log10(s_r/s_1) = {best[0][0]:.2f} at {np.round(best[0][1], 4)}{flag}  [{time.time() - t0:.0f}s]",
                  flush=True)
elif group == "stack_free":
    sp = cx.Spec("stack", "free", xiA=0.2, bA=3.0)
    faces = {k: v[:4] + [("f" if True else None)] + v[4:] for k, v in STACK_FACES.items()}
    # tn free in all faces above; add tn faces
    faces.update({
        "tn=0": [F, F, F, F, 0.0, F, F],
        "tn=0&R2=0": [F, F, F, F, 0.0, F, 0.0],
        "Delta&tn=0": [F, F, F, ("=", 2), 0.0, F, F],
        "Delta&tn=0&R2=0": [F, F, F, ("=", 2), 0.0, F, 0.0],
        "tn=tR&R2=0": [F, F, F, F, ("=", 2), F, 0.0],
        "tn=tl&R2=0": [F, F, F, F, ("=", 3), F, 0.0],
        "tn=tl&tR=0&R2=0": [F, F, 0.0, F, ("=", 3), F, 0.0],
        "Pi&tn=0&R2=0": [F, F, 0.0, 0.0, 0.0, F, 0.0],
        "Pi&tn=0&R2=0&b=bs": [F, 1.0, 0.0, 0.0, 0.0, F, 0.0],
    })
    claimed = {"Delta": 6, "Delta&R2=0": 6, "Pi": 6, "Pi&R2=0": 6, "D2": 6, "D2&R2=0": 6, "Delta&tn=0": 6,
               "Delta&tn=0&R2=0": 5, "Pi&tn=0&R2=0": 5, "Pi&tn=0&R2=0&b=bs": 4}
    for face, rules in faces.items():
        r = claimed.get(face, 7)
        best = search(sp, rules, r)
        key = f"stack free c_n | {face} | claimed rank {r}"
        results[key] = best
        flag = "  <-- LOW" if best[0][0] < -8 else ""
        print(f"{key}: min log10(s_r/s_1) = {best[0][0]:.2f} at {np.round(best[0][1], 4)}{flag}  [{time.time() - t0:.0f}s]",
              flush=True)
elif group == "bare":
    # bare film (xi, b, tR, tl[, R]); read-out factor eps0 (0: T(0); 1: theta_e), bottom eps
    BF = {
        "interior": [F, F, F, F, F], "R=0": [F, F, F, F, 0.0], "Delta": [F, F, F, ("=", 2), F],
        "Delta&R=0": [F, F, F, ("=", 2), 0.0], "b=bs&R=0": [F, 1.0, F, F, 0.0], "D2": [F, 1.0, F, ("=", 2), F],
        "D2&R=0": [F, 1.0, F, ("=", 2), 0.0], "Pi": [F, F, 0.0, 0.0, F], "Pi&R=0": [F, F, 0.0, 0.0, 0.0],
        "tl=0&R=0": [F, F, F, 0.0, 0.0], "tR=0&R=0": [F, F, 0.0, F, 0.0], "tl=0": [F, F, F, 0.0, F],
    }
    for (e0, e) in [(0.0, 0.0), (0.0, 1.0), (0.0, 0.084), (0.0, -0.5), (1.0, 0.0), (1.0, -1.0), (1.0, 0.5),
                    (1.0, 1.0), (0.5, 0.0), (0.3, -0.6)]:
        sp = cx.Spec("bare", "tied", e0, e, R_known=None)
        for face, rules in BF.items():
            # claimed generic rank (Theorem C2 / Proposition C4, R unknown)
            if e0 == 0 and e == 0:
                r = {"Delta": 4, "Delta&R=0": 3, "D2": 4, "D2&R=0": 3, "Pi": 4, "Pi&R=0": 3}.get(face, 5)
            elif e0 == 0:
                r = {"Pi&R=0": 4}.get(face, 5)
            elif e0 == 1 and e == 0:
                r = {"Delta&R=0": 4, "D2&R=0": 4, "Pi&R=0": 4}.get(face, 5)
            else:
                r = {"Pi&R=0": 4}.get(face, 5)
                if e == 0 and face in ("Delta&R=0", "D2&R=0"):
                    r = 4
            best = search(sp, rules, r)
            key = f"bare eps0={e0:g} eps={e:g} R unknown | {face} | claimed rank {r}"
            results[key] = best
            flag = "  <-- LOW" if best[0][0] < -8 else ""
            print(f"{key}: min log10(s_r/s_1) = {best[0][0]:.2f} at {np.round(best[0][1], 4)}{flag}  [{time.time() - t0:.0f}s]",
                  flush=True)

with open(DIR + f"cx_07_adversarial_{group}.json", "w") as f:
    json.dump(results, f, indent=1)
