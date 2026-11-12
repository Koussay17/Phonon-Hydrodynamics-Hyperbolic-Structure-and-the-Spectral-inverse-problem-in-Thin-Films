"""cx_07b: search for exact function-level rank drops by solving  G(th) c = 0, |c| = 1  (G: relative real Jacobian
at 25 frequencies in [1e-3, 1e3], full parameter directions c in R^n) with th restricted to a face of the domain.
Batched Levenberg-Marquardt in (log th_free, c); dG/dth by central differences of the exact (forward-mode) Jacobian.
Near-degenerate corners stall at residual ~1e-6; exact degenerate sets converge to ~1e-13.  Each converged point
is classified against the known strata, and the claimed rank of the face is used to remove its generic kernel
(we require c orthogonal to the claimed kernel directions of the face, given as a callable).

Usage: python cx_07b_kernel_solve.py <group> [walkers]
"""
import json
import sys
import time
if len(sys.argv) > 1 and sys.argv[1] in ("tied2", "mixed", "free"):
    sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
    import throttle
    throttle.wait_slot(2)          # light wait before any numerical import
import numpy as np

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
group = sys.argv[1] if len(sys.argv) > 1 else "validate"
MW = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
OM = np.logspace(-3, 3, 25)
rng = np.random.default_rng(77)
F = "f"


def make_face(n, rules):
    free = [i for i, r in enumerate(rules) if r == "f"]

    def th_of(Y):
        m = Y.shape[0]
        TH = np.zeros((m, n))
        for k, i in enumerate(free):
            TH[:, i] = np.exp(Y[:, k])
        for i, r in enumerate(rules):
            if isinstance(r, float):
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


def Gmat(spec, TH):
    Z, J = cx.ZJ_np_batch(spec, TH, OM)
    G = J / Z[:, None, :]
    return np.concatenate([G.real, G.imag], axis=2).transpose(0, 2, 1)     # (m, 2nf, n)


def solve(spec, rules, known_kernel, M=1500, iters=120, h=1e-6):
    """known_kernel(TH) -> (m, q, n) array of the claimed generic kernel basis on the face (q may be 0)."""
    n = spec.n
    free, th_of = make_face(n, rules)
    k = len(free)
    LO, HI = box(spec, free)
    Y = LO + (HI - LO) * rng.uniform(0, 1, (M, k))
    TH = th_of(Y)
    G0 = Gmat(spec, TH)
    S = 1.0 / np.maximum(np.linalg.norm(G0, axis=1), 1e-300)          # fixed column scaling (m, n)
    # start c: smallest right singular vector of the scaled Jacobian, projected off the known kernel
    _, _, Vt = np.linalg.svd(G0 * S[:, None, :])
    C = Vt[:, -1, :].copy()
    lam = np.full(M, 1e-3)

    def proj_known(TH, C):
        Kn = known_kernel(TH)
        if Kn is None or Kn.shape[1] == 0:
            return C, None
        # scaled coordinates: direction d = S*c in raw parameters; known kernel K in raw -> K/S in scaled
        Ks = Kn / S[:, None, :]
        Q, _ = np.linalg.qr(Ks.transpose(0, 2, 1))           # (m, n, q)
        C = C - np.einsum("mnq,mq->mn", Q, np.einsum("mnq,mn->mq", Q, C))
        return C, Q

    C, _Q0 = proj_known(TH, C)
    C /= np.linalg.norm(C, axis=1, keepdims=True)

    def resid(Y, C):
        TH = th_of(np.clip(Y, LO, HI))
        G = Gmat(spec, TH) * S[:, None, :]
        r = np.einsum("mij,mj->mi", G, C)
        Kn = known_kernel(TH)
        if Kn is None or Kn.shape[1] == 0:
            return r, G, None
        Ks = Kn / S[:, None, :]
        Qn, _ = np.linalg.qr(Ks.transpose(0, 2, 1))
        return np.concatenate([r, np.einsum("mnq,mn->mq", Qn, C)], axis=1), G, Qn

    with np.errstate(all="ignore"):
        r, G, Q = resid(Y, C)
        cost = np.sum(r * r, axis=1)
        for it in range(iters):
            # Jacobian wrt Y by central differences (includes the y-dependence of the projector)
            dR = []
            for j in range(k):
                e = np.zeros(k); e[j] = h
                rp, _, _ = resid(Y + e, C)
                rm, _, _ = resid(Y - e, C)
                dR.append((rp - rm) / (2 * h))
            Jy = np.stack(dR, axis=2)                          # (m, nrows, k)
            Jc = G if Q is None else np.concatenate([G, Q.transpose(0, 2, 1)], axis=1)
            # constraint |c|^2 = 1
            rows_c = 2 * C[:, None, :]
            A = np.concatenate([np.concatenate([Jy, Jc], axis=2),
                                np.concatenate([np.zeros((M, 1, k)), rows_c], axis=2)], axis=1)
            b = np.concatenate([r, (np.sum(C * C, axis=1) - 1)[:, None]], axis=1)
            AtA = np.einsum("mik,mil->mkl", A, A)
            Atb = np.einsum("mik,mi->mk", A, b)
            dia = np.einsum("mii->mi", AtA)
            step = -np.linalg.solve(AtA + lam[:, None, None] * (np.eye(k + n)[None] * (dia[:, :, None] + 1e-14)),
                                    Atb[:, :, None])[:, :, 0]
            Yn = np.clip(Y + np.clip(step[:, :k], -0.5, 0.5), LO, HI)
            Cn = C + step[:, k:]
            Cn /= np.linalg.norm(Cn, axis=1, keepdims=True)
            rn, Gn, Qn = resid(Yn, Cn)
            cn = np.sum(rn * rn, axis=1)
            cn[~np.isfinite(cn)] = np.inf
            acc = cn < cost
            Y[acc], C[acc], r[acc], G[acc], cost[acc] = Yn[acc], Cn[acc], rn[acc], Gn[acc], cn[acc]
            if Q is not None:
                Q[acc] = Qn[acc]
            lam[acc] = np.maximum(lam[acc] / 3, 1e-12)
            lam[~acc] = np.minimum(lam[~acc] * 4, 1e8)
    res = np.sqrt(cost)
    onb = np.any(np.abs(Y - LO) < 1e-4, axis=1) | np.any(np.abs(Y - HI) < 1e-4, axis=1)
    return th_of(Y), C * S, res, onb


def none_kernel(TH):
    return None


def delta_kernel(n, iR=2, il=3):
    def f(TH):
        K = np.zeros((TH.shape[0], 1, n)); K[:, 0, iR] = 1; K[:, 0, il] = 1
        return K
    return f


def kvecs(*fns):
    """known kernel built from callables v(TH) -> (m, n) vectors"""
    def f(TH):
        if not fns:
            return None
        return np.stack([fn(TH) for fn in fns], axis=1)
    return f


def unit(n, *idx_val):
    def v(TH):
        out = np.zeros((TH.shape[0], n))
        for i, val in idx_val:
            out[:, i] = val
        return out
    return v


def thin_slab(n, i_xi=0, i_b=1, i_R=None):
    def v(TH):
        out = np.zeros((TH.shape[0], n))
        b = TH[:, i_b]
        out[:, i_xi] = b / (b * b - 1.0)
        out[:, i_R] = 1.0
        return out
    return v


def stack_faces(n=6, free_cn=False):
    """faces of the stack domain; coordinates (xi, b, tR, tl, [tn,] R1, R2)"""
    i2 = 6 if free_cn else 5
    base = [F] * n
    def mk(**fix):
        r = list(base)
        for k, v in fix.items():
            r[int(k[1:])] = v
        return r
    faces = {
        "interior": mk(), "R2=0": mk(**{f"i{i2}": 0.0}), "tR=0": mk(i2=0.0), "tl=0": mk(i3=0.0),
        "Delta": mk(i3=("=", 2)), "Delta&R2=0": mk(i3=("=", 2), **{f"i{i2}": 0.0}), "b=bs": mk(i1=1.0),
        "b=bs&R2=0": mk(i1=1.0, **{f"i{i2}": 0.0}), "tl=0&R2=0": mk(i3=0.0, **{f"i{i2}": 0.0}),
        "tR=0&R2=0": mk(i2=0.0, **{f"i{i2}": 0.0}) if not free_cn else mk(i2=0.0, i6=0.0),
        "Pi": mk(i2=0.0, i3=0.0), "Pi&R2=0": mk(i2=0.0, i3=0.0, **{f"i{i2}": 0.0}),
        "D2": mk(i1=1.0, i3=("=", 2)), "D2&R2=0": mk(i1=1.0, i3=("=", 2), **{f"i{i2}": 0.0}),
    }
    # fix the special index-2 key clash (tR) for tR=0 faces
    faces["tR=0"] = [F, F, 0.0] + [F] * (n - 3)
    faces["tR=0&R2=0"] = [F, F, 0.0] + [F] * (n - 3)
    faces["tR=0&R2=0"][i2] = 0.0
    faces["Pi"] = [F, F, 0.0, 0.0] + [F] * (n - 4)
    faces["Pi&R2=0"] = [F, F, 0.0, 0.0] + [F] * (n - 4)
    faces["Pi&R2=0"][i2] = 0.0
    if free_cn:
        faces.update({
            "tn=0": [F, F, F, F, 0.0, F, F], "tn=0&R2=0": [F, F, F, F, 0.0, F, 0.0],
            "Delta&tn=0": [F, F, F, ("=", 2), 0.0, F, F], "Delta&tn=0&R2=0": [F, F, F, ("=", 2), 0.0, F, 0.0],
            "tn=tR&R2=0": [F, F, F, F, ("=", 2), F, 0.0], "tn=tl&R2=0": [F, F, F, F, ("=", 3), F, 0.0],
            "tn=tl&tR=0&R2=0": [F, F, 0.0, F, ("=", 3), F, 0.0], "Pi&tn=0&R2=0": [F, F, 0.0, 0.0, 0.0, F, 0.0],
            "D2&tn=0&R2=0": [F, 1.0, F, ("=", 2), 0.0, F, 0.0]})
    return faces


def known_stack(e0, e, face, n=6, free_cn=False):
    i2 = 6 if free_cn else 5
    res = unit(n, (2, 1.0), (3, 1.0))
    if free_cn:
        fns = []
        if "Delta" in face or "Pi" in face or "D2" in face:
            fns.append(res)
        if face in ("Delta&tn=0&R2=0", "Pi&tn=0&R2=0"):
            fns.append(thin_slab(n, i_R=i2))
        if face == "D2&tn=0&R2=0":
            fns.append(unit(n, (0, 1.0)))
        return kvecs(*fns)
    if e0 == 0 and e == 0:
        if face in ("Delta", "Pi", "D2"):
            return kvecs(res)
        if face in ("Delta&R2=0", "Pi&R2=0"):
            return kvecs(res, thin_slab(n, i_R=i2))
        if face == "D2&R2=0":
            return kvecs(res, unit(n, (0, 1.0)))
        return kvecs()
    fns = []
    if face == "Pi&R2=0":
        fns.append(thin_slab(n, i_R=i2))
    if e == 0 and face == "Delta&R2=0":
        fns.append(thin_slab(n, i_R=i2))
    if e == 0 and face == "D2&R2=0":
        fns.append(unit(n, (0, 1.0)))
    return kvecs(*fns)


def claimed_kernel_point(e0, e, th, free_cn=False):
    """claimed kernel (Tables 0.3, Proposition C4 / mixed paragraph) at the stratum on which th lies"""
    n = len(th)
    xi, b, tR, tl = th[:4]
    R2 = th[-1]
    tn = th[4] if free_cn else None
    D = abs(tR - tl) <= 1e-6 * max(tR, tl, 1e-300)
    P = tR < 1e-12 and tl < 1e-12
    R0 = R2 < 1e-12
    B1 = abs(b - 1) < 1e-9

    def ev(*iv):
        v = np.zeros(n)
        for i, x in iv:
            v[i] = x
        return v
    res = ev((2, 1.0), (3, 1.0))
    exi = ev((0, 1.0))
    thin = ev((0, b / (b * b - 1.0)), (n - 1, 1.0)) if not B1 else exi
    V = []
    if free_cn:
        T0 = abs(tn) < 1e-12
        if D:
            V.append(res)
        if D and T0 and R0:
            V.append(thin)
            if P and B1:
                V.append(ev((2, -2.0), (4, 1.0)))
        return V
    if e0 == 0 and e == 0:
        if D:
            V.append(res)
        if D and R0:
            V.append(thin)
        return V
    if e0 == e:
        if P and R0:
            V.append(thin)
            if B1:
                V.append(ev((2, 1 - 2 * e), (3, 1.0)))
        return V
    if P and R0:
        V.append(thin)
    elif D and R0 and e == 0:
        V.append(thin)
    if D and (not R0) and e0 == 1 and e <= 0:
        t = tR
        if e == 0 and b < 1 and abs(t / (4 * xi * xi * (1 - b * b)) - 1) < 1e-5 and abs(R2 / (2 * xi * (1 - b * b) / b) - 1) < 1e-5:
            V.append(ev((0, 1.0), (2, -t / xi), (3, t / xi), (5, (b * b - 1) / b)))
        if (e < 0 and abs(b / np.sqrt(1 - e) - 1) < 1e-5 and abs(t / (4 * abs(e) * xi * xi / (1 + abs(e))) - 1) < 1e-5
                and abs(R2 / (2 * abs(e) * xi / b) - 1) < 1e-5):
            V.append(ev((0, 1.0), (2, -t / xi), (3, t / xi), (5, -e / b)))
    return V


def explained(e0, e, th, c, free_cn=False):
    V = claimed_kernel_point(e0, e, th, free_cn)
    if not V:
        return False, "no claimed kernel at this stratum"
    A = np.array(V).T
    coef, *_ = np.linalg.lstsq(A, c, rcond=None)
    rel = np.linalg.norm(A @ coef - c) / np.linalg.norm(c)
    return rel < 1e-5, f"in span of {len(V)} claimed vector(s), rel. residual {rel:.1e}"


results = {}
t0 = time.time()
if group in ("tied", "tied2", "mixed", "free"):
    if group in ("tied", "tied2"):
        convs = [(0.0, 0.0), (1.0, 1.0), (0.084, 0.084), (0.5, 0.5), (2.5, 2.5), (-0.5, -0.5), (1 / 3, 1 / 3), (-2.0, -2.0)]
    elif group == "mixed":
        convs = [(1.0, 0.0), (1.0, -1.0), (0.3, 0.7), (0.0, 1.0), (0.5, 0.0), (1.0, 0.5), (2.0, -1.0)]
    else:
        convs = [None]
    tests = []
    for cv in convs:
        if cv is None:
            for fname, rules in stack_faces(7, True).items():
                tests.append(("stackfree", (None, None), fname, rules, known_stack(None, None, fname, 7, True)))
        else:
            for fname, rules in stack_faces(6, False).items():
                tests.append(("stack", cv, fname, rules, known_stack(cv[0], cv[1], fname)))
    for kind, (e0, e), fname, rules, kk in tests:
        if group in ("tied2", "mixed", "free"):
            throttle.wait_slot(2)      # pause between faces if two other jobs are running
        sp = cx.Spec("stack", "free", xiA=0.2, bA=3.0) if kind == "stackfree" else cx.Spec("stack", "tied", e0, e, xiA=0.2, bA=3.0)
        TH, Cs, res, onb = solve(sp, rules, kk, M=MW, iters=100)
        good = res < 1e-10
        key = f"{kind} eps0={e0} eps={e} | {fname}"
        pts = []
        for i in np.where(good)[0]:
            th = TH[i]; c = Cs[i] / np.max(np.abs(Cs[i]))
            ok, why = explained(e0, e, th, c, free_cn=(kind == "stackfree"))
            pts.append(dict(theta=th.tolist(), c=c.tolist(), res=float(res[i]), known=why if ok else "", why=why))
        unknown = [p for p in pts if not p["known"]]
        results[key] = dict(walkers=MW, n_exact=len(pts), n_unexplained=len(unknown), best_res=float(np.min(res)),
                            unexplained=unknown[:20], known_examples=[p for p in pts if p["known"]][:3])
        flag = "   <-- UNEXPLAINED" if unknown else ""
        print(f"{key}: exact {len(pts)} (unexplained {len(unknown)}), best residual {np.min(res):.1e} [{time.time() - t0:.0f}s]{flag}",
              flush=True)
        for p in unknown[:3]:
            print("      th", np.round(p["theta"], 6), "c", np.round(p["c"], 4), "res %.1e" % p["res"], flush=True)
        with open(DIR + f"cx_07b_{group}.json", "w") as f:
            json.dump(results, f, indent=1)
    raise SystemExit
if group == "validate":
    tests = [("stack", (1.0, 0.0), "Delta", [F, F, F, ("=", 2), F, F], none_kernel),
             ("stack", (1.0, -1.0), "Delta", [F, F, F, ("=", 2), F, F], none_kernel),
             ("stack", (0.0, 0.0), "Delta (local; generic kernel removed)", [F, F, F, ("=", 2), F, F], delta_kernel(6))]
else:
    raise SystemExit
for kind, (e0, e), fname, rules, kk in tests:
    sp = cx.Spec(kind, "tied", e0, e, xiA=0.2, bA=3.0)
    TH, Cs, res, onb = solve(sp, rules, kk, M=MW)
    good = (res < 1e-10)
    print(f"{kind} eps0={e0} eps={e} {fname}: walkers {MW}, exact solutions {int(good.sum())}, "
          f"best residual {np.min(res):.1e}, boundary-stalled {int(onb.sum())} [{time.time() - t0:.0f}s]", flush=True)
    for i in np.where(good)[0][:5]:
        th = TH[i]; c = Cs[i] / np.max(np.abs(Cs[i]))
        xi, b, t = th[0], th[1], th[2]
        extra = ""
        if e0 == 1 and e == 0 and b < 1:
            extra = f" E0 test t/(4xi^2(1-b^2)) = {t / (4 * xi * xi * (1 - b * b)):.8f}, R2/(2xi(1-b^2)/b) = {th[5] / (2 * xi * (1 - b * b) / b):.8f}"
        if e0 == 1 and e == -1:
            extra = f" E_-1 test b/sqrt2 = {b / np.sqrt(2):.8f}, t/(2xi^2) = {t / (2 * xi * xi):.8f}, R2/(sqrt2 xi) = {th[5] / (np.sqrt(2) * xi):.8f}"
        print("    th =", np.round(th, 6), " c =", np.round(c, 4), " res %.1e" % res[i], "(box edge)" if onb[i] else "", extra, flush=True)
    results[f"{kind} eps0={e0} eps={e} {fname}"] = dict(n_exact=int(good.sum()),
                                                        points=[TH[i].tolist() for i in np.where(good)[0][:50]])
with open(DIR + f"cx_07b_{group}.json", "w") as f:
    json.dump(results, f, indent=1)
