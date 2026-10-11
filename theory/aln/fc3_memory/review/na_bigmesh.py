"""Numerical audit, part 5: reference memory time on additional meshes (mesh refinement for S8), memory-lean.

Reference event rates as exported by aln_events.py (ev_g = W pp_ref), orbit-tied by the group average over the 24
operations, operator assembled in place (dense, chunked np.add.at), C + gamma e e^T (gamma = tr C / n) factorised
in place (Cholesky), b projected onto H = e^perp before the solve. The lowest eigenpairs of C on H are obtained by
shift-invert Lanczos with the same factor (no second dense matrix). Reported: tau, kappa, lambda_min, the share of N
carried by the lowest eigenpairs and by low-frequency modes, the RTA analogue, and diag C of the lowest modes.
Option 'block' (4th argument): instead of the orbit (group) average, the gauge-invariant degenerate-block average of
na_blockavg.py, computed from the exported rates: pp_eff = g / (g0w W_shape) with g0w the Gaussian of the detuning and
W_shape = 1/(4 s_p s_a s_b) or 1/(8 s_p s_a^2); block mean over ordered band triples; g = g0w W_shape x mean.
Option sigma' (5th argument, <= sigma of the export): exact Gaussian re-weighting of every event to the smaller width
sigma' with the cut at 4 sigma' (the export holds all events with |detuning| <= 4 sigma).
Option pp file (6th argument): rates of another fc3 tensor on the same event list, g = conv pp g0w / (4 s_p s_a s_b)
(or / 8 s_p s_a^2), from a p3_events.py output; otherwise the exported reference rates are used.
usage: python -B na_bigmesh.py <events.npz> <out.json> [nsig_max|none] [orbit|block] [sigma'|none] [pp.npz]
"""
import json
import os
import sys
import time

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "4")
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

FOURPI = 4 * np.pi


def main():
    t0 = time.perf_counter()
    path, out = sys.argv[1], sys.argv[2]
    nsig_max = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "none" else None
    tying = sys.argv[4] if len(sys.argv) > 4 else "orbit"
    Z = np.load(path)
    f = Z["freqs"]
    Nq, nb = f.shape
    mesh = Z["mesh"].tolist()
    sigma = float(Z["sigma"])
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(int(valid.sum()))
    fv = f.ravel()[valid]
    n = len(fv)
    T = float(Z["T"])
    kT = float(Z["KB"]) * T / float(Z["THzToEv"])
    x = fv / kT
    nocc = 1 / np.expm1(x)
    d = np.sqrt(nocc * (1 + nocc))
    v = Z["gv"].reshape(Nq * nb, 3)[valid]
    e = d * fv
    ehat = e / np.linalg.norm(e)
    b = (d * fv)[:, None] * v
    kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (2 * Nq * float(Z["KB"]) * T ** 2)
    ev_p, ev_a, ev_b = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    g = Z["ev_g"].copy()
    dl = Z["ev_dl"]
    m = len(g)
    s2 = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != "none" else None
    ppfile = sys.argv[6] if len(sys.argv) > 6 else None
    if ppfile is not None:
        PP = np.load(ppfile)
        nidx = newidx
        Pv_, Av_, Bv_ = nidx[ev_p], nidx[ev_a], nidx[ev_b]
        sh_ = np.sinh(fv / kT / 2)
        rep_ = Av_ == Bv_
        conv = float(PP["conv"])
        ppw = PP["pp"] * PP["g0w"] * conv
        g_pp = np.where(rep_, ppw / (8 * sh_[Pv_] * sh_[Av_] ** 2), ppw / (4 * sh_[Pv_] * sh_[Av_] * sh_[Bv_]))
        del ppw
        rel_ref = float(np.max(np.abs(g_pp - g) / g))
        g = g_pp
        del Pv_, Av_, Bv_, sh_, rep_, g_pp, PP
    if s2 is not None:
        assert s2 <= sigma
        g = np.where(np.abs(dl) <= 4 * s2, g * (sigma / s2) * np.exp(-dl ** 2 / 2 * (1 / s2 ** 2 - 1 / sigma ** 2)), 0.0)
    R = {"mesh": mesh, "sigma": sigma, "sigma_prime": s2, "n_modes": n, "n_events": m, "nsig_max": nsig_max, "tying": tying,
         "pp_file": ppfile, "pp_rates_vs_exported_maxrel": (rel_ref if ppfile is not None else None)}
    if tying == "block":
        Pv, Av, Bv = newidx[ev_p], newidx[ev_a], newidx[ev_b]
        qv = np.nonzero(valid)[0] // nb
        ordm = np.lexsort((fv, qv))
        grp = np.empty(n, dtype=np.int64)
        gid, pq, pf = -1, -1, None
        for i in ordm:
            if qv[i] != pq or fv[i] - pf > 1e-6:
                gid += 1
            grp[i] = gid
            pq, pf = qv[i], fv[i]
        sizes = np.bincount(grp)
        sh = np.sinh(fv / kT / 2)
        rep = Av == Bv
        gw = np.exp(-dl ** 2 / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma)
        Wsh = gw * np.where(rep, 1 / (8 * sh[Pv] * sh[Av] ** 2), 1 / (4 * sh[Pv] * sh[Av] * sh[Bv]))
        ga, gb, gp = grp[Av], grp[Bv], grp[Pv]
        G = int(grp.max()) + 1
        keyb = (gp * G + np.minimum(ga, gb)) * G + np.maximum(ga, gb)
        uniq, inv = np.unique(keyb, return_inverse=True)
        same = (ga == gb) & ~rep
        wsum = np.bincount(inv, weights=np.where(same, 2.0, 1.0) * g / Wsh, minlength=len(uniq))
        g = Wsh * wsum[inv] / (sizes[gp] * sizes[ga] * sizes[gb])
        R["n_blocks"] = int(len(uniq))
        R["n_degenerate_modes"] = int(np.sum(sizes[grp] > 1))
        del Pv, Av, Bv, keyb, uniq, inv, same, wsum, Wsh, gw, ga, gb, gp
    # ---- orbit tying by group average (memory-lean)
    M = Nq * nb
    key = (ev_p * M + ev_a) * M + ev_b
    order = np.argsort(key)
    skey = key[order]
    del key
    acc = np.zeros(m)
    n_missing_images = 0
    for Rm in Z["rot_maps"]:
        p2 = Rm[ev_p // nb] * nb + ev_p % nb
        a2 = Rm[ev_a // nb] * nb + ev_a % nb
        b2 = Rm[ev_b // nb] * nb + ev_b % nb
        k2 = (p2 * M + np.minimum(a2, b2)) * M + np.maximum(a2, b2)
        del p2, a2, b2
        pos = np.minimum(np.searchsorted(skey, k2), m - 1)
        found = skey[pos] == k2
        n_missing_images += int((~found).sum())
        acc += np.where(found, g[order[pos]], 0.0)   # a missing image (selection-rule zero, not exported) counts as 0
        del k2, pos, found
    R["n_missing_symmetry_images"] = n_missing_images
    R["tying_max_change_rel_to_max_rate"] = float(np.max(np.abs(acc / len(Z["rot_maps"]) - g)) / g.max())
    if tying == "orbit":
        g = acc / len(Z["rot_maps"])
    del acc, order, skey
    if nsig_max is not None:
        g = np.where(np.abs(dl) <= nsig_max * sigma, g, 0.0)
    R["t_tie"] = time.perf_counter() - t0
    # ---- assembly
    P, A, B = newidx[ev_p], newidx[ev_a], newidx[ev_b]
    del ev_p, ev_a, ev_b, dl
    rep = A == B
    lp = np.sqrt((fv[A] + fv[B]) / fv[P])
    sl = np.stack([P, A, B], axis=1)
    ac = np.stack([-lp, np.where(rep, 2 / lp, 1 / lp), np.where(rep, 0.0, 1 / lp)], axis=1) / d[sl]
    del P, A, B, lp, rep
    C = np.zeros((n, n))
    Cf = C.reshape(-1)
    ch = 1_000_000
    for s in range(0, m, ch):
        sls, acs, gs = sl[s:s + ch], ac[s:s + ch], g[s:s + ch]
        for k in range(3):
            for l in range(3):
                np.add.at(Cf, sls[:, k] * n + sls[:, l], gs * acs[:, k] * acs[:, l])
    del sl, ac
    R["t_assemble"] = time.perf_counter() - t0
    r = np.diag(C).copy()
    gam = float(r.sum() / n)
    Ce = C @ ehat
    R["energy_residual_rel"] = float(np.linalg.norm(Ce) / r.max())
    R["symmetry_check_rel"] = float(np.max(np.abs(C[:200, :] - C[:, :200].T)) / r.max())
    # RTA analogue and low-mode data
    for c in (0, 2):
        bb = b[:, c] ** 2
        KR, NR = float(np.sum(bb / r)), float(np.sum(bb / r ** 2))
        R[f"RTA_{'xyz'[c]}"] = {"tau_ps": NR / KR / FOURPI, "kappa": kf * KR}
    lo = np.argsort(fv)[:8]
    R["lowest_mode_freqs"] = [float(fv[i]) for i in lo]
    R["diagC_lowest_modes"] = [float(r[i]) for i in lo]
    # ---- regularise and factorise in place (C is exactly symmetric; C.T is its Fortran-ordered view)
    for s in range(0, n, 500):
        C[s:s + 500] += gam * ehat[s:s + 500, None] * ehat[None, :]
    cf = sla.cho_factor(C.T, lower=True, overwrite_a=True, check_finite=False)
    R["t_factor"] = time.perf_counter() - t0
    bins = [0, 1, 2, 3, 4, 6, 8, 12, 30]
    R["freq_bins_THz"] = bins
    R["mode_counts_per_bin"] = np.histogram(fv, bins=bins)[0].tolist()
    xs = {}
    for c in (0, 2):
        a = "xyz"[c]
        bH = b[:, c] - (b[:, c] @ ehat) * ehat
        xx = sla.cho_solve(cf, bH, check_finite=False)
        K, N = float(bH @ xx), float(xx @ xx)
        xs[a] = (bH, xx, K, N)
        R[a] = {"K": K, "N": N, "tau_ps": N / K / FOURPI, "kappa": kf * K,
                "b_dot_ehat_rel": float(abs(b[:, c] @ ehat) / np.linalg.norm(b[:, c])),
                "N_share_by_freq": [float(v_) for v_ in np.histogram(fv, bins=bins, weights=xx * xx)[0] / N],
                "K_share_by_freq": [float(v_) for v_ in np.histogram(fv, bins=bins, weights=bH * xx)[0] / K]}
        print(a, R[a]["tau_ps"], R[a]["kappa"], flush=True)
    # ---- lowest eigenpairs on H by shift-invert Lanczos with the Cholesky factor
    OPinv = spla.LinearOperator((n, n), matvec=lambda y: sla.cho_solve(cf, y, check_finite=False), dtype=float)

    def nomatvec(y):
        raise RuntimeError("A matvec not expected in shift-invert mode")

    Aop = spla.LinearOperator((n, n), matvec=nomatvec, dtype=float)
    try:
        lam, V = spla.eigsh(Aop, k=24, sigma=0.0, which="LM", OPinv=OPinv, tol=1e-10, maxiter=5000)
        o = np.argsort(lam)
        lam, V = lam[o], V[:, o]
        R["lam_lowest_H"] = [float(t) for t in lam[:12]]
        R["lam_min_H"] = float(lam[0])
        R["lowest_vec_null_overlap_max"] = float(np.max(np.abs(V.T @ ehat)))
        for a in ("x", "z"):
            bH, xx, K, N = xs[a]
            beta = V.T @ bH
            cN = beta ** 2 / lam ** 2 / N
            R[a]["N_frac_lowest"] = {str(k): float(cN[:k].sum()) for k in (1, 5, 20)}
            R[a]["N_top_pair_in_lowest24"] = {"lambda": float(lam[np.argmax(cN)]), "frac_N": float(cN.max())}
        lowvec = []
        for k in range(4):
            w = V[:, k] ** 2
            t = np.argsort(-w)[:4]
            lowvec.append({"lambda": float(lam[k]), "PR": float(1 / np.sum(w ** 2)), "mean_f": float(w @ fv),
                           "top": [(float(fv[i]), int(np.nonzero(valid)[0][i] % nb), Z["addr"][np.nonzero(valid)[0][i] // nb].tolist(),
                                    float(w[i])) for i in t]})
        R["lowest_eigvecs"] = lowvec
    except Exception as ex:  # noqa: BLE001
        R["eigsh_error"] = repr(ex)
    R["t_total"] = time.perf_counter() - t0
    with open(out, "w") as fh:
        json.dump(R, fh, indent=1)
    print(json.dumps({k: R[k] for k in ("mesh", "n_modes", "n_events", "t_total")}), flush=True)


if __name__ == "__main__":
    main()
