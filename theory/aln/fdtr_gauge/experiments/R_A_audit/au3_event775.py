"""AU3: independent recomputation of tau_R and tau_l of the film-regularised conserving event operator on 7x7x5
for REF and the twelve three-mesh optima, with two energy-conserving completions of the |v_z|/L regulator.

Independent of FC3Class/FastOps/memgap (only the raw event list, the stored amplitudes and the symmetry maps are read):
  * modes: full grid (q, j) with f > cutoff; e = sqrt(n(n+1)) eps, b = e v_z (aln_geometry conventions re-derived);
  * event vectors a_alpha: -lp/dscale at P, ld (2 ld if A = B)/dscale at A, ld/dscale at B, lp = sqrt((fa+fb)/fp);
  * rates g = W pp, W = conv g0w / (4 s_P s_A s_B) (repeated daughter: / (8 s_P s_A^2)), s = sinh(x/2);
  * orbit tying: event images under the 24 stored grid maps located by sorted integer keys, orbits = connected
    components of the image graph (scipy.sparse.csgraph), g averaged per orbit;
  * C = sum g a a^T assembled by a COO sparse sum (duplicates summed), dense;
  * regulator d = |v_z| / (4 pi L[A]) in operator units, two completions on the even sector:
      E-QDQ : C + Q D Q                     (campaign, ra1_film_event_ext.py)
      E-phys: C + D - (D e)(D e)^T/(e^T D e)  (relaxation of the wall process to its own local temperature, i.e.
              diffuse re-emission with the flux-weighted temperature; identical to E-QDQ on the odd sector)
  * solves by the KKT bordered system [[Chat, e], [e^T, 0]] (dense LU): x = Chat^+ b, y = Chat^+ u, u = Q V x;
    tau_R = |x|^2/(4 pi K), tau_l = |e|^2 (u.y)/(4 pi K^2) (operator units: physical rate = 4 pi g).
Output: out/au3_event775.json
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import scipy.sparse as sps  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
RA1 = HERE.parent / "R_A" / "out" / "ra1_film_m775.json"
EV = FC3 / "results_c/m775x/events_m775_s0.1.npz"
PPD = FC3 / "results_c/m775x/film2_pp"
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
LS = (30.0, 100.0, 300.0, 500.0, 1000.0)
FOURPI = 4 * np.pi


def geometry():
    Z = np.load(EV)
    f = Z["freqs"]
    Nq, nb = f.shape
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(valid.sum())
    eps = f.ravel()[valid]
    kT = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    vel = Z["gv"].reshape(Nq * nb, 3)[valid]
    p_full, a_full, b_full = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    P, A, B = newidx[p_full], newidx[a_full], newidx[b_full]
    assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
    x = eps / kT
    nocc = 1.0 / np.expm1(x)
    dsc = np.sqrt(nocc * (1 + nocc))
    e = dsc * eps
    rep = A == B
    lp = np.sqrt((eps[A] + eps[B]) / eps[P])
    ld = 1.0 / lp
    vP = -lp / dsc[P]
    vA = np.where(rep, 2 * ld, ld) / dsc[A]
    vB = np.where(rep, 0.0, ld / dsc[B])
    s = np.sinh(x / 2)
    kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (2 * Nq * float(Z["KB"]) * float(Z["T"]) ** 2)
    # orbits from the stored grid maps (full-grid event keys)
    NF = Nq * nb
    lo, hi = np.minimum(a_full, b_full), np.maximum(a_full, b_full)
    key = (p_full * NF + lo) * NF + hi
    order = np.argsort(key)
    skey = key[order]
    m = key.size
    rows, cols = [np.arange(m)], [np.arange(m)]
    n_missing = 0
    for R in Z["rot_maps"]:
        def img(full):
            q, j = full // nb, full % nb
            return R[q] * nb + j
        ip, ia, ib = img(p_full), img(a_full), img(b_full)
        k2 = (ip * NF + np.minimum(ia, ib)) * NF + np.maximum(ia, ib)
        pos = np.searchsorted(skey, k2)
        pos = np.clip(pos, 0, m - 1)
        found = skey[pos] == k2
        n_missing += int((~found).sum())
        rows.append(np.arange(m)[found])
        cols.append(order[pos[found]])
    G = sps.coo_matrix((np.ones(sum(r.size for r in rows)), (np.concatenate(rows), np.concatenate(cols))), shape=(m, m))
    n_orb, orb = connected_components(G, directed=False)
    return {"eps": eps, "vz": vel[:, 2], "e": e, "P": P, "A": A, "B": B, "rep": rep, "vP": vP, "vA": vA, "vB": vB,
            "s": s, "kf": kf, "orb": orb, "n_orb": int(n_orb), "n_missing_images": n_missing, "n": int(eps.size),
            "m": int(m), "umk": Z["ev_umk"]}


def rates(G, pp, conv, g0w):
    P, A, B, rep, s = G["P"], G["A"], G["B"], G["rep"], G["s"]
    W = np.where(rep, conv * g0w / (8 * s[P] * s[A] * s[A]), conv * g0w / (4 * s[P] * s[A] * s[np.where(rep, A, B)]))
    g = W * pp
    sums = np.bincount(G["orb"], weights=g, minlength=G["n_orb"])
    cnt = np.bincount(G["orb"], minlength=G["n_orb"])
    return (sums / cnt)[G["orb"]]


def assemble(G, g):
    n = G["n"]
    idx = [G["P"], G["A"], G["B"]]
    val = [G["vP"], G["vA"], G["vB"]]
    r_, c_, w_ = [], [], []
    for i in range(3):
        for j in range(3):
            r_.append(idx[i]); c_.append(idx[j]); w_.append(g * val[i] * val[j])
    Cs = sps.coo_matrix((np.concatenate(w_), (np.concatenate(r_), np.concatenate(c_))), shape=(n, n))
    C = Cs.toarray()          # COO -> dense sums duplicates
    return 0.5 * (C + C.T)


def kinetic(C, e, vz, d, variant):
    n = e.size
    ehat = e / np.linalg.norm(e)
    if variant == "E-QDQ":
        h = d * ehat
        Ch = C + np.diag(d) - np.outer(ehat, h) - np.outer(h, ehat) + (ehat @ h) * np.outer(ehat, ehat)
    elif variant == "E-phys":
        g = d * e
        Ch = C + np.diag(d) - np.outer(g, g) / (e @ g)
    else:
        raise KeyError(variant)
    Kkt = np.zeros((n + 1, n + 1))
    Kkt[:n, :n] = Ch
    Kkt[:n, n] = e
    Kkt[n, :n] = e
    lu = sla.lu_factor(Kkt, check_finite=False)
    b = vz * e
    x = sla.lu_solve(lu, np.append(b, 0.0), check_finite=False)[:n]
    K = b @ x
    u = vz * x
    u = u - ehat * (ehat @ u)
    y = sla.lu_solve(lu, np.append(u, 0.0), check_finite=False)[:n]
    M = u @ y
    return {"K": float(K), "tau_ps": float(x @ x / K / FOURPI), "taul_ps": float((e @ e) * M / (FOURPI * K * K)),
            "residual_energy": float(abs(ehat @ x) / np.linalg.norm(x))}


def main():
    t0 = time.perf_counter()
    G = geometry()
    print("geometry", G["n"], "modes", G["m"], "events", G["n_orb"], "orbits", "missing images", G["n_missing_images"],
          f"{time.perf_counter() - t0:.0f}s", flush=True)
    ref = np.load(PPD / "pp_r01_REF.npz")
    conv, g0w = float(ref["conv"]), ref["g0w"]
    ra1 = json.loads(RA1.read_text())
    res = {"n_modes": G["n"], "n_events": G["m"], "n_orbits": G["n_orb"], "n_missing_images": G["n_missing_images"],
           "orbits_match_campaign": G["n_orb"] == ra1["n_event_orbits"], "tensors": {}}
    for nm in ["REF"] + WIT:
        Pz = np.load(PPD / f"pp_r01_{nm}.npz")
        g0w_same = bool(np.array_equal(Pz["g0w"], g0w))
        g = rates(G, Pz["pp"], conv, g0w)
        C = assemble(G, g)
        e, vz = G["e"], G["vz"]
        en_res = float(np.linalg.norm(C @ e) / np.linalg.norm(np.diag(C) * e))
        ent = {"g0w_same_as_REF": g0w_same, "energy_residual": en_res}
        for L in LS:
            d = np.abs(vz) / (FOURPI * L * 10.0)
            row = {}
            for var in ("E-QDQ", "E-phys"):
                row[var] = kinetic(C, e, vz, d, var)
            r1 = ra1[nm][f"{L:.6g}"]
            row["campaign"] = {"tau_z": r1["tau_z"], "taul_kin": r1["taul_kin"], "kappa_z": r1["kappa_z"]}
            row["rel_diff_tau_vs_campaign"] = abs(row["E-QDQ"]["tau_ps"] / r1["tau_z"] - 1)
            row["rel_diff_taul_vs_campaign"] = abs(row["E-QDQ"]["taul_ps"] / r1["taul_kin"] - 1)
            row["rel_diff_kappa_vs_campaign"] = abs(G["kf"] * row["E-QDQ"]["K"] / r1["kappa_z"] - 1)
            row["phys_over_QDQ_taul"] = row["E-phys"]["taul_ps"] / row["E-QDQ"]["taul_ps"]
            ent[f"{L:.6g}"] = row
        res["tensors"][nm] = ent
        r5 = ent["500"]
        print(nm, f"en_res {en_res:.1e}", f"500nm tau {r5['E-QDQ']['tau_ps']:.4f} taul QDQ {r5['E-QDQ']['taul_ps']:.3f} "
              f"phys {r5['E-phys']['taul_ps']:.3f} | rel diff vs campaign tau {r5['rel_diff_tau_vs_campaign']:.1e} "
              f"taul {r5['rel_diff_taul_vs_campaign']:.1e} kappa {r5['rel_diff_kappa_vs_campaign']:.1e}",
              f"({time.perf_counter() - t0:.0f}s)", flush=True)
        (OUT / "au3_event775.json").write_text(json.dumps(res, indent=1))
    # excursions relative to REF (both completions)
    exc = {}
    for nm in WIT:
        exc[nm] = {}
        for L in LS:
            k = f"{L:.6g}"
            a, r = res["tensors"][nm][k], res["tensors"]["REF"][k]
            exc[nm][k] = {var: {"tau": a[var]["tau_ps"] / r[var]["tau_ps"] - 1,
                                "taul": a[var]["taul_ps"] / r[var]["taul_ps"] - 1,
                                "delta": (a[var]["tau_ps"] - a[var]["taul_ps"]) / (r[var]["tau_ps"] - r[var]["taul_ps"]) - 1}
                          for var in ("E-QDQ", "E-phys")}
    res["excursions"] = exc
    res["max_rel_diff_vs_campaign"] = {q: max(res["tensors"][nm][f"{L:.6g}"][q] for nm in res["tensors"] for L in LS)
                                       for q in ("rel_diff_tau_vs_campaign", "rel_diff_taul_vs_campaign",
                                                 "rel_diff_kappa_vs_campaign")}
    res["seconds"] = time.perf_counter() - t0
    (OUT / "au3_event775.json").write_text(json.dumps(res, indent=1))
    print("max rel diff vs campaign:", res["max_rel_diff_vs_campaign"], flush=True)


if __name__ == "__main__":
    main()
