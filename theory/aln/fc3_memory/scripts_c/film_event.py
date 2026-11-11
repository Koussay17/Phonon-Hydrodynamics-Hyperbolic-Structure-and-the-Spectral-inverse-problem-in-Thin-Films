"""E30: film-regularised memory time of the conserving event operator.

C_L = C(g) + diag(d), d_alpha = |v_z,alpha| / (4 pi L) (THz units of the operator: physical rate 4 pi r [1/ps];
v in THz*A, L in A). C_L is positive definite (e^T D e > 0), so no energy-mode regulariser is added.
K_L = b^T C_L^-1 b, N_L = |C_L^-1 b|^2, tau_mem(L) = N_L / K_L / (4 pi) [ps]; kappa_L = kf K_L (film-suppressed).
L = inf reproduces FC3Class.metrics (bulk, regularised on e). Same Matthiessen regulator as E26 (ir_film.py), now on
the full operator instead of RTA. Rates orbit-tied as in FC3Class.
usage: python -B film_event.py <events.npz> <pp_dir> <out.json> name [name ...]   (REF is always evaluated)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI  # noqa: E402

LS_NM = (30, 100, 300, 1000, 3000, 10000, None)


def _C_padded(M, g):
    """C(g) on the (n+1)^2 bincount buffer, without a sliced copy: the dummy slot n is decoupled (row/col 0, diagonal 1),
    so the factor of the padded matrix restricted to the first n entries is the factor of C (memory: one n^2 buffer)."""
    ops = M.ops
    n1 = ops.n1
    Cf = np.bincount(ops.pair_flat.ravel(), weights=(g[:, None] * ops.pair_val).ravel(), minlength=n1 * n1)
    Cf = Cf.reshape(n1, n1)
    Cf[-1, :] = 0.0
    Cf[:, -1] = 0.0
    Cf[-1, -1] = 1.0
    return Cf


def film_metrics(M, g, cols=(0, 2)):
    n = M.geom.n
    vz = np.abs(M.geom.vel[:, 2])
    dC = M.ops.diag_full(g)
    trC = float(dC.sum())
    di = np.diag_indices(n)
    out = {}
    for L in LS_NM:
        CL = _C_padded(M, g)                  # rebuilt per L and factorised in place
        if L is None:
            e = M.geom.ehat
            a = trC / n
            for i0 in range(0, n, 1024):      # rank-one update by row blocks (no n^2 temporary)
                CL[i0:min(i0 + 1024, n), :n] += (a * e[i0:i0 + 1024])[:, None] * e[None, :]
            key = "bulk"
            dg = dC
        else:
            CL[di] += vz / (FOURPI * L * 10.0)
            key = str(L)
            dg = dC + vz / (FOURPI * L * 10.0)
        cf = sla.cho_factor(CL, lower=True, overwrite_a=True, check_finite=False)
        del CL
        r = {}
        for c in cols:
            b = M.geom.b[:, c]
            x = sla.cho_solve(cf, np.append(b, 0.0), check_finite=False)[:n]
            K, N = float(b @ x), float(x @ x)
            nm = "xyz"[c]
            r[f"kappa_{nm}"], r[f"tau_{nm}"] = M.kf * K, N / K / FOURPI
            # RTA on the same diagonal (diag C + d), for comparison with E26
            r[f"tauRTA_{nm}"] = float(np.sum(b * b / dg ** 2) / np.sum(b * b / dg)) / FOURPI
        del cf
        out[key] = r
    return out


def main():
    ev, d, outp, names = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4:]
    t0 = time.perf_counter()
    M = FC3Class(ev, d / "pp_r01_REF.npz", allow_missing=True)
    res = {"events": ev, "mesh": M.geom.labels["mesh"], "n_modes": int(M.geom.n), "L_nm": [str(x) for x in LS_NM],
           "n_missing_symmetry_images": M.tie["n_missing_images"]}
    if M.tie["n_missing_images"]:
        gr = M.rates(np.load(d / "pp_r01_REF.npz")["pp"])
        res["missing_images_max_rate_rel"] = float(gr[M.tie["events_with_missing_images"]].max() / gr.max())
    for nm in ["REF"] + [n for n in names if n != "REF"]:
        g = M.rates(np.load(d / f"pp_r01_{nm}.npz")["pp"])
        res[nm] = film_metrics(M, g)
        print(nm, json.dumps({k: [round(v["tau_x"], 3), round(v["tau_z"], 3)] for k, v in res[nm].items()}),
              flush=True)
    ref = res["REF"]
    for nm in names:
        if nm == "REF":
            continue
        res[nm + "_ratio"] = {k: {q: res[nm][k][q] / ref[k][q] for q in ref[k]} for k in ref}
    res["seconds"] = time.perf_counter() - t0
    outp.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
