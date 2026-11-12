"""RS0: exploration of the 7x7x5 / 5x5x3 mode sets for the full-operator slab: zero-velocity modes, parity map
(q_x, q_y, q_z) -> (q_x, q_y, -q_z) among the stored grid maps, symmetry of eps and v_z under it, velocity and
mean-free-path ranges (RTA diagonal of REF). Output: out/rs0_explore.json"""
import json, os, sys, time
from pathlib import Path
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np
FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c")
HERE = Path(__file__).resolve().parent
EV = {"m775": FC3 / "m775x/events_m775_s0.1.npz", "m553": FC3 / "m553/events_m553_s0.1.npz"}
res = {}
for mesh, path in EV.items():
    Z = np.load(path)
    f = Z["freqs"]; Nq, nb = f.shape
    addr = Z["addr"]; msh = Z["mesh"]; maps = Z["rot_maps"]
    gv = Z["gv"]
    # parity map: addr (a,b,c) -> (a,b,-c) mod mesh
    target = np.stack([addr[:, 0] % msh[0], addr[:, 1] % msh[1], (-addr[:, 2]) % msh[2]], 1)
    key = lambda A: (A[:, 0] % msh[0]) * msh[1] * msh[2] + (A[:, 1] % msh[1]) * msh[2] + (A[:, 2] % msh[2])
    kt = key(target)
    found = []
    for r, R in enumerate(maps):
        if np.array_equal(key(addr[R]), kt):
            found.append(r)
    out = {"Nq": int(Nq), "nb": int(nb), "parity_map_indices": found}
    R = maps[found[0]]
    valid = f > float(Z["cutoff"])
    # symmetry of eps and v under (q, j) -> (Rq, j)
    dfe = np.abs(f[R] - f)[valid[R] & valid]
    vz = gv[:, :, 2] * 100.0      # THz*A -> m/s
    dv = np.abs(vz[R] + vz)[valid]
    dvx = np.abs(gv[R][:, :, :2] - gv[:, :, :2]).max()
    out.update({"max_abs_dfreq_THz": float(dfe.max()), "max_rel_dfreq": float((dfe / f[valid]).max()),
                "max_abs_vz_plus_vzP_mps": float(dv.max()), "max_abs_dvxy_THzA": float(dvx)})
    # mismatches beyond tolerance and whether they sit in degenerate subspaces
    bad = np.argwhere((np.abs(vz[R] + vz) > 1e-6 * (np.abs(vz) + 1.0)) & valid)
    out["n_vz_mismatch"] = int(len(bad))
    degen = 0
    for (q, j) in bad:
        if np.any(np.abs(np.delete(f[q], j) - f[q, j]) < 1e-6):
            degen += 1
    out["n_vz_mismatch_in_degenerate"] = degen
    avz = np.abs(vz[valid])
    zero = avz < 1e-12
    out.update({"n_valid": int(valid.sum()), "n_zero_vz": int(zero.sum()), "n_moving": int((~zero).sum()),
                "min_abs_vz_moving_mps": float(avz[~zero].min()), "max_abs_vz_mps": float(avz.max()),
                "quantiles_abs_vz_moving": [float(x) for x in np.quantile(avz[~zero], [0.001, 0.01, 0.1, 0.5, 0.9])],
                "max_abs_vz_on_zero_set": float(avz[zero].max()) if zero.any() else None})
    # q_z = 0 plane check
    plane = (addr[:, 2] % msh[2] == 0)
    out["zero_vz_equals_qz0_plane"] = bool(np.array_equal(zero, np.repeat(plane[:, None], nb, 1)[valid]))
    res[mesh] = out
    print(mesh, json.dumps(out), flush=True)
(HERE / "out" / "rs0_explore.json").write_text(json.dumps(res, indent=1))
