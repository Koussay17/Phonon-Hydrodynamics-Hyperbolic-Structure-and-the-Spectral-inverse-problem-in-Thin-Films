"""T2 for the film quantity: smearing dependence of the film memory time tau_z(L) (E30 operator) at fixed mesh.

Exact Gaussian re-weighting of the stored events from sigma = 0.1 THz to sigma' <= 0.1 THz, as review/na_spectral.py
(audit 2.4): g' = g (sigma / sigma') exp(-dl^2/2 (1/sigma'^2 - 1/sigma^2)), events kept for |dl| <= 4 sigma'.
For each sigma' the film metrics of film_event.py are evaluated for REF and for the named tensors, and their ratios to
REF at the same sigma' are recorded.
usage: python -B film_smearing.py <events.npz> <out.json> REFDIR [name=pp_dir ...]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class  # noqa: E402
from film_event import film_metrics  # noqa: E402

SIGMAS = (0.05, 0.06, 0.08, 0.10)


def main():
    ev, outp, refdir = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    named = dict(a.split("=", 1) for a in sys.argv[4:])
    t0 = time.perf_counter()
    M = FC3Class(ev, refdir / "pp_r01_REF.npz", allow_missing=True)
    sigma, dl = float(M.geom.labels["sigma"]), np.asarray(M.geom.labels["delta"])
    g = {"REF": M.rates(np.load(refdir / "pp_r01_REF.npz")["pp"])}
    for nm, d in named.items():
        g[nm] = M.rates(np.load(Path(d) / f"pp_r01_{nm}.npz")["pp"])
    res = {"events": ev, "mesh": M.geom.labels["mesh"], "sigma_stored": sigma, "sigmas": list(SIGMAS),
           "named": named, "n_missing_symmetry_images": M.tie["n_missing_images"]}
    for s2 in SIGMAS:
        assert s2 <= sigma + 1e-12
        wgt = (sigma / s2) * np.exp(-dl ** 2 / 2 * (1 / s2 ** 2 - 1 / sigma ** 2))
        keep = np.abs(dl) <= 4 * s2
        r = {nm: film_metrics(M, np.where(keep, gg * wgt, 0.0)) for nm, gg in g.items()}
        ref = r["REF"]
        out = {"REF": ref}
        for nm in named:
            out[nm + "_ratio"] = {k: {q: r[nm][k][q] / ref[k][q] for q in ref[k]} for k in ref}
        res[f"{s2:g}"] = out
        print(s2, json.dumps({k: round(v["tau_z"], 3) for k, v in ref.items()}),
              {nm: round(out[nm + "_ratio"]["100"]["tau_z"], 4) for nm in named}, flush=True)
    res["seconds"] = time.perf_counter() - t0
    outp.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
