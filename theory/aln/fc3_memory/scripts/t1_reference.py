"""Task 1 + 2a (phono3py route): reproduce the memory-gap reference from fc3, and the r0-convention variant.

Checks: W * pp_ref == memory-gap event rates ev_g (bitwise up to rounding); orbit-averaged rates give
kappa_xx 214.30, kappa_zz 201.97 W/(m K), tau_mem 34.25 / 27.45 ps; then make_r0_average = False.
usage: python -B t1_reference.py <events.npz> <resultdir> <label>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, summary  # noqa: E402


def main():
    ev_path, rdir, label = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    geom = M.geom
    out = {"n_modes": geom.n, "n_events": geom.m, "n_event_orbits": M.n_orb,
           "n_mode_orbits": int(M.tie["n_mode_orb"])}
    pp1 = np.load(rdir / "pp_r01_REF.npz")["pp"]
    g_raw = M.W * pp1
    rel = np.abs(g_raw - geom.gphys) / np.maximum(np.abs(geom.gphys), 1e-300)
    out["W_pp_vs_evg_relmax"] = float(rel.max())
    out["W_pp_vs_evg_relmedian"] = float(np.median(rel))
    # reference (memory-gap convention: orbit average of ev_g)
    m_ref = M.metrics(M.orbit_avg(geom.gphys))
    out["memgap_reference"] = summary(m_ref)
    m1 = M.metrics(M.rates(pp1))
    out["fc3_route_r0avgTrue"] = summary(m1)
    m1raw = M.metrics(M.rates(pp1, tie=False))
    out["fc3_route_r0avgTrue_untied"] = summary(m1raw)
    pp0 = np.load(rdir / "pp_r00_REF.npz")["pp"]
    m0 = M.metrics(M.rates(pp0))
    out["fc3_route_r0avgFalse"] = summary(m0)
    # per-mode lifetime change between conventions
    r1, r0 = m1["r"], m0["r"]
    d = r0 / r1 - 1
    out["r0False_vs_True_lifetime_relchange"] = {"min": float(d.min()), "max": float(d.max()),
                                                  "median_abs": float(np.median(np.abs(d))),
                                                  "rms": float(np.sqrt(np.mean(d * d)))}
    # per-event rate ratio statistics (allowed events only)
    ratio = (M.rates(pp0) / np.maximum(M.rates(pp1), 1e-300))
    allowed = M.rates(pp1) > 1e-9 * M.rates(pp1).max()
    q = np.quantile(ratio[allowed], [0.001, 0.01, 0.5, 0.99, 0.999])
    out["r0False_vs_True_event_rate_ratio_quantiles(0.001,0.01,0.5,0.99,0.999)"] = q.tolist()
    out["n_allowed_events"] = int(allowed.sum())
    Path(rdir / f"t1_reference_{label}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
