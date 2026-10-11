"""Task 5 (part 1): independent amplitude engine (amp_engine.py, NumPy) versus phono3py's Rust kernel.

For every event of the fixed list, pp_engine = |Phi|^2 * unit_conversion is compared with phono3py's pp
exported by p3_events.py at the same representative triplet and band triple. Both phase conventions.
Optional extra fc3 files (compact .npy) with their phono3py pp files: name=fc3.npy:pp.npz.
usage: python -B t5_engine_check.py <resultdir> <fc3.hdf5> [name=fc3.npy:pp.npz ...]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from amp_engine import AmpEngine  # noqa: E402


def compare(eng, fc3c, P):
    t0 = time.perf_counter()
    phi = eng.phi_events(fc3c, P["tp"], P["bands"])
    dt = time.perf_counter() - t0
    pp_e = np.abs(phi) ** 2 * eng.uc
    pp_p = P["pp"]
    big = pp_p > 1e-9 * pp_p.max()
    rel = np.abs(pp_e[big] - pp_p[big]) / pp_p[big]
    return {"seconds": dt, "n_compared_rel": int(big.sum()), "rel_max": float(rel.max()),
            "rel_median": float(np.median(rel)),
            "abs_max_over_ppmax_all": float(np.abs(pp_e - pp_p).max() / pp_p.max())}


def main():
    rdir = Path(sys.argv[1])
    with h5py.File(sys.argv[2], "r") as h:
        fc3 = h["fc3"][:]
    S = np.load(rdir / "static_r01.npz")
    p2s = S["p2s"]
    fc3c = fc3[p2s] if fc3.shape[0] == fc3.shape[1] else fc3
    out = {}
    gam = {"freq_maxdiff": float(np.abs(S["gamma_f_parentGamma"] - S["gamma_f_parent1"]).max()),
           "eigv_maxdiff": float(np.abs(S["gamma_e_parentGamma"] - S["gamma_e_parent1"]).max())}
    out["gamma_context"] = gam
    for r0 in (1, 0):
        P = np.load(rdir / f"pp_r0{r0}_REF.npz")
        eng = AmpEngine(rdir / f"static_r0{r0}.npz", r0avg=bool(r0), unit_conversion=float(P["unit_conversion"]))
        out[f"REF_r0avg{r0}"] = compare(eng, fc3c, P)
        print(json.dumps({f"REF_r0avg{r0}": out[f"REF_r0avg{r0}"]}), flush=True)
    for spec in sys.argv[3:]:
        name, rest = spec.split("=", 1)
        fpath, ppath = rest.split(":", 1)
        P = np.load(ppath)
        r0 = 1 if "_r01_" in Path(ppath).name else 0
        eng = AmpEngine(rdir / f"static_r0{r0}.npz", r0avg=bool(r0), unit_conversion=float(P["unit_conversion"]))
        out[name] = compare(eng, np.load(fpath), P)
        print(json.dumps({name: out[name]}), flush=True)
    (rdir / "t5_engine_check.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
