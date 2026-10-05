"""Infrared scaling of K and tau_mem under mesh refinement (no optimisation).

For each (d, alpha, N): Debye-type event set, prescribed diagonal r = A eps^alpha (A fixed by the
total rate of the Klemens reference), reference rates g = KL projection of the Klemens rates onto
{W g = r} (exact diagonal), K0 = K(g), tau_ref = tau(g); RTA sums with the same r.
Reported K's are normalised by N^d (continuum normalisation of the mode sum).
Prediction (RTA, b^2 ~ const, r ~ q^alpha): K finite iff alpha < d; tau finite iff 2 alpha < d;
tau_RTA(N) ~ N^(2 alpha - d) for d/2 < alpha < d, ~ log N at 2 alpha = d, ~ N^alpha for alpha > d.
usage: python ir_scaling.py <d> <alpha,alpha,...> <N,N,...> <outfile.json>
"""
from __future__ import annotations

import json
import os
import sys
import time

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "4")
sys.dont_write_bytecode = True

import numpy as np

from debye_events import build
from memgap import rta
from optim import Slice


def one(d, alpha, N):
    t0 = time.time()
    geom = build(d, N)
    g_phys = geom.gphys
    r_phys = geom.W @ g_phys
    if alpha is None:
        r, g = r_phys, g_phys
    else:
        r = geom.eps ** alpha
        r = r * (r_phys.sum() / r.sum())
        S = Slice(geom, r, 1.0)
        g, ok = S.project_diag(g_phys)
        if not ok:
            return {"d": d, "alpha": alpha, "N": N, "error": "diag projection failed"}
    R = geom.response(g)
    KR, NR, tR = rta(geom, r)
    b = geom.b[:, 0]
    nk, w = geom.check_kernel(g)
    rec = {"d": d, "alpha": alpha, "N": N, "n": geom.n, "m": geom.m,
           "K0_norm": R["K"] / N ** d, "tau_ref": R["tau"], "K_RTA_norm": KR / N ** d, "tau_RTA": tR,
           "N_ref_norm": R["N"] / N ** d, "N_RTA_norm": NR / N ** d,
           "CS_lower": R["K"] / float(b @ b), "kernel_dim": nk, "lambda_min": float(np.sort(w)[1]),
           "seconds": time.time() - t0}
    return rec


def main():
    d = int(sys.argv[1])
    alphas = [None if a == "phys" else float(a) for a in sys.argv[2].split(",")]
    Ns = [int(x) for x in sys.argv[3].split(",")]
    out = sys.argv[4]
    recs = []
    if os.path.exists(out):
        recs = json.load(open(out))
    done = {(r["alpha"], r["N"]) for r in recs}
    for a in alphas:
        for N in Ns:
            if (a, N) in done:
                continue
            rec = one(d, a, N)
            recs.append(rec)
            print(json.dumps(rec), flush=True)
            json.dump(recs, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
