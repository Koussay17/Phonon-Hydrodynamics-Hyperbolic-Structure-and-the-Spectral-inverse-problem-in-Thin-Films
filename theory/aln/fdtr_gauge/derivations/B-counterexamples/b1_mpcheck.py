"""Exact versus practical rank deficiency: recompute the equilibrated singular values in mpmath.

For a parameter point (beta, tR, tl) (xi1 = bs = 1) and a log grid, compute the column-equilibrated
relative Jacobian at dps digits and compare its singular values with the double-precision ones.
An exact rank deficiency shows s_4 ~ 10^-dps; a practical one gives the same s_4 at every precision.
Usage: python b1_mpcheck.py beta tR tl [pad ppd dps variant]
"""
import sys, os, json
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


def mp_equilibrated_svals(beta, tR, tl, pad=2.0, ppd=12, dps=50, variant="code", tscale="log",
                          wlist=None, with_R=False, R=0, equilibrate=True):
    mp.mp.dps = dps
    if wlist is None:
        w, wts = zm.wide_grid(1.0, float(beta), float(tR), float(tl), 1.0, pad=pad, ppd=ppd)
    else:
        w = np.asarray(wlist, float)
        wts = np.ones_like(w)
    ws = [mp.mpf(x) for x in w]
    wt = [mp.mpf(x) for x in wts]
    M = zm.real_rel_jac_mp(ws, mp.mpf(1), mp.mpf(beta), mp.mpf(tR), mp.mpf(tl), 1, mp.mpf(R),
                           variant, with_R, tscale, wt)
    if equilibrate:
        for j in range(M.cols):
            nrm = mp.sqrt(mp.fsum(M[i, j] ** 2 for i in range(M.rows)))
            for i in range(M.rows):
                M[i, j] /= nrm
    sv = zm.svals_mp(M)
    # double precision on the same grid
    Md = zm.real_rel_jac_np(w, 1.0, float(beta), float(tR), float(tl), 1.0, float(R), variant,
                            with_R, tscale, wts)
    if equilibrate:
        Md = Md / np.linalg.norm(Md, axis=0)
    svd = np.linalg.svd(Md, compute_uv=False)
    return [float(x) for x in sv], svd.tolist(), len(w)


if __name__ == "__main__":
    beta, tR, tl = (float(x) for x in sys.argv[1:4])
    pad = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
    ppd = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    dps = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    variant = sys.argv[7] if len(sys.argv) > 7 else "code"
    sv, svd, n = mp_equilibrated_svals(beta, tR, tl, pad, ppd, dps, variant)
    print(json.dumps({"beta": beta, "tR": tR, "tl": tl, "pad": pad, "ppd": ppd, "dps": dps,
                      "n_omega": n, "svals_mp": sv, "svals_double": svd}))
