"""Finite-beam (twin gaussian_response, w_p = 10 um, w_r = 8 um, isotropic Cattaneo AlN) vs 1D tau_R phase signature."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, r"C:/Users/Koussay/these/src")
import twin_fdtr as tw
from ab_common import BASE, dlnZ

HERE = os.path.dirname(os.path.abspath(__file__))
f = np.logspace(5, np.log10(2e8), 16)
t0 = time.time()
def H(tau):
    aln = tw.Material(BASE["lam"], BASE["C"], tau_r=tau)
    st = tw.au_aln_sapphire(aln, aln_thickness=BASE["e"], au_thickness=BASE["e_au"], r_au_aln=BASE["R1"], r_aln_sapphire=BASE["R2"])
    return tw.gaussian_response(f, st, 10e-6, 8e-6, passivity="ignore")
h = 1e-3
s3d = np.angle(H(27.5e-12 * (1 + h)) / H(27.5e-12 * (1 - h))) / (2 * h)
s1d = dlnZ(2 * np.pi * f, BASE, "tau_r").imag
for fi, a, b in zip(f, s3d, s1d):
    print("f = %9.3e Hz  3D %.4e  1D %.4e deg/w  ratio %.4f" % (fi, np.degrees(a), np.degrees(b), a / b))
json.dump({"f": f.tolist(), "s3d_rad": s3d.tolist(), "s1d_rad": s1d.tolist(), "runtime_s": time.time() - t0},
          open(os.path.join(HERE, "results", "x2_beam_vs_1d.json"), "w"), indent=1)
print("runtime %.1f s" % (time.time() - t0))
