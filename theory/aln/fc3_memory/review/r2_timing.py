import os, sys, time
nth = sys.argv[1]
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = nth
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261005-fc3-memory\review")
import numpy as np
from r2_opt import ScaledModel
R = r"C:\Users\Koussay\ResearchLab\orchestration\campaigns\20261005-fc3-memory\results\m553"
t0 = time.time()
m = ScaledModel(R + r"\events_m553_s0.1.npz", R, R, "")
print("init", time.time() - t0, flush=True)
u = np.load(R + r"\t3_orb\orb0.1_c0_max_y.npy")
t0 = time.time()
for _ in range(5):
    ev = m.evaluate(u, grad=True)
    g = m.chain(ev["dK"][0], ev["Phi"])
    rel2, By = m.orbit_rel2(u)
    gg = m.Bmat.T @ By
print("eval+grad", (time.time() - t0) / 5)
print("tau", ev["tau"][0] / (4 * np.pi), ev["tau"][2] / (4 * np.pi), "orbmax", np.sqrt(rel2.max()))
