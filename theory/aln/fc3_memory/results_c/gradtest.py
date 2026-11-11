import sys, numpy as np
sys.path.insert(0, r"C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/scripts")
import joint_mesh as J
from optimize_fc3 import FC3ModelScaled
m = J.MeshModel(*J.MESH["m553"], "cal1nz", J.RES / "m553")
ref = FC3ModelScaled(str(J.MESH["m553"][0]), J.MESH["m553"][1], tag="cal1nz")
rng = np.random.default_rng(0); u = 0.02 * rng.standard_normal(m.n_var)
a, b = m.evaluate(u), ref.evaluate(u)
print("tau model vs FC3ModelScaled", a["tau"][0] / b["tau"][0] - 1, a["tau"][2] / b["tau"][2] - 1)
W = np.load(J.RES / "m553" / "scaled_cal1nz_W.npy"); cb = np.load(J.RES / "m553" / "cluster_orbits.npz")["col_block"]
import types
models = {"m553": m}
# directional derivative test of the AL function at u with nonzero multipliers via one outer step mu
best, hist = None, None
f = None
def make(mu=50.0):
    out = {}
    def capture(models, obj, cm, Wm, cbk, eps, col, sense, **kw):
        pass
    return out
# reuse internals: build fun by running run_joint with maxouter=0 is not possible; finite-difference check of the objective part
ev = m.evaluate(u, grad=True)
D = -(ev["dN"][0] / ev["N"][0] - ev["dK"][0] / ev["K"][0])
g = m.chain(D, ev["Phi"])
d = rng.standard_normal(m.n_var); h = 1e-6
fp = -np.log(m.evaluate(u + h * d)["tau"][0]); fm = -np.log(m.evaluate(u - h * d)["tau"][0])
print("grad check", (fp - fm) / (2 * h), g @ d)
