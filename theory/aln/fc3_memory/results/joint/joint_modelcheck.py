import sys, json, numpy as np
sys.path.insert(0, r"C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/scripts")
import joint_mesh as J
m = J.MeshModel(*J.MESH["m775"], "cal1", J.RES / "m553")
T = json.load(open(J.RES / "m775x" / "mesh_transfer.json"))
for nm in ("cal1_c0_max", "cal1_c2_max", "cal1_c0_min"):
    u = np.load(J.RES / "m553" / "t3_cal" / f"{nm}_y.npy")
    ev = m.evaluate(u)
    tr = [ev["tau"][c] / m.base["tau"][c] for c in (0, 2)]
    kr = [ev["K"][c] / m.K0[c] for c in (0, 2)]
    print(nm, "linear model 7x7x5 block-avg tau", tr, "kappa", kr, "| phono3py route orbit-avg", T[nm]["tau_ratio"], T[nm]["kappa_ratio"])
