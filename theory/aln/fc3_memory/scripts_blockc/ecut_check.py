"""Pre-registered cutoff choice for Block C (rule fixed before the test runs).
For ecut in (60/480, 80/480): dF = F(disp-00001) - F(disp-00000) at that cutoff, compared with dF at 120/600.
Accept the lowest cutoff with max|dF - dF_120| <= 0.005 * max|dF_120| (max over atoms and components).
Writes ecut_test/result.json and ecut.choice ("ECUTW ECUTR")."""
import json, re, numpy as np
def forces(f):
    t = open(f).read(); assert "JOB DONE" in t, f
    blk = t[t.index("Forces acting on atoms"):t.index("Total force")]
    return np.array([[float(x) for x in m] for m in re.findall(r"force =\s+(\S+)\s+(\S+)\s+(\S+)", blk)])
ref = forces("pbe/disp-00001.out") - forces("pbe/disp-00000.out")
res, choice = {}, "120 600"
for e in ("80_480", "60_480"):
    d = forces("ecut_test/%s/disp-00001.out" % e) - forces("ecut_test/%s/disp-00000.out" % e)
    rel = float(np.abs(d - ref).max() / np.abs(ref).max())
    res[e] = {"max_rel_dF_error": rel, "max_abs_dF_error_Ry_bohr": float(np.abs(d - ref).max())}
res["max_abs_dF_120_Ry_bohr"] = float(np.abs(ref).max())
for e in ("60_480", "80_480"):
    if res[e]["max_rel_dF_error"] <= 0.005:
        choice = e.replace("_", " "); break
res["choice"] = choice; res["rule"] = __doc__
json.dump(res, open("ecut_test/result.json", "w"), indent=1); open("ecut.choice", "w").write(choice + "\n")
print(json.dumps({k: v for k, v in res.items() if k != "rule"}))
