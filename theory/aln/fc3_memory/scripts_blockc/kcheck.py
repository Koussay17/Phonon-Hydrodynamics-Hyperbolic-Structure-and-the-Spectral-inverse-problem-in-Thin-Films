"""Code-path check at the production cutoff: disp-00001 computed with K_POINTS gamma (Gamma tricks, production path)
and with K_POINTS automatic 1 1 1 0 0 0 (complex path) must agree: |dE| < 1e-6 Ry and max|dF| < 1e-5 Ry/bohr.
Motivation: in QE 6.7 (Ubuntu build) the Gamma-point 72-atom supercell gave energies below the variational minimum at
80/480 and 120/600 Ry (folding test against 18 x primitive cell on a 3x3x2 k-mesh), while 60/480 passed to 3.5e-7 Ry."""
import re, json, numpy as np
def read(f):
    t = open(f).read(); assert "JOB DONE" in t, f
    E = float(re.findall(r"^!\s+total energy\s+=\s+(\S+)", t, re.M)[-1])
    i = t.index("Forces acting on atoms"); j = t.index("Total force", i)
    F = np.array(re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", t[i:j]), float)
    return E, F
Eg, Fg = read("pbe/disp-00001.out"); Ek, Fk = read("kcheck/disp-00001.out")
r = {"dE_Ry": abs(Eg - Ek), "max_dF_Ry_bohr": float(np.abs(Fg - Fk).max())}
r["pass"] = r["dE_Ry"] < 1e-6 and r["max_dF_Ry_bohr"] < 1e-5
json.dump(r, open("kcheck/result.json", "w"), indent=1); print("kcheck", r)
if r["pass"]: open("kcheck.ok", "w").write("ok\n")
raise SystemExit(0 if r["pass"] else 1)
