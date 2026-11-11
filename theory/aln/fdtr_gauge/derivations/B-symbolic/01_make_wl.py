"""Writes 01_rank.wl: Wolfram check of the singular values at a subset of the points of
results/01_points.json (quadrupole product, symbolic D, 50-digit SVD).  The output of the
Wolfram run is stored in results/01_rank_wl.txt and compared with results/01_rank_mp.txt."""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
pts = json.load(open(HERE / "results" / "01_points.json"))
SEL = [0, 3, 7, 24, 28, 32, 36, 38]
SUB = {1: [3], 2: [1, 5], 4: [0, 2, 5, 7], 8: list(range(8))}
def wl(fr):
    n, d = fr.split("/")
    return f"({n}/{d})"
lines = []
lines.append("zq[p_, x1_, x2_, x3_, x4_, bs_] := Module[{aa, bb, zs, m},")
lines.append("  aa = x1 Sqrt[p (1 + x3 p)/(1 + x4 p)]; bb = x2 Sqrt[p (1 + x4 p)/(1 + x3 p)];")
lines.append("  zs = 1/(bs Sqrt[p]); m = {{Cosh[aa], Sinh[aa]/bb}, {bb Sinh[aa], Cosh[aa]}};")
lines.append("  (m[[1, 1]] zs + m[[1, 2]])/(m[[2, 1]] zs + m[[2, 2]])];")
lines.append("grad = D[zq[p, x1, x2, x3, x4, bs], {{x1, x2, x3, x4}}];")
lines.append("eq[mat_] := Transpose[#/Norm[#] & /@ Transpose[mat]];")
lines.append("sv[mat_] := Module[{s = SingularValueList[eq[mat], Tolerance -> 0]}, N[s/First[s], 4]];")
lines.append("rowsAt[pt_, oms_] := Table[N[grad /. Join[pt, {p -> I om}], 60], {om, oms}];")
lines.append("realJ[r_] := Flatten[{Re[#], Im[#]} & /@ r, 1];")
lines.append("res = {};")
for i in SEL:
    pt = pts[i]
    rule = (f"{{x1 -> {wl(pt['xi1'])}, x2 -> {wl(pt['b'])}, x3 -> {wl(pt['tau_R'])}, "
            f"x4 -> {wl(pt['tau_l'])}, bs -> {wl(pt['bs'])}}}")
    oms = "{" + ", ".join(wl(x) for x in pt["omegas"]) + "}"
    lines.append(f"pt = {rule}; oms = {oms}; r = rowsAt[pt, oms];")
    parts = []
    for N, idx in SUB.items():
        ii = "{" + ", ".join(str(k + 1) for k in idx) + "}"
        parts.append(f'"N={N} real" -> sv[realJ[r[[{ii}]]]], "N={N} complex" -> sv[r[[{ii}]]]')
    lines.append(f'AppendTo[res, {{"point {i:02d}", ' + ", ".join(parts) + "}];")
lines.append("res")
(HERE / "01_rank.wl").write_text("\n".join(lines) + "\n", encoding="utf-8")
print((HERE / "01_rank.wl").read_text())
