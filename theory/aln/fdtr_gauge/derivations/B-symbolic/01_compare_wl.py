"""Compares the Wolfram singular values (results/01_rank_wl.txt) with an mpmath recomputation
(symbolic SymPy derivatives, 50 digits) at the same points and frequency subsets."""
import json, re, sys
from fractions import Fraction
from pathlib import Path
import mpmath as mp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M
mp.mp.dps = 50
txt = (HERE / "results" / "01_rank_wl.txt").read_text()
pts = json.load(open(HERE / "results" / "01_points.json"))
F = lambda s: Fraction(*map(int, s.split("/")))
SUB = {1: [3], 2: [1, 5], 4: [0, 2, 5, 7], 8: list(range(8))}
jsym = M.sympy_jacobian_mp(with_R=False)
blocks = re.findall(r'\{"point (\d+)", (.*?)\}\}', txt.replace("\n", " "))
worst = 0
for pid, body in blocks:
    pt = pts[int(pid)]
    rows = [jsym(mp.mpc(0, M._mpf(F(o))), F(pt["xi1"]), F(pt["b"]), F(pt["tau_R"]), F(pt["tau_l"]), F(pt["bs"]), 0)
            for o in pt["omegas"]]
    for m in re.finditer(r'"N=(\d) (real|complex)" -> \{([^}]*)', body + "}"):
        N, kind, vals = int(m.group(1)), m.group(2), m.group(3)
        wv = [float(v.split("`")[0]) for v in vals.split(",")]
        sel = [rows[i] for i in SUB[N]]
        if kind == "real":
            s = M.singular_values(M.equilibrate(M.real_jacobian(sel)))
        else:
            s = M.singular_values(M.equilibrate(M.complex_jacobian(sel)), complex_=True)
        s = [x / s[0] for x in s]
        d = max(abs(float(x) - y) for x, y in zip(s, wv))
        worst = max(worst, d)
print(f"compared {len(blocks)} points x 8 matrices; max |s_mpmath - s_wolfram| (normalised) = {worst:.2e}")
(HERE / "results" / "01_compare_wl.txt").write_text(
    f"compared {len(blocks)} points x 8 matrices; max |s_mpmath - s_wolfram| (normalised) = {worst:.2e}\n"
    "(Wolfram values printed with ~19 significant digits; agreement at that level means the two\n"
    " independent symbolic differentiations and SVDs agree.)\n")
