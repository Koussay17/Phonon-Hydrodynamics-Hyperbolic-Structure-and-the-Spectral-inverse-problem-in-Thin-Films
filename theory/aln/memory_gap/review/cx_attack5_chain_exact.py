"""Attack 5c: exact-arithmetic counterexample to the 'self-energy linewidth' reading of L3.

Chain network (energies 1, 2, 4, 8, 16; events 2->1+1, 4->2+2, 8->4+4, 16->8+8; s_max = 2),
entropy scale D = I (any positive diagonal is admissible in the abstract class), rates g = 1.
r_diag  = diag C (L3's r);  r_SE = self-energy-type linewidth in which a repeated-daughter slot
contributes coefficient 2 (instead of 2^2 = 4) -- verified on the AlN data: (diagC - Gamma_phono3py)/Gamma
= 0.503 x (repeated-daughter share), correlation 0.992 (cx_attack5_rse_aln.json).
We exhibit an explicit rational b orthogonal to e with K/K_RTA(r_SE) < 1/s_max = 1/2, and check
K/K_RTA(r_diag) >= 1/2 for the same b (L3 itself holds). Exact rationals throughout (sympy).
Output: review/cx_attack5_chain_exact.json
"""
import json

import numpy as np
import scipy.linalg as sla
import sympy as spy

eps = [1, 2, 4, 8, 16]
S = [[2, -1, 0, 0, 0], [0, 2, -1, 0, 0], [0, 0, 2, -1, 0], [0, 0, 0, 2, -1]]
n = 5
A = spy.Matrix(S).T                       # n x m
C = A * A.T                               # g = 1
e = spy.Matrix(eps)
r_diag = [C[i, i] for i in range(n)]
r_se = [sum((abs(S[k][i]) if abs(S[k][i]) != 2 else 2) for k in range(4) if S[k][i] != 0) for i in range(n)]
# float adversarial b for r_SE: minimise b^T C^+ b / b^T R_SE^{-1} b over b orthogonal to e
Cf = np.array(C.tolist(), float)
ef = np.array(eps, float) / np.linalg.norm(eps)
Q = sla.null_space(ef[None, :])
CH = Q.T @ Cf @ Q
Rinv = Q.T @ np.diag(1 / np.array(r_se, float)) @ Q
w, V = sla.eigh(np.linalg.inv(CH), Rinv)
b_f = Q @ V[:, 0]
b_f /= np.max(np.abs(b_f))
# rationalise and project exactly onto e-perp
b = spy.Matrix([spy.Rational(round(float(x) * 1000), 1000) for x in b_f])
b = b - e * (e.dot(b) / e.dot(e))
assert e.dot(b) == 0
Creg = C + e * e.T
x = Creg.LUsolve(b)
assert (C * x - b).applyfunc(spy.simplify) == spy.zeros(n, 1)
K = spy.Rational(b.dot(x))
KR_diag = sum(b[i] ** 2 / r_diag[i] for i in range(n))
KR_se = sum(b[i] ** 2 / r_se[i] for i in range(n))
out = {"r_diag": [int(v) for v in r_diag], "r_selfenergy": r_se, "b": [str(v) for v in b],
       "K": str(K), "K_RTA_diag": str(KR_diag), "K_RTA_selfenergy": str(KR_se),
       "K_over_KRTA_diag": str(K / KR_diag), "K_over_KRTA_diag_float": float(K / KR_diag),
       "K_over_KRTA_selfenergy": str(K / KR_se), "K_over_KRTA_selfenergy_float": float(K / KR_se),
       "s_max": 2, "min_generalised_eig_float_selfenergy": float(w[0])}
out["L3_holds_with_diag"] = bool(K / KR_diag >= spy.Rational(1, 2))
out["selfenergy_reading_violated"] = bool(K / KR_se < spy.Rational(1, 2))
print(json.dumps(out, indent=1))
json.dump(out, open(r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review\cx_attack5_chain_exact.json", "w"), indent=1)
