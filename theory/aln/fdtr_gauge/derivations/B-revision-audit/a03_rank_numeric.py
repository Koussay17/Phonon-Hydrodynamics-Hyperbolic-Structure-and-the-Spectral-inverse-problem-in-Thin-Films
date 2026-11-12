"""Audit of the rank tables (Theorems C2, T3; Propositions C4, C5): 60-digit function-level Jacobians at selected
points, with an implementation independent of models_mp.py (single-valued C, X, Y form of the film quadrupole,
B-analytic Lemma 2, so that one-sided derivatives at tau = 0 are derivatives of the analytic extension).

Film block  Phi = N0 (C L + N X)/(Y L + N C),  L = R + 1/(bs sqrt p),  N = 1 + tau_n p,  N0 = 1 + tau_o p,
C = cosh(sqrt z), X = (xi/b) g S(z), Y = xi b p S(z), S(z) = sinh(sqrt z)/sqrt z, z = xi^2 p g.
Stack       Z_tot = Mob(T_A; R1 + Phi(R -> R2)),  T_A Fourier layer (xiA, bA) = (0.2, 3).
Units b_s = 1.  Derivatives: central differences, h = 1e-22 at 60 digits (truncation ~1e-44).
The main target is the point where the audit's hand analysis disagrees with Proposition C4 (a)/(b):
theta_e read-out with bottom coefficient eps = 1 at the Fourier corner, b = b_s, R = 0.
Output: a03_rank_numeric.json
"""
import os
import json
import mpmath as mp

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
H = mp.mpf(10) ** -22
WS = [mp.mpf(10) ** (mp.mpf(k) / 6 - 2) for k in range(25)]      # 1e-2 .. 1e2


def film_block(par, p, conv):
    xi, b, tR, tl = par["xi"], par["b"], par["tR"], par["tl"]
    R = par.get("R", 0)
    eps = conv.get("eps", 0)
    tn = par["tn"] if conv.get("free") else eps * tl
    obs = conv["obs"]
    to = {"T": 0, "theta_e": tl, "theta_q": tn}[obs]
    g = (1 + tR * p) / (1 + tl * p)
    z = xi ** 2 * p * g
    rz = mp.sqrt(z)
    C = mp.cosh(rz)
    Sz = mp.sinh(rz) / rz if z != 0 else mp.mpf(1)
    X = xi / b * g * Sz
    Y = xi * b * p * Sz
    L = R + 1 / mp.sqrt(p)
    N, N0 = 1 + tn * p, 1 + to * p
    return N0 * (C * L + N * X) / (Y * L + N * C)


def response(par, w, conv):
    p = mp.mpc(0, w)
    if conv.get("stack"):
        q = dict(par)
        q["R"] = par["R2"]
        Win = par["R1"] + film_block(q, p, conv)
        A = mp.mpf("0.2") * mp.sqrt(p)
        B = 3 * mp.sqrt(p)
        return (mp.cosh(A) * Win + mp.sinh(A) / B) / (B * mp.sinh(A) * Win + mp.cosh(A))
    return film_block(par, p, conv)


def jac(par, names, conv):
    cols = []
    for nm in names:
        a, c = dict(par), dict(par)
        a[nm] = par[nm] + H
        c[nm] = par[nm] - H
        cols.append([(response(a, w, conv) - response(c, w, conv)) / (2 * H) for w in WS])
    J = mp.matrix(2 * len(WS), len(names))
    for j, col in enumerate(cols):
        for i, d in enumerate(col):
            J[2 * i, j], J[2 * i + 1, j] = mp.re(d), mp.im(d)
    return J


def analyse(par, names, conv, kernels):
    J = jac(par, names, conv)
    norms = [mp.sqrt(sum(J[i, j] ** 2 for i in range(J.rows))) for j in range(J.cols)]
    Jn = J.copy()
    big = max(norms)
    zero_cols = [j for j in range(J.cols) if norms[j] <= mp.mpf(10) ** -30 * big]   # identically zero to working accuracy
    for j in range(J.cols):
        for i in range(J.rows):
            Jn[i, j] = J[i, j] / norms[j] if j not in zero_cols else 0
    sv = sorted([abs(x) for x in mp.svd_r(Jn, compute_uv=False)], reverse=True)
    sv = [x / sv[0] for x in sv]
    rank = sum(1 for x in sv if x > mp.mpf(10) ** -15)
    nJ = mp.sqrt(sum(J[i, j] ** 2 for i in range(J.rows) for j in range(J.cols)))
    kres = []
    for c in kernels:
        v = J * mp.matrix([mp.mpf(x) if not isinstance(x, mp.mpf) else x for x in c])
        kres.append(float(mp.sqrt(sum(v[i] ** 2 for i in range(v.rows))) / (nJ * mp.sqrt(sum(mp.mpf(x) ** 2 for x in c)))))
    return rank, [mp.nstr(x, 5) for x in sv], kres, zero_cols


one = mp.mpf(1)
s2 = mp.sqrt(2)
CASES = [
    # ---- Proposition C4 at eps = 1 (theta_e read-out, bottom coefficient tau_l): Fourier corner, b = b_s, R = 0
    ("C4 eps=1 corner b=bs R=0, R known", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R=0), ["xi", "b", "tR", "tl"],
     dict(obs="theta_e", eps=1), [[1, 0, 0, 0], [0, 0, -1, 1]], "text: rank 3, kernel e_xi; audit: rank 2"),
    ("C4 eps=1 corner b=bs R=0, R unknown", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R=0), ["xi", "b", "tR", "tl", "R"],
     dict(obs="theta_e", eps=1), [[1, 0, 0, 0, 0], [0, 0, -1, 1, 0]], "text: rank 4, kernel e_xi; audit: rank 3"),
    ("C4 eps=0.5 corner b=bs R=0, R unknown (control)", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R=0),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_e", eps=mp.mpf("0.5")), [[1, 0, 0, 0, 0]], "text and audit: rank 4"),
    ("C4 eps=1 corner b=2 R=0, R unknown (thin slab)", dict(xi=mp.mpf("0.7"), b=mp.mpf(2), tR=0, tl=0, R=0),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_e", eps=1), [[1, 0, 0, 0, mp.mpf(3) / 2]], "rank 4, kernel (1,0,0,0,(b^2-1)/b)"),
    # ---- Proposition C4 exotic sets
    ("C4 E_0 (eps=0)", dict(xi=one, b=mp.mpf("0.6"), tR=mp.mpf("2.56"), tl=mp.mpf("2.56"), R=mp.mpf(32) / 15),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_e", eps=0),
     [[1, 0, -mp.mpf("2.56"), mp.mpf("2.56"), (mp.mpf("0.36") - 1) / mp.mpf("0.6")]], "rank 4"),
    ("C4 E_eps (eps=-1)", dict(xi=one, b=s2, tR=mp.mpf(2), tl=mp.mpf(2), R=s2), ["xi", "b", "tR", "tl", "R"],
     dict(obs="theta_e", eps=-1), [[1, 0, -2, 2, 1 / s2]], "rank 4"),
    ("C4 Delta R>0 off E_0 (eps=0, control)", dict(xi=one, b=mp.mpf("0.6"), tR=mp.mpf("2.56"), tl=mp.mpf("2.56"), R=mp.mpf(2)),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_e", eps=0), [], "rank 5"),
    # ---- Proposition C5 (theta at both faces, eps = 0.3) corner
    ("C5 eps=0.3 corner b=bs R=0, R unknown", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R=0), ["xi", "b", "tR", "tl", "R"],
     dict(obs="theta_q", eps=mp.mpf("0.3")), [[1, 0, 0, 0, 0], [0, 0, mp.mpf("0.4"), 1, 0]], "rank 3"),
    ("C5 eps=0.5 tl=0 tR>0 R=0, R unknown", dict(xi=mp.mpf("0.7"), b=mp.mpf("1.7"), tR=mp.mpf("0.4"), tl=0, R=0),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_q", eps=mp.mpf("0.5")), [], "rank 5"),
    ("C5 eps=1 Delta t>0 R=0, R unknown", dict(xi=mp.mpf("0.7"), b=mp.mpf("1.7"), tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R=0),
     ["xi", "b", "tR", "tl", "R"], dict(obs="theta_q", eps=1), [], "rank 5"),
    # ---- Theorem C2: T(0) read-out
    ("C2 tied eps=0.084 Delta t>0 R=0, R unknown", dict(xi=mp.mpf("0.7"), b=mp.mpf("1.7"), tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R=0),
     ["xi", "b", "tR", "tl", "R"], dict(obs="T", eps=mp.mpf("0.084")), [], "rank 5"),
    ("C2 free Delta tn=0 b=bs R=0, R known", dict(xi=mp.mpf("0.7"), b=one, tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R=0, tn=0),
     ["xi", "b", "tR", "tl", "tn"], dict(obs="T", free=True), [[1, 0, 0, 0, 0], [0, 0, 1, 1, 0]], "rank 3"),
    ("C2 free Delta tn=-0.3 R=0, R unknown", dict(xi=mp.mpf("0.7"), b=mp.mpf("1.7"), tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R=0, tn=mp.mpf("-0.3")),
     ["xi", "b", "tR", "tl", "tn", "R"], dict(obs="T", free=True), [[0, 0, 1, 1, 0, 0]], "rank 5"),
    # ---- Theorem T3: stack
    ("T3 entropic Pi R2=0 b=bs", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R1=mp.mpf("0.3"), R2=0),
     ["xi", "b", "tR", "tl", "R1", "R2"], dict(obs="theta_q", eps=1, stack=True),
     [[1, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 0]], "rank 4"),
    ("T3 entropic Delta R2=0 t>0 b=bs", dict(xi=mp.mpf("0.7"), b=one, tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R1=mp.mpf("0.3"), R2=0),
     ["xi", "b", "tR", "tl", "R1", "R2"], dict(obs="theta_q", eps=1, stack=True), [], "rank 6"),
    ("T3 local Delta R2=0 b=2", dict(xi=mp.mpf("0.7"), b=mp.mpf(2), tR=mp.mpf("0.4"), tl=mp.mpf("0.4"), R1=mp.mpf("0.3"), R2=0),
     ["xi", "b", "tR", "tl", "R1", "R2"], dict(obs="theta_q", eps=0, stack=True),
     [[0, 0, 1, 1, 0, 0], [mp.mpf(2) / 3, 0, 0, 0, 0, 1]], "rank 4"),
    ("T3 eps=-0.5 Pi R2=0 b=bs", dict(xi=mp.mpf("0.7"), b=one, tR=0, tl=0, R1=mp.mpf("0.3"), R2=0),
     ["xi", "b", "tR", "tl", "R1", "R2"], dict(obs="theta_q", eps=mp.mpf("-0.5"), stack=True),
     [[1, 0, 0, 0, 0, 0], [0, 0, 2, 1, 0, 0]], "rank 4"),
]

OUT = []
for name, par, names, conv, kers, claim in CASES:
    par = {k: mp.mpf(v) for k, v in par.items()}
    rk, sv, kres, zc = analyse(par, names, conv, kers)
    rec = {"case": name, "rank": rk, "n": len(names), "normalized_singular_values": sv,
           "kernel_residuals": kres, "zero_columns": [names[j] for j in zc], "claim": claim}
    OUT.append(rec)
    print(f"{name:52s} rank {rk}/{len(names)}  sv {sv[-3:]}  zero cols {[names[j] for j in zc]}  ker_res {['%.1e' % x for x in kres]}  [{claim}]", flush=True)

with open(os.path.join(HERE, "a03_rank_numeric.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
