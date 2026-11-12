"""Scan of kernel dimensions at special values of the interface coefficient (audit of C2, C4, C5, mixed faces).

The exact strata check of b2_02 keeps epsilon symbolic, i.e. it gives the rank for generic epsilon on each stratum;
the real quantifier elimination of b2_03 decides only {rank < n}.  Neither sees a larger kernel at a special epsilon
inside a degenerate stratum.  This script builds the Lemma K system itself (from (i), (O), (E), verified in
a02_lemmaK_symbolic.py), substitutes epsilon in {generic, 0, 1/2, 1, -1, 2, 1/3} and reports every (model, stratum,
epsilon) whose exact rank differs from the generic-epsilon rank on that stratum.
Output: a06_special_eps_scan.out
"""
import os
import itertools
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


p = sp.Symbol("p")
xi, b, bs = sp.symbols("xi b b_s", positive=True)
tR, tl, R, t = sp.symbols("tau_R tau_l R t", nonnegative=True)
eps, eps0 = sp.symbols("epsilon epsilon_0", real=True)
cxi, cb, ctR, ctl, cR = sp.symbols("cxi cb ctR ctl cR")


def system(obs, R_unknown, eps_bot, eps_top):
    """Rows of the linear system (i), (O), (E) in the direction components; tau_n = eps_bot tl, tau_o = eps_top tl
    (obs='T': eps_top = 0; 'theta_e': eps_top = 1; 'theta_q': eps_top = eps_bot; 'mixed': eps_top given)."""
    comps = [cxi, cb, ctR, ctl] + ([cR] if R_unknown else [])
    c_R = cR if R_unknown else 0
    tn, to = eps_bot * tl, eps_top * tl
    dtn, dto = eps_bot * ctl, eps_top * ctl
    N, N0 = 1 + tn * p, 1 + to * p
    gn, gd = 1 + tR * p, 1 + tl * p
    # (i) times gn gd N0:  2 c_b/b gn gd N0 - [c_tR p gd - c_tl p gn] N0 - 2 dto p gn gd
    I1 = 2 * cb / b * gn * gd * N0 - (ctR * p * gd - ctl * p * gn) * N0 - 2 * dto * p * gn * gd
    kappa = cxi * N0 - xi * dto * p
    I2 = N * dto - N0 * dtn - 2 * b * R * kappa
    I3 = c_R * N * N0 * b * gd - kappa * (b ** 2 * gd / bs ** 2 - b ** 2 * R ** 2 * p * gd - gn * N ** 2)
    rows = []
    for I in (I1, I2, I3):
        P = sp.Poly(sp.expand(I), p)
        for co in P.all_coeffs():
            row = [sp.expand(sp.diff(co, c)) for c in comps]
            if any(r != 0 for r in row):
                rows.append(row)
    return sp.Matrix(rows), len(comps)


STRATA = {
    "generic": {},
    "Delta": {tR: t, tl: t},
    "Delta,R=0": {tR: t, tl: t, R: 0},
    "Delta,b=bs": {tR: t, tl: t, b: bs},
    "Delta,b=bs,R=0": {tR: t, tl: t, b: bs, R: 0},
    "Pi": {tR: 0, tl: 0},
    "Pi,R=0": {tR: 0, tl: 0, R: 0},
    "Pi,R=0,b=bs": {tR: 0, tl: 0, R: 0, b: bs},
    "Pi,b=bs": {tR: 0, tl: 0, b: bs},
    "tl=0": {tl: 0},
    "tl=0,R=0": {tl: 0, R: 0},
    "tl=0,R=0,b=bs": {tl: 0, R: 0, b: bs},
    "tR=0": {tR: 0},
    "tR=0,R=0": {tR: 0, R: 0},
    "tR=0,R=0,b=bs": {tR: 0, R: 0, b: bs},
    "b=bs": {b: bs},
    "b=bs,R=0": {b: bs, R: 0},
    "R=0": {R: 0},
}
SPECIAL = [sp.Integer(0), sp.Rational(1, 2), sp.Integer(1), sp.Integer(-1), sp.Integer(2), sp.Rational(1, 3)]
MODELS = [("T", None), ("theta_e", 1), ("theta_q", "same")]
findings = []
for obs, top in MODELS:
    for R_unknown in (False, True):
        for sname, sub in STRATA.items():
            def rank_for(ev, e0=None):
                et = 0 if obs == "T" else (1 if obs == "theta_e" else ev)
                if e0 is not None:
                    et = e0
                M, n = system(obs, R_unknown, ev, et)
                Ms = M.subs(sub).applyfunc(sp.expand)
                return Ms.rank(simplify=True), n, Ms
            rg, n, Mg = rank_for(eps)
            line = f"{obs:8s} R_unknown={R_unknown!s:5s} {sname:16s} generic-eps rank {rg}/{n}"
            extra = []
            for ev in SPECIAL:
                rk, _, Ms = rank_for(ev)
                if rk != rg:
                    ker = [list(v.T) for v in Ms.nullspace(simplify=True)]
                    extra.append(f"eps={ev}: rank {rk} ker {[[sp.simplify(e) for e in k] for k in ker]}")
                    findings.append((obs, R_unknown, sname, str(ev), rk, rg, n))
            log(line + ("   SPECIAL: " + "; ".join(extra) if extra else ""))

# mixed faces: top eps0, bottom eps, R unknown, on the corner and Delta strata (eps0 != eps)
log("---- mixed faces (top eps0 tl, bottom eps tl), R unknown")
for sname in ("Delta", "Delta,R=0", "Pi,R=0", "Pi,R=0,b=bs", "tl=0", "tl=0,R=0"):
    sub = STRATA[sname]
    for e0, e1 in itertools.product([sp.Integer(0), sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 10)],
                                    [sp.Integer(0), sp.Rational(7, 10), sp.Integer(-1), sp.Integer(1)]):
        if e0 == e1:
            continue
        M, n = system("mixed", True, e1, e0)
        Ms = M.subs(sub).applyfunc(sp.expand)
        rk = Ms.rank(simplify=True)
        if rk < n:
            ker = [list(v.T) for v in Ms.nullspace(simplify=True)]
            log(f"  {sname:12s} (eps0, eps) = ({e0}, {e1}): rank {rk}/{n} ker {[[sp.simplify(e) for e in k] for k in ker]}")
log("special-epsilon findings (obs, R_unknown, stratum, eps, rank, generic rank, n):")
for f_ in findings:
    log("  ", f_)
with open(os.path.join(HERE, "a06_special_eps_scan.out"), "w") as fh:
    fh.write("\n".join(out) + "\n")
