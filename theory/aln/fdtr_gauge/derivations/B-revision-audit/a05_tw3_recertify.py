"""Independent re-certification of the three-frequency twins of Proposition TW3 (local convention).

Independent of ad_arb.py / models_mp.py: reflection form for the film, tanh form for the transducer,
    Z_film = (1/B)(1 + r E)/(1 - r E),  r = (B L - 1)/(B L + 1),  E = exp(-2A),  L = R2 + 1/sqrt(p),
    W = R1 + Z_film,  Z_tot = (W + th/BA)/(1 + BA th W),  th = (1 - EA)/(1 + EA),  EA = exp(-2 xiA sqrt(p)),
    BA = bA sqrt(p),  (xiA, bA) = (1/5, 3),  b_s = 1,  p = i omega,  A = xi sqrt(p g), B = b sqrt(p/g).
Own forward-mode jets over Arb complex balls.  Exact rational target and frequencies.  For each twin:
 (1) Newton refinement in 400-bit ball arithmetic from the 40-digit values of b4_07_certify_twin.json;
 (2) Krawczyk test on the box X = x~ + [-r, r]^6 with r = 1e-35 (theirs: 1e-40), 300 bits:
        K(X) = x~ - Y F(x~) + (I - Y J(X)) (X - x~),  Y ~ J(x~)^-1;  K(X) in int(X)  =>  unique zero in X, J(X) regular;
 (3) certified bounds: tau_R - tau_l > 0 on X (off Delta), the twin is separated from the target, and the response
     mismatch at omega = 9/10 is certified nonzero.
Output: a05_tw3_recertify.json
"""
import os
import json
from flint import acb, arb, arb_mat, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
NP = 6


class J:
    """Jet: value and gradient (6 components), entries acb."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = acb(v)
        self.d = d if d is not None else [acb(0)] * NP

    def __add__(self, o):
        o = o if isinstance(o, J) else J(o)
        return J(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
    __radd__ = __add__

    def __neg__(self):
        return J(-self.v, [-a for a in self.d])

    def __sub__(self, o):
        o = o if isinstance(o, J) else J(o)
        return self + (-o)

    def __rsub__(self, o):
        return J(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, J) else J(o)
        return J(self.v * o.v, [a * o.v + self.v * b for a, b in zip(self.d, o.d)])
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, J) else J(o)
        inv = 1 / o.v
        q = self.v * inv
        return J(q, [(a - q * b) * inv for a, b in zip(self.d, o.d)])

    def __rtruediv__(self, o):
        return J(o) / self


def jexp(x):
    e = x.v.exp()
    return J(e, [e * a for a in x.d])


def jsqrt(x):
    r = x.v.sqrt()
    h = 1 / (2 * r)
    return J(r, [h * a for a in x.d])


def ztot(theta, w, grad=True):
    """theta: 6 arb balls (xi, b, tR, tl, R1, R2); w: arb.  Returns a jet."""
    V = []
    for i, t in enumerate(theta):
        d = [acb(0)] * NP
        if grad:
            d[i] = acb(1)
        V.append(J(acb(t), d))
    xi, b, tR, tl, R1, R2 = V
    p = J(acb(0, w))
    sq = jsqrt(p)
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * jsqrt(p * g)
    B = b * jsqrt(p / g)
    L = R2 + 1 / sq
    BL = B * L
    r = (BL - 1) / (BL + 1)
    E = jexp(-2 * A)
    Zf = (1 + r * E) / (B * (1 - r * E))
    W = R1 + Zf
    xiA, bA = arb(fmpq(1, 5)), arb(3)
    EA = jexp(-2 * xiA * sq)
    th = (1 - EA) / (1 + EA)
    BA = bA * sq
    return (W + th / BA) / (1 + BA * th * W)


def residual_and_jac(theta, target, ws, grad=True):
    F, Jm = [], []
    for w, z0 in zip(ws, target):
        z = ztot(theta, w, grad)
        dz = z.v - z0
        F += [dz.real, dz.imag]
        if grad:
            Jm.append([a.real for a in z.d])
            Jm.append([a.imag for a in z.d])
    return F, Jm


def mid(x):
    return arb(x.mid())


CASES = [
    ("local:generic:N3a", [fmpq(4, 5), fmpq(2), fmpq(1, 2), fmpq(1, 10), fmpq(3, 10), fmpq(1, 5)],
     [fmpq(1, 2), fmpq(4), fmpq(17, 10)]),
    ("local:generic:N3b", [fmpq(4, 5), fmpq(2), fmpq(1, 2), fmpq(1, 10), fmpq(3, 10), fmpq(1, 5)],
     [fmpq(1, 5), fmpq(11, 10), fmpq(6)]),
    ("local:b=bs:N3b", [fmpq(4, 5), fmpq(1), fmpq(1, 2), fmpq(1, 10), fmpq(3, 10), fmpq(1, 5)],
     [fmpq(1, 5), fmpq(11, 10), fmpq(6)]),
]
prev = json.load(open(os.path.join(HERE, "..", "B-revision", "b4_07_certify_twin.json")))
OUT = {}
for key, tgt_q, ws_q in CASES:
    start = prev[key + "#0"]["twin_40_digits"]
    # ---------------- (1) Newton refinement at 400 bits
    ctx.prec = 400
    target = [ztot([arb(q) for q in tgt_q], arb(w), grad=False).v for w in ws_q]
    ws = [arb(w) for w in ws_q]
    x = [arb(sv) for sv in start]
    for it in range(8):
        F, Jm = residual_and_jac([mid(t) for t in x], target, ws)
        Jmat = arb_mat([[mid(a) for a in row] for row in Jm])
        step = Jmat.solve(arb_mat([[mid(f)] for f in F]))
        x = [mid(x[i] - step[i, 0]) for i in range(NP)]
    Fn, _ = residual_and_jac(x, target, ws, grad=False)
    res_newton = max(float(abs(f).mid()) for f in Fn)
    # ---------------- (2) Krawczyk at 300 bits, r = 1e-35
    ctx.prec = 300
    target = [ztot([arb(q) for q in tgt_q], arb(w), grad=False).v for w in ws_q]
    ws = [arb(w) for w in ws_q]
    xm = [mid(t) for t in x]
    rad = arb("1e-35").mid()
    Fm, Jmid = residual_and_jac(xm, target, ws)
    Y = arb_mat([[mid(a) for a in row] for row in Jmid]).inv()
    Y = arb_mat([[mid(Y[i, j]) for j in range(NP)] for i in range(NP)])
    box = [arb(xm[i], rad) for i in range(NP)]
    _, JX = residual_and_jac(box, target, ws)
    JX = arb_mat(JX)
    I6 = arb_mat([[1 if i == j else 0 for j in range(NP)] for i in range(NP)])
    dX = arb_mat([[box[i] - xm[i]] for i in range(NP)])
    K = arb_mat([[t] for t in xm]) - Y * arb_mat([[f] for f in Fm]) + (I6 - Y * JX) * dX
    dev = [(K[i, 0] - xm[i]).abs_upper() for i in range(NP)]
    inside = all(d < rad for d in dev)
    # ---------------- (3) certified side conditions
    off_delta = (box[2] - box[3]) > 0
    sep = max(float(abs(box[i] - arb(tgt_q[i])).lower()) for i in range(NP))
    z1 = ztot(box, arb(fmpq(9, 10)), grad=False).v
    z0 = ztot([arb(q) for q in tgt_q], arb(fmpq(9, 10)), grad=False).v
    mism = abs(z1 - z0) / abs(z0)
    detJX = JX.det()
    OUT[key] = {"twin_35_digits": [xm[i].str(36, radius=False) for i in range(NP)],
                "newton_residual_400bit": res_newton,
                "box_radius": 1e-35, "krawczyk_inside": bool(inside),
                "K_minus_mid_upper": [float(d) for d in dev],
                "det_J_over_box": detJX.str(10),
                "tau_R_minus_tau_l_positive_on_box": bool(off_delta),
                "max_component_separation_from_target_lower_bound": sep,
                "relative_mismatch_at_omega_9_10": mism.str(10),
                "mismatch_certified_nonzero": bool(mism > 0),
                "agreement_with_b4_07_40_digits": max(abs(float((xm[i] - arb(start[i])).mid())) for i in range(NP))}
    print(key, json.dumps(OUT[key], indent=1), flush=True)

with open(os.path.join(HERE, "a05_tw3_recertify.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
