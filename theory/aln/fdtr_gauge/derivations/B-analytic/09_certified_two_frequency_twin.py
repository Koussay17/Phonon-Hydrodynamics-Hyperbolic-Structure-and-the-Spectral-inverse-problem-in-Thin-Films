"""Workstream B: certified example that two frequencies do not give global uniqueness.

Target theta* = (1, 1, 0.7, 0.3) (b = b_s), frequencies (0.5, 4.0) (units xi_1^2 = 1, b_s = 1).
05_injectivity_search found a second exact root theta' ~ (0.3821, 0.9946, 1.5832, 0.7732).
Here: (1) refine theta' to 60 digits (mpmath Newton); (2) prove existence and uniqueness of a root of
f(x) = [Re, Im](Z(i w_k; x) - Z(i w_k; theta*)), k = 1, 2, in a small box around theta' by the Krawczyk
test in Arb ball arithmetic; (3) certify that the responses at a third frequency differ (consistent
with Theorem 1: the full responses are distinct).
Run:  python 09_certified_two_frequency_twin.py > 09_certified_two_frequency_twin.out
"""
import mpmath as mp
from flint import arb, acb, arb_mat, ctx

mp.mp.dps = 70          # set before any constant is created (0.7 and 0.3 must be the exact decimals)
W = [mp.mpf('0.5'), mp.mpf(4)]
TH0 = [mp.mpf(1), mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.3')]
GUESS = [mp.mpf('0.38213304'), mp.mpf('0.99455365'), mp.mpf('1.58320216'), mp.mpf('0.77319049')]


def Z_mp(w, xi, b, tR, tl, bs=1):
    p = mp.mpc(0, w)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); S = bs*mp.sqrt(p)
    r = (B - S)/(B + S); E = mp.exp(-2*A)
    return (1/B)*(1 + r*E)/(1 - r*E)


data = [Z_mp(w, *TH0) for w in W]


def f_mp(*x):
    out = []
    for w, d in zip(W, data):
        z = Z_mp(w, *x) - d
        out += [mp.re(z), mp.im(z)]
    return out


root = mp.findroot(f_mp, GUESS, tol=mp.mpf(10)**-60, maxsteps=100)
root = [root[i] for i in range(4)]
print("refined second root theta' =", [mp.nstr(x, 30) for x in root])
print("residual |f(theta')| =", mp.nstr(max(abs(v) for v in f_mp(*root)), 5))

# ---------------- Krawczyk test in ball arithmetic ----------------
ctx.prec = 300


def Zb(w, xi, b, tR, tl):
    p = acb(0, w); one = acb(1)
    g = (one + tR*p)/(one + tl*p)
    A = xi*(p*g).sqrt(); B = b*(p/g).sqrt(); S = p.sqrt()
    r = (B - S)/(B + S); E = (-2*A).exp()
    return (one/B)*(one + r*E)/(one - r*E), A, B, S, r, E


def Jb(w, xi, b, tR, tl):
    """Unscaled dZ/d(xi, b, tR, tl) from the log-derivative formulas."""
    Z, A, B, S, r, E = Zb(w, xi, b, tR, tl)
    one = acb(1)
    den = one - r*r*E*E
    u = Z*(-4*A*r*E/den)
    v = Z*(-one + 4*B*S*E/((B + S)**2*den))
    p = acb(0, w)
    hR, hl = p/(one + tR*p), p/(one + tl*p)
    Wd = u - v
    return [u/xi, v/b, Wd*hR/2, -Wd*hl/2]


def to_arb(x):
    return arb(mp.nstr(x, 85))


th0b = [arb(1), arb(1), arb(7)/10, arb(3)/10]
datab = [Zb(arb(w), *th0b)[0] for w in (arb(1)/2, arb(4))]
wb = [arb(1)/2, arb(4)]


def fb(x):
    out = []
    for w, d in zip(wb, datab):
        z = Zb(w, *x)[0] - d
        out += [z.real, z.imag]
    return out


def Jmat(x):
    rows = []
    for w in wb:
        J = Jb(w, *x)
        rows.append([c.real for c in J])
        rows.append([c.imag for c in J])
    return arb_mat(rows)


xt = [to_arb(x) for x in root]                      # exact midpoint (85 digits)
rad = arb('1e-40')
X = [arb(xm.mid(), rad) for xm in xt]              # box
fx = arb_mat([[v] for v in fb(xt)])
Jmid = Jmat([arb(x.mid()) for x in xt])
# approximate inverse Y (midpoints only)
Jm_mid = arb_mat([[Jmid[i, j].mid() for j in range(4)] for i in range(4)])
Y = Jm_mid.inv()
Y = arb_mat([[Y[i, j].mid() for j in range(4)] for i in range(4)])
JX = Jmat(X)
I4 = arb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
dX = arb_mat([[X[i] - xt[i]] for i in range(4)])
K = arb_mat([[xt[i]] for i in range(4)]) - Y*fx + (I4 - Y*JX)*dX
inside = all(X[i].contains_interior(K[i, 0]) for i in range(4))
print("Krawczyk box radius 1e-40; K(X) strictly inside X:", inside)
for i in range(4):
    print(f"   component {i}: K midpoint - box midpoint = {float((K[i, 0] - xt[i]).mid()):.3e}, K radius = {float(K[i, 0].rad()):.3e}")
if inside:
    print("=> CERTIFIED: a unique root of the two-frequency equations lies in the box around theta',"
          " and theta' != theta* (the boxes are disjoint).")

# third frequency: responses differ
w3 = arb(17)/10
z0 = Zb(w3, *th0b)[0]
z1 = Zb(w3, *X)[0]
diff = z1 - z0
print("Z(i*1.7; theta') - Z(i*1.7; theta*) =", diff.str(15),
      "; excludes 0:", not diff.real.contains(0) or not diff.imag.contains(0))
