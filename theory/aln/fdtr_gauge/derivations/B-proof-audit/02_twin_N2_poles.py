"""Audit of proofs/B-analytic.md: independent checks of Proposition TW, Theorem N2 and of the
absence of poles of Z on the closed right half-plane (used implicitly by Lemma F / Theorem 2).

(A) Proposition TW. Z from the cosh/sinh closed form Z = (B cosh A + S sinh A)/(B (B sinh A + S cosh A));
    Jacobian from the quotient rule on that form (not the log-derivative formulas of 09_*.py), validated
    against numerical differentiation of the quadrupole matrix product. Krawczyk test (python-flint, 300 bits)
    on a box of radius 1e-15 centred on the 20-digit values PRINTED in the proof.
(B) Theorem N2. det J_2(theta; sigma^2, lam^2 sigma^2) with derivatives by mpmath numerical differentiation
    of Z (tanh form), compared with sigma^|K| t_K(lam) G_K, G_K from 01_riccati_moments.py, for lam = 3 and 1/2.
(C) |r E| < 1 on Re p >= 0, p != 0 (random sampling incl. boundary strata and R > 0).
Run:  python 02_twin_N2_poles.py > 02_twin_N2_poles.out
"""
import random
import mpmath as mp
import sympy as sp
from flint import arb, acb, arb_mat, ctx

random.seed(11112026)

# ------------------------------------------------------------------------------------------- (A)
print("=== (A) Proposition TW: independent Krawczyk certificate ===")
ctx.prec = 300
ONE = acb(1)


def parts_b(w, th):
    xi, b, tR, tl = th
    p = acb(0, w)
    g = (ONE + tR*p)/(ONE + tl*p)
    A = xi*(p*g).sqrt()
    B = b*(p/g).sqrt()
    S = p.sqrt()                                        # b_s = 1
    return p, g, A, B, S


def Z_b(w, th):
    p, g, A, B, S = parts_b(w, th)
    ch, sh = A.cosh(), A.sinh()
    return (B*ch + S*sh)/(B*(B*sh + S*ch))


def J_b(w, th):
    """dZ/d(xi, b, tR, tl) by the quotient rule on N/D, N = B ch + S sh, D = B (B sh + S ch)."""
    xi, b, tR, tl = th
    p, g, A, B, S = parts_b(w, th)
    ch, sh = A.cosh(), A.sinh()
    N = B*ch + S*sh
    D = B*(B*sh + S*ch)
    NA = B*sh + S*ch
    NB = ch
    DA = B*(B*ch + S*sh)
    DB = 2*B*sh + S*ch
    ZA = (NA*D - N*DA)/(D*D)
    ZB = (NB*D - N*DB)/(D*D)
    hR, hl = p/(ONE + tR*p), p/(ONE + tl*p)
    dA = [A/xi, acb(0), A*hR/2, -A*hl/2]
    dB = [acb(0), B/b, -B*hR/2, B*hl/2]
    return [ZA*dA[k] + ZB*dB[k] for k in range(4)]


# validation of the quotient-rule Jacobian against numerical differentiation of the matrix product (mpmath)
mp.mp.dps = 50


def Z_mat(w, xi, b, tR, tl, bs=1):
    p = mp.mpc(0, w)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); S = bs*mp.sqrt(p)
    M = mp.matrix([[mp.cosh(A), mp.sinh(A)/B], [B*mp.sinh(A), mp.cosh(A)]])
    v = M*mp.matrix([1/S, 1])
    return v[0]/v[1]


thp_str = ['0.38213304403876405621', '0.99455364692717951823', '1.58320216151666891086', '0.77319048590574559013']
thp_mp = [mp.mpf(x) for x in thp_str]
worst = 0
for w in [mp.mpf('0.5'), mp.mpf(4), mp.mpf('1.7')]:
    Jq = J_b(arb(str(w)), [arb(x) for x in thp_str])
    for kk in range(4):
        def f(x, kk=kk):
            a = list(thp_mp); a[kk] = x
            return Z_mat(w, *a)
        jn = mp.diff(f, thp_mp[kk])
        ja = mp.mpc(float(Jq[kk].real.mid()), float(Jq[kk].imag.mid()))
        worst = max(worst, abs(ja - jn)/abs(jn))
print("quotient-rule Jacobian vs numerical differentiation of the matrix product: max rel. diff =",
      mp.nstr(worst, 3), "(double-precision transfer of the ball midpoints)")

W = [arb(1)/2, arb(4)]
TH0 = [arb(1), arb(1), arb(7)/10, arb(3)/10]
DATA = [Z_b(w, TH0) for w in W]


def f_b(th):
    out = []
    for w, d in zip(W, DATA):
        z = Z_b(w, th) - d
        out += [z.real, z.imag]
    return out


def Jm_b(th):
    rows = []
    for w in W:
        J = J_b(w, th)
        rows.append([x.real for x in J]); rows.append([x.imag for x in J])
    return arb_mat(rows)


xt = [arb(x) for x in thp_str]                         # midpoint = printed digits (ball of tiny radius)
rad = arb('1e-15')
X = [arb(x.mid(), rad) for x in xt]
fx = arb_mat([[v] for v in f_b(xt)])
Jmid = Jm_b([arb(x.mid()) for x in xt])
Ym = arb_mat([[Jmid[i, j].mid() for j in range(4)] for i in range(4)]).inv()
Ym = arb_mat([[Ym[i, j].mid() for j in range(4)] for i in range(4)])
JX = Jm_b(X)
I4 = arb_mat([[1 if i == j else 0 for j in range(4)] for i in range(4)])
dX = arb_mat([[X[i] - xt[i]] for i in range(4)])
K = arb_mat([[xt[i]] for i in range(4)]) - Ym*fx + (I4 - Ym*JX)*dX
inside = all(X[i].contains_interior(K[i, 0]) for i in range(4))
print("Krawczyk: K(X) strictly inside X (radius 1e-15 around the printed theta'):", inside)
for i in range(4):
    print(f"   comp {i}: |K mid - X mid| = {float(abs((K[i, 0] - xt[i]).mid())):.3e}, rad K = {float(K[i, 0].rad()):.3e}")
detJ = Jmid.det()
print("det J_2(theta') =", detJ.str(10), "(regular root)")
print("|f(printed theta')| components:", [v.str(5) for v in f_b(xt)])
# third frequency difference with the certified box
z3 = Z_b(arb(17)/10, X) - Z_b(arb(17)/10, TH0)
print("Z(1.7 i; box) - Z(1.7 i; theta*) =", z3.str(12))

# ------------------------------------------------------------------------------------------- (B)
print("\n=== (B) Theorem N2: leading term by an independent implementation ===")
mp.mp.dps = 60


def Z_tanh(p, xi, b, tR, tl, bs=1):
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); S = bs*mp.sqrt(p)
    t = mp.tanh(A)
    return (B + S*t)/(B*(S + B*t))


def J2_det(th, sig, lamv):
    rows = []
    for sg in (sig, lamv*sig):
        p = mp.mpc(0, sg**2)
        J = []
        for kk in range(4):
            def f(x, kk=kk):
                a = list(th); a[kk] = x
                return Z_tanh(p, *a)
            J.append(mp.diff(f, th[kk]))
        rows.append([mp.re(x) for x in J]); rows.append([mp.im(x) for x in J])
    return mp.det(mp.matrix(rows))


xi_, b_, bs_, tR_, tl_, lam_ = sp.symbols('xi b b_s tau_R tau_l lam', positive=True)
G0124 = -2*xi_**4*(b_ - bs_)*(b_ + bs_)*(tR_ - tl_)/(b_**4*bs_**3)
G0234 = 2*xi_**4*(tR_ - tl_)*(2*b_**4*xi_**2 - 5*b_**2*bs_**2*xi_**2 - 3*bs_**4*tR_ + 3*bs_**4*tl_ + 3*bs_**4*xi_**2)/(3*b_**4*bs_**5)
t0124 = sp.sqrt(2)*lam_*(lam_ - 1)**2*(lam_ + 1)*(lam_**2 + 1)/2
t0234 = sp.sqrt(2)*lam_**2*(lam_ - 1)**2*(lam_ + 1)*(lam_**2 + 1)/2
for lamv in [mp.mpf(3), mp.mpf(1)/2]:
    for label, th in [("interior", (1, 2, '0.7', '0.3')), ("b = b_s", (1, 1, '0.7', '0.3')),
                      ("tau_l = 0", ('1.3', '0.6', '0.4', 0)), ("tau_R = 0", ('0.8', '1.7', 0, '0.9')),
                      ("b = b_s, tau_R = 0", ('1.2', 1, 0, '0.4')), ("tau_l > tau_R", ('0.8', '0.5', '0.2', '0.9'))]:
        th = [mp.mpf(x) for x in th]
        sub = {xi_: th[0], b_: th[1], bs_: 1, tR_: th[2], tl_: th[3], lam_: lamv}
        if th[1] != 1:
            m0, lead = 7, sp.N((t0124*G0124).subs(sub), 70)
        else:
            m0, lead = 9, sp.N((t0234*G0234).subs(sub), 70)
        lead = mp.mpf(str(lead))
        rat = [J2_det(th, sg, lamv)/(sg**m0*lead) for sg in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'))]
        print(f"  lam = {mp.nstr(lamv, 3)}, {label}: m0 = {m0}, ratio det/(sigma^m0 t G) at sigma = 1e-2,1e-3,1e-4:",
              [mp.nstr(x, 8) for x in rat])

# ------------------------------------------------------------------------------------------- (C)
print("\n=== (C) |r E| < 1 on Re p >= 0, p != 0 ===")
mp.mp.dps = 30
worst = mp.mpf(0)
worst_re_BL = mp.mpf('inf')
n = 0
for _ in range(20000):
    xi = mp.mpf(10)**random.uniform(-2, 2); b = mp.mpf(10)**random.uniform(-2, 2)
    tR = 0 if random.random() < 0.15 else mp.mpf(10)**random.uniform(-3, 3)
    tl = 0 if random.random() < 0.15 else mp.mpf(10)**random.uniform(-3, 3)
    R = 0 if random.random() < 0.3 else mp.mpf(10)**random.uniform(-3, 2)
    bs = mp.mpf(10)**random.uniform(-2, 2)
    phi = random.choice([mp.pi/2, -mp.pi/2, random.uniform(-1, 1)*mp.pi/2])
    p = mp.mpf(10)**random.uniform(-6, 6)*mp.expj(phi)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); S = bs*mp.sqrt(p)
    L = R + 1/S
    r = (B*L - 1)/(B*L + 1)
    E = mp.exp(-2*A)
    worst = max(worst, abs(r*E))
    worst_re_BL = min(worst_re_BL, mp.re(B*L)/abs(B*L))
    n += 1
print(f"  {n} random points: max |r E| = {mp.nstr(worst, 8)};  min Re(BL)/|BL| = {mp.nstr(worst_re_BL, 5)} (> 0 => |r| < 1)")
