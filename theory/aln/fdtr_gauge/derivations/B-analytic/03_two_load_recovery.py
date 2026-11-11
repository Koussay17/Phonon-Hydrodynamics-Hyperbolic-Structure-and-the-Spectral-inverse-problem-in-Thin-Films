"""Workstream B: numerical check of the two-load (substrate second sheet) recovery.

(a) Monodromy: continuing all square roots continuously around closed loops in the p-plane,
    the film quadrupole M returns to itself (single-valued), while S = b_s sqrt(p) changes sign
    after a loop around p = 0. Loops around -1/tau_R only, around 0 only, and around both are tested.
(b) Recovery: from F(s) and F(-s) alone, the 2x3 linear system gives (C, X, Y) up to a factor,
    fixed by C^2 - X Y = 1 and continuity (C -> 1 as p -> 0); compared with the true entries.
(c) Parameters: from B^2 = Y/X and z = (arccosh C)^2 the parameters (xi, b, tau_R, tau_l) are recovered.
Run:  python 03_two_load_recovery.py > 03_two_load_recovery.out
"""
import mpmath as mp

mp.mp.dps = 50
bs = mp.mpf(1)


def film_entire(p, xi, b, tR, tl):
    """C, X, Y as single-valued functions of p (entire functions of z = xi^2 p g)."""
    g = (1 + tR*p)/(1 + tl*p)
    z = xi**2*p*g
    # cosh(sqrt z) and sinh(sqrt z)/sqrt z via series-safe forms
    sq = mp.sqrt(z)
    C = mp.cosh(sq)
    Sc = mp.sinh(sq)/sq if z != 0 else mp.mpf(1)
    return C, xi/b*g*Sc, xi*b*p*Sc


def F_of_s(s, th):
    xi, b, tR, tl = th
    p = s**2
    C, X, Y = film_entire(p, xi, b, tR, tl)
    S = bs*s
    return (C + X*S)/(Y + C*S)


def monodromy(th, center, radius, n=4000):
    """Continue sqrt(p), sqrt(p g), sqrt(p/g) continuously along a circle; return relative change of M and sign of S."""
    xi, b, tR, tl = th
    g = lambda p: (1 + tR*p)/(1 + tl*p)
    p0 = center + radius
    sp_, sa, sb = mp.sqrt(p0), mp.sqrt(p0*g(p0)), mp.sqrt(p0/g(p0))

    def M(sa, sb):
        A, B = xi*sa, b*sb
        return mp.matrix([[mp.cosh(A), mp.sinh(A)/B], [B*mp.sinh(A), mp.cosh(A)]])
    M0 = M(sa, sb)
    S0 = sp_
    for k in range(1, n + 1):
        p = center + radius*mp.exp(2j*mp.pi*k/n)
        cands = [(mp.sqrt(p), sp_), (mp.sqrt(p*g(p)), sa), (mp.sqrt(p/g(p)), sb)]
        new = []
        for val, prev in cands:
            new.append(val if abs(val - prev) < abs(val + prev) else -val)
        sp_, sa, sb = new
    M1 = M(sa, sb)
    rel = mp.norm(M1 - M0)/mp.norm(M0)
    return rel, mp.nstr(sp_/S0, 5), mp.nstr(sa/mp.sqrt(p0*g(p0)), 5), mp.nstr(sb/mp.sqrt(p0/g(p0)), 5)


th = [mp.mpf(1), mp.mpf(2), mp.mpf('0.7'), mp.mpf('0.3')]   # -1/tR = -1.43, -1/tl = -3.33
print("=== (a) monodromy of the film quadrupole and of S ===")
for label, c, rad in [("around 0 only", mp.mpf(0), mp.mpf('0.5')),
                      ("around 0 and -1/tau_R", mp.mpf(0), mp.mpf(2)),
                      ("around -1/tau_R only", mp.mpf(-1)/th[2], mp.mpf('0.3'))]:
    rel, sS, sA, sB = monodromy(th, c, rad)
    print(f"  loop {label}: |M_end - M_start|/|M| = {mp.nstr(rel, 3)};  S_end/S_start = {sS};"
          f"  sqrt(pg)_end/start = {sA};  sqrt(p/g)_end/start = {sB}")

print("\n=== (b) recovery of (C, X, Y) from F(s), F(-s) ===")


def recover(s, th):
    Fp, Fm = F_of_s(s, th), F_of_s(-s, th)
    S = bs*s
    L = mp.matrix([[Fp*S - 1, -S, Fp], [-Fm*S - 1, S, Fm]])
    # kernel of the 2x3 matrix: cross product of the rows
    r1, r2 = [L[0, j] for j in range(3)], [L[1, j] for j in range(3)]
    k = [r1[1]*r2[2] - r1[2]*r2[1], r1[2]*r2[0] - r1[0]*r2[2], r1[0]*r2[1] - r1[1]*r2[0]]
    lam2 = 1/(k[0]**2 - k[1]*k[2])          # C^2 - X Y = 1
    lam = mp.sqrt(lam2)
    C, X, Y = [lam*x for x in k]
    if mp.re(C) < 0:                          # sign fixed by continuity from C(0) = 1 (small |p| here)
        C, X, Y = -C, -X, -Y
    return C, X, Y


for thv in [th, [mp.mpf(1), mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.3')],      # b = b_s, off resonance
            [mp.mpf('1.3'), mp.mpf('0.6'), mp.mpf(0), mp.mpf('0.5')],         # tau_R = 0
            [mp.mpf('0.8'), mp.mpf(2), mp.mpf('0.4'), mp.mpf(0)]]:            # tau_l = 0
    worst = 0
    for s in [mp.mpc('0.2', '0.2'), mp.mpc('0.3', '-0.1'), mp.mpc('0.05', '0.4'), mp.mpc('0.4', '0.25')]:
        C, X, Y = recover(s, thv)
        Ct, Xt, Yt = film_entire(s**2, *thv)
        worst = max(worst, abs(C - Ct)/abs(Ct), abs(X - Xt)/abs(Xt), abs(Y - Yt)/abs(Yt))
    print(f"  theta = {[mp.nstr(x, 3) for x in thv]}: max relative error of recovered (C,X,Y) = {mp.nstr(worst, 3)}")

print("\n=== (c) parameters from the recovered quadrupole ===")
for thv in [th, [mp.mpf(1), mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.3')], [mp.mpf('1.3'), mp.mpf('0.6'), mp.mpf(0), mp.mpf('0.5')]]:
    est = []
    for (s1, s2) in [(mp.mpc('0.2', '0.2'), mp.mpc('0.35', '0.1'))]:
        out = []
        for s in (s1, s2):
            p = s**2
            C, X, Y = recover(s, thv)
            B2 = Y/X
            z = mp.acosh(C)**2
            kap = mp.sqrt(z*B2)/p              # xi b, up to the sign of the square root
            kap = kap if mp.re(kap) > 0 else -kap   # sign fixed by kappa = xi b > 0 (branch-cut safe)
            rho = mp.sqrt(z/B2)                # (xi/b) g
            rho = rho if mp.re(rho) > 0 else -rho
            out.append((p, kap, rho))
        # g(p) = rho(p)/(xi/b), with xi/b = lim_{p->0} rho; use g's Mobius form: rho(p)(1 + tl p) = (xi/b)(1 + tR p)
        (p1, k1, r1), (p2, k2, r2) = out
        # unknowns q = xi/b, a = q tR, c = tl :  r1 + r1 c p1 = q + a p1 ;  r2 + r2 c p2 = q + a p2 ; need one more: use kappa consistency
        # third point
        s3 = mp.mpc('0.1', '0.3'); p3 = s3**2; C3, X3, Y3 = recover(s3, thv); z3 = mp.acosh(C3)**2; r3 = mp.sqrt(z3/(Y3/X3)); r3 = r3 if mp.re(r3) > 0 else -r3
        Mx = mp.matrix([[1, p1, -r1*p1], [1, p2, -r2*p2], [1, p3, -r3*p3]])
        q, a, c = mp.lu_solve(Mx, mp.matrix([r1, r2, r3]))
        xi_ = mp.sqrt(k1*q); b_ = mp.sqrt(k1/q)
        est = [xi_, b_, a/q, c]
    print(f"  true {[mp.nstr(x, 6) for x in thv]}  recovered {[mp.nstr(mp.re(x), 12) for x in est]}")
