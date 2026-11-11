"""Audit: function-level rank of the 5-parameter model on Delta ∩ {R = 0, b != b_s} (Theorem 2', third row).

Independent implementation: Z from the tanh form with load L = R + 1/(b_s sqrt p), all derivatives by
mpmath numerical differentiation (two-sided; the formula is analytic in R and tau near R = 0, tau = 0 for
p = i omega). Singular values of the stacked real Jacobian over 41 log-spaced frequencies, columns
unscaled. Predicted kernel: (0,0,1,1,0) and (b b_s^2/(b^2 - b_s^2), 0, 0, 0, 1).
Run:  python 04_resonance_R0_stratum.py > 04_resonance_R0_stratum.out
"""
import mpmath as mp

mp.mp.dps = 50
BS = mp.mpf(1)


def Z(p, xi, b, tR, tl, R):
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); S = BS*mp.sqrt(p)
    L = R + 1/S
    t = mp.tanh(A)
    return (L + t/B)/(1 + B*L*t)


def jac(p, th):
    out = []
    for k in range(5):
        def f(x, k=k):
            a = list(th); a[k] = x
            return Z(p, *a)
        out.append(mp.diff(f, th[k]))
    return out


omegas = [mp.mpf(10)**(mp.mpf(k)/10) for k in range(-20, 21)]
for label, th in [("Delta, R = 0, b = 2", [mp.mpf(1), mp.mpf(2), mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf(0)]),
                  ("Delta, R = 0, b = 0.4", [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('0.2'), mp.mpf('0.2'), mp.mpf(0)]),
                  ("Delta, R = 0.3, b = 2 (control: rank 4)", [mp.mpf(1), mp.mpf(2), mp.mpf('0.5'), mp.mpf('0.5'), mp.mpf('0.3')]),
                  ("off Delta, R = 0, b = 2 (control: rank 5)", [mp.mpf(1), mp.mpf(2), mp.mpf('0.7'), mp.mpf('0.3'), mp.mpf(0)])]:
    rows = []
    for w in omegas:
        J = jac(mp.mpc(0, w), th)
        rows.append([mp.re(x) for x in J]); rows.append([mp.im(x) for x in J])
    M = mp.matrix(rows)
    U, Sv, V = mp.svd_r(M)
    sv = [Sv[i]/Sv[0] for i in range(5)]
    b = th[1]
    k = b*BS**2/(b**2 - BS**2)
    res = []
    for kv in ([0, 0, 1, 1, 0], [k, 0, 0, 0, 1]):
        v = mp.matrix(kv)
        res.append(mp.norm(M*v)/(mp.norm(v)*Sv[0]))
    print(f"[{label}] singular values / max: {[mp.nstr(x, 3) for x in sv]}")
    print(f"     |J k|/(|k| s_max) for (0,0,1,1,0): {mp.nstr(res[0], 3)};  for (b b_s^2/(b^2-b_s^2),0,0,0,1): {mp.nstr(res[1], 3)}")
