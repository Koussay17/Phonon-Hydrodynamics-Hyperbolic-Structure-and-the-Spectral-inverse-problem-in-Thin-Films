"""Workstream B: numerical rank of the parameter Jacobian on each stratum (mpmath, 60 digits).

Two independent derivative computations (analytic log-derivative formulas vs numerical
differentiation of the quadrupole product) are compared first. Then the singular values of the
stacked real Jacobian over 60 log-spaced frequencies are printed for every stratum, with the
numerical kernel compared with the predicted kernel.
Run:  python 02_numerical_rank.py > 02_numerical_rank.out
"""
import mpmath as mp
from model import jac_analytic, jac_numeric, real_jacobian, Z_quad, Z_refl

mp.mp.dps = 60
bs = mp.mpf(1)
omegas = [mp.mpf(10)**(mp.mpf(k)/10) for k in range(-30, 31)]   # 1e-3 .. 1e3

print("=== derivative cross-check (analytic vs numerical differentiation) ===")
worst = 0
for th, R in [((1, 2, 0.7, 0.3), None), ((1, 1, 0.7, 0.3), None), ((1.3, 0.6, 0, 0.5), None),
              ((0.8, 2, 0.4, 0), None), ((1, 2, 0.7, 0.3), 0.3), ((1.1, 0.7, 0, 0.6), 0.25)]:
    th = [mp.mpf(x) for x in th]
    for w in [mp.mpf('0.01'), mp.mpf('0.5'), mp.mpf(3), mp.mpf(40)]:
        p = mp.mpc(0, w)
        z1, z2 = Z_quad(p, *th, bs, 0 if R is None else R), Z_refl(p, *th, bs, 0 if R is None else R)
        Ja = jac_analytic(p, *th, bs, R=R)
        Jn = jac_numeric(p, *th, bs, R=None if R is None else mp.mpf(R))
        err = max(abs(a - n)/max(abs(n), mp.mpf(10)**-40) for a, n in zip(Ja, Jn))
        worst = max(worst, err, abs(z1 - z2)/abs(z1))
print("max relative discrepancy (Jacobian and Z forms):", mp.nstr(worst, 5))


def report(label, th, R=None, predicted_kernel=None):
    th = [mp.mpf(x) for x in th]
    Rm = None if R is None else mp.mpf(R)
    J = real_jacobian(omegas, th, bs, R=Rm)
    U, Sv, V = mp.svd_r(J)
    sv = [Sv[i] for i in range(len(Sv))]
    rel = [x/sv[0] for x in sv]
    n = len(sv)
    print(f"\n[{label}] theta = {[mp.nstr(x, 4) for x in th]}" + ("" if R is None else f", R = {mp.nstr(Rm, 4)}"))
    print("   singular values / max:", [mp.nstr(x, 3) for x in rel])
    rank = sum(1 for x in rel if x > mp.mpf(10)**-30)
    print("   numerical rank (threshold 1e-30):", rank, "of", n)
    if predicted_kernel:
        # columns were scaled by parameter values (when > 0); convert predicted kernel to scaled coordinates
        thR = th + ([] if R is None else [Rm])
        ker = []
        for kv in predicted_kernel:
            kv = [mp.mpf(x) for x in kv]
            ker.append([kv[i]/(thR[i] if thR[i] != 0 else 1) for i in range(n)])
        for kv in ker:
            res = J*mp.matrix(kv)
            print("   |J k|/|k| for predicted kernel vector", [mp.nstr(x, 4) for x in kv], ":",
                  mp.nstr(mp.norm(res)/mp.norm(mp.matrix(kv))/sv[0], 3))
    return rank


b0 = 2
print("\n=== four-parameter model ===")
report("generic interior", (1, 2, 0.7, 0.3))
report("b = b_s, tau_R != tau_l", (1, 1, 0.7, 0.3))
report("high-frequency matched b sqrt(tl/tR) = b_s", (1, 2, 0.8, 0.2))
report("Cattaneo boundary tau_l = 0", (1, 2, 0.7, 0))
report("pure nonlocal boundary tau_R = 0", (1, 2, 0, 0.5))
report("b = b_s and tau_l = 0", (1, 1, 0.7, 0))
report("b = b_s and tau_R = 0", (1, 1, 0, 0.5))
report("near resonance tau_l = tau_R(1 + 1e-4)", (1, 2, 0.5, 0.50005))
report("resonance tau_R = tau_l, b != b_s", (1, 2, 0.5, 0.5), predicted_kernel=[(0, 0, 1, 1)])
report("Fourier point tau_R = tau_l = 0, b != b_s", (1, 2, 0, 0), predicted_kernel=[(0, 0, 1, 1)])
report("D2: resonance and b = b_s", (1, 1, 0.5, 0.5), predicted_kernel=[(1, 0, 0, 0), (0, 0, 1, 1)])
report("D2 at the Fourier point", (1, 1, 0, 0), predicted_kernel=[(1, 0, 0, 0), (0, 0, 1, 1)])

print("\n=== five-parameter model (interface resistance R) ===")
report("generic, R > 0", (1, 2, 0.7, 0.3), R=0.3)
report("R = 0, tau_R != tau_l", (1, 2, 0.7, 0.3), R=0)
report("R = 0, b = b_s, tau_R != tau_l", (1, 1, 0.7, 0.3), R=0)
report("pure nonlocal tau_R = 0, R > 0", (1, 2, 0, 0.5), R=0.3)
report("pure nonlocal tau_R = 0, R > 0, b = b_s", (1, 1, 0, 0.5), R=0.3)
report("Cattaneo tau_l = 0, R > 0", (1, 2, 0.7, 0), R=0.3)
report("resonance, R > 0", (1, 2, 0.5, 0.5), R=0.3, predicted_kernel=[(0, 0, 1, 1, 0)])
report("resonance, R > 0, b = b_s", (1, 1, 0.5, 0.5), R=0.3, predicted_kernel=[(0, 0, 1, 1, 0)])
kx = mp.mpf(b0)*bs**2/(mp.mpf(b0)**2 - bs**2)
report("resonance, R = 0, b != b_s", (1, 2, 0.5, 0.5), R=0,
       predicted_kernel=[(0, 0, 1, 1, 0), (kx, 0, 0, 0, 1)])
report("Fourier point, R = 0, b != b_s", (1, 2, 0, 0), R=0,
       predicted_kernel=[(0, 0, 1, 1, 0), (kx, 0, 0, 0, 1)])
report("D2, R = 0", (1, 1, 0.5, 0.5), R=0, predicted_kernel=[(1, 0, 0, 0, 0), (0, 0, 1, 1, 0)])
