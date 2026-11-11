"""Workstream B: certified (ball-arithmetic) non-vanishing of finite-frequency Jacobian minors.

For the 4-parameter model with N = 2 frequencies the real 4x4 matrix
    [Re J(w1); Im J(w1); Re J(w2); Im J(w2)],  J = d log Z / d log theta (d/d tau when tau = 0),
is evaluated in Arb ball arithmetic (python-flint, 256-bit). If the determinant ball excludes 0,
the determinant is proved non-zero at that point; by real-analyticity in (theta, w1, w2) it is then
non-zero on an open dense, full-measure subset of the connected stratum containing the point
(Lemma F in proofs/B-analytic.md). For the 5-parameter model (interface resistance) the same is
done with N = 3 frequencies and a 5x5 minor of the 6x5 matrix.
Units: time in xi_1^2-free units, effusivity in units of b_s (b_s = 1); see the .md file.
Run:  python 04_certified_finite_frequency.py > 04_certified_finite_frequency.out
"""
from flint import arb, acb, arb_mat, ctx

ctx.prec = 256
bs = arb(1)


def jac_rel(w, xi, b, tR, tl, R=None):
    p = acb(0, w)
    one = acb(1)
    g = (one + tR*p)/(one + tl*p)
    A = xi*(p*g).sqrt()
    B = b*(p/g).sqrt()
    S = bs*p.sqrt()
    L = (arb(0) if R is None else R) + one/S
    r = (B*L - one)/(B*L + one)
    E = (-2*A).exp()
    den = one - r*r*E*E
    ulog = -4*A*r*E/den                              # A dlogZ/dA
    vlog = -one + 2*E*(2*B*L/((B*L + one)**2))/den   # B dlogZ/dB
    hR, hl = p/(one + tR*p), p/(one + tl*p)
    W = ulog - vlog
    th = [xi, b, tR, tl] + ([] if R is None else [R])
    J = [ulog/xi, vlog/b, W*hR/2, -W*hl/2]
    if R is not None:
        J.append(2*E*(2*B/((B*L + one)**2))/den)
    # log-scaling where the parameter is positive
    out = []
    for k, Jk in enumerate(J):
        sc = th[k] if not (th[k] == 0) else arb(1)
        out.append(Jk*sc)
    return out


def real_matrix(ws, th, R=None):
    rows = []
    for w in ws:
        J = jac_rel(arb(w), *th, R=R)
        rows.append([x.real for x in J])
        rows.append([x.imag for x in J])
    return rows


def certify(label, th, ws, R=None, drop_row=None):
    th = [arb(x) for x in th]
    Rb = None if R is None else arb(R)
    rows = real_matrix(ws, th, Rb)
    if drop_row is not None:
        rows = [rw for i, rw in enumerate(rows) if i != drop_row]
    M = arb_mat(rows)
    d = M.det()
    ok = not d.contains(0)
    print(f"[{'CERTIFIED NONZERO' if ok else 'NOT CERTIFIED'}] {label}: theta={[float(x.mid()) for x in th]}"
          + ("" if R is None else f", R={float(Rb.mid())}") + f", omegas={ws}" + (f", minor without row {drop_row}" if drop_row is not None else ""))
    print(f"      det midpoint = {float(d.mid()):.12e}, rigorous radius = {float(d.rad()):.3e}")
    return ok


print("=== 4-parameter model, N = 2 frequencies ===")
pts = [("interior", (1, 2, 0.7, 0.3)),
       ("b = b_s, off resonance", (1, 1, 0.7, 0.3)),
       ("Cattaneo tau_l = 0", (1, 2, 0.7, 0)),
       ("pure nonlocal tau_R = 0", (1, 2, 0, 0.5)),
       ("b = b_s, tau_l = 0", (1, 1, 0.7, 0)),
       ("b = b_s, tau_R = 0", (1, 1, 0, 0.5)),
       ("AlN-like 500 nm, units xi1^2 = 1.88 ns, b_s = 10.3e3", (1, 2.70, 0.532, 0.266))]
for label, th in pts:
    certify(label, th, [0.5, 4.0])
# the AlN-like point at two physical frequencies (10 MHz, 100 MHz -> 0.118, 1.18 in these units)
certify("AlN-like 500 nm at 10 and 100 MHz", (1, 2.70, 0.532, 0.266), [0.118, 1.18])
certify("AlN-like Cattaneo 500 nm at 10 and 100 MHz", (1, 2.70, 0.532, 0), [0.118, 1.18])

print("\n=== on the resonance (expected: every 4x4 determinant vanishes; ball must contain 0) ===")
certify("resonance tau_R = tau_l", (1, 2, 0.5, 0.5), [0.5, 4.0])

print("\n=== 5-parameter model (interface resistance), N = 3 frequencies, 5x5 minors ===")
pts5 = [("interior, R > 0", (1, 2, 0.7, 0.3), 0.3),
        ("R = 0 boundary", (1, 2, 0.7, 0.3), 0),
        ("b = b_s, R > 0", (1, 1, 0.7, 0.3), 0.3),
        ("tau_R = 0, R > 0", (1, 2, 0, 0.5), 0.3),
        ("tau_l = 0, R > 0", (1, 2, 0.7, 0), 0.3),
        ("AlN-like, R = 1e-8 m^2K/W (2.4 in units)", (1, 2.70, 0.532, 0.266), 2.4)]
for label, th, R in pts5:
    certify(label, th, [0.3, 1.5, 6.0], R=R, drop_row=5)
