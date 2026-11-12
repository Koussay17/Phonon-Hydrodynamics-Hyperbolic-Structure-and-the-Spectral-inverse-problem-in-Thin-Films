"""cx_10: poles of the response in the closed right half-plane Re p >= 0 (instability / non-passivity)
for the face conventions tau_n = eps tau_l (bare film read-out T(0); stack with the same coefficient at
both film faces, as in the twin).

In s = sqrt(p) the right half-plane is the sector |arg s| < pi/4, where Re A > 0 and cosh A != 0.
Poles of Z are the zeros of the analytic functions
   bare : D(s) = B L tanh A + N
   stack: D(s) = b_A s tanh(A_A) (R1 (B L2 tanh A + N) + N0 (L2 + N tanh A / B)) + (B L2 tanh A + N)
(obtained by clearing the film and transducer denominators and dividing by cosh A cosh A_A).
Zeros are counted by the argument principle on the boundary of {rho1 < |s| < rho2, |arg s| < pi/4 - d}
and located by Newton from grid minima; every reported root is re-verified in mpmath.
"""
import json
import sys
import numpy as np
import mpmath as mp

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"


def D_np(s, prm, kind, xiA=0.2, bA=3.0, bs=1.0):
    xi, b, tR, tl, tn, to, R1, R2 = prm
    p = s * s
    g = (1 + tR * p) / (1 + tl * p)
    sg = np.sqrt(g)
    A = xi * s * sg
    B = b * s / sg
    N = 1 + tn * p
    N0 = 1 + to * p
    L = R2 + 1 / (bs * s)
    T = np.tanh(A)
    if kind == "bare":
        return B * L * T + N
    TA = np.tanh(xiA * s)
    return bA * s * TA * (R1 * (B * L * T + N) + N0 * (L + N * T / B)) + (B * L * T + N)


def winding(prm, kind, rho1=1e-3, rho2=3e3, d=1e-3, n=40000):
    a = np.pi / 4 - d
    r = np.logspace(np.log10(rho1), np.log10(rho2), n)
    th = np.linspace(-a, a, n // 4)
    path = np.concatenate([r * np.exp(-1j * a), rho2 * np.exp(1j * th), r[::-1] * np.exp(1j * a),
                           rho1 * np.exp(1j * th[::-1])])
    with np.errstate(all="ignore"):
        f = D_np(path, prm, kind)
    ph = np.unwrap(np.angle(f))
    jumps = np.max(np.abs(np.diff(ph)))
    return (ph[-1] - ph[0]) / (2 * np.pi), jumps


def locate(prm, kind, rho1=1e-3, rho2=3e3):
    """grid minima of |D| in the sector, refined by mpmath findroot; returns roots with Re p >= 0."""
    lr = np.linspace(np.log(rho1), np.log(rho2), 500)
    ang = np.linspace(-np.pi / 4, np.pi / 4, 121)
    LR, ANG = np.meshgrid(lr, ang)
    S = np.exp(LR + 1j * ANG)
    with np.errstate(all="ignore"):
        F = np.abs(D_np(S, prm, kind)) / (1 + np.abs(S) ** 2)
    roots = []
    for i in range(1, F.shape[0] - 1):
        for j in range(1, F.shape[1] - 1):
            v = F[i, j]
            if v < F[i - 1, j] and v < F[i + 1, j] and v < F[i, j - 1] and v < F[i, j + 1]:
                s0 = S[i, j]
                try:
                    mp.mp.dps = 30
                    xi, b, tR, tl, tn, to, R1, R2 = [mp.mpf(float(x)) for x in prm]

                    def fm(s):
                        p = s * s
                        g = (1 + tR * p) / (1 + tl * p)
                        sg = mp.sqrt(g)
                        A = xi * s * sg
                        B = b * s / sg
                        N = 1 + tn * p
                        N0 = 1 + to * p
                        L = R2 + 1 / s
                        T = mp.tanh(A)
                        if kind == "bare":
                            return B * L * T + N
                        TA = mp.tanh(mp.mpf("0.2") * s)
                        return 3 * s * TA * (R1 * (B * L * T + N) + N0 * (L + N * T / B)) + (B * L * T + N)
                    sr = mp.findroot(fm, mp.mpc(s0.real, s0.imag))
                    if abs(fm(sr)) < 1e-20 and abs(mp.arg(sr)) <= mp.pi / 4 + 1e-12 and abs(sr) > rho1:
                        if not any(abs(sr - q) < 1e-10 * abs(sr) for q in roots):
                            roots.append(sr)
                except Exception:
                    pass
    return [complex(r * r) for r in roots]   # return p


rng = np.random.default_rng(11)
out = {}
EPS = [-2.0, -1.0, -0.5, -0.2, -0.05, 0.0, 0.084, 0.3, 0.5, 1.0, 1.5, 2.5, 5.0]
for kind in ("bare", "stack"):
    for eps in EPS:
        nunst, nposw, ex = 0, 0, []
        ntr = 60
        for t in range(ntr):
            xi = float(np.exp(rng.uniform(np.log(0.2), np.log(5))))
            b = float(np.exp(rng.uniform(np.log(0.2), np.log(5))))
            tR = float(np.exp(rng.uniform(np.log(1e-3), np.log(3))))
            tl = float(np.exp(rng.uniform(np.log(1e-3), np.log(3))))
            if t % 4 == 0:
                tl = tR
            R1 = float(np.exp(rng.uniform(np.log(1e-2), np.log(2))))
            R2 = 0.0 if t % 3 == 0 else float(np.exp(rng.uniform(np.log(1e-2), np.log(2))))
            tn = eps * tl
            to = 0.0 if kind == "bare" else tn
            prm = (xi, b, tR, tl, tn, to, R1, R2)
            w, jmp = winding(prm, kind)
            nz = int(round(w))
            if nz != 0:
                nunst += 1
                if len(ex) < 4:
                    roots = locate(prm, kind)
                    ex.append(dict(params=dict(xi=xi, b=b, tR=tR, tl=tl, R1=R1, R2=R2), winding=w,
                                   roots_p=[[r.real, r.imag] for r in roots]))
        out[f"{kind} eps={eps:g}"] = dict(trials=ntr, n_with_RHP_poles=nunst, examples=ex)
        print(f"{kind:5s} eps={eps:6g}: RHP poles in {nunst}/{ntr} random parameter sets", flush=True)
        for e in ex[:2]:
            print("      e.g.", {k: round(v, 4) for k, v in e["params"].items()}, "winding %.3f" % e["winding"],
                  "poles p =", [complex(round(a, 4), round(b_, 4)) for a, b_ in e["roots_p"]][:4], flush=True)

# AlN note-16 stack (nondimensional) with kinetic coefficients of either sign
sys.path.insert(0, DIR)
import cxmodel as cx
th, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=10e-12)
for eps in [-1.0, -0.5, -0.1, 0.0, 0.084, 1.0]:
    xi, b, tR, tl, R1, R2 = th
    prm = (xi, b, tR, tl, eps * tl, eps * tl, R1, R2)
    w, _ = winding(prm, "stack", rho1=1e-3, rho2=3e4)
    out[f"AlN note16 stack eps={eps:g}"] = dict(winding=w)
    print(f"AlN note-16 stack eps={eps:g}: winding (RHP poles) {w:.3f}; expected growth rate scale "
          f"1/(|eps| tau_l) = {1 / (abs(eps) * tl * t0) if eps else float('inf'):.3e} 1/s", flush=True)
with open(DIR + "cx_10_stability.json", "w") as f:
    json.dump(out, f, indent=1)
