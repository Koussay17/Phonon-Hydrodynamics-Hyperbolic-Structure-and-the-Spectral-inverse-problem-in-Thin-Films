"""cx_11: practical identifiability for the Au/AlN/sapphire note-16 stack, 1D, 80 frequencies 10 kHz-200 MHz.

Misfit m(th') = RMS_j |Z(omega_j; th')/Z(omega_j; th*) - 1| (complex relative; 1e-4 ~ 0.006 deg phase).
(1) profile in tau_l': min over the other five parameters with tau_l' fixed (0 .. 100 ps), target tau_l = 10 ps,
    entropic and local conventions;
(2) profile in tau_R' (same);
(3) convention confusion: data from the entropic model, fitted with the local model (and conversely),
    all six parameters free.
No instrument nuisances (they can only lower the misfits).
"""
import sys as _s
_s.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle as _t
_t.wait_slot(2)          # light wait before any numerical import
import json
import sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
th_t, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=10e-12)
f = np.logspace(4, np.log10(2e8), 80)
om = 2 * np.pi * f * t0
ps = 1e-12 / t0     # one picosecond in model time units
out = dict(t0=t0, target=list(th_t), transducer=(xiA, bA))


def fit(spec, Zt, x0, fixed):
    """least squares in log variables for the free coordinates; fixed: dict index -> value."""
    free = [i for i in range(6) if i not in fixed]

    def th_of(y):
        th = np.zeros(6)
        for k, i in enumerate(free):
            th[i] = np.exp(y[k])
        for i, v in fixed.items():
            th[i] = v
        return th

    def res(y):
        z = cx.Z_np(spec, th_of(y), om)
        r = z / Zt - 1
        return np.concatenate([r.real, r.imag])

    def jac(y):
        th = th_of(y)
        z, J = cx.ZJ_np(spec, th, om)
        G = (J / Zt[None, :])[free] * th[free][:, None]
        return np.concatenate([G.real.T, G.imag.T], axis=0)

    best = None
    for y0 in x0:
        y0 = np.log(np.asarray(y0)[free])
        r = least_squares(res, y0, jac=jac, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        m = np.sqrt(np.mean(np.abs(r.fun[:len(om)] + 1j * r.fun[len(om):]) ** 2))
        if best is None or m < best[0]:
            best = (m, th_of(r.x))
    return best


rng = np.random.default_rng(2)
for conv, eps in (("entropic", 1.0), ("local", 0.0)):
    sp = cx.Spec("stack", "tied", eps, eps, xiA=xiA, bA=bA)
    Zt = cx.Z_np(sp, th_t, om)
    starts = [th_t] + [np.array(th_t) * np.exp(rng.uniform(-0.3, 0.3, 6)) for _ in range(4)]
    prof_l = []
    for tl_ps in [0.0, 1.0, 5.0, 20.0, 40.0, 70.0, 100.0]:
        if tl_ps == 0.0:
            fixed = {3: 0.0}
        else:
            fixed = {3: tl_ps * ps}
        st = [np.where(np.arange(6) == 3, max(tl_ps * ps, 1e-12), s) for s in starts]
        m, th = fit(sp, Zt, st, fixed)
        prof_l.append(dict(tau_l_ps=tl_ps, misfit=float(m), theta=[float(x) for x in th],
                           tau_R_ps=float(th[2] / ps)))
        print(f"{conv}: tau_l' = {tl_ps:5.1f} ps fixed -> min misfit {m:.2e}, tau_R' = {th[2]/ps:8.2f} ps, "
              f"R1' {th[4]/th_t[4]:.4f}x, R2' {th[5]/th_t[5]:.4f}x", flush=True)
    prof_R = []
    for tR_ps in [0.0, 5.0, 15.0, 40.0, 80.0]:
        fixed = {2: tR_ps * ps}
        st = [np.where(np.arange(6) == 2, max(tR_ps * ps, 1e-12), s) for s in starts]
        m, th = fit(sp, Zt, st, fixed)
        prof_R.append(dict(tau_R_ps=tR_ps, misfit=float(m), theta=[float(x) for x in th], tau_l_ps=float(th[3] / ps)))
        print(f"{conv}: tau_R' = {tR_ps:5.1f} ps fixed -> min misfit {m:.2e}, tau_l' = {th[3]/ps:8.2f} ps", flush=True)
    out[conv] = dict(profile_tau_l=prof_l, profile_tau_R=prof_R)

# convention confusion
for data_conv, fit_conv in (("entropic", "local"), ("local", "entropic")):
    e_d = 1.0 if data_conv == "entropic" else 0.0
    e_f = 1.0 if fit_conv == "entropic" else 0.0
    spd = cx.Spec("stack", "tied", e_d, e_d, xiA=xiA, bA=bA)
    spf = cx.Spec("stack", "tied", e_f, e_f, xiA=xiA, bA=bA)
    Zt = cx.Z_np(spd, th_t, om)
    starts = [th_t] + [np.array(th_t) * np.exp(rng.uniform(-0.3, 0.3, 6)) for _ in range(6)]
    m, th = fit(spf, Zt, starts, {})
    out[f"confusion data={data_conv} fit={fit_conv}"] = dict(misfit=float(m), theta=[float(x) for x in th],
                                                              tau_R_ps=float(th[2] / ps), tau_l_ps=float(th[3] / ps))
    print(f"data {data_conv}, fitted with {fit_conv}: min misfit {m:.2e} at tau_R' = {th[2]/ps:.2f} ps, tau_l' = {th[3]/ps:.2f} ps",
          flush=True)
with open(DIR + "cx_11_practical_profiles.json", "w") as f:
    json.dump(out, f, indent=1)
