"""Cross-check of the 1D stack (ab_common.stack_Z, impedance recursion) against the project's digital twin
(src/twin_fdtr.py, kernel at k = 0: modal linear system), Cattaneo and GK (local temperature) AlN layer."""
import json, os, sys
import numpy as np
sys.path.insert(0, r"C:/Users/Koussay/these/src")
import twin_fdtr as tw
from ab_common import BASE, with_, stack_Z

HERE = os.path.dirname(os.path.abspath(__file__))
f = np.logspace(4, np.log10(2e8), 80)
out = {}
for label, tau_r, ell in (("Cattaneo 27.5 ps", 27.5e-12, 0.0), ("GK tau_l = tau_R/2", 27.5e-12, None), ("Fourier", 0.0, 0.0)):
    lam, C = BASE["lam"], BASE["C"]
    a = lam / C
    alpha = 1.0 / 3.0
    if ell is None:   # choose ell so that tau_l = (1+alpha) ell^2 / a = tau_R/2
        ell = np.sqrt(0.5 * tau_r * a / (1 + alpha))
    aln = tw.Material(lam, C, tau_r=tau_r, ell=ell, alpha=alpha)
    st = tw.au_aln_sapphire(aln, aln_thickness=BASE["e"], au_thickness=BASE["e_au"],
                            r_au_aln=BASE["R1"], r_aln_sapphire=BASE["R2"])
    Ztw = tw.kernel(f, 0.0, st)
    P = with_(BASE, tau_r=tau_r, tau_l=aln.tau_l)
    Zab = stack_Z(2 * np.pi * f, P)
    err = float(np.max(np.abs(Ztw / Zab - 1)))
    dph = float(np.max(np.abs(np.angle(Ztw / Zab))) * 180 / np.pi)
    out[label] = {"max_rel_err": err, "max_phase_diff_deg": dph, "tau_l": aln.tau_l}
    print("%-20s tau_l = %.3e s  max |Z_twin/Z_1D - 1| = %.2e, max phase diff = %.2e deg" % (label, aln.tau_l, err, dph))
# phase sensitivity to tau_R via the twin (finite difference w = 1e-3) vs dual-number 1D
aln1 = tw.Material(BASE["lam"], BASE["C"], tau_r=27.5e-12 * 1.001)
aln0 = tw.Material(BASE["lam"], BASE["C"], tau_r=27.5e-12 * 0.999)
mk = lambda m: tw.kernel(f, 0.0, tw.au_aln_sapphire(m, aln_thickness=BASE["e"], au_thickness=BASE["e_au"], r_au_aln=BASE["R1"], r_aln_sapphire=BASE["R2"]))
s_tw = np.angle(mk(aln1) / mk(aln0)) / 0.002
from ab_common import dlnZ
s_ab = dlnZ(2 * np.pi * f, BASE, "tau_r").imag
out["dphi_dlntau_max_abs_diff_deg"] = float(np.max(np.abs(s_tw - s_ab)) * 180 / np.pi)
out["dphi_dlntau_max_deg"] = float(np.max(np.abs(s_ab)) * 180 / np.pi)
print("d phi/d ln tau_R: twin FD vs 1D dual: max diff %.2e deg (signal max %.3e deg)" % (out["dphi_dlntau_max_abs_diff_deg"], out["dphi_dlntau_max_deg"]))
json.dump(out, open(os.path.join(HERE, "results", "x1_twin_crosscheck.json"), "w"), indent=1)
