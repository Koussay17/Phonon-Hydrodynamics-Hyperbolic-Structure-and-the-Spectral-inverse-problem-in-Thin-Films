"""Audit check: origin of the non-zero b.e_hat at 9x9x6 (zone-boundary plane q_z = 1/2 of an even c-mesh).
b.e = sum_mu (d f)^2 v_mu cancels between time-reversal partners (q, -q); on an even mesh the plane q_z = 1/2 maps
onto itself, so the cancellation there relies on v_z(q_z = 1/2) = 0 being reproduced by the group-velocity solver."""
import sys
import numpy as np
for mk in ("m553", "m775", "m996"):
    p = {"m553": "m553/events_m553_s0.1.npz", "m775": "m775x/events_m775_s0.1.npz", "m996": "m996/events_m996_s0.1.npz"}[mk]
    Z = np.load("C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results/" + p)
    f = Z["freqs"]; gv = Z["gv"]; addr = Z["addr"]; mesh = Z["mesh"]
    kT = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    valid = f > float(Z["cutoff"])
    x = np.where(valid, f / kT, 1.0); n = 1 / np.expm1(x); d = np.sqrt(n * (1 + n))
    e = np.where(valid, d * f, 0.0)
    w = e ** 2
    tot = [float(np.sum(w[..., None] * gv, axis=(0, 1))[k]) for k in range(3)]
    bnd = np.abs(addr[:, 2]) * 2 == mesh[2]
    onb = [float(np.sum((w[bnd])[..., None] * gv[bnd], axis=(0, 1))[k]) for k in range(3)]
    vz_bnd = float(np.abs(gv[bnd][..., 2]).max()) if bnd.any() else None
    vz_all = float(np.abs(gv[..., 2]).max())
    print(mk, "sum (d f)^2 v =", ["%.3e" % t for t in tot], " from q_z=1/2 plane:", ["%.3e" % t for t in onb],
          " n_boundary_q:", int(bnd.sum()), " max|v_z| on plane:", vz_bnd, " max|v_z| overall:", "%.3e" % vz_all)
