"""Referee spot-check: minimal re-implementation of the event operator from the written definitions
(02-assumptions A/B, memgap docstring), applied to the campaign's stored phono3py pp files.
No campaign module is imported. No orbit averaging (raw per-event rates): differences of O(1e-4)
to the campaign's orbit-averaged convention are expected (E1, 09 item 2).

Outputs, per mesh and tensor: kappa_x, kappa_z, tau_x, tau_z (ps), degeneracy-averaged r = diag C,
and for the reference: frequency-binned shares of K and N (N = |C^+ b|^2) and RTA proxies.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

ROOT = Path(r"C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results")
FOURPI = 4.0 * np.pi


class Mesh:
    def __init__(self, ev_path):
        Z = np.load(ev_path)
        f = Z["freqs"]
        self.Nq, self.nb = f.shape
        cutoff = float(Z["cutoff"])
        valid = (f > cutoff).ravel()
        self.valid = valid
        newidx = -np.ones(self.Nq * self.nb, dtype=np.int64)
        newidx[valid] = np.arange(valid.sum())
        eps = f.ravel()[valid]
        self.eps = eps
        self.n = n = len(eps)
        T, KB, THzToEv = float(Z["T"]), float(Z["KB"]), float(Z["THzToEv"])
        self.kT = kT = KB * T / THzToEv
        vel = Z["gv"].reshape(self.Nq * self.nb, 3)[valid]
        self.vel = vel
        P, A, B = newidx[Z["ev_p"]], newidx[Z["ev_a"]], newidx[Z["ev_b"]]
        assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
        rep = A == B
        fp, fa, fb = eps[P], eps[A], eps[B]
        lp = np.sqrt((fa + fb) / fp)
        ld = 1.0 / lp
        nocc = 1.0 / np.expm1(eps / kT)
        d = np.sqrt(nocc * (1.0 + nocc))
        self.d = d
        idx = np.stack([P, A, np.where(rep, n, B)], axis=1)          # slot n = dummy
        coef = np.stack([-lp, np.where(rep, 2 * ld, ld), np.where(rep, 0.0, ld)], axis=1)
        dd = np.append(d, 1.0)
        vals = coef / dd[idx]
        n1 = n + 1
        self.n1 = n1
        self.pair_flat = np.stack([idx[:, k] * n1 + idx[:, l] for k in range(3) for l in range(3)], axis=1).ravel()
        self.pair_val = np.stack([vals[:, k] * vals[:, l] for k in range(3) for l in range(3)], axis=1)
        s = np.sinh(eps / kT / 2)
        self.sp, self.sa = s[P], s[A]
        self.sb = np.where(rep, self.sa, s[np.where(rep, 0, B)])
        self.rep = rep
        self.e = d * eps
        self.ehat = self.e / np.linalg.norm(self.e)
        self.b = (d * eps)[:, None] * vel
        self.kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * THzToEv ** 2 / (2 * self.Nq * KB * T ** 2)
        # sanity: energy conservation of every event vector
        self.qidx = (np.arange(self.Nq * self.nb) // self.nb)[valid]
        self.grp = self._degroups()
        self.m = len(P)

    def _degroups(self):
        f = self.eps
        grp = -np.ones(self.n, int)
        g = 0
        for q in np.unique(self.qidx):
            ii = np.nonzero(self.qidx == q)[0]
            ii = ii[np.argsort(f[ii])]
            st = 0
            for k in range(1, len(ii) + 1):
                if k == len(ii) or f[ii[k]] - f[ii[k - 1]] > 1e-6:
                    grp[ii[st:k]] = g
                    g += 1
                    st = k
        return grp

    def degavg(self, r):
        return (np.bincount(self.grp, weights=r) / np.bincount(self.grp))[self.grp]

    def rates(self, ppfile, g0w_ref=None):
        Pp = np.load(ppfile)
        conv, g0w, pp = float(Pp["conv"]), Pp["g0w"], Pp["pp"]
        if g0w_ref is not None:
            assert np.allclose(g0w, g0w_ref, rtol=1e-12, atol=0), "g0w differs between pp files"
        W = np.where(self.rep, conv * g0w / (8 * self.sp * self.sa * self.sa),
                     conv * g0w / (4 * self.sp * self.sa * self.sb))
        return W * pp, g0w

    def C(self, g):
        w = (g[:, None] * self.pair_val).ravel()
        Cf = np.bincount(self.pair_flat, weights=w, minlength=self.n1 * self.n1)
        return Cf.reshape(self.n1, self.n1)[: self.n, : self.n].copy()

    def metrics(self, g, want_x=False):
        C = self.C(g)
        r = np.diag(C).copy()
        # energy invariant check
        en = float(np.linalg.norm(C @ self.e) / np.linalg.norm(r * self.e))
        Creg = C + (np.trace(C) / self.n) * np.outer(self.ehat, self.ehat)
        del C
        cf = sla.cho_factor(Creg, lower=True, overwrite_a=True, check_finite=False)
        out = {"r": r, "energy_residual": en}
        for c, nm in ((0, "x"), (2, "z")):
            b = self.b[:, c]
            x = sla.cho_solve(cf, b, check_finite=False)
            K, N = float(b @ x), float(x @ x)
            out[f"kappa_{nm}"] = self.kf * K
            out[f"tau_{nm}"] = N / K / FOURPI
            out[f"tauRTA_{nm}"] = float(np.sum(b * b / r ** 2) / np.sum(b * b / r)) / FOURPI
            out[f"kappaRTA_{nm}"] = self.kf * float(np.sum(b * b / r))
            if want_x:
                out[f"x_{nm}"] = x
                out[f"K_{nm}"], out[f"N_{nm}"] = K, N
        return out


def lifetime_stats(mesh, r, r0):
    dev = mesh.degavg(r) / mesh.degavg(r0) - 1.0
    # one value per degenerate group
    _, first = np.unique(mesh.grp, return_index=True)
    dg = dev[first]
    return {"min": float(dg.min()), "max": float(dg.max()), "rms": float(np.sqrt(np.mean(dg ** 2))),
            "frac_outside_10pct": float(np.mean(np.abs(dg) > 0.10)),
            "frac_outside_20pct": float(np.mean(np.abs(dg) > 0.20)),
            "n_groups": int(len(dg))}


def lowfreq_profile(mesh, m0):
    f = mesh.eps
    bins = [0, 2, 4, 6, 8, 10, 15, 40]
    prof = {}
    for nm, c in (("x", 0), ("z", 2)):
        x = m0[f"x_{nm}"]
        b = mesh.b[:, c]
        Kc, Nc = b * x, x * x
        rows = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            sel = (f >= lo) & (f < hi)
            rows.append({"f_THz": [lo, hi], "n_modes": int(sel.sum()),
                         "K_share": float(Kc[sel].sum() / Kc.sum()), "N_share": float(Nc[sel].sum() / Nc.sum())})
        # share of N from the 2% lowest-frequency modes
        order = np.argsort(f)
        k = max(1, int(0.02 * len(f)))
        prof[nm] = {"bins": rows, "N_share_lowest2pct_modes": float(Nc[order[:k]].sum() / Nc.sum()),
                    "K_share_lowest2pct_modes": float(Kc[order[:k]].sum() / Kc.sum()),
                    "f_lowest2pct_max": float(f[order[k - 1]])}
    prof["f_min_THz"] = float(f.min())
    return prof


def low_freq_rate_exponent(mesh, r):
    """Fit r ~ f^n on acoustic-like modes below 4 THz (degeneracy averaged), per mesh (diagnostic only)."""
    f = mesh.eps
    ra = mesh.degavg(r)
    sel = f < 4.0
    if sel.sum() < 5:
        return None
    A = np.vstack([np.log(f[sel]), np.ones(sel.sum())]).T
    coef, *_ = np.linalg.lstsq(A, np.log(ra[sel]), rcond=None)
    return {"n_fit_modes": int(sel.sum()), "exponent": float(coef[0])}


def run(meshname, ev_path, ref_pp, others, out):
    t0 = time.perf_counter()
    M = Mesh(ev_path)
    g0, g0w = M.rates(ref_pp)
    m0 = M.metrics(g0, want_x=True)
    res = {"mesh": meshname, "n_modes": M.n, "n_events": M.m,
           "REF": {k: m0[k] for k in ("kappa_x", "kappa_z", "tau_x", "tau_z", "tauRTA_x", "tauRTA_z",
                                       "kappaRTA_x", "kappaRTA_z", "energy_residual")},
           "lowfreq_profile": lowfreq_profile(M, m0),
           "rate_exponent_below_4THz": low_freq_rate_exponent(M, m0["r"])}
    print(meshname, "REF", json.dumps(res["REF"]), flush=True)
    for nm, pth in others.items():
        g, _ = M.rates(pth, g0w_ref=g0w)
        m = M.metrics(g)
        res[nm] = {"kappa_ratio": [m["kappa_x"] / m0["kappa_x"], m["kappa_z"] / m0["kappa_z"]],
                   "tau_ratio": [m["tau_x"] / m0["tau_x"], m["tau_z"] / m0["tau_z"]],
                   "lifetimes_degavg": lifetime_stats(M, m["r"], m0["r"])}
        print(meshname, nm, json.dumps(res[nm]), flush=True)
    res["seconds"] = time.perf_counter() - t0
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    which = sys.argv[1]
    outdir = Path(sys.argv[2])
    if which == "m553":
        d = ROOT / "m553"
        others = {"alamode": d / "alamode" / "pp_alm_r01_alamode.npz",
                  "orb0.1_c0_max(E15)": d / "witness" / "pp_orb_r01_orb0.1_c0_max.npz",
                  "orb0.1_c2_max(E15)": d / "witness" / "pp_orb_r01_orb0.1_c2_max.npz",
                  "cal1_c0_max(E18)": d / "witness" / "pp_cal_r01_cal1_c0_max.npz",
                  "cal1_c2_max(E18)": d / "witness" / "pp_cal_r01_cal1_c2_max.npz",
                  "cal1nz_J_c0_max(E20)": d / "joint_pp" / "pp_r01_cal1nz_J_c0_max.npz",
                  "cal1nz_J_c2_max(E20)": d / "joint_pp" / "pp_r01_cal1nz_J_c2_max.npz",
                  "orb0.1_J_c0_max(E20)": d / "joint_pp" / "pp_r01_orb0.1_J_c0_max.npz",
                  "orb0.1_J_c2_max(E20)": d / "joint_pp" / "pp_r01_orb0.1_J_c2_max.npz"}
        run("m553", d / "events_m553_s0.1.npz", d / "pp_r01_REF.npz", others, outdir / "ref_m553.json")
    elif which == "m775":
        d = ROOT / "m775x"
        others = {"orb0.1_c0_max(E15)": d / "pp_r01_orb0.1_c0_max.npz",
                  "orb0.1_c2_max(E15)": d / "pp_r01_orb0.1_c2_max.npz",
                  "cal1_c0_max(E18)": d / "pp_r01_cal1_c0_max.npz",
                  "cal1_c2_max(E18)": d / "pp_r01_cal1_c2_max.npz",
                  "cal1nz_J_c0_max(E20)": d / "pp_r01_cal1nz_J_c0_max.npz",
                  "cal1nz_J_c2_max(E20)": d / "pp_r01_cal1nz_J_c2_max.npz",
                  "orb0.1_J_c0_max(E20)": d / "pp_r01_orb0.1_J_c0_max.npz",
                  "orb0.1_J_c2_max(E20)": d / "pp_r01_orb0.1_J_c2_max.npz"}
        run("m775", d / "events_m775_s0.1.npz", d / "pp_r01_REF.npz", others, outdir / "ref_m775.json")
    elif which == "m996":
        d = ROOT / "m996"
        others = {"alamode": d / "pp_r01_alamode.npz",
                  "cal1_c0_max(E18)": d / "pp_r01_cal1_c0_max.npz",
                  "cal1nz_J_c0_max(E20)": d / "joint" / "pp_r01_cal1nz_J_c0_max.npz",
                  "cal1nz_J_c2_max(E20)": d / "joint" / "pp_r01_cal1nz_J_c2_max.npz",
                  "orb0.1_J_c0_max(E20)": d / "joint" / "pp_r01_orb0.1_J_c0_max.npz",
                  "orb0.1_J_c2_max(E20)": d / "joint" / "pp_r01_orb0.1_J_c2_max.npz"}
        run("m996", d / "events_m996_s0.1.npz", d / "pp_r01_REF.npz", others, outdir / "ref_m996.json")
