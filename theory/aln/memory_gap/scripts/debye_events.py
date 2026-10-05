"""Debye-type event sets: isotropic linear dispersion on periodic meshes.

Model (declared):
* simple-cubic reciprocal lattice, mesh N^d (N odd: unique Brillouin-zone representatives
  kappa in {-(N-1)/2..(N-1)/2}^d, so q(-k) = -q(k) and no zone-boundary self-conjugate points);
* branches with linear isotropic dispersion eps_b(k) = c_b |kappa(k)| (Euclidean norm of the
  zone representative), group velocity v = c_b kappa/|kappa|; the zero-frequency modes at k = 0
  are excluded;
* three-phonon events p -> a + b' with k_p = k_a + k_b' (mod N); normal if no folding is needed,
  umklapp otherwise.

Resonance:
* exact (integer arithmetic) in d = 1, where energies are integers c_b |kappa|;
* in d >= 2 a declared tolerance |eps_p - eps_a - eps_b'| <= tol * (2 pi / N) * c_min, with
  conserving stoichiometry s = -l_p e_p + l_d (e_a + e_b'), l_p = sqrt((eps_a+eps_b')/eps_p),
  l_d = 1/l_p, so s^T eps = 0 exactly (rank-one, PSD, energy-conserving event).
Reference amplitude (declared, Klemens-type): Gamma = eps_p eps_a eps_b' / eps_max^3 times the
resonance weight (1 for exact events in d = 1 after 1/N normalization; box 1/(2 tol) otherwise),
divided by 2 for repeated daughters; g = Gamma * sqrt(N_p(1+N_a)(1+N_b) (1+N_p) N_a N_b) / N^d.
Units: energies in units of the zone-edge transverse energy scale; rates are model units.
"""
from __future__ import annotations

import itertools

import numpy as np

from memgap import EventGeometry


def mesh_modes(d, N, speeds):
    assert N % 2 == 1, "use odd N"
    h = (N - 1) // 2
    ks = np.array(list(itertools.product(range(-h, h + 1), repeat=d)), dtype=int)
    ks = ks[np.any(ks != 0, axis=1)]                     # drop Gamma
    nq = len(ks)
    modes_k = np.vstack([ks for _ in speeds])
    modes_b = np.repeat(np.arange(len(speeds)), nq)
    unit = 2 * np.pi / N
    q = modes_k * unit
    qn = np.linalg.norm(q, axis=1)
    c = np.asarray(speeds, dtype=float)[modes_b]
    eps = c * qn
    vel = c[:, None] * q / qn[:, None]
    return modes_k, modes_b, eps, vel, nq


def fold(k, N):
    h = (N - 1) // 2
    return ((k + h) % N) - h


def build(d, N, speeds=(1.0, 2.0), kT=None, tol=None, exact=None, include_repeated=True,
          amp_power=1.0, max_events=5_000_000):
    """Return EventGeometry for the Debye-type model.

    exact: if None, exact for d == 1 (requires integer speeds), tolerance otherwise.
    tol: tolerance in units of c_min * (2 pi / N) (d >= 2), default 0.5.
    kT: temperature in energy units; default 0.5 * max energy.
    amp_power: Gamma proportional to (eps_p eps_a eps_b)^amp_power.
    """
    modes_k, modes_b, eps, vel, nq = mesh_modes(d, N, speeds)
    n = len(eps)
    if exact is None:
        exact = (d == 1)
    if kT is None:
        kT = 0.5 * eps.max()
    unit = 2 * np.pi / N
    if exact:
        # integer energies: c_b |kappa| with integer speeds
        assert d == 1 and all(float(c).is_integer() for c in speeds)
        eint = (np.asarray(speeds, dtype=int)[modes_b] * np.abs(modes_k[:, 0])).astype(np.int64)
    if tol is None:
        tol = 0.5
    tol_abs = tol * min(speeds) * unit
    # index lookup: (branch, k-vector) -> mode index
    h = (N - 1) // 2
    def kid(kv):
        z = np.zeros(len(kv), dtype=np.int64)
        for j in range(d):
            z = z * N + (kv[:, j] + h)
        return z
    lut = -np.ones((len(speeds), N ** d), dtype=np.int64)
    lut[modes_b, kid(modes_k)] = np.arange(n)
    ev, co, kind, dlt = [], [], [], []
    for a in range(n):
        bs = np.arange(a, n) if include_repeated else np.arange(a + 1, n)
        ksum = modes_k[a][None, :] + modes_k[bs]
        kp = fold(ksum, N)
        umk = np.any(kp != ksum, axis=1)
        nonzero = np.any(kp != 0, axis=1)
        for bp in range(len(speeds)):
            p = lut[bp, kid(kp)]
            ok = nonzero & (p >= 0)
            if exact:
                good = ok & (eint[p] == eint[a] + eint[bs])
                delta = np.zeros(len(bs))
            else:
                delta = eps[p] - eps[a] - eps[bs]
                good = ok & (np.abs(delta) <= tol_abs)
            idx = np.nonzero(good)[0]
            for i in idx:
                bb = bs[i]
                pp = p[i]
                if pp == a or pp == bb:
                    continue
                ev.append((pp, a, bb if bb != a else -1))
                kind.append(int(umk[i]))
                dlt.append(delta[i])
        if len(ev) > max_events:
            raise MemoryError("too many events")
    ev = np.array(ev, dtype=np.int64)
    kind = np.array(kind, dtype=np.int64)
    dlt = np.array(dlt, dtype=float)
    m = len(ev)
    rep = ev[:, 2] < 0
    ep = eps[ev[:, 0]]
    ea = eps[ev[:, 1]]
    eb = np.where(rep, ea, eps[np.maximum(ev[:, 2], 0)])
    sd = ea + eb
    if exact:
        lp = np.ones(m); ld = np.ones(m)
    else:
        lp = np.sqrt(sd / ep); ld = 1.0 / lp
    coef = np.zeros((m, 3))
    coef[:, 0] = -lp
    coef[:, 1] = np.where(rep, 2.0 * ld, ld)
    coef[:, 2] = np.where(rep, 0.0, ld)
    # reference rates
    x = eps / kT
    nocc = 1.0 / np.expm1(x)
    Np, Na = nocc[ev[:, 0]], nocc[ev[:, 1]]
    Nb = np.where(rep, Na, nocc[np.maximum(ev[:, 2], 0)])
    fwd = Np * (1 + Na) * (1 + Nb)
    bwd = (1 + Np) * Na * Nb
    emax = eps.max()
    Gam = ((ep * ea * eb) / emax ** 3) ** amp_power
    Gam = np.where(rep, 0.5 * Gam, Gam)
    w = 1.0 if exact else 1.0 / (2 * tol_abs)
    gphys = Gam * w * np.sqrt(fwd * bwd) / (N ** d)
    geom = EventGeometry(eps=eps, kT=kT, vel=vel, events=ev, coef=coef, kind=kind, gphys=gphys,
                         labels={"d": d, "N": N, "speeds": tuple(speeds), "exact": exact, "tol": tol,
                                 "modes_k": modes_k, "modes_b": modes_b, "delta": dlt,
                                 "amp_power": amp_power})
    return geom


if __name__ == "__main__":
    import sys
    d = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    g = build(d, N)
    print("n", g.n, "m", g.m, "umklapp", int(g.kind.sum()), "repeated", int((g.events[:, 2] < 0).sum()))
    nk, w = g.check_kernel(g.gphys)
    print("kernel dim", nk, "smallest eig", w[:4])
    print("max |s^T eps| / eps_p", np.max(np.abs((g.coef[:, 0] * g.eps[g.events[:, 0]] + g.coef[:, 1] * g.eps[g.events[:, 1]] + np.where(g.events[:, 2] >= 0, g.coef[:, 2] * g.eps[np.maximum(g.events[:, 2], 0)], 0)) / g.eps[g.events[:, 0]])))
    print("b.e", float(g.b[:, 0] @ g.e) / np.linalg.norm(g.b[:, 0]) / np.linalg.norm(g.e))
