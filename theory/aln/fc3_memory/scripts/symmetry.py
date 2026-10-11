"""Symmetry tying of event rates.

For a group of mode permutations (each preserving energies, entropy scales and the event set),
event orbits get one common rate gamma_o (g = P^T gamma), and the diagonal constraint keeps one
row per mode orbit (the other rows are identical for symmetric g).
"""
from __future__ import annotations

import itertools

import numpy as np
import scipy.sparse as sp


def signed_permutations(d):
    ops = []
    for perm in itertools.permutations(range(d)):
        for signs in itertools.product((1, -1), repeat=d):
            M = np.zeros((d, d), dtype=int)
            for i, (p, s) in enumerate(zip(perm, signs)):
                M[i, p] = s
            ops.append(M)
    return ops


def debye_mode_maps(geom):
    """Mode permutations for the cubic point group (signed permutations) of a Debye mesh."""
    d = geom.labels["d"]
    N = geom.labels["N"]
    mk = geom.labels["modes_k"]
    mb = geom.labels["modes_b"]
    h = (N - 1) // 2
    nb = int(mb.max()) + 1
    key = {}
    for i, (b, k) in enumerate(zip(mb, map(tuple, mk))):
        key[(int(b), k)] = i
    maps = []
    for M in signed_permutations(d):
        kk = mk @ M.T
        perm = np.array([key[(int(b), tuple(k))] for b, k in zip(mb, kk)], dtype=np.int64)
        maps.append(perm)
    return maps


def canonical_event(p, a, b):
    if b < 0:
        return (p, a, -1)
    return (p, min(a, b), max(a, b))


def orbits(geom, maps, check=True):
    ev = geom.events
    m = len(ev)
    index = {canonical_event(*map(int, e)): i for i, e in enumerate(ev)}
    parent = np.arange(m)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm in maps:
        for i, (p, a, b) in enumerate(ev):
            img = canonical_event(int(perm[p]), int(perm[a]), int(perm[b]) if b >= 0 else -1)
            j = index.get(img)
            if j is None:
                raise ValueError("event set not invariant under symmetry")
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    roots = np.array([find(i) for i in range(m)])
    uniq, ev_orb = np.unique(roots, return_inverse=True)
    # mode orbits
    n = geom.n
    mpar = np.arange(n)

    def mfind(i):
        while mpar[i] != i:
            mpar[i] = mpar[mpar[i]]
            i = mpar[i]
        return i

    for perm in maps:
        for i in range(n):
            ri, rj = mfind(i), mfind(int(perm[i]))
            if ri != rj:
                mpar[ri] = rj
    mroots = np.array([mfind(i) for i in range(n)])
    muniq, mode_orb = np.unique(mroots, return_inverse=True)
    reps = np.array([np.nonzero(mode_orb == o)[0][0] for o in range(len(muniq))])
    if check:
        # invariance of energies and entropy scales
        for perm in maps:
            assert np.allclose(geom.eps[perm], geom.eps)
    P = sp.coo_matrix((np.ones(m), (ev_orb, np.arange(m))), shape=(len(uniq), m)).tocsr()
    return {"P": P, "ev_orb": ev_orb, "mode_orb": mode_orb, "reps": reps,
            "n_ev_orb": len(uniq), "n_mode_orb": len(muniq)}


def time_reversal_maps(geom):
    """Only k -> -k (time reversal / inversion)."""
    d = geom.labels["d"]
    N = geom.labels["N"]
    mk = geom.labels["modes_k"]
    mb = geom.labels["modes_b"]
    key = {(int(b), tuple(k)): i for i, (b, k) in enumerate(zip(mb, map(tuple, mk)))}
    perm = np.array([key[(int(b), tuple(-k))] for b, k in zip(mb, mk)], dtype=np.int64)
    return [np.arange(geom.n), perm]
