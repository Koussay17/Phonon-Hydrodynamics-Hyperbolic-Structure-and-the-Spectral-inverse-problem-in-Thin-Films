"""Kinetic nonlocal time tau_l of a linearised, symmetrised collision operator (revision R-A, objection A-2).

Symmetrised linear BTE (one cross-plane coordinate z): (p + i k V + Chat) f = 0, V = diag(v_z), Chat symmetric
positive semidefinite with the energy invariant e in its kernel, current b = V e (b orthogonal to e by parity).
Eliminating the temperature theta = e.f/|e|^2 (energy balance p |e|^2 theta + i k b.g = 0, f = theta e + g, g _|_ e):
    g = -i k theta [p + Chat + i k Q V Q]^{-1}_{e-perp} b,    Q = 1 - ehat ehat^T,
so the heat current is q = -i k theta kappa(k, p) with the exact kernel kappa(k, p) = b.[p + Chat + i k QVQ]^{-1} b.
Its expansion at small (k, p) (the odd term vanishes by parity):
    kappa = K [1 - tau_R p - L^2 k^2 + ...],  K = b.x,  x = Chat^{-1} b,
    tau_R = x.x / K                                   (memory time, = -K'(0)/K(0)),
    L^2   = u.Chat^{-1} u / K,   u = Q V x            (nonlocal coefficient of tau_R q_t + q = -lam T_z + L^2 q_zz).
With the diffusivity a = K / |e|^2 of the same reduction (p |e|^2 + k^2 kappa = 0), the Guyer-Krumhansl time is
    tau_l = L^2 / a = |e|^2 (u.Chat^{-1} u) / K^2 ,   and   delta = tau_R - tau_l.
For the gray isotropic Callaway operator (one speed, N conserving energy and momentum, R conserving energy) this
gives tau_l = (4/5) tau_c, tau_c = (1/tau_N + 1/tau_R)^-1, i.e. the conserving closure of note 14; without the
projection Q (energy mode not separated) it would give (9/5) tau_c, the historical coefficient.
Kinematic factor gamma_g = |e|^2 |Q V b|^2 / |b|^4 (= 4/5 for the isotropic one-speed model); tau_l = gamma_g tau_c,eff
defines an effective viscous relaxation time tau_c,eff of a non-gray operator.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla


def kinetic_times_dense(Chat, e, vz):
    """tau_R, tau_l, L2, a, gamma_g for a dense symmetric operator with Chat e = 0 (any consistent units).

    Solves on e-perp through an orthonormal basis of e-perp (Householder), so no regularisation constant enters.
    """
    e = np.asarray(e, float)
    n = e.size
    ehat = e / np.linalg.norm(e)
    # Householder reflector H with H ehat = +-e_0; columns 1..n-1 of H span e-perp
    v = ehat.copy()
    v[0] += np.sign(v[0]) if v[0] != 0 else 1.0
    v /= np.linalg.norm(v)
    H = np.eye(n) - 2.0 * np.outer(v, v)
    B = H[:, 1:]
    Cr = B.T @ Chat @ B
    b = vz * e
    br = B.T @ b
    xr = sla.solve(Cr, br, assume_a="sym")
    x = B @ xr
    K = float(b @ x)
    tauR = float(x @ x) / K
    u = vz * x
    u = u - ehat * (ehat @ u)
    yr = sla.solve(Cr, B.T @ u, assume_a="sym")
    M = float((B.T @ u) @ yr)
    e2 = float(e @ e)
    taul = e2 * M / K ** 2
    QVb = vz * b - ehat * (ehat @ (vz * b))
    gam = e2 * float(QVb @ QVb) / float(b @ b) ** 2
    return {"tauR": tauR, "taul": taul, "L2": M / K, "a": K / e2, "K": K, "gamma_g": gam}


def kernel_k2_numeric(Chat, e, vz, ks):
    """Exact static kernel kappa(k, 0) = b.[Chat + i k QVQ]^{-1}_{e-perp} b at the given k (complex solves)."""
    e = np.asarray(e, float)
    n = e.size
    ehat = e / np.linalg.norm(e)
    v = ehat.copy()
    v[0] += np.sign(v[0]) if v[0] != 0 else 1.0
    v /= np.linalg.norm(v)
    H = np.eye(n) - 2.0 * np.outer(v, v)
    B = H[:, 1:]
    Cr = B.T @ Chat @ B
    Vr = B.T @ (vz[:, None] * B)
    br = B.T @ (vz * e)
    out = []
    for k in ks:
        A = Cr + 1j * k * Vr
        out.append(complex(br @ np.linalg.solve(A, br.astype(complex))))
    return np.array(out)


def gray_isotropic_operator(n_mu=8, n_phi=16, v=1.0, tauR=1.0, tauN=1.0):
    """Symmetrised gray Callaway operator on a product quadrature of the unit sphere.

    Modes = directions (one speed v), capacity weights w_i (sum 1); e_i = sqrt(w_i), b = v n_z e.
    Chat = (1/tauR)(1 - ehat ehat^T) + (1/tauN)(1 - ehat ehat^T - sum_j m_j m_j^T), m_j the normalised drift modes
    (m_j,i ~ sqrt(w_i) n_j,i): N conserves energy and momentum, R conserves energy only.
    """
    mu, wmu = np.polynomial.legendre.leggauss(n_mu)
    phi = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    MU, PHI = np.meshgrid(mu, phi, indexing="ij")
    W = (wmu[:, None] * np.ones(n_phi)[None, :] / n_phi / 2.0).ravel()
    st = np.sqrt(1 - MU ** 2).ravel()
    nvec = np.stack([st * np.cos(PHI.ravel()), st * np.sin(PHI.ravel()), MU.ravel()], axis=1)
    e = np.sqrt(W)
    ehat = e / np.linalg.norm(e)
    n = e.size
    P_e = np.outer(ehat, ehat)
    P_m = np.zeros((n, n))
    for j in range(3):
        m = e * nvec[:, j]
        m /= np.linalg.norm(m)
        P_m += np.outer(m, m)
    I = np.eye(n)
    Chat = (I - P_e) / tauR + (I - P_e - P_m) / tauN
    return Chat, e, v * nvec[:, 2]
