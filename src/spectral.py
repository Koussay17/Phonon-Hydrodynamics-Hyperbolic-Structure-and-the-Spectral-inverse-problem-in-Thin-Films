"""Mode-resolved thermal diagnostics in SI units; no inferred GK times.

The conserving collision helper is a relaxation-model construction, not an
ab initio collision matrix. RTA lifetimes are not resistive/normal lifetimes.
"""
import numpy as np
from scipy.constants import Boltzmann as KB, hbar


def mode_heat_capacity(omega_rad_s, temperature):
    """Heat capacity per oscillator (J/K), including its zero-frequency limit."""
    w = np.asarray(omega_rad_s, dtype=float)
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("angular frequencies must be finite and nonnegative")
    x = hbar * w / (KB * temperature)
    # x/(2 sinh(x/2)) squared; avoid overflow and cancellation.
    out = np.zeros_like(x)
    small = x < 1e-4
    out[small] = 1 - x[small]**2/12 + x[small]**4/240
    mid = (~small) & (x < 700)
    out[mid] = x[mid]**2 * np.exp(-x[mid]) / np.expm1(-x[mid])**2
    return KB * out


def rta_moments(capacity, velocity_sq, rate):
    """DC kappa and two distinct directional times for diagonal RTA.

    capacity includes quadrature weight and volume, in J/(m^3 K);
    velocity_sq is the squared directional velocity or basal-plane average.
    tau_static matches kappa with the exact velocity moment.
    tau_memory matches the derivative of kappa(p) at p=0.
    """
    c, v2, r = np.broadcast_arrays(
        np.asarray(capacity, float), np.asarray(velocity_sq, float),
        np.asarray(rate, float))
    if any(np.any(~np.isfinite(a)) for a in (c, v2, r)):
        raise ValueError("finite mode inputs required")
    if np.any(c < 0) or np.any(v2 < 0) or np.any(r < 0):
        raise ValueError("nonnegative mode inputs required")
    a = c * v2
    active = a > 0
    if not np.any(active) or np.any(r[active] == 0):
        raise ValueError("finite relaxation required for heat-carrying modes")
    tau = np.divide(1., r, out=np.zeros_like(r), where=r > 0)
    weights = a * tau
    kappa = np.sum(weights)
    return {"kappa": float(kappa),
            "velocity_moment": float(np.sum(a)),
            "tau_static_s": float(kappa / np.sum(a)),
            "tau_memory_s": float(np.sum(weights * tau) / kappa),
            "contributions": weights}


def inplane_suppression(kn_normal, specularity=0.):
    """Steady film RTA with identical planar reflecting boundaries.

    kn_normal = |v_normal| tau_bulk / thickness. Applies to heat flow PARALLEL
    to the film, with specularity in [0,1], no transmission and uniform
    in-plane gradient. This is not a cross-plane FDTR suppression law.
    """
    kn = np.asarray(kn_normal, float)
    p = float(specularity)
    if np.any(~np.isfinite(kn)) or np.any(kn < 0):
        raise ValueError("finite nonnegative normal Knudsen numbers required")
    if not np.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("specularity must be between zero and one")
    if p == 1:
        return np.ones_like(kn)
    t = np.divide(1., kn, out=np.full_like(kn, np.inf), where=kn > 0)
    # S = [(1-p)(1-(1-exp(-t))/t)+p(1-exp(-t))] /
    #     [(1-p)+p(1-exp(-t))], stably at t -> 0.
    u = -np.expm1(-t)
    a = np.ones_like(t)
    small = t < 1e-3
    a[small] = t[small]/2 - t[small]**2/6 + t[small]**3/24 - t[small]**4/120
    normal = (~small) & np.isfinite(t)
    a[normal] = 1 - u[normal]/t[normal]
    return ((1-p)*a + p*u) / ((1-p) + p*u)


def conserving_collision_action(rate, invariants, perturbation):
    """Apply D - D H (H.T D H)^+ H.T D without a dense mode matrix.

    H columns span conserved moments in entropy-normalized coordinates.
    For N collisions use energy + crystal momentum; for R use energy only.
    The mode quadrature weights must already be included in H and y.
    Input rates are bare model parameters: projection changes the resulting
    collision diagonal. Positive rates required. SVD avoids an ill-conditioned
    Gram inversion.
    """
    r = np.asarray(rate, float)
    h = np.asarray(invariants, float)
    y = np.asarray(perturbation, float)
    if r.ndim != 1 or y.shape != r.shape or h.ndim != 2 or h.shape[0] != len(r):
        raise ValueError("incompatible rate, invariant and perturbation shapes")
    if (np.any(~np.isfinite(r)) or np.any(r <= 0)
            or np.any(~np.isfinite(h)) or np.any(~np.isfinite(y))):
        raise ValueError("finite values and strictly positive rates required")
    if h.shape[1] == 0:
        return r*y
    scale = np.linalg.norm(h, axis=0)
    if np.any(scale == 0):
        raise ValueError("invariants must have nonzero columns")
    root = np.sqrt(r)
    u, singular, _ = np.linalg.svd(root[:, None]*(h/scale), full_matrices=False)
    rank = np.sum(singular > singular[0]*max(h.shape)*np.finfo(float).eps)
    u = u[:, :rank]
    z = root*y
    return root*(z-u@(u.T@z))


def collision_order_fractions(normal_rate, resistive_rate, boundary_rate,
                              weights, separation=10.):
    """Weighted *ordering* diagnostics, not a proof of hydrodynamics.

    Rates refer to the SAME modes. Boundary rate is a separate diagnostic;
    resistive_rate must be bulk only (no boundary double counting).
    """
    rn, rr, rb, w = np.broadcast_arrays(
        *[np.asarray(a, float) for a in
          (normal_rate, resistive_rate, boundary_rate, weights)])
    if (any(np.any(~np.isfinite(a)) or np.any(a < 0) for a in (rn, rr, rb, w))
            or w.sum() <= 0 or not np.isfinite(separation) or separation <= 1):
        raise ValueError("nonnegative finite rates/weights and separation > 1 required")
    candidate = (rn >= separation*rb) & (rb >= separation*rr) & (rb > 0)
    boundary_first = rb >= separation*(rn+rr)
    return {"candidate_weight_fraction": float(w[candidate].sum()/w.sum()),
            "boundary_first_weight_fraction": float(w[boundary_first & (rb > 0)].sum()/w.sum())}
