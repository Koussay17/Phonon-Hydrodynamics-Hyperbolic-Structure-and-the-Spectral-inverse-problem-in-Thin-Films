"""Finite exactly resonant, distinct-mode three-phonon event reference.

This is a validation implementation, not an AlN exporter. Modes have equal
full-grid weights. Each row lists (parent, daughter, daughter) once, with
the reverse reaction already included. Wave-vector reversal is a separate
event. Prefactors include all counting/volume conventions supplied by the
caller and are NOT inferred from diagonal lifetimes.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator


class DecayEventOperator:
    """Positive entropy-coordinate generator C, with dy/dt = -C y.

    energy_over_kbt: positive dimensionless mode energies x.
    events: integer (n_events, 3) array with x[parent] = x[j] + x[k].
    prefactors: nonnegative Gamma in inverse seconds for the explicitly
        defined population flux Gamma*[n_i(1+n_j)(1+n_k)-(1+n_i)n_j*n_k].
    resonance_rtol: roundoff tolerance on relative energy mismatch AND
        on absolute mismatch in energy/(k_B T). Not spectral broadening.

    Coordinates y = delta_n/sqrt(n0*(1+n0)); uniform mode weights only.
    Repeated indices, irreducible grids and detuned physical events are
    deliberately outside this contract.
    """

    def __init__(self, energy_over_kbt, events, prefactors, resonance_rtol=1e-12):
        x = np.asarray(energy_over_kbt, dtype=float)
        raw = np.asarray(events)
        gamma = np.asarray(prefactors, dtype=float)
        if x.ndim != 1 or x.size == 0 or np.any(~np.isfinite(x)) or np.any(x <= 0):
            raise ValueError("finite strictly positive mode energies required")
        if raw.ndim != 2 or raw.shape[1] != 3 or raw.dtype.kind not in "iu":
            raise ValueError("events must be an integer array of shape (m, 3)")
        if np.any(raw < 0) or np.any(raw >= x.size):
            raise ValueError("event mode index out of range")
        event = raw.astype(np.intp, copy=True)
        if np.any(np.diff(np.sort(event, axis=1), axis=1) == 0):
            raise ValueError("repeated indices require a separately validated counting convention")
        if gamma.shape != (len(event),) or np.any(~np.isfinite(gamma)) or np.any(gamma < 0):
            raise ValueError("one finite nonnegative prefactor per event required")
        if not np.isfinite(resonance_rtol) or not 0 <= resonance_rtol <= 1e-8:
            raise ValueError("resonance tolerance is a roundoff tolerance, at most 1e-8")
        # Divide first to avoid overflow when large dimensionless energies add.
        energies = x[event]
        scale = np.max(energies, axis=1) if len(event) else np.empty(0)
        residual = (energies[:, 0]/scale
                    - energies[:, 1]/scale - energies[:, 2]/scale)
        detuning = (energies[:, 0] - energies[:, 1]) - energies[:, 2]
        if (np.any(np.abs(residual) > resonance_rtol)
                or np.any(np.abs(detuning) > resonance_rtol)):
            raise ValueError("events must conserve energy; broadening is not supported")
        # log(1-exp(-x)) remains accurate at small x and finite at large x.
        log_one_minus = np.log(-np.expm1(-x))
        log_gamma = np.full(len(event), -np.inf)
        active = gamma > 0
        log_gamma[active] = np.log(gamma[active])
        # Form energy differences BEFORE adding log prefactors: otherwise
        # cancellation at huge x can erase a finite log(Gamma).
        exponent = (0.5*log_gamma[:, None]
                    - 0.5*(energies[:, :1] - energies)
                    + log_one_minus[event]
                    - 0.5*np.sum(log_one_minus[event], axis=1)[:, None])
        with np.errstate(over="raise", invalid="raise"):
            values = np.exp(exponent) * np.array([-1., 1., 1.])
        if np.any(~np.isfinite(values)):
            raise ValueError("event factor is outside representable floating-point range")
        self.energy_over_kbt = x.copy()
        self.events = event
        self.prefactors = gamma.copy()
        self.relative_energy_residual = residual
        self.detailed_balance_log_residual = detuning
        self.factor = coo_matrix(
            (values.ravel(), (np.repeat(np.arange(len(event)), 3), event.ravel())),
            shape=(len(event), x.size)).tocsr()
        with np.errstate(over="ignore"):
            self.entropy_scale = np.exp(-x/2)/(-np.expm1(-x))
        if np.any(~np.isfinite(self.entropy_scale)):
            raise ValueError("entropy scale exceeds floating-point range")
        # A finite factor need not give a representable collision diagonal.
        with np.errstate(over="ignore"):
            diagonal = np.asarray(self.factor.power(2).sum(axis=0)).ravel()
        if np.any(~np.isfinite(diagonal)):
            raise ValueError("collision diagonal exceeds floating-point range")

    def action(self, perturbation):
        """Apply C = R.T R without allocating a dense mode-by-mode matrix."""
        y = np.asarray(perturbation)
        if y.shape != self.energy_over_kbt.shape or np.any(~np.isfinite(y)):
            raise ValueError("finite mode vector with the operator shape required")
        return self.factor.T @ (self.factor @ y)

    def as_linear_operator(self):
        """Real self-adjoint operator, also applicable to complex vectors."""
        n = len(self.energy_over_kbt)
        def matvec(y):
            return self.action(np.asarray(y).reshape(n))
        def matmat(y):
            return self.factor.T @ (self.factor @ y)
        return LinearOperator((n, n), matvec=matvec, rmatvec=matvec,
                              matmat=matmat, rmatmat=matmat, dtype=np.dtype(float))

    def population_rhs(self, occupation):
        """Independent nonlinear population flux; intended for moderate inputs.

        Direct products can overflow or cancel at extreme occupations.
        Use arbitrary precision for such reference checks.
        """
        n = np.asarray(occupation, dtype=float)
        if n.shape != self.energy_over_kbt.shape or np.any(~np.isfinite(n)) or np.any(n < 0):
            raise ValueError("finite nonnegative mode occupations required")
        i, j, k = self.events.T
        flux = self.prefactors * (n[i]*(1+n[j])*(1+n[k]) - (1+n[i])*n[j]*n[k])
        out = np.zeros_like(n)
        for index, sign in ((i, -1.), (j, 1.), (k, 1.)):
            np.add.at(out, index, sign*flux)
        return out
