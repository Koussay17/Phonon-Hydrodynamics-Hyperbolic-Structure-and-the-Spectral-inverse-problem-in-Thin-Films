"""Finite-matrix Schur elimination audit. Run: py second_generation_schur.py.

Writes only second_generation_results.json and second_generation_console.txt
beside this script. All rates are dimensionless. No AlN data enter the tests.
"""
from __future__ import annotations

import itertools
import json
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la

OUT = Path(__file__).resolve().parent
SEED = 20260917
MACHEPS = np.finfo(float).eps


def norm(a):
    return float(la.norm(a, 2))


def minherm(a):
    return float(la.eigvalsh((a + a.conj().T) / 2)[0])


def model(kscale):
    """PSD Schur construction, exact invariant, diagonal mode velocities."""
    rng = np.random.default_rng(SEED)
    n, p = 8, 3
    e = np.ones(n) / np.sqrt(n)
    trial = rng.normal(size=(n, n))
    trial[:, 0] = e
    t, _ = la.qr(trial)
    if t[:, 0] @ e < 0:
        t[:, 0] *= -1
    od, _ = la.qr(rng.normal(size=(n - p, n - p)))
    d = od @ np.diag([0.5, 0.7, 1.1, 1.6, 2.0]) @ od.T
    coupling = kscale * rng.normal(size=(n - p, p))
    coupling[:, 0] = 0
    a = np.diag([0, 0.025, 0.08]) + coupling.T @ la.solve(d, coupling)
    cb = np.block([[a, coupling.T], [coupling, d]])
    c = t @ cb @ t.T
    c = (c + c.T) / 2
    v = np.diag([-1.7, -1.1, -0.6, -0.2, 0.2, 0.6, 1.1, 1.7])
    # Re-extract from the matrices used by the direct full solve.
    cb, vb = t.T @ c @ t, t.T @ v @ t
    return dict(n=n, p=p, T=t, C=c, V=v, e=e, A=cb[:p, :p], K=cb[p:, :p],
                D=cb[p:, p:], U1=vb[:p, :p], B1=vb[p:, :p], W1=vb[p:, p:], kscale=kscale)


def evaluate(m, z, k):
    n, p, t = m['n'], m['p'], m['T']
    c, v, a, coupling, d = m['C'], m['V'], m['A'], m['K'], m['D']
    u, b, w = k*m['U1'], k*m['B1'], k*m['W1']
    gamma = float(la.eigvalsh(d)[0])
    di = la.solve(d, np.eye(n-p))
    rq = la.solve(z*np.eye(n-p) + d + 1j*w, np.eye(n-p))
    x, y = coupling + 1j*b, coupling.T + 1j*b.T
    s = z*np.eye(p) + a + 1j*u - y@rq@x
    s0 = z*np.eye(p) + a + 1j*u - y@di@x
    a0 = a - coupling.T@di@coupling
    u0 = u - coupling.T@di@b - b.T@di@coupling
    s0_expanded = z*np.eye(p) + a0 + 1j*u0 + b.T@di@b
    # Resolvent identity avoids subtractive cancellation in the small error.
    zp = z*np.eye(n-p) + 1j*w
    delta = y@rq@zp@di@x
    b5 = (norm(coupling)+norm(b))**2 * (abs(z)+norm(w)) / (gamma*(gamma+z.real))
    b3 = (abs(z)+norm(w)) / (gamma*(gamma+z.real))
    main = z*np.eye(n) + c + 1j*k*v
    direct = la.solve(main, np.eye(n))
    slow = la.solve(s, np.concatenate([np.eye(p), -y@rq], axis=1))
    fast_id = np.concatenate([np.zeros((n-p, p)), np.eye(n-p)], axis=1)
    fast = rq@(fast_id-x@slow)
    reduced = t@np.concatenate([slow, fast], axis=0)@t.T
    pdirect, preduced = t[:, :p].T@direct, slow@t.T
    si, s0i = la.solve(s, np.eye(p)), la.solve(s0, np.eye(p))
    invdiff, sdiff = norm(si-s0i), norm(delta)
    condfull, conds = float(np.linalg.cond(main)), float(np.linalg.cond(s))
    err = norm(reduced-direct)/norm(direct)
    source_error = y@rq@zp@di
    return dict(kscale=m['kscale'], z=[z.real, z.imag], k=k, gamma=gamma,
        condition_full=condfull, condition_schur=conds,
        full_reconstruction_relative_error=err,
        projected_resolvent_relative_error=norm(pdirect-preduced)/norm(pdirect),
        full_backward_residual=norm(main@reduced-np.eye(n))/(norm(main)*norm(reduced)+1),
        direct_full_backward_residual=norm(main@direct-np.eye(n))/(norm(main)*norm(direct)+1),
        forward_error_over_eps_condition=err/(MACHEPS*max(condfull, conds)),
        S0_expansion_discrepancy=norm(s0-s0_expanded),
        delta_difference_roundoff=norm(delta-(s-s0)), delta_norm=sdiff,
        A5_bound=b5, A5_ratio=sdiff/b5 if b5 else 0,
        A3_ratio=norm(rq@zp@di)/b3,
        A6_actual_inverse_difference=invdiff, A6_bound=sdiff/z.real**2,
        A6_ratio=invdiff*z.real**2/sdiff if sdiff else 0,
        S_inverse_times_sigma=norm(si)*z.real, S0_inverse_times_sigma=norm(s0i)*z.real,
        source_correction_ratio=norm(source_error)/(norm(y)*b3),
        S0_generator_min_hermitian_eigenvalue=minherm(s0-z*np.eye(p)),
        S_min_hermitian_eigenvalue_minus_sigma=minherm(s)-z.real,
        perturbation_eta_actual=norm(s0i)*sdiff,
        local_inverse_relative_error=invdiff/norm(si))


def random_suite():
    rows, structure = [], []
    for scale in [0, 0.01, 0.2, 1, 5]:
        m = model(scale)
        c, v = m['C'], m['V']
        structure.append(dict(kscale=scale, dimension=m['n'], slow_dimension=m['p'],
            C_eigenvalues=la.eigvalsh(c).tolist(),
            invariant_relative_residual=norm((c@m['e']).reshape(-1, 1))/norm(c),
            streaming_commutator_norm=norm(c@v-v@c), K_norm=norm(m['K']),
            D_min_eigenvalue=float(la.eigvalsh(m['D'])[0])))
        for sigma, omega, k in itertools.product(
            [1e-7, 1e-3, 0.1, 1], [-2, -0.3, 0, 0.4, 1.5], [0, 1e-4, 0.02, 0.2, 1, 3]):
            rows.append(evaluate(m, complex(sigma, omega), k))
    keys = ['full_reconstruction_relative_error', 'projected_resolvent_relative_error',
        'full_backward_residual', 'direct_full_backward_residual', 'condition_full',
        'condition_schur', 'forward_error_over_eps_condition', 'A5_ratio', 'A3_ratio',
        'A6_ratio', 'S_inverse_times_sigma', 'S0_inverse_times_sigma', 'source_correction_ratio',
        'S0_expansion_discrepancy', 'delta_difference_roundoff']
    summary = {'count': len(rows), **{'max_'+key: max(r[key] for r in rows) for key in keys}}
    summary['minimum_static_hermitian_part'] = min(r['S0_generator_min_hermitian_eigenvalue'] for r in rows)
    summary['minimum_exact_schur_margin'] = min(r['S_min_hermitian_eigenvalue_minus_sigma'] for r in rows)
    moderate = [r for r in rows if r['condition_full'] < 1e6]
    summary['moderate_condition_count'] = len(moderate)
    summary['moderate_max_projected_error'] = max(r['projected_resolvent_relative_error'] for r in moderate)
    summary['worst_projected_case'] = max(rows, key=lambda r: r['projected_resolvent_relative_error'])
    return dict(summary=summary, models=structure, rows=rows)


def cubic_suite():
    m, rows = model(0), []
    for k in np.logspace(-1, -6, 11):
        r = evaluate(m, complex(0.3, 0.4)*k, float(k))
        row = dict(k=float(k), delta_norm=r['delta_norm'], A5_bound=r['A5_bound'],
                   scaled_error_delta_over_k3=r['delta_norm']/k**3)
        if rows:
            row['observed_power'] = float(np.log(row['delta_norm']/rows[-1]['delta_norm']) / np.log(k/rows[-1]['k']))
        rows.append(row)
    return rows


def poor_suite():
    mp.mp.dps = 90
    rows = []
    for etext in ['1e-2', '1e-4', '1e-6', '1e-8', '1e-10', '1e-12']:
        eps = mp.mpf(etext)
        for k in [mp.mpf(0), eps]:
            mat = mp.matrix([[1+eps+mp.j*k, 1], [1, 1-mp.j*k]])
            determinant, trace = eps+k*k-mp.j*eps*k, 2+eps
            slow = determinant / ((trace+mp.sqrt(trace**2-4*determinant))/2)
            direct_slow = min(mp.eig(mat, left=False, right=False), key=abs)
            z = eps/10
            full = mp.lu_solve(mat+z*mp.eye(2), mp.matrix([1, 0]))
            s = z+1+eps+mp.j*k-1/(z+1-mp.j*k)
            s0 = z+eps+mp.j*k
            b5, delta = (abs(z)+abs(k))/(1+z), abs(s-s0)
            rows.append(dict(epsilon=etext, k_over_epsilon=float(k/eps),
                full_slow_eigenvalue=[mp.nstr(mp.re(slow),25), mp.nstr(mp.im(slow),25)],
                static_eigenvalue=[mp.nstr(eps,25), mp.nstr(k,25)],
                static_to_exact_decay_rate=float(eps/mp.re(slow)),
                static_complex_eigenvalue_relative_error=float(abs(eps+mp.j*k-slow)/abs(slow)),
                mass_corrected_eigenvalue_relative_error=float(abs(eps/2-slow)/abs(slow)),
                mp_dense_eigenvalue_relative_error=float(abs(direct_slow-slow)/abs(slow)),
                mp_dense_vs_schur_relative_error=float(abs(full[0]-1/s)/abs(full[0])),
                causal_z_over_epsilon=0.1, delta_norm=float(delta), A5_bound=float(b5),
                A5_ratio=float(delta/b5), causal_static_response_relative_error=float(abs(1/s0-full[0])/abs(full[0])),
                perturbation_eta_actual=float(delta/abs(s0))))
    z = mp.mpf('0.3')
    return dict(rows=rows, A5_sharp_calibration=dict(z=0.3, k=0, actual_error=float(z/(1+z)), bound=float(z/(1+z))),
                streaming_commutator_norm_for_k1=2.0)


def near_pole_suite():
    mp.mp.dps = 110
    kf = 0.01
    km = mp.mpf(kf)  # Reference the same exact binary input as float64.
    rate = 2*km**2/(1+mp.sqrt(1-4*km**2))
    rows = []
    for dtext in ['1e-3', '1e-6', '1e-8', '1e-10', '1e-12', '1e-14', '1e-16', '1e-18', '1e-20']:
        zf = float(-rate+mp.mpf(dtext))
        zm = mp.mpf(zf)
        matmp = mp.matrix([[zm, mp.j*km], [mp.j*km, zm+1]])
        fullmp = mp.lu_solve(matmp, mp.matrix([1, 0]))
        smp, s0mp = zm+km**2/(1+zm), zm+km**2
        matnp = np.array([[zf, 1j*kf], [1j*kf, zf+1]], dtype=complex)
        fullnp = la.solve(matnp, np.eye(2))
        snp = zf+kf**2/(1+zf)
        full00, schur00 = mp.mpc(complex(fullnp[0,0])), mp.mpc(complex(1/snp if snp != 0 else np.inf))
        ref = fullmp[0]
        delta = abs(smp-s0mp)
        rows.append(dict(intended_distance_to_exact_pole=dtext,
            actual_float_input_distance_to_exact_pole=mp.nstr(zm+rate,25), z_float=zf,
            inside_stated_A5_half_plane=zf>=0, condition_full=float(np.linalg.cond(matnp)),
            full_float64_response_relative_error=float(abs(full00-ref)/abs(ref)),
            schur_float64_response_relative_error=float(abs(schur00-ref)/abs(ref)),
            mp_dense_schur_relative_error=float(abs(ref-1/smp)/abs(ref)),
            full_float64_backward_residual=norm(matnp@fullnp-np.eye(2))/(norm(matnp)*norm(fullnp)+1),
            static_response_relative_error=float(abs(1/s0mp-ref)/abs(ref)),
            operator_difference=float(delta), scalar_extended_bound=float(km**2*abs(zm)/(1+zm)),
            perturbation_eta_actual=float(delta/abs(s0mp))))
    return dict(k_float=kf, exact_slow_rate=mp.nstr(rate,35), static_slow_rate=kf*kf,
        pole_shift=mp.nstr(rate-km**2,35), wrong_Y_equals_X_adjoint_static_rate=-kf*kf, rows=rows)


def aligned_causal_suite():
    rows = []
    for k, alpha in itertools.product([0.1, 0.01, 0.001, 0.0001], [0, 0.1, 1]):
        z = alpha*k*k
        exact, static = z+k*k/(1+z), z+k*k
        rows.append(dict(k=k, z_over_k2=alpha, response_relative_error=abs(exact/static-1),
            analytic_relative_error=alpha*k*k/((1+alpha)*(1+alpha*k*k))))
    return rows


def main():
    results = dict(metadata=dict(seed=SEED, python=platform.python_version(), numpy=np.__version__,
        scipy=scipy.__version__, mpmath=mp.__version__, float_epsilon=MACHEPS,
        mpmath_digits_poor=90, mpmath_digits_near_pole=110,
        domain='One spatial Fourier mode of a finite accretive kinetic matrix',
        physical_status='Abstract operator tests; no AlN coefficients'),
        random_sweep=random_suite(), K_zero_scaling=cubic_suite(), poor_coordinate=poor_suite(),
        near_pole=near_pole_suite(), causal_aligned_control=aligned_causal_suite())
    (OUT/'second_generation_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    brief = {key: (val['summary'] if key=='random_sweep' else val) for key,val in results.items()}
    console = json.dumps(brief,indent=2)+'\n'
    (OUT/'second_generation_console.txt').write_text(console,encoding='utf-8')
    print(console)


if __name__ == '__main__':
    main()
