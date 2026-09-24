"""Bounded independent import check against the frozen DecayEventOperator.
Run with py -B second_import_check.py. Writes second_import_results.json only.
No implementation or first-pass file is edited; bytecode writes are disabled.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
sys.dont_write_bytecode = True
import numpy as np
import scipy
import scipy.linalg as la

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent
FROZEN = CAMPAIGN / 'review-inputs' / 'collision_events.py'
REFERENCE = HERE / 'C_parity_events.py'
EXPECTED_HASH = '8daa09d8f851da3b7823bdd5cc1ad86f98ced5892c814a380dcb58ce50f020c6'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def norm(a):
    return float(la.norm(a, 2))

def cpair(z):
    return [float(np.real(z)), float(np.imag(z))]

def matrix_action(op, n):
    eye = np.eye(n)
    return np.column_stack([op.action(eye[:, i]) for i in range(n)])

def rejected(thunk):
    try:
        thunk()
    except ValueError as error:
        return dict(rejected=True, exception_type=type(error).__name__, message=str(error))
    return dict(rejected=False, exception_type=None, message=None)

def main():
    before = dict(frozen=digest(FROZEN), independent_reference=digest(REFERENCE))
    assert before['frozen'] == EXPECTED_HASH
    repo = load('_second_frozen_collision_events', FROZEN)
    ref = load('_second_independent_c_reference', REFERENCE)
    adapted = np.array([(parent, a, b) for a, b, parent in ref.EVENTS], dtype=int)
    gamma = np.ones(len(adapted))
    op = repo.DecayEventOperator(ref.ENERGY, adapted, gamma)
    got = matrix_action(op, ref.N)
    expected = ref.C
    even_got, even_ref = ref.te.T@got@ref.te, ref.te.T@expected@ref.te
    odd_got, odd_ref = ref.to.T@got@ref.to, ref.to.T@expected@ref.to
    probe = np.eye(ref.N)+1j*np.arange(ref.N*ref.N).reshape(ref.N, ref.N)/31
    wrapper_result = op.as_linear_operator().matmat(probe)
    action_checks = dict(
        full_relative_error=norm(got-expected)/norm(expected),
        even_relative_error=norm(even_got-even_ref)/norm(even_ref),
        odd_relative_error=norm(odd_got-odd_ref)/norm(odd_ref),
        mixed_collision_block_norm=norm(ref.te.T@got@ref.to),
        entropy_scale_relative_error=norm(op.entropy_scale-ref.sw)/norm(ref.sw),
        complex_action_relative_error=norm(op.action(probe[:, 2])-expected@probe[:, 2])/norm(expected@probe[:, 2]),
        complex_linear_operator_relative_error=norm(wrapper_result-expected@probe)/norm(expected@probe),
        energy_relative_residual=norm(got@ref.e)/(norm(got)*norm(ref.e)),
        momentum_relative_residual=norm(got@ref.p)/(norm(got)*norm(ref.p)),
        eigenvalues=la.eigvalsh(got).tolist())

    finite_difference_rows = []
    directions = dict(even=ref.D[-1], odd=ref.to[:, 0]+.3*ref.to[:, 2],
                      generic=np.array([1., -2., .4, .7, -1.2, .6, .9]))
    h = 1e-3
    for name, y in directions.items():
        y = y/norm(y)
        plus, minus = ref.n0+h*ref.sw*y, ref.n0-h*ref.sw*y
        assert np.min(plus)>0 and np.min(minus)>0
        independent_rhs = ref.nonlinear(plus)
        frozen_rhs = op.population_rhs(plus)
        derivative = -(op.population_rhs(plus)-op.population_rhs(minus))/(2*h*ref.sw)
        finite_difference_rows.append(dict(direction=name, h=h,
            nonlinear_rhs_relative_error=norm(frozen_rhs-independent_rhs)/norm(independent_rhs),
            entropy_derivative_relative_error=norm(derivative-expected@y)/norm(expected@y)))

    response_rows = []
    samples = [(complex(.25), k) for k in [0., .03, .1, .3, 1.]]
    samples.append((complex(.25, .4), .3))
    for z, k in samples:
        full = z*np.eye(ref.N)+got+1j*k*ref.v
        base = z*np.eye(ref.N)+expected+1j*k*ref.v
        rr, rb = la.solve(full, np.eye(ref.N)), la.solve(base, np.eye(ref.N))
        er, eb = ref.te.T@rr@ref.te, ref.te.T@rb@ref.te
        ore, ob = ref.to.T@rr@ref.to, ref.to.T@rb@ref.to
        jr, jb = ref.j@rr@ref.j, ref.j@rb@ref.j
        tr, tb = ref.e@rr@ref.e, ref.e@rb@ref.e
        response_rows.append(dict(z=cpair(z), k=k,
            full_relative_error=norm(rr-rb)/norm(rb),
            even_block_relative_error=norm(er-eb)/norm(eb),
            odd_block_relative_error=norm(ore-ob)/norm(ob),
            parity_cross_absolute_error=norm(ref.te.T@(rr-rb)@ref.to),
            current_response=cpair(jr), independent_current_response=cpair(jb),
            current_response_relative_error=float(abs(jr-jb)/abs(jb)),
            energy_response=cpair(tr), independent_energy_response=cpair(tb),
            energy_response_relative_error=float(abs(tr-tb)/abs(tb)),
            normalized_equation_residual=norm(full@rr-np.eye(ref.N))/(norm(full)*norm(rr)+1),
            condition_full=float(np.linalg.cond(full))))

    # Deliberately add the self-reciprocal partner with daughters exchanged.
    duplicate_rows = np.vstack([adapted, adapted[-1, [0, 2, 1]]])
    duplicate_op = repo.DecayEventOperator(ref.ENERGY, duplicate_rows, np.ones(4))
    duplicate_matrix = matrix_action(duplicate_op, ref.N)
    difference = duplicate_matrix-got
    z, k = .25, 1.
    good_res = la.solve(z*np.eye(ref.N)+got+1j*k*ref.v, np.eye(ref.N))
    dup_res = la.solve(z*np.eye(ref.N)+duplicate_matrix+1j*k*ref.v, np.eye(ref.N))
    even_y = ref.D[-1]/norm(ref.D[-1])
    duplicate_control = dict(rows=duplicate_rows.tolist(),
        accepted_by_constructor=True,
        independent_double_model_relative_error=norm(duplicate_matrix-ref.CDOUBLE)/norm(ref.CDOUBLE),
        full_action_relative_change=norm(difference)/norm(got),
        even_block_relative_change=norm(ref.te.T@difference@ref.te)/norm(even_got),
        odd_block_absolute_change=norm(ref.to.T@difference@ref.to),
        even_direction_relative_change=norm(difference@even_y)/norm(got@even_y),
        z=z, k=k,
        energy_response_relative_change=float(abs(ref.e@(dup_res-good_res)@ref.e)/abs(ref.e@good_res@ref.e)),
        current_response_relative_change=float(abs(ref.j@(dup_res-good_res)@ref.j)/abs(ref.j@good_res@ref.j)))

    detuned = ref.ENERGY.copy()
    detuned[[2, 5]] += 1e-6
    negative_controls = dict(
        original_daughters_first_order=rejected(lambda: repo.DecayEventOperator(ref.ENERGY, np.array(ref.EVENTS), gamma)),
        detuned_reciprocal_pair=rejected(lambda: repo.DecayEventOperator(detuned, adapted, gamma)),
        repeated_daughter=rejected(lambda: repo.DecayEventOperator(np.array([1., 2.]), np.array([[1, 0, 0]]), np.ones(1))))
    negative_controls['detuning_added'] = 1e-6
    negative_controls['constructor_default_resonance_tolerance'] = 1e-12
    after = dict(frozen=digest(FROZEN), independent_reference=digest(REFERENCE))
    checks = dict(
        frozen_hash_matches_manifest=before['frozen']==EXPECTED_HASH,
        sources_unchanged=before==after,
        full_operator=action_checks['full_relative_error']<1e-13,
        both_parity_operators=max(action_checks['even_relative_error'],action_checks['odd_relative_error'])<1e-13,
        nonlinear_population=max(r['nonlinear_rhs_relative_error'] for r in finite_difference_rows)<1e-10,
        nonlinear_derivative=max(r['entropy_derivative_relative_error'] for r in finite_difference_rows)<1e-10,
        complex_action=max(action_checks['complex_action_relative_error'],action_checks['complex_linear_operator_relative_error'])<1e-13,
        full_and_both_parity_responses=max(max(r['full_relative_error'], r['even_block_relative_error'], r['odd_block_relative_error']) for r in response_rows)<1e-12,
        detuned_event_rejected=negative_controls['detuned_reciprocal_pair']['rejected'],
        repeated_index_explicitly_unsupported=negative_controls['repeated_daughter']['rejected'],
        wrong_adapter_order_rejected=negative_controls['original_daughters_first_order']['rejected'],
        duplicate_matches_independent_full_failure=duplicate_control['independent_double_model_relative_error']<1e-13,
        odd_only_check_misses_duplicate=duplicate_control['odd_block_absolute_change']<1e-13 and duplicate_control['even_direction_relative_change']>.9)
    result = dict(metadata=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        source_frozen=str(FROZEN), source_reference=str(REFERENCE), source_hash_before=before,
        source_hash_after=after, script_sha256=digest(Path(__file__)),
        scope='Frozen snapshot, distinct positive-frequency exact-resonance events, equal full-grid weights',
        bytecode_writes_disabled=sys.dont_write_bytecode),
        adapter=dict(original_daughters_daughters_parent=[list(t) for t in ref.EVENTS],
                     adapted_parent_daughters=adapted.tolist(), prefactors=gamma.tolist()),
        operator=action_checks, nonlinear=finite_difference_rows, responses=response_rows,
        duplicate_self_reciprocal_control=duplicate_control, negative_controls=negative_controls,
        checks=checks)
    (HERE/'second_import_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    assert all(checks.values()), checks

if __name__=='__main__':
    main()

