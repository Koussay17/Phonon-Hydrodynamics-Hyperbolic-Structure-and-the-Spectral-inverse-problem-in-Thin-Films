"""Independent population-dynamics and invariant checks for event actions."""
import mpmath as mp
import numpy as np
import pytest
from src.collision_events import DecayEventOperator


def pair():
    return DecayEventOperator([3., 1., 2., 3., 1., 2.],
                             [[0, 1, 2], [3, 4, 5]], [2., 2.])


def dense(op):
    return np.column_stack([op.action(v) for v in np.eye(len(op.energy_over_kbt))])


def test_nonlinear_population_derivative_matches_entropy_action():
    op = pair()
    n = 1/np.expm1(op.energy_over_kbt)
    d = op.entropy_scale
    direction = np.array([.3, -.2, .4, .1, -.3, .2])
    for h in [1e-2, 1e-3, 1e-4]:
        numerical = (op.population_rhs(n+h*d*direction)
                     - op.population_rhs(n-h*d*direction))/(2*h*d)
        # Cubic Bose terms cancel. The flux is quadratic: a central
        # derivative has only roundoff error, not a second-order truncation.
        np.testing.assert_allclose(numerical, -op.action(direction), atol=2e-10)
    np.testing.assert_allclose(op.population_rhs(n), 0, atol=2e-16)


@pytest.mark.parametrize("scale", [1e-6, 1., 100., 1000.])
def test_factor_against_high_precision_nonlinear_jacobian(scale):
    x = np.array([3., 1., 2.])*scale
    op = DecayEventOperator(x, [[0, 1, 2]], [2.])
    with mp.workdps(100):
        xm = list(map(mp.mpf, x))
        n = [1/mp.expm1(a) for a in xm]
        d = [mp.sqrt(a*(1+a)) for a in n]
        def rhs(y):
            v = [n[i]+d[i]*y[i] for i in range(3)]
            f = 2*(v[0]*(1+v[1])*(1+v[2])-(1+v[0])*v[1]*v[2])
            return [-f/d[0], f/d[1], f/d[2]]
        jac = mp.matrix([[mp.diff(lambda z: rhs(
            [z if k==j else mp.mpf(0) for k in range(3)])[i], 0)
                          for j in range(3)] for i in range(3)])
        expected = np.array(jac.tolist(), dtype=float)
    np.testing.assert_allclose(dense(op), -expected, rtol=2e-12, atol=1e-300)


def test_energy_normal_momentum_psd_and_reversal():
    op = pair()
    c = dense(op)
    e = op.energy_over_kbt*op.entropy_scale
    q = np.array([.3, .1, .2, -.3, -.1, -.2])
    np.testing.assert_allclose(c@e, 0, atol=8e-16)
    np.testing.assert_allclose(c@(q*op.entropy_scale), 0, atol=2e-16)
    np.testing.assert_allclose(c, c.T, atol=1e-16)
    assert np.linalg.eigvalsh(c).min() > -1e-15
    reversal = np.array([3, 4, 5, 0, 1, 2])
    np.testing.assert_allclose(c[np.ix_(reversal, reversal)], c, atol=0)
    y = np.arange(6.) + 1j*np.arange(6.)[::-1]
    np.testing.assert_allclose(op.as_linear_operator()@y, c@y)
    np.testing.assert_allclose(np.vdot(y, c@y).real,
                               np.linalg.norm(op.factor@y)**2)


def test_umklapp_retains_energy_but_not_crystal_momentum():
    op = pair()
    # Reciprocal coordinate unit G=1: q_parent-q_j-q_k=1.
    q = np.array([.4, -.3, -.3, -.4, .3, .3])
    assert np.linalg.norm(op.action(q*op.entropy_scale)) > .1
    np.testing.assert_allclose(op.action(op.energy_over_kbt*op.entropy_scale),
                               0, atol=8e-16)


def test_unequal_reciprocal_weights_break_parity_without_breaking_energy():
    op = DecayEventOperator([3., 1., 2., 3., 1., 2.],
                            [[0, 1, 2], [3, 4, 5]], [2., 3.])
    c = dense(op)
    j = np.eye(6)[[3, 4, 5, 0, 1, 2]]
    assert np.linalg.norm(c@j-j@c) > 1
    np.testing.assert_allclose(c@(op.energy_over_kbt*op.entropy_scale), 0, atol=1e-15)


def test_duplicate_listing_doubles_rate():
    a = DecayEventOperator([3., 1., 2.], [[0, 1, 2]], [2.])
    # Deliberate counting error: the reverse reaction is already included.
    b = DecayEventOperator([3., 1., 2.], [[0, 1, 2], [0, 2, 1]], [2., 2.])
    np.testing.assert_allclose(dense(b), 2*dense(a), atol=1e-15)


def test_zero_events_and_zero_weight():
    a = DecayEventOperator([1.], np.empty((0, 3), dtype=int), [])
    np.testing.assert_array_equal(a.action([4.]), [0.])
    b = DecayEventOperator([3., 1., 2.], [[0, 1, 2]], [0.])
    np.testing.assert_array_equal(b.action([1., 2., 3.]), [0., 0., 0.])


@pytest.mark.parametrize("x,events,gamma", [
    ([3.1, 1., 2.], [[0, 1, 2]], [1.]),
    ([2., 1.], [[0, 1, 1]], [1.]),
    ([0., 1., 1.], [[0, 1, 2]], [1.]),
    ([3., 1., 2.], [[0, 1, 3]], [1.]),
    ([3., 1., 2.], [[0., 1., 2.]], [1.]),
    ([3., 1., 2.], [[0, 1, 2]], [-1.]),
])
def test_out_of_contract_events_rejected(x, events, gamma):
    with pytest.raises(ValueError):
        DecayEventOperator(x, events, gamma)


def test_relative_resonance_alone_does_not_certify_detailed_balance():
    # Fractional mismatch is tiny but the Boltzmann factor ratio is exp(-100).
    with pytest.raises(ValueError, match="conserve energy"):
        DecayEventOperator([3e15+100., 1e15, 2e15], [[0, 1, 2]], [1.])


def test_linear_operator_column_and_multiple_right_hand_sides():
    op = pair()
    c = dense(op)
    linear = op.as_linear_operator()
    y = np.arange(6.)[:, None]
    np.testing.assert_allclose(linear@y, c@y)
    np.testing.assert_allclose(linear@np.eye(6), c)
    np.testing.assert_allclose(linear.H@np.eye(6), c)


def test_unrepresentable_entropy_scale_is_rejected():
    with pytest.raises(ValueError, match="entropy scale"):
        DecayEventOperator([3e-320, 1e-320, 2e-320], [[0, 1, 2]], [1.])


def test_extreme_low_occupation_limit_retains_prefactor():
    op = DecayEventOperator([3e16, 1e16, 2e16], [[0, 1, 2]], [2.])
    # Spontaneous decay survives as thermal daughter occupations vanish.
    np.testing.assert_allclose(dense(op), np.diag([2., 0., 0.]), rtol=2e-15)
