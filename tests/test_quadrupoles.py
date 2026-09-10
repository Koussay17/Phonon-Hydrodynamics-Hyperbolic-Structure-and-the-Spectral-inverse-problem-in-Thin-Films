"""Unit tests for src/quadrupoles.py.

Run from the repository root:

    python -m pytest tests/ -v

Test 4 is the most important one: it validates the matrix and the composition
rule at the same time. Test 6 validates the matrix, the impedance and the
front-face formula in a single analytical identity.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import quadrupoles as q  # noqa: E402


# Reference material: dense alumina, order of magnitude only.
LAM = 30.0          # W / (m K)
RHO_C = 3.0e6       # J / (m^3 K)
E = 5.0e-6          # m

# Harmonic regime over a decade range that stays clear of overflow.
OMEGA = 2.0 * np.pi * np.logspace(1, 5, 40)
P = 1j * OMEGA


def test_property_round_trip():
    """(lambda, rho c) -> (a, b) -> (lambda, rho c) must be the identity."""
    a = q.diffusivity(LAM, RHO_C)
    b = q.effusivity(LAM, RHO_C)
    lam2, rho_c2 = q.from_diffusivity_effusivity(a, b)
    assert lam2 == pytest.approx(LAM, rel=1e-12)
    assert rho_c2 == pytest.approx(RHO_C, rel=1e-12)


def test_lambda_k_equals_b_sqrt_p():
    """The identity lambda * k = b sqrt(p) underpins both parameterisations."""
    a = q.diffusivity(LAM, RHO_C)
    b = q.effusivity(LAM, RHO_C)
    k = np.sqrt(P / a)
    assert np.allclose(LAM * k, b * np.sqrt(P), rtol=1e-12)


def test_determinant_is_one():
    """Necessary but not sufficient: the determinant does not detect a
    convention mismatch on the cross-section factor."""
    m = q.homogeneous_wall(P, LAM, RHO_C, E)
    det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
    assert np.allclose(det, 1.0, rtol=1e-10, atol=1e-12)


def test_entries_match_analytical_expression():
    """Term-by-term comparison against the closed form, computed independently."""
    a = q.diffusivity(LAM, RHO_C)
    b = q.effusivity(LAM, RHO_C)
    ke = np.sqrt(P / a) * E
    bsp = b * np.sqrt(P)

    m = q.homogeneous_wall(P, LAM, RHO_C, E)
    assert np.allclose(m[..., 0, 0], np.cosh(ke), rtol=1e-12)
    assert np.allclose(m[..., 1, 1], np.cosh(ke), rtol=1e-12)
    assert np.allclose(m[..., 0, 1], np.sinh(ke) / bsp, rtol=1e-12)
    assert np.allclose(m[..., 1, 0], bsp * np.sinh(ke), rtol=1e-12)


def test_composition_of_two_identical_layers():
    """Two slabs of thickness e in series must equal one slab of thickness 2e.

    This validates the matrix and the composition rule together.
    """
    m1 = q.homogeneous_wall(P, LAM, RHO_C, E)
    m_double = q.homogeneous_wall(P, LAM, RHO_C, 2.0 * E)
    assert np.allclose(m1 @ m1, m_double, rtol=1e-10, atol=1e-12)


def test_zero_thickness_is_identity():
    """A slab of vanishing thickness must not alter the state vector."""
    m = q.homogeneous_wall(P, LAM, RHO_C, 0.0)
    eye = np.broadcast_to(np.eye(2, dtype=complex), m.shape)
    assert np.allclose(m, eye, rtol=1e-12, atol=1e-14)


def test_semi_infinite_closure_is_exact():
    """A slab backed by a substrate of the same effusivity is a semi-infinite
    medium, whose front-face temperature is exactly power / (b sqrt(p)).

    Validates the matrix, the impedance and the front-face formula at once.
    """
    b = q.effusivity(LAM, RHO_C)
    m = q.homogeneous_wall(P, LAM, RHO_C, E)
    z = q.semi_infinite_impedance(P, b)

    theta = q.front_face_temperature(m, z, power=1.0, h=0.0)
    expected = 1.0 / (b * np.sqrt(P))
    assert np.allclose(theta, expected, rtol=1e-10)


def test_semi_infinite_phase_is_minus_45_degrees():
    """The reference result of photothermal metrology."""
    b = q.effusivity(LAM, RHO_C)
    theta = 1.0 / (b * np.sqrt(P))
    phase = np.degrees(np.angle(theta))
    assert np.allclose(phase, -45.0, atol=1e-10)


def test_overflow_guard_triggers():
    """A layer far too thick for the frequency range must raise, not silently
    return inf or nan."""
    with pytest.raises(OverflowError):
        q.homogeneous_wall(1j * 1e12, LAM, RHO_C, 1.0)


def test_scalar_input_returns_2x2():
    m = q.homogeneous_wall(1j * 100.0, LAM, RHO_C, E)
    assert m.shape == (2, 2)


# ==========================================================================
# Graded layer with a linear metaproperty profile
# ==========================================================================

B0 = q.effusivity(LAM, RHO_C)
B1 = 3.0 * B0                      # effusivity triples across the layer
XI1 = q.xi_from_thickness(E, q.diffusivity(LAM, RHO_C))


def test_graded_reduces_to_homogeneous():
    """With equal effusivities on both faces, equation (27) must give back the
    homogeneous wall exactly. This is the unit test of reference."""
    m_graded = q.graded_linear_layer(P, B0, B0, XI1, form="T")
    m_wall = q.homogeneous_wall(P, LAM, RHO_C, E)
    assert np.allclose(m_graded, m_wall, rtol=1e-10, atol=1e-14)


def test_graded_determinant_is_one():
    for form in ("T", "phi"):
        m = q.graded_linear_layer(P, B0, B1, XI1, form=form)
        det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
        assert np.allclose(det, 1.0, rtol=1e-9, atol=1e-11), form


def test_graded_layer_splits_consistently():
    """A linear profile split in two remains linear on each half, so the
    product of the two sub-layers must equal the whole layer.

    This validates the matrix and the composition rule together.
    """
    s0, s1 = np.sqrt(B0), np.sqrt(B1)
    b_mid = ((s0 + s1) / 2.0) ** 2

    left = q.graded_linear_layer(P, B0, b_mid, XI1 / 2.0, form="T")
    right = q.graded_linear_layer(P, b_mid, B1, XI1 / 2.0, form="T")
    whole = q.graded_linear_layer(P, B0, B1, XI1, form="T")

    assert np.allclose(left @ right, whole, rtol=1e-9, atol=1e-11)


def test_two_forms_are_distinct():
    """Krapez section 4: two distinct dual profiles can never give the same
    thermal response. The matrices must therefore differ."""
    m_t = q.graded_linear_layer(P, B0, B1, XI1, form="T")
    m_phi = q.graded_linear_layer(P, B0, B1, XI1, form="phi")
    assert not np.allclose(m_t, m_phi, rtol=1e-6)


def test_profile_endpoints():
    """The profile must reach the prescribed effusivities at both faces."""
    for form in ("T", "phi"):
        assert q.linear_effusivity_profile(0.0, B0, B1, XI1, form) == pytest.approx(B0)
        assert q.linear_effusivity_profile(XI1, B0, B1, XI1, form) == pytest.approx(B1)


def test_profile_is_monotonic():
    xi = np.linspace(0.0, XI1, 200)
    for form in ("T", "phi"):
        b = q.linear_effusivity_profile(xi, B0, B1, XI1, form)
        assert np.all(np.diff(b) > 0.0), form


def test_graded_invalid_form_raises():
    with pytest.raises(ValueError):
        q.graded_linear_layer(P, B0, B1, XI1, form="x")
