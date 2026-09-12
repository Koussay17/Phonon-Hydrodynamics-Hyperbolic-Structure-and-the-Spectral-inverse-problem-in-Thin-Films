"""Unit tests for src/multilayer.py.

The decisive tests are the consistency ones: a stack must reproduce the
single-film model when it has one layer, and splitting a layer must not
change the result. Together they check the matrix, the ordering and the
interleaving of the contact resistances.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import multilayer as ml  # noqa: E402
import quadrupoles as q  # noqa: E402

FREQ = np.logspace(3, 10, 40)
P = 2j * np.pi * FREQ


def film():
    return ml.Layer(lam=60.0, rho_c=2.41e6, thickness=500e-9)


# --------------------------------------------------------------------------
# Consistency with the single-film model
# --------------------------------------------------------------------------

def test_single_layer_matches_sample():
    st = ml.Stack([film()], sub_lam=35.0, sub_rho_c=3.03e6)
    sa = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9,
                   sub_lam=35.0, sub_rho_c=3.03e6)
    assert np.allclose(ml.response(P, st), fm.response(P, sa), rtol=1e-12)


def test_single_layer_with_resistance_matches_sample():
    r = 1e-8
    st = ml.Stack([film()], contact_resistances=[r])
    sa = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9,
                   contact_resistance=r)
    assert np.allclose(ml.response(P, st), fm.response(P, sa), rtol=1e-12)


def test_graded_layer_matches_sample():
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 500e-9,
                            lam_back=120.0, rho_c_back=2.41e6)])
    sa = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9,
                   film_lam_back=120.0, film_rho_c_back=2.41e6)
    assert np.allclose(ml.response(P, st), fm.response(P, sa), rtol=1e-12)


# --------------------------------------------------------------------------
# Internal consistency
# --------------------------------------------------------------------------

def test_splitting_a_layer_changes_nothing():
    """One layer of thickness e, or two of thickness e/2 in perfect contact,
    must give the same matrix."""
    whole = ml.Stack([film()])
    halves = ml.Stack([ml.Layer(60.0, 2.41e6, 250e-9),
                       ml.Layer(60.0, 2.41e6, 250e-9)])
    assert np.allclose(ml.stack_matrix(P, whole), ml.stack_matrix(P, halves),
                       rtol=1e-10)


def test_ten_identical_layers_equal_one_thick_layer():
    thin = ml.Stack([ml.Layer(60.0, 2.41e6, 50e-9) for _ in range(10)])
    thick = ml.Stack([film()])
    assert np.allclose(ml.response(P, thin), ml.response(P, thick), rtol=1e-9)


def test_uniform_stack_is_semi_infinite():
    """Every layer and the substrate sharing the same properties: the whole
    stack is one homogeneous semi-infinite medium."""
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 200e-9),
                   ml.Layer(60.0, 2.41e6, 300e-9)],
                  sub_lam=60.0, sub_rho_c=2.41e6)
    b = q.effusivity(60.0, 2.41e6)
    assert np.allclose(ml.response(P, st), 1.0 / (b * np.sqrt(P)), rtol=1e-9)


def test_determinant_of_stack_is_one():
    """Unimodularity is preserved by the product, as long as the entries stay
    within reach of double precision."""
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 200e-9),
                   ml.Layer(150.0, 1.66e6, 100e-9),
                   ml.Layer(30.0, 3.0e6, 400e-9)],
                  contact_resistances=[1e-9, 2e-9, 5e-9])
    p = 2j * np.pi * np.logspace(3, 8, 40)
    m = ml.stack_matrix(p, st)
    det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
    assert np.allclose(det, 1.0, rtol=1e-8)


def test_determinant_degrades_at_very_high_frequency():
    """Documented limitation of the transfer-matrix formalism.

    The entries grow like cosh(sqrt(p) xi), so a thick stack at high frequency
    produces numbers around 1e40. The determinant, being a difference of such
    numbers, then loses all its significant digits to cancellation.

    The consequence is practical: unimodularity cannot be used as a health
    check outside the range where the entries stay moderate, and a scaled
    formulation becomes necessary if that range must be exceeded.
    """
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 200e-9),
                   ml.Layer(150.0, 1.66e6, 100e-9),
                   ml.Layer(30.0, 3.0e6, 400e-9)])
    p = 2j * np.pi * 1e10
    m = ml.stack_matrix(p, st)
    det = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]

    assert np.max(np.abs(m)) > 1e20          # the entries have exploded
    assert abs(det - 1.0) > 1e-6             # the determinant is no longer usable


def test_layer_order_matters():
    """Two different materials do not commute: swapping them changes the
    front-face response, which is what makes depth profiling possible at all.
    """
    a = ml.Layer(60.0, 2.41e6, 300e-9)
    b = ml.Layer(150.0, 1.66e6, 300e-9)
    r_ab = ml.response(P, ml.Stack([a, b]))
    r_ba = ml.response(P, ml.Stack([b, a]))
    # atol must be zero: the responses themselves are around 1e-10, far below
    # the default absolute tolerance, which would declare anything equal.
    assert not np.allclose(r_ab, r_ba, rtol=1e-3, atol=0.0)
    assert np.max(np.abs(r_ab - r_ba) / np.abs(r_ab)) > 0.05


# --------------------------------------------------------------------------
# Derived quantities and validation
# --------------------------------------------------------------------------

def test_liouville_thicknesses_add():
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 200e-9),
                   ml.Layer(30.0, 3.0e6, 300e-9)])
    expected = sum(l.xi for l in st.layers) ** 2
    assert st.total_diffusion_time == pytest.approx(expected)
    assert st.total_thickness == pytest.approx(500e-9)


def test_empty_stack_rejected():
    with pytest.raises(ValueError):
        ml.Stack([])


def test_wrong_number_of_resistances_rejected():
    with pytest.raises(ValueError):
        ml.Stack([film(), film()], contact_resistances=[0.0])


def test_pulsed_uniform_stack_matches_closed_form():
    st = ml.Stack([ml.Layer(60.0, 2.41e6, 200e-9),
                   ml.Layer(60.0, 2.41e6, 300e-9)],
                  sub_lam=60.0, sub_rho_c=2.41e6)
    b = q.effusivity(60.0, 2.41e6)
    t = np.logspace(-10, -6, 20)
    got = ml.pulsed_response(t, st)
    assert np.allclose(got, 1.0 / (b * np.sqrt(np.pi * t)), rtol=1e-5)
