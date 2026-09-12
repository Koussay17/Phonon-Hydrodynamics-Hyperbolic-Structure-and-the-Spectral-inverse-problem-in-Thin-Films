"""Unit tests for src/forward_model.py.

The strongest test is `test_uniform_sample_is_semi_infinite`: when film and
substrate share the same properties and the contact is perfect, the whole
stack is a homogeneous semi-infinite medium, whose response is known exactly.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import quadrupoles as q  # noqa: E402


FREQ = np.logspace(0, 7, 40)


def default_sample(**kw):
    return fm.Sample(**kw)


# --------------------------------------------------------------------------
# Sample description
# --------------------------------------------------------------------------

def test_derived_quantities_are_consistent():
    s = default_sample()
    assert s.film_a == pytest.approx(s.film_lam / s.film_rho_c)
    assert s.film_b == pytest.approx(np.sqrt(s.film_lam * s.film_rho_c))
    assert s.xi1 == pytest.approx(s.thickness / np.sqrt(s.film_a))
    assert s.xi1 ** 2 == pytest.approx(s.diffusion_time)
    assert not s.is_graded


def test_graded_flag():
    s = default_sample(film_lam_back=120.0, film_rho_c_back=2.41e6)
    assert s.is_graded
    assert s.film_b_back == pytest.approx(np.sqrt(120.0 * 2.41e6))


# --------------------------------------------------------------------------
# The reference case
# --------------------------------------------------------------------------

def test_uniform_sample_is_semi_infinite():
    """Film and substrate identical, perfect contact: the stack is a single
    homogeneous semi-infinite medium, so theta = 1 / (b sqrt(p)) exactly."""
    s = default_sample(sub_lam=60.0, sub_rho_c=2.41e6)
    p = 2j * np.pi * FREQ
    got = fm.response(p, s)
    expected = 1.0 / (s.film_b * np.sqrt(p))
    assert np.allclose(got, expected, rtol=1e-9)


def test_uniform_sample_phase_is_minus_45():
    s = default_sample(sub_lam=60.0, sub_rho_c=2.41e6)
    _, phase = fm.modulated_response(FREQ, s)
    assert np.allclose(phase, -45.0, atol=1e-8)


# --------------------------------------------------------------------------
# Interface resistance
# --------------------------------------------------------------------------

def test_zero_resistance_changes_nothing():
    p = 2j * np.pi * FREQ
    a = fm.response(p, default_sample(contact_resistance=0.0))
    b = fm.response(p, default_sample())
    assert np.allclose(a, b, rtol=1e-12)


def test_resistance_matrix_is_unimodular():
    m = fm.interface_resistance(2j * np.pi * FREQ, 1e-8)
    det = m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0]
    assert np.allclose(det, 1.0)


def test_resistance_raises_front_temperature():
    """An interface resistance impedes evacuation towards the substrate, so
    the front face runs hotter, never cooler."""
    p = 2j * np.pi * FREQ
    cold = np.abs(fm.response(p, default_sample(contact_resistance=0.0)))
    hot = np.abs(fm.response(p, default_sample(contact_resistance=1e-7)))
    assert np.all(hot >= cold * (1.0 - 1e-12))
    assert hot[0] > cold[0] * 1.001


def test_resistance_becomes_invisible_once_the_wave_is_confined():
    """The interface can only be seen while the thermal wave reaches it.

    For the default 500 nm film the characteristic frequency is about 16 MHz,
    so confinement requires approaching the gigahertz. This sets a hard
    requirement on the measurement bandwidth: below it, the interface
    resistance and the substrate cannot be separated from the film.
    """
    s = default_sample()
    f_hi = 200.0 / (2.0 * np.pi * s.diffusion_time)      # |sqrt(p) xi1| ~ 14
    p = 2j * np.pi * f_hi
    cold = np.abs(fm.response(p, default_sample(contact_resistance=0.0)))
    hot = np.abs(fm.response(p, default_sample(contact_resistance=1e-7)))
    assert hot == pytest.approx(cold, rel=1e-6)


# --------------------------------------------------------------------------
# Asymptotic behaviour
# --------------------------------------------------------------------------

def test_apparent_effusivity_limits():
    """Tends to the film effusivity at high frequency, to the substrate
    effusivity at low frequency."""
    s = default_sample()
    wide = np.logspace(-4, 10, 200)
    b_app = fm.apparent_effusivity(wide, s)
    assert b_app[-1] == pytest.approx(s.film_b, rel=1e-3)
    assert b_app[0] == pytest.approx(s.sub_b, rel=1e-2)


def test_apparent_effusivity_follows_the_contrast():
    """The direction of variation is set by the sign of the contrast. With the
    default values the sapphire substrate is less effusive than the film, so
    the apparent effusivity rises with frequency.
    """
    s = default_sample()
    assert s.sub_b < s.film_b
    b_app = fm.apparent_effusivity(np.logspace(-2, 9, 300), s)
    assert b_app[-1] > b_app[0]

    s_rev = default_sample(sub_lam=150.0, sub_rho_c=3.03e6)
    assert s_rev.sub_b > s_rev.film_b
    b_rev = fm.apparent_effusivity(np.logspace(-2, 9, 300), s_rev)
    assert b_rev[-1] < b_rev[0]


def test_apparent_effusivity_shows_undulations():
    """Near the transition the curve is not monotonic.

    These undulations are the interference of the thermal wave between the
    front face and the interface, reported by Krapez in the range of
    normalised frequency 0.1 to 1. They are a physical feature, not a
    numerical artefact.
    """
    s = default_sample()
    f = np.logspace(-2, 10, 600)
    b_app = fm.apparent_effusivity(f, s)
    d = np.diff(b_app)
    assert np.any(d > 0) and np.any(d < 0)


# --------------------------------------------------------------------------
# Pulsed regime
# --------------------------------------------------------------------------

def test_pulsed_uniform_sample_matches_closed_form():
    """Uniform stack under a pulse: theta(t) = energy / (b sqrt(pi t))."""
    s = default_sample(sub_lam=60.0, sub_rho_c=2.41e6)
    t = np.logspace(-9, -5, 20)
    got = fm.pulsed_response(t, s)
    expected = 1.0 / (s.film_b * np.sqrt(np.pi * t))
    assert np.allclose(got, expected, rtol=1e-5)


def test_pulsed_response_decreases():
    s = default_sample()
    t = np.logspace(-10, -6, 30)
    got = fm.pulsed_response(t, s)
    assert np.all(np.diff(got) < 0.0)


# --------------------------------------------------------------------------
# Noise
# --------------------------------------------------------------------------

def test_noise_is_reproducible():
    amp, ph = fm.modulated_response(FREQ, default_sample())
    a1, p1 = fm.add_noise(amp, ph, rng=42)
    a2, p2 = fm.add_noise(amp, ph, rng=42)
    assert np.array_equal(a1, a2) and np.array_equal(p1, p2)


def test_noise_has_requested_magnitude():
    amp, ph = fm.modulated_response(np.logspace(0, 7, 4000), default_sample())
    na, npz = fm.add_noise(amp, ph, relative_amplitude=0.02,
                           absolute_phase_deg=0.3, rng=0)
    assert np.std((na - amp) / amp) == pytest.approx(0.02, rel=0.1)
    assert np.std(npz - ph) == pytest.approx(0.3, rel=0.1)


def test_noise_amplitude_only():
    amp, _ = fm.modulated_response(FREQ, default_sample())
    out = fm.add_noise(amp, rng=1)
    assert out.shape == amp.shape


# --------------------------------------------------------------------------
# Effusivity contrast
# --------------------------------------------------------------------------

def test_reflection_coefficient_definition():
    s = default_sample()
    r = s.sub_b / s.film_b
    assert s.effusivity_ratio == pytest.approx(r)
    assert s.reflection_coefficient == pytest.approx((1.0 - r) / (1.0 + r))
    assert s.blind_film_effusivity == pytest.approx(s.sub_b)


def test_matched_effusivity_makes_the_interface_invisible():
    """When film and substrate share the same effusivity, the front-face
    response is exactly that of a semi-infinite medium, whatever the thickness
    and whatever the diffusivity contrast.

    This is the numerical counterpart of Krapez's result that the response
    depends on the effusivity alone once expressed in the Liouville
    coordinate. Here the diffusivities differ by a factor of 1.6 and the
    thickness spans two decades, yet nothing of either is visible.
    """
    b = 10298.0
    rc_f, rc_s = 2.41e6, 3.03e6
    lam_f, lam_s = b ** 2 / rc_f, b ** 2 / rc_s

    p = 2j * np.pi * np.logspace(3, 10, 20)
    ref = 1.0 / (b * np.sqrt(p))

    for e in (50e-9, 500e-9, 5e-6):
        s = default_sample(film_lam=lam_f, film_rho_c=rc_f, thickness=e,
                           sub_lam=lam_s, sub_rho_c=rc_s)
        assert s.reflection_coefficient == pytest.approx(0.0, abs=1e-12)
        assert np.allclose(fm.response(p, s), ref, rtol=1e-12)


def test_uncertainty_diverges_at_the_blind_point():
    """The blind point is approached continuously: the problem is merely
    ill-conditioned nearby and exactly non-identifiable at the point itself.

    Measured on AlN-like values over sapphire, the blind conductivity sits at
    44 W/(m K), squarely inside the plausible range for a thin film.
    """
    import inversion as inv

    rc, sub_l, sub_rc = 2.41e6, 35.0, 3.03e6
    blind = q.effusivity(sub_l, sub_rc) ** 2 / rc

    def sigma(lam):
        s = default_sample(film_lam=lam, film_rho_c=rc,
                           sub_lam=sub_l, sub_rho_c=sub_rc)
        fc = np.log10(s.characteristic_frequency)
        f = np.logspace(fc, fc + 1, 60)
        _, cov = inv.fisher_analysis(f, s, ["film_lam", "film_rho_c"], 0.01, 0.1)
        return np.sqrt(np.diag(cov)).max()

    assert blind == pytest.approx(44.0, rel=1e-3)
    assert sigma(blind) > 100.0 * sigma(blind * 1.25)
    assert sigma(blind * 1.01) > sigma(blind * 1.10)


def test_contrast_report_warns_when_blind():
    rc, sub_l, sub_rc = 2.41e6, 35.0, 3.03e6
    blind = q.effusivity(sub_l, sub_rc) ** 2 / rc
    near = fm.contrast_report(default_sample(film_lam=blind, film_rho_c=rc,
                                             sub_lam=sub_l, sub_rho_c=sub_rc))
    far = fm.contrast_report(default_sample(film_lam=250.0, film_rho_c=rc,
                                            sub_lam=sub_l, sub_rho_c=sub_rc))
    assert "WARNING" in near
    assert "WARNING" not in far
