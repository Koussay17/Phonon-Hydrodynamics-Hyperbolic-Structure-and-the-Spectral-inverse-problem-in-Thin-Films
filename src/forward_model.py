"""Forward model of a film-on-substrate sample, and synthetic data generation.

Everything is parameterised. The default values are placeholders of the right
order of magnitude for an AlN film on sapphire; they carry no authority and
are meant to be replaced by measured values.

The stack is, from the heated front face inwards:

    front surface (optional heat losses h)
    film            homogeneous, or graded with a linear metaproperty profile
    interface       optional thermal contact resistance
    substrate       semi-infinite, entering through its impedance

Both regimes are available. The modulated regime returns amplitude and phase
against frequency; the pulsed regime returns temperature against time, through
the Stehfest inversion.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import laplace as lp
import quadrupoles as q

__all__ = [
    "Sample",
    "contrast_report",
    "interface_resistance",
    "stack_matrix",
    "response",
    "modulated_response",
    "pulsed_response",
    "apparent_effusivity",
    "add_noise",
]


# --------------------------------------------------------------------------
# Sample description
# --------------------------------------------------------------------------

@dataclass
class Sample:
    """Film on a semi-infinite substrate.

    All values are placeholders and must be replaced by laboratory data.

    Attributes
    ----------
    film_lam, film_rho_c : float
        Film conductivity and volumetric heat capacity. The defaults are a
        thin-film value, well below bulk AlN, since interfaces and
        microstructure reduce conductivity substantially.
    thickness : float
        Film thickness in metres.
    film_lam_back, film_rho_c_back : float or None
        If given, the film is treated as graded: properties vary between the
        front values and these back values. If None, the film is homogeneous.
    sub_lam, sub_rho_c : float
        Substrate properties. Only the effusivity enters, through the
        impedance of the semi-infinite medium.
    contact_resistance : float
        Thermal contact resistance at the film-substrate interface, in
        m^2 K / W. Zero means perfect contact.
    front_losses : float
        Heat exchange coefficient on the free front face, in W / (m^2 K).
        Zero is the adiabatic case.
    relaxation_time, sub_relaxation_time : float
        Relaxation times of the Cattaneo constitutive law, in seconds, for the
        film and the substrate. Zero recovers Fourier conduction.

        The law replaces the spectral parameter p by p (1 + tau p) in the
        layer concerned. Nothing else changes: the Liouville coordinate, the
        effusivity and the potential keep their Fourier definitions.
    """

    film_lam: float = 60.0
    film_rho_c: float = 2.41e6
    thickness: float = 500e-9

    film_lam_back: float | None = None
    film_rho_c_back: float | None = None

    sub_lam: float = 35.0
    sub_rho_c: float = 3.03e6

    contact_resistance: float = 0.0
    front_losses: float = 0.0

    relaxation_time: float = 0.0
    sub_relaxation_time: float = 0.0

    # ---- derived quantities -------------------------------------------

    @property
    def film_a(self) -> float:
        return q.diffusivity(self.film_lam, self.film_rho_c)

    @property
    def film_b(self) -> float:
        return q.effusivity(self.film_lam, self.film_rho_c)

    @property
    def film_b_back(self) -> float | None:
        if self.film_lam_back is None or self.film_rho_c_back is None:
            return None
        return q.effusivity(self.film_lam_back, self.film_rho_c_back)

    @property
    def sub_b(self) -> float:
        return q.effusivity(self.sub_lam, self.sub_rho_c)

    @property
    def xi1(self) -> float:
        """Liouville thickness of the film."""
        return q.xi_from_thickness(self.thickness, self.film_a)

    @property
    def diffusion_time(self) -> float:
        """Characteristic diffusion time through the film, thickness^2 / a."""
        return self.thickness ** 2 / self.film_a

    @property
    def relaxation_frequency(self) -> float:
        """Frequency at which omega tau reaches one, in hertz.

        The relaxation time only becomes measurable once the highest measured
        frequency reaches this value; below it the uncertainty grows like the
        inverse of the gap.
        """
        if self.relaxation_time <= 0.0:
            return np.inf
        return 1.0 / (2.0 * np.pi * self.relaxation_time)

    @property
    def effusivity_ratio(self) -> float:
        """Substrate effusivity divided by film effusivity.

        Written b32 in the photothermal literature. Together with the transit
        time it exhausts what the front-face measurement can determine about
        the coating.
        """
        return self.sub_b / self.film_b

    @property
    def reflection_coefficient(self) -> float:
        """Amplitude reflection coefficient of the thermal wave at the interface.

            Gamma = (1 - b32) / (1 + b32)

        The response can be written as a series of waves reflected between the
        front face and the interface, with Gamma as the reflection factor. It
        vanishes when film and substrate share the same effusivity, and the
        interface then becomes strictly invisible: the front-face response
        equals that of a semi-infinite medium, for any thickness and any
        diffusivity contrast.

        Krapez and Rigollet, arXiv:1708.07362 (2017), section on the
        coating-substrate response.
        """
        r = self.effusivity_ratio
        return (1.0 - r) / (1.0 + r)

    @property
    def blind_film_effusivity(self) -> float:
        """Film effusivity at which the interface becomes invisible.

        Equal to the substrate effusivity. A sample whose film sits near this
        value cannot yield its transit time, hence neither its thickness nor
        its diffusivity, whatever the inversion method.
        """
        return self.sub_b

    @property
    def characteristic_frequency(self) -> float:
        """Frequency at which the thermal wave just spans the film.

        Equal to 1 / (2 pi * diffusion_time). Below it the response is
        dominated by the substrate; well above it the wave is confined to the
        film and the interface becomes invisible. Confinement requires
        roughly two decades above this value.

        For a 500 nm film of the default properties this sits near 16 MHz,
        which is why thin-film thermal metrology relies on pump-probe methods
        rather than low-frequency modulation.
        """
        return 1.0 / (2.0 * np.pi * self.diffusion_time)

    @property
    def is_graded(self) -> bool:
        return self.film_b_back is not None


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------

def effective_p(p, tau: float):
    """Spectral parameter under the Cattaneo constitutive law.

        P = p (1 + tau p)

    With tau = 0 this is the identity, so Fourier conduction is the special
    case rather than a separate code path.

    In the modulated regime p = i omega, hence P = i omega - tau omega^2: a
    real part appears, growing as the square of the frequency. The dimension-
    less group governing the departure from Fourier is omega tau.
    """
    p = np.asarray(p, dtype=complex)
    return p * (1.0 + tau * p)


def interface_resistance(p, r: float) -> np.ndarray:
    """Quadrupole of a thermal contact resistance.

        [[1, r], [0, 1]]

    Temperature jumps by r times the flux; the flux itself is continuous.
    """
    p = np.asarray(p, dtype=complex)
    m = np.zeros(p.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = 1.0
    m[..., 0, 1] = r
    m[..., 1, 1] = 1.0
    return m if p.shape else m.reshape(2, 2)


def stack_matrix(p, s: Sample) -> np.ndarray:
    """Quadrupole of film plus interface, front face to substrate.

    The film is evaluated at its own effective spectral parameter, so a
    relaxation time in the film does not affect the substrate and conversely.
    """
    pf = effective_p(p, s.relaxation_time)

    if s.is_graded:
        m = q.graded_linear_layer(pf, s.film_b, s.film_b_back, s.xi1, form="T")
    else:
        m = q.homogeneous_wall(pf, s.film_lam, s.film_rho_c, s.thickness)

    if s.contact_resistance != 0.0:
        m = m @ interface_resistance(pf, s.contact_resistance)
    return m


def response(p, s: Sample, power=1.0):
    """Front-face temperature in the Laplace domain, per unit area."""
    m = stack_matrix(p, s)
    z = q.semi_infinite_impedance(effective_p(p, s.sub_relaxation_time), s.sub_b)
    return q.front_face_temperature(m, z, power=power, h=s.front_losses)


# --------------------------------------------------------------------------
# The two regimes
# --------------------------------------------------------------------------

def modulated_response(freq, s: Sample, power=1.0):
    """Amplitude and phase of the front-face temperature.

    Returns
    -------
    amplitude : ndarray
        In kelvin per unit power density.
    phase_deg : ndarray
        In degrees. The reference of a homogeneous semi-infinite medium is
        -45 degrees at every frequency.
    """
    freq = np.asarray(freq, dtype=float)
    theta = response(2j * np.pi * freq, s, power=power)
    return np.abs(theta), np.degrees(np.angle(theta))


def pulsed_response(t, s: Sample, energy=1.0, n_stehfest: int = 12):
    """Front-face temperature after a Dirac pulse of surface energy density.

    The inversion holds over roughly two decades of decay; see the note in
    `laplace.py`. Times far beyond that return numerical noise.
    """
    return lp.stehfest_inverse(
        lambda pp: np.real(response(pp, s, power=energy)), t, n=n_stehfest
    )


def apparent_effusivity(freq, s: Sample):
    """Frequency-dependent apparent effusivity, in absolute units.

    Defined from the inverse amplitude contrast with respect to a homogeneous
    semi-infinite medium having the film effusivity. Tends to the film
    effusivity at high frequency and to the substrate effusivity at low
    frequency, provided the contact is perfect.
    """
    freq = np.asarray(freq, dtype=float)
    amp, _ = modulated_response(freq, s)
    ref = 1.0 / (s.film_b * np.sqrt(2.0 * np.pi * freq))
    return s.film_b * ref / amp


# --------------------------------------------------------------------------
# Noise
# --------------------------------------------------------------------------

def add_noise(amplitude, phase_deg=None, relative_amplitude=0.01,
              absolute_phase_deg=0.1, rng=None):
    """Add Gaussian noise to a synthetic measurement.

    Amplitude noise is multiplicative, phase noise additive, which is what a
    lock-in measurement produces. Pass an integer or a Generator as `rng` to
    make the draw reproducible.

    With `phase_deg` omitted, only the noisy amplitude is returned, which
    suits a pulsed measurement.
    """
    rng = np.random.default_rng(rng)
    amplitude = np.asarray(amplitude, dtype=float)

    noisy_amp = amplitude * (
        1.0 + relative_amplitude * rng.standard_normal(amplitude.shape)
    )
    if phase_deg is None:
        return noisy_amp

    phase_deg = np.asarray(phase_deg, dtype=float)
    noisy_phase = phase_deg + absolute_phase_deg * rng.standard_normal(phase_deg.shape)
    return noisy_amp, noisy_phase


def contrast_report(s: Sample, warn_below: float = 0.05) -> str:
    """Readable assessment of the effusivity contrast of a sample.

    A small reflection coefficient means the interface is nearly invisible and
    the transit time nearly unmeasurable. The threshold is indicative: the
    uncertainty grows continuously as the coefficient approaches zero, it does
    not jump at any particular value.
    """
    g = s.reflection_coefficient
    lines = [
        f"film effusivity      : {s.film_b:10.0f}",
        f"substrate effusivity : {s.sub_b:10.0f}",
        f"ratio b32            : {s.effusivity_ratio:10.4f}",
        f"reflection Gamma     : {g:+10.5f}",
        f"blind at effusivity  : {s.blind_film_effusivity:10.0f}",
    ]
    if abs(g) < warn_below:
        lines.append("")
        lines.append(
            "WARNING: the contrast is weak. The interface returns little "
            "signal, so the transit time, and with it the thickness and the "
            "diffusivity, are poorly determined. Raising the contrast requires "
            "a different substrate, a transducer layer, or a measurement "
            "geometry that brings in a length other than the thickness."
        )
    return "\n".join(lines)


def analogy_energy(freq, tau: float):
    """Energy of the Schrodinger analogy under the Cattaneo law.

        E = -P = tau omega^2 - i omega

    Its argument, -arctan(1 / (omega tau)), depends on the product omega tau
    alone. It equals -90 degrees in the diffusive limit, exactly -45 degrees
    at omega tau = 1, and tends to zero as the wave limit is approached.
    """
    omega = 2.0 * np.pi * np.asarray(freq, dtype=float)
    return -effective_p(1j * omega, tau)
