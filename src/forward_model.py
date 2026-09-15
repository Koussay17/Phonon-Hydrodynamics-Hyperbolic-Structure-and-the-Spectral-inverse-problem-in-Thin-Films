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

_RELAX_OVERFLOW = 700.0

__all__ = [
    "Sample",
    "contrast_report",
    "regime_report",
    "relaxation_pair",
    "relaxation_wall",
    "relaxation_impedance",
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
        Resistive relaxation time of the constitutive law, in seconds. Zero
        recovers Fourier conduction.
    nonlocal_time, sub_nonlocal_time : float
        Nonlocal time of the Guyer-Krumhansl law, tau_l = 3 ell^2 / a, in
        seconds. Zero reduces the law to Cattaneo.

        Under Guyer-Krumhansl the effective conductivity becomes frequency
        dependent,

            lambda_eff(p) = (lambda + 3 ell^2 p rho c) / (1 + tau_R p),

        and the response depends on the layer through two combinations only:

            sigma e      = xi_1 sqrt[ p (1 + tau_R p) / (1 + tau_l p) ]
            lambda_eff s = b    sqrt[ p (1 + tau_l p) / (1 + tau_R p) ]

        Substituting p by p (1 + tau_R p) throughout the Fourier quadrupole is
        NOT equivalent: it gives the right argument but the wrong flux
        coefficient, by a factor (1 + tau_R p). The correct form reproduces
        the closed-form phase of Camacho de la Rosa et al. (2025), which runs
        from -45 to 0 degrees; the naive substitution runs the other way.
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
    nonlocal_time: float = 0.0
    sub_nonlocal_time: float = 0.0

    phonon_velocity: float = 6000.0

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
    def mean_free_path(self) -> float:
        """Total phonon mean free path, from kinetic theory.

            Lambda = 3 lambda / (rho c v)

        This is the length that decides whether transport is ballistic. It
        must not be confused with the mean free path of the normal processes
        alone, v tau_N, which decides whether phonon hydrodynamics is possible
        at all and is typically three orders of magnitude larger.
        """
        return 3.0 * self.film_lam / (self.film_rho_c * self.phonon_velocity)

    @property
    def knudsen_number(self) -> float:
        """Mean free path divided by thickness.

        Above one the carriers cross the film without scattering and transport
        is ballistic. Below about a tenth it is diffusive. In between the
        extracted conductivity is an apparent value, reduced by boundary
        scattering and dependent on thickness.
        """
        return self.mean_free_path / self.thickness

    @property
    def transport_regime(self) -> str:
        kn = self.knudsen_number
        if kn > 1.0:
            return "ballistic"
        if kn > 0.1:
            return "transitional"
        return "diffusive"

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
    """Spectral parameter of the propagation constant, P = p (1 + tau p).

    It governs the argument of the hyperbolic functions, not the flux
    coefficient. See `relaxation_pair`.
    """
    p = np.asarray(p, dtype=complex)
    return p * (1.0 + tau * p)


def relaxation_pair(p, b: float, xi1: float, tau_r: float, tau_l: float):
    """The two combinations through which a layer enters the response.

    Returns (sigma e, lambda_eff sigma). With both times zero these reduce to
    xi_1 sqrt(p) and b sqrt(p), the Fourier case, exactly.

    When the two times coincide the square roots cancel and the pair is the
    Fourier one whatever their common value: such a medium is thermally
    indistinguishable from a Fourier medium.
    """
    p = np.asarray(p, dtype=complex)
    num = 1.0 + tau_r * p
    den = 1.0 + tau_l * p
    arg = xi1 * np.sqrt(p * num / den)
    coeff = b * np.sqrt(p * den / num)
    return arg, coeff


def relaxation_wall(p, b: float, xi1: float, tau_r: float = 0.0,
                    tau_l: float = 0.0) -> np.ndarray:
    """Quadrupole of a homogeneous layer under Guyer-Krumhansl, per unit area."""
    p = np.asarray(p, dtype=complex)
    arg, coeff = relaxation_pair(p, b, xi1, tau_r, tau_l)

    if np.any(np.abs(arg) > _RELAX_OVERFLOW):
        raise OverflowError(
            f"|sigma e| exceeds {_RELAX_OVERFLOW:.0f}; cosh and sinh overflow."
        )

    ch, sh = np.cosh(arg), np.sinh(arg)
    m = np.empty(p.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = ch
    m[..., 0, 1] = sh / coeff
    m[..., 1, 0] = coeff * sh
    m[..., 1, 1] = ch
    return m if p.shape else m.reshape(2, 2)


def relaxation_impedance(p, b: float, tau_r: float = 0.0, tau_l: float = 0.0):
    """Impedance of a semi-infinite medium under Guyer-Krumhansl.

    Equal to the reciprocal of the flux coefficient. Reduces to 1/(b sqrt(p))
    when both times vanish.
    """
    _, coeff = relaxation_pair(p, b, 1.0, tau_r, tau_l)
    return 1.0 / coeff


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
    if s.is_graded:
        if s.relaxation_time != 0.0 or s.nonlocal_time != 0.0:
            raise NotImplementedError(
                "a graded layer with a finite relaxation time is not treated; "
                "the Liouville transformation would have to be redone with a "
                "frequency-dependent effective conductivity"
            )
        m = q.graded_linear_layer(p, s.film_b, s.film_b_back, s.xi1, form="T")
    else:
        m = relaxation_wall(p, s.film_b, s.xi1,
                            s.relaxation_time, s.nonlocal_time)

    if s.contact_resistance != 0.0:
        m = m @ interface_resistance(p, s.contact_resistance)
    return m


def response(p, s: Sample, power=1.0):
    """Front-face temperature in the Laplace domain, per unit area."""
    m = stack_matrix(p, s)
    z = relaxation_impedance(p, s.sub_b,
                             s.sub_relaxation_time, s.sub_nonlocal_time)
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


def regime_report(s: Sample) -> str:
    """Readable assessment of the transport regime of a sample.

    Two lengths are reported and must not be confused. The total mean free
    path decides whether transport is ballistic; the mean free path of the
    normal processes alone decides whether phonon hydrodynamics is possible,
    and is typically a thousand times larger.
    """
    kn = s.knudsen_number
    lines = [
        f"thickness            : {s.thickness * 1e9:10.1f} nm",
        f"mean free path       : {s.mean_free_path * 1e9:10.1f} nm",
        f"Knudsen number       : {kn:10.3f}",
        f"transport regime     : {s.transport_regime}",
    ]
    if kn > 0.1:
        lines.append("")
        lines.append(
            "WARNING: the mean free path is not negligible against the thickness. "
            "A conductivity extracted here is an apparent value, reduced by "
            "boundary scattering and dependent on thickness; it is not an "
            "intrinsic property of the material. Reporting it as such, or "
            "comparing it with a bulk value, is not meaningful."
        )
    return "\n".join(lines)
