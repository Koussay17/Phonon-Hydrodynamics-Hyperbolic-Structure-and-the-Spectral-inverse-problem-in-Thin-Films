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
        Film conductivity and volumetric heat capacity. The defaults are
        illustrative placeholders, not measurements of a particular AlN film.
    thickness : float
        Film thickness in metres.
    film_lam_back, film_rho_c_back : float or None
        Supply both or neither. A variable diffusivity also requires
        graded_xi1, the integral of dz/sqrt(a). Otherwise a constant
        diffusivity throughout the layer is assumed, requiring equal endpoint
        diffusivities. Endpoint values alone do not determine the integral.
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
        Nonlocal time tau_l = L^2 / a, in seconds. Historical GK uses
        L^2 = 3 ell^2; the conserving grey closure uses L^2 = 4 ell^2 / 3.
        Zero reduces the law to Cattaneo. See note 14 for the conversion.

        Under Guyer-Krumhansl the effective conductivity becomes frequency
        dependent,

            lambda_eff(p) = lambda * (1 + tau_l p) / (1 + tau_R p),

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
    graded_xi1: float | None = None

    def __post_init__(self):
        for name in ("film_lam", "film_rho_c", "thickness", "sub_lam",
                     "sub_rho_c", "phonon_velocity"):
            q._positive(**{name: getattr(self, name)})
        for name in ("contact_resistance", "front_losses", "relaxation_time",
                     "sub_relaxation_time", "nonlocal_time", "sub_nonlocal_time"):
            value = getattr(self, name)
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if (self.film_lam_back is None) != (self.film_rho_c_back is None):
            raise ValueError("both back-face properties are required for a graded layer")
        if self.film_lam_back is not None:
            q._positive(film_lam_back=self.film_lam_back,
                        film_rho_c_back=self.film_rho_c_back)
            if self.graded_xi1 is None:
                a_back = self.film_lam_back / self.film_rho_c_back
                if not np.isclose(a_back, self.film_a, rtol=1e-10, atol=0):
                    raise ValueError("a graded layer with variable diffusivity requires graded_xi1")
        if self.graded_xi1 is not None:
            q._positive(graded_xi1=self.graded_xi1)
            if self.film_lam_back is None:
                raise ValueError("graded_xi1 requires graded back-face properties")


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
        """Liouville thickness; explicit integral for variable graded diffusivity."""
        if self.graded_xi1 is not None:
            return self.graded_xi1
        return q.xi_from_thickness(self.thickness, self.film_a)

    @property
    def diffusion_time(self) -> float:
        """Square of Liouville thickness; thickness^2/a for a homogeneous film."""
        return self.xi1 ** 2

    @property
    def relaxation_frequency(self) -> float:
        """Frequency at which omega tau reaches one, in hertz.

        This is a crossover, not a detectability cutoff. Identifiability also
        depends on noise, calibration, nuisance parameters and the sampled band.
        """
        if self.relaxation_time <= 0.0:
            return np.inf
        return 1.0 / (2.0 * np.pi * self.relaxation_time)

    @property
    def mean_free_path(self) -> float:
        """Grey conductivity-derived transport length; not an all-collision mean free path.

            Lambda = 3 lambda / (rho c v)

        A single-velocity kinetic estimate. It does not resolve the phonon
        spectrum and cannot establish the transport regime of a real film.
        Mode-resolved normal-process lengths v*tau_N are different quantities.
        """
        return 3.0 * self.film_lam / (self.film_rho_c * self.phonon_velocity)

    @property
    def knudsen_number(self) -> float:
        """Mean free path divided by thickness.

        The conventional cutoffs 0.1 and 1 are heuristic grey-model labels.
        A spectrum of mean free paths requires a spectral transport calculation.
        """
        return self.mean_free_path / self.thickness

    @property
    def transport_regime(self) -> str:
        """Heuristic grey-model label; not a measured material classification."""
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

        For Fourier media at perfect contact (or matched non-Fourier times),
        the response can be written as a series of waves reflected between the
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
        film and interface sensitivity decays. There is no exact visibility
        cutoff; evaluate sensitivities and noise on the proposed band.

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
    Fourier one whatever their common value: this homogeneous source-free, zero-initial-perturbation layer has the
    Fourier quadrupole (Fourier resonance; Kovacs 2018).
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
    q._positive(b=b)
    if xi1 < 0 or tau_r < 0 or tau_l < 0:
        raise ValueError("thickness and relaxation times must be nonnegative")
    arg, coeff = relaxation_pair(p, b, xi1, tau_r, tau_l)

    if np.any(np.abs(arg.real) > _RELAX_OVERFLOW):
        raise OverflowError(
            f"|Re(sigma e)| exceeds {_RELAX_OVERFLOW:.0f}; cosh and sinh overflow."
        )

    ch, sh = np.cosh(arg), np.sinh(arg)
    m = np.empty(p.shape + (2, 2), dtype=complex)
    m[..., 0, 0] = ch
    m[..., 0, 1] = xi1 / b * (1 + tau_r*p)/(1 + tau_l*p) * q.sinhc(arg)
    m[..., 1, 0] = b * xi1 * p * q.sinhc(arg)
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
    p = np.asarray(p, dtype=complex)
    if np.any(~np.isfinite(p)) or np.any(p == 0):
        raise ValueError("response requires finite nonzero spectral parameters")
    z = relaxation_impedance(p, s.sub_b,
                             s.sub_relaxation_time, s.sub_nonlocal_time)
    if s.is_graded:
        return q.front_face_temperature(stack_matrix(p, s), z,
                                        power=power, h=s.front_losses)
    # Divide the transfer relation by cosh before evaluating it. This has a
    # finite thick-layer limit, unlike a product of exponentially large matrices.
    z = z + s.contact_resistance
    arg, coeff = relaxation_pair(p, s.film_b, s.xi1,
                                 s.relaxation_time, s.nonlocal_time)
    t = q.stable_tanh(arg)
    zin = (z + t / coeff) / (1 + coeff * z * t)
    return power * zin / (1 + s.front_losses * zin)


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
        -45 degrees at every frequency for a Fourier half-space.
    """
    freq = np.asarray(freq, dtype=float)
    if np.any(~np.isfinite(freq)) or np.any(freq <= 0):
        raise ValueError("frequencies must be finite and positive")
    theta = response(2j * np.pi * freq, s, power=power)
    return np.abs(theta), np.degrees(np.angle(theta))


def pulsed_response(t, s: Sample, energy=1.0, n_stehfest: int = 12):
    """Front-face temperature after a Dirac pulse of surface energy density.

    Validate Stehfest on the signal and its derivatives; the exponential
    benchmark in laplace.py does not give a universal time window.
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
        "Fourier/perfect-contact contrast diagnostic; not a full identifiability test.",
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

    Reports only a grey conductivity-derived length. Mode-resolved normal
    scattering times are not supplied and hydrodynamics cannot be diagnosed.
    """
    kn = s.knudsen_number
    lines = [
        f"thickness            : {s.thickness * 1e9:10.1f} nm",
        f"grey transport length: {s.mean_free_path * 1e9:10.1f} nm",
        f"Knudsen number       : {kn:10.3f}",
        f"grey-model regime    : {s.transport_regime} (heuristic)",
    ]
    if kn > 0.1:
        lines.append("")
        lines.append(
            "WARNING: the grey length is not negligible against the thickness. "
            "Boundary scattering may matter; model adequacy requires checking. "
            "Report an apparent film conductivity with thickness and assumptions. "
            "Comparison with bulk values is meaningful as a suppression comparison."
        )
    return "\n".join(lines)
