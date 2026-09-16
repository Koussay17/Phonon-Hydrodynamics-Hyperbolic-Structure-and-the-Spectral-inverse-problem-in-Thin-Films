"""Axisymmetric Fourier FDTR baseline, with Gaussian pump/probe averaging.

Surface absorption, constant properties, coaxial beams, semi-infinite substrate.
No optical penetration, electron/phonon nonequilibrium or hydrodynamic boundary
conditions are implied. Frequencies use exp(+i omega t); radii are 1/e^2.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad_vec

@dataclass(frozen=True)
class Medium:
    lam_z: float
    rho_c: float
    lam_r: float | None = None

    def __post_init__(self):
        if self.lam_r is None:
            object.__setattr__(self, "lam_r", self.lam_z)
        for value in (self.lam_z,self.lam_r,self.rho_c):
            if not np.isfinite(value) or value <= 0:
                raise ValueError("conductivities and volumetric capacity must be positive")

@dataclass(frozen=True)
class Film:
    medium: Medium
    thickness: float
    resistance_below: float = 0.

    def __post_init__(self):
        if not isinstance(self.medium, Medium):
            raise TypeError("film medium must be a Medium")
        if not np.isfinite(self.thickness) or self.thickness <= 0:
            raise ValueError("thickness must be positive")
        if not np.isfinite(self.resistance_below) or self.resistance_below < 0:
            raise ValueError("interface resistance must be nonnegative")

def _frequencies(freq):
    f=np.atleast_1d(np.asarray(freq,dtype=float))
    if f.ndim!=1 or not f.size or not np.all(np.isfinite(f)) or np.any(f<=0):
        raise ValueError("frequencies must be a nonempty positive vector")
    return f

def radial_impedance(freq, radial_wavenumber, films, substrate):
    """Surface T/normal flux for one radial wavenumber, in m^2 K/W.

    Films are ordered top-to-bottom; each carries its bottom contact.
    """
    f=_frequencies(freq)
    k=float(radial_wavenumber)
    if not np.isfinite(k) or k<0:
        raise ValueError("radial wavenumber must be finite and nonnegative")
    p=2j*np.pi*f
    def admittance(medium):
        gamma=np.sqrt((medium.lam_r*k*k+medium.rho_c*p)/medium.lam_z)
        return gamma,medium.lam_z*gamma
    _, ys=admittance(substrate)
    z=1/ys
    for film in reversed(tuple(films)):
        z=z+film.resistance_below
        gamma,y=admittance(film.medium)
        argument=gamma*film.thickness
        t=np.ones_like(argument)
        use=argument.real<20
        t[use]=np.tanh(argument[use])
        z=(z+t/y)/(1+y*z*t)
    return z

def gaussian_response(freq, films, substrate, pump_radius, probe_radius,
                      epsrel=1e-8):
    """Probe-averaged complex temperature per absorbed watt, K/W.

    Pump: q(r)=2 P/(pi wp^2) exp(-2r^2/wp^2).
    Probe: normalized weight 2/(pi wr^2) exp(-2r^2/wr^2).
    Adaptive quadrature integrates the radial Fourier response, not its phase.
    """
    f=_frequencies(freq)
    films=tuple(films)
    for value in (pump_radius,probe_radius,epsrel):
        if not np.isfinite(value) or value<=0:
            raise ValueError("radii and quadrature tolerance must be positive")
    alpha=(pump_radius**2+probe_radius**2)/8
    def integrand(u):
        z=radial_impedance(f,u/np.sqrt(alpha),films,substrate)
        return u*np.exp(-u*u)*z/(2*np.pi*alpha)
    # Resolve crossover near the origin, especially as frequency approaches DC.
    points=np.sqrt(alpha*substrate.rho_c*2*np.pi*f/substrate.lam_r)
    points=np.unique(points[(points>0)&(points<10)])
    value,error,info=quad_vec(integrand,0.,10.,epsabs=1e-12,epsrel=epsrel,
                              points=points,full_output=True)
    if not info.success or not np.all(np.isfinite(value)):
        raise RuntimeError(f"radial integration failed: {info.message}")
    return value

def plane_response(freq, films, substrate, pump_radius, probe_radius):
    """Local 1D approximation to the same Gaussian-weighted observable."""
    for value in (pump_radius,probe_radius):
        if not np.isfinite(value) or value<=0:
            raise ValueError("radii must be positive")
    return radial_impedance(freq,0.,films,substrate)*2/(np.pi*(pump_radius**2+probe_radius**2))
