"""Physical normalization and independent ODE benchmarks for the FDTR baseline."""
import sys
from pathlib import Path
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
from scipy.special import erfcx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from fdtr import Medium, Film, radial_impedance, gaussian_response, plane_response

@pytest.mark.parametrize("freq",[[1e-4,1e2,1e6,2e8],[1e4,1e7]])
def test_gaussian_halfspace_against_closed_form(freq):
    lam,C=60.,2.41e6
    wp,wr=10e-6,8e-6
    alpha=(wp*wp+wr*wr)/8
    beta=np.sqrt(2j*np.pi*np.array(freq)*C/lam)
    exact=erfcx(beta*np.sqrt(alpha))/(4*lam*np.sqrt(np.pi*alpha))
    actual=gaussian_response(freq,[],Medium(lam,C),wp,wr)
    assert_allclose(actual,exact,rtol=1e-8)

def test_radial_layer_against_physical_energy_balance():
    f=np.array([1e5,2e6])
    k=2e5
    sub=Medium(35.,3.03e6)
    film=Film(Medium(60.,2.41e6,120.),500e-9,2e-8)
    for frequency in f:
        p=2j*np.pi*frequency
        ysub=np.sqrt(sub.lam_z*(sub.lam_r*k*k+sub.rho_c*p))
        # Physical z-coordinate ODE, backward from interface T jump.
        def rhs(z,state):
            T,q=state
            return [-q/film.medium.lam_z,
                    -(film.medium.rho_c*p+film.medium.lam_r*k*k)*T]
        initial=np.array([1/ysub+film.resistance_below,1.],complex)
        sol=solve_ivp(rhs,(film.thickness,0.),initial,rtol=1e-11,atol=1e-13)
        assert sol.success
        assert_allclose(radial_impedance([frequency],k,[film],sub)[0],
                         sol.y[0,-1]/sol.y[1,-1],rtol=1e-9)

def test_split_layer_and_large_spot_limit():
    sub=Medium(35.,3.03e6)
    m=Medium(60.,2.41e6)
    f=np.array([1e6,2e8])
    whole=[Film(m,500e-9,2e-8)]
    split=[Film(m,200e-9),Film(m,300e-9,2e-8)]
    assert_allclose(gaussian_response(f,whole,sub,10e-6,8e-6),
                    gaussian_response(f,split,sub,10e-6,8e-6),rtol=1e-10)
    assert_allclose(gaussian_response(f,whole,sub,1e-2,1e-2),
                    plane_response(f,whole,sub,1e-2,1e-2),rtol=1e-5)

@pytest.mark.parametrize("args",[(0,2.),(1.,-2.),(1.,2.,np.nan)])
def test_nonphysical_medium_rejected(args):
    with pytest.raises(ValueError):
        Medium(*args)

@pytest.mark.parametrize("freq",[[0.],[-1.],[np.nan],[]])
def test_nonphysical_frequency_rejected(freq):
    with pytest.raises(ValueError):
        gaussian_response(freq,[],Medium(1.,2.),1.,1.)
