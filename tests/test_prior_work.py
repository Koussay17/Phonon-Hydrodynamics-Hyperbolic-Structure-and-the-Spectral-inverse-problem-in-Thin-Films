"""Independent equation checks used in the prior-work comparison."""
import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from forward_model import relaxation_impedance

def test_hennessy_longitudinal_limit_and_factor_three():
    lam,C,tau,ell=150.,1.692e6,42e-12,185e-9
    a=lam/C
    p=2j*np.pi*np.geomspace(1e5,1e10,50)
    gamma=np.sqrt(p*(1+p*tau)/(a+3*ell**2*p))
    lam_eff=(lam+3*C*ell**2*p)/(1+p*tau)
    independent=1/(lam_eff*gamma)
    assert_allclose(independent,relaxation_impedance(p,np.sqrt(lam*C),tau,3*ell**2/a),rtol=2e-14)
    zf=relaxation_impedance(p,np.sqrt(lam*C))
    assert_allclose(relaxation_impedance(p,np.sqrt(lam*C),3*ell**2/a,3*ell**2/a),zf,rtol=2e-14)
    assert np.max(abs(relaxation_impedance(p,np.sqrt(lam*C),ell**2/a,3*ell**2/a)/zf-1))>.1

@pytest.mark.parametrize("p",[.01+.1j, .2+1j, 2+3j])
@pytest.mark.parametrize("u0",[.4,1.1])
def test_darboux_profile_against_energy_flux_ode(p,u0):
    # Dimensionless xi/c, A=c=1. At xi=20, b differs from 1 below roundoff.
    root=np.sqrt(p)
    def ode(x,y):
        b=np.tanh(x+u0)**2
        return [-y[1]/b,-p*b*y[0]]
    sol=solve_ivp(ode,(20.,0.),np.array([1.,root],complex),rtol=2e-11,atol=1e-12)
    assert sol.success
    numerical=sol.y[0,-1]/sol.y[1,-1]
    t=np.tanh(u0)
    exact=(root+t)/(t*root*(1+t*root))
    assert_allclose(numerical,exact,rtol=2e-9)
