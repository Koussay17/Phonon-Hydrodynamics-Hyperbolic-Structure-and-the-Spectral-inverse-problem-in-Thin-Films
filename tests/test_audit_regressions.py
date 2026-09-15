"""Independent physical balances and counterexamples from the September audit."""
from pathlib import Path
import sys
from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import forward_model as fm
import quadrupoles as q
import inversion as inv
import admissibility as adm

@pytest.mark.parametrize("tr,tl",[(0.,0.),(1e-9,0.),(1e-9,2e-10),(2e-10,1e-9)])
def test_homogeneous_matrix_against_independent_energy_and_flux_odes(tr,tl):
    s=fm.Sample(film_lam=321)
    for p in (1e6,2j*np.pi*1e6,2j*np.pi*1e8):
        keff=s.film_lam*(1+tl*p)/(1+tr*p)
        # Physical balances dT/dz=-q/keff, dq/dz=-p C T.
        a=np.array([[0,-s.thickness*s.film_b/keff],
                    [-s.thickness*p*s.film_rho_c/s.film_b,0]],complex)
        sol=solve_ivp(lambda z,u:(a@u.reshape(2,2)).ravel(),(1,0),
                      np.eye(2,dtype=complex).ravel(),method="DOP853",
                      rtol=1e-11,atol=1e-12)
        assert sol.success
        m=fm.relaxation_wall(p,s.film_b,s.xi1,tr,tl)
        got=np.diag([1,1/s.film_b])@m@np.diag([1,s.film_b])
        ref=sol.y[:,-1].reshape(2,2)
        assert np.linalg.norm(got-ref)/np.linalg.norm(ref)<1e-9

@pytest.mark.parametrize("form",["T","phi"])
def test_graded_matrix_against_independent_ode_including_dc(form):
    b0,b1,xi=12000.,24000.,1e-4
    sign=1 if form=="T" else -1
    for p in (0,1e6,2j*np.pi*1e8):
        def rhs(x,u):
            # Construct the profile directly; do not call the profile helper.
            meta=(1-x)*b0**(sign/2)+x*b1**(sign/2)
            b=meta**(2/sign)
            a=np.array([[0,-xi*b0/b],[-xi*p*b/b0,0]],complex)
            return (a@u.reshape(2,2)).ravel()
        sol=solve_ivp(rhs,(1,0),np.eye(2,dtype=complex).ravel(),
                      method="DOP853",rtol=1e-11,atol=1e-12)
        assert sol.success
        m=q.graded_linear_layer(p,b0,b1,xi,form)
        got=np.diag([1,1/b0])@m@np.diag([1,b0])
        ref=sol.y[:,-1].reshape(2,2)
        assert np.linalg.norm(got-ref)/np.linalg.norm(ref)<1e-9

def test_dc_matrix_is_resistance_and_zero_thickness_is_identity():
    expected=np.array([[1,5e-7/60],[0,1]])
    assert np.allclose(q.homogeneous_wall(0,60,2.41e6,5e-7),expected,atol=1e-20)
    assert np.array_equal(q.homogeneous_wall(0,60,2.41e6,0),np.eye(2))
    s=fm.Sample()
    assert np.allclose(fm.relaxation_wall(0,s.film_b,s.xi1,1e-9,2e-9),
                       expected,atol=1e-20)

def test_graded_geometry_cannot_be_inferred_from_incomplete_data():
    with pytest.raises(ValueError):
        fm.Sample(film_lam_back=120)
    with pytest.raises(ValueError):
        fm.Sample(film_lam_back=120,film_rho_c_back=2.41e6)
    s=fm.Sample(film_lam_back=120,film_rho_c_back=2.41e6,graded_xi1=2e-4)
    assert s.xi1==2e-4 and s.diffusion_time==pytest.approx(4e-8)

def test_scale_gauge_does_not_produce_a_finite_covariance():
    s=fm.Sample()
    with pytest.raises(inv.NonIdentifiableError):
        inv.fisher_analysis(np.logspace(4,8,40),s,["film_lam","film_rho_c","thickness"])
    with pytest.raises(inv.NonIdentifiableError):
        inv.information_from_jacobian(np.array([[1,2],[2,4],[3,6]]))
    # Tiny independent sensitivity means large uncertainty, not rank loss.
    _,cov=inv.information_from_jacobian(np.diag([1,1e-12]))
    assert cov[1,1]==pytest.approx(1e24)

def test_low_frequency_film_sensitivity_is_not_zero():
    s=fm.Sample(film_lam=321)
    p=2j*np.pi*np.logspace(3,6,30)
    change=np.max(abs(fm.response(p,replace(s,thickness=1e-6))/fm.response(p,s)-1))
    assert .15<change<.18

def test_equal_effusivity_is_not_blind_with_contact_resistance():
    s=fm.Sample(film_lam=35*3.03e6/2.41e6,contact_resistance=1e-8)
    p=2j*np.pi*1e7
    assert abs(fm.response(p,s)/fm.relaxation_impedance(p,s.sub_b)-1)>.05

def test_sub_crossover_relaxation_can_have_resolvable_ideal_phase():
    x=2*np.pi*2e8*1e-11
    phase=np.degrees(np.arctan(x))/2
    relative_error=np.radians(.01)/(x/(2*(1+x*x)))
    assert phase==pytest.approx(.35998,rel=1e-4)
    assert relative_error==pytest.approx(.02778,rel=1e-3)

def test_wavenumber_branch_decays_into_the_sample():
    k=adm.wavenumber(np.array([1e6,1e8]),60/2.41e6,1e-9,2e-10)
    assert np.all(k.imag>0)
