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


# --- 2 October 2026 audit: graded-film parametrization and silent fit success ---

def _constant_diffusivity_graded(**kw):
    # 30/1.205e6 equals 60/2.41e6, so the default graded layer has constant diffusivity.
    return fm.Sample(film_lam_back=30.0, film_rho_c_back=1.205e6, **kw)

@pytest.mark.parametrize("name",["film_lam","film_rho_c","film_lam_back","film_rho_c_back"])
def test_constant_diffusivity_graded_property_fit_is_rejected_explicitly(name):
    f=np.logspace(4,8,20); a,p=fm.modulated_response(f,_constant_diffusivity_graded())
    with pytest.raises(ValueError,match="graded_xi1"):
        inv.fit_modulated(f,a,p,_constant_diffusivity_graded(),[name])
    with pytest.raises(ValueError,match="graded_xi1"):
        inv.fisher_analysis(f,_constant_diffusivity_graded(),[name])

def test_constant_diffusivity_graded_thickness_still_fits():
    f=np.logspace(4,8,40); a,p=fm.modulated_response(f,_constant_diffusivity_graded())
    r=inv.fit_modulated(f,a,p,_constant_diffusivity_graded(thickness=450e-9),["thickness"])
    assert r.success and r.values[0]==pytest.approx(500e-9,rel=1e-8)

def test_variable_diffusivity_graded_fit_with_explicit_xi1():
    xi1=1.1*fm.Sample().xi1
    truth=fm.Sample(film_lam_back=30.0,film_rho_c_back=2.0e6,graded_xi1=xi1)
    f=np.logspace(4,8,40); a,p=fm.modulated_response(f,truth)
    r=inv.fit_modulated(f,a,p,replace(truth,film_lam=50.0),["film_lam"])
    assert r.success and r.values[0]==pytest.approx(60.0,rel=1e-8)

def test_fit_reports_failure_when_solution_is_a_sentinel():
    fun=np.full(4,inv._INVALID); jac=np.zeros((4,1))
    assert not inv._solution_is_valid(fun,jac)
    assert not inv._solution_is_valid(np.zeros(4),np.zeros((4,1)))
    assert inv._solution_is_valid(np.zeros(4),np.ones((4,1)))

def test_none_parameter_raises_value_error_not_type_error():
    f=np.logspace(4,8,10); a,p=fm.modulated_response(f,fm.Sample())
    with pytest.raises(ValueError):
        inv.fit_modulated(f,a,p,fm.Sample(),["film_lam_back"])
    with pytest.raises(ValueError):
        inv.fit_pulsed(np.logspace(-9,-6,10),np.ones(10),fm.Sample(),["film_lam_back"])
