"""Independent limits and published-file checks for spectral processing."""
from pathlib import Path
import json
import numpy as np
import pytest
from scipy.constants import Boltzmann as KB, hbar
from scipy.integrate import quad, solve_ivp
from src.spectral import (mode_heat_capacity, rta_moments, inplane_suppression,
                          conserving_collision_action, collision_order_fractions)

ROOT=Path(__file__).resolve().parents[1]


def test_bose_capacity_matches_energy_derivative_and_limits():
    w=np.array([0.,1e11,1e13,1e14,1e16])
    assert mode_heat_capacity(w,300)[0]==KB
    def energy(t):
        x=hbar*w[1:]/(KB*t)
        return hbar*w[1:]*np.exp(-x)/(-np.expm1(-x))
    numerical=(energy(300.001)-energy(299.999))/.002
    np.testing.assert_allclose(mode_heat_capacity(w,300)[1:],numerical,rtol=1e-7,atol=1e-40)
    np.testing.assert_allclose(mode_heat_capacity(w[:3],1e9),KB,rtol=1e-12)


def test_debye_low_temperature_capacity_integral():
    # Single acoustic branch, linear dispersion and infinite low-T cutoff:
    # C = 2 pi^2 kB^4 T^3 / (15 hbar^3 v^3).
    t=10.;v=5000.
    actual=quad(lambda x: float(mode_heat_capacity(x*KB*t/hbar,t))*x*x,0,100,epsabs=1e-30)[0]
    actual*= (KB*t/hbar)**3/(2*np.pi**2*v**3)
    expected=2*np.pi**2*KB**4*t**3/(15*hbar**3*v**3)
    np.testing.assert_allclose(actual,expected,rtol=1e-10)


def test_static_and_memory_means_have_different_weights():
    c=np.array([2.,3.]);v2=np.array([4.,9.]);tau=np.array([1.,10.])
    m=rta_moments(c,v2,1/tau)
    np.testing.assert_allclose(m["tau_static_s"],278/35)
    np.testing.assert_allclose(m["tau_memory_s"],2708/278)
    # Derivative of independent modal Laplace response.
    h=1e-7
    response=lambda p: np.sum(c*v2*tau/(1+p*tau))
    derivative=(response(h)-response(-h))/(2*h)
    np.testing.assert_allclose(-derivative/response(0),m["tau_memory_s"],rtol=1e-8)
    assert m["tau_memory_s"]>=m["tau_static_s"]


@pytest.mark.parametrize("kn,p",[(.03,0),(.8,0),(8.,.5),(1e3,.9)])
def test_film_suppression_against_transport_ode(kn,p):
    # h+ satisfies kn dh+/dz + h+ = 1, h+(0)=p h-(0).
    # h-(0)=h+(1) for a symmetric film. Integrate with numerical shooting.
    from scipy.optimize import brentq
    def solve(initial):
        return solve_ivp(lambda z,y:[(1-y[0])/kn,y[0]],(0,1),[initial,0],
                         rtol=2e-11,atol=2e-12)
    initial=brentq(lambda a:a-p*solve(a).y[0,-1],0,1,xtol=1e-13)
    np.testing.assert_allclose(inplane_suppression(kn,p),solve(initial).y[1,-1],rtol=1e-8)


def test_fuchs_sondheimer_thick_isotropic_limit_and_specular_limit():
    # Angular average of v_x^2 in an isotropic solid.
    kn=.001
    avg=quad(lambda mu:1.5*(1-mu*mu)*float(inplane_suppression(kn*mu)),0,1)[0]
    np.testing.assert_allclose(avg,1-3*kn/8,rtol=1e-12)
    np.testing.assert_array_equal(inplane_suppression([0,1,1e15],1),[1,1,1])
    s=inplane_suppression([1e6,1e12],0)
    np.testing.assert_allclose(s,[.5e-6,.5e-12],rtol=1e-6)


def test_projected_collision_conserves_moments_is_positive_and_basis_invariant():
    rng=np.random.default_rng(178)
    r=np.geomspace(.01,100,20);h=rng.normal(size=(20,4));y=rng.normal(size=20)
    out=conserving_collision_action(r,h,y)
    np.testing.assert_allclose(h.T@out,0,atol=2e-12)
    assert y@out>=0
    np.testing.assert_allclose(conserving_collision_action(r,h,h@np.arange(4)),0,atol=2e-12)
    np.testing.assert_allclose(out,conserving_collision_action(r,h@np.diag([1e-20,1,1e10,3]),y),rtol=1e-10,atol=1e-12)
    # Construct its matrix and verify reciprocity/nonnegative eigenvalues.
    matrix=np.column_stack([conserving_collision_action(r,h,v) for v in np.eye(20)])
    np.testing.assert_allclose(matrix,matrix.T,atol=2e-12)
    assert np.linalg.eigvalsh(matrix).min()>-1e-12


def test_ordering_requires_same_mode_and_strong_separation():
    out=collision_order_fractions([1000,1,1],[1,1,1000],[20,1000,1],[2,3,5])
    assert out["candidate_weight_fraction"]==.2
    assert out["boundary_first_weight_fraction"]==.3


def test_published_aln_capacity_conductivity_and_scattering_sum():
    with np.load(ROOT/"theory/aln/rao_300K_modes.npz") as z:
        w,v,r=z["omega_rad_s"],z["velocity_m_s"],z["total_rate_s"]
        c=mode_heat_capacity(w,300)*z["degeneracy"][:,None]/z["degeneracy"].sum()/z["cell_volume_m3"]
        np.testing.assert_allclose(r,z["anharmonic_rate_s"]+z["isotope_rate_s"],rtol=2e-8,atol=1e3)
        assert w.shape==(793,12) and z["degeneracy"].sum()==24**3
    ref=json.loads((ROOT/"theory/aln/sources.json").read_text())["rao"]["references"]
    np.testing.assert_allclose(c.sum(),ref["C_J_m3K"],rtol=1e-6)
    for v2,idx in (((v[:,:,0]**2+v[:,:,1]**2)/2,0),(v[:,:,2]**2,2)):
        np.testing.assert_allclose(rta_moments(c,v2,r)["kappa"],
                                  ref["kappa_RTA_W_mK"][idx][idx],rtol=3e-5)


@pytest.mark.parametrize("func,args",[
    (mode_heat_capacity,([1],0)),(mode_heat_capacity,([-1],300)),
    (inplane_suppression,([-1],)),(inplane_suppression,([1],1.1)),
    (rta_moments,([1],[1],[0])),
    (conserving_collision_action,([0],[[1]],[1])),
    (collision_order_fractions,([1],[1],[1],[0]))])
def test_invalid_physical_inputs_rejected(func,args):
    with pytest.raises(ValueError):func(*args)
