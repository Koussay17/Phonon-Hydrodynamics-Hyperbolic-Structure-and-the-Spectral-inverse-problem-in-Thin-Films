"""Independent angular checks of the grey conserving closure (note 14)."""
import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from forward_model import relaxation_impedance

def sphere():
    mu, w = np.polynomial.legendre.leggauss(12)
    phi = np.arange(32) * (2*np.pi/32)
    m, p = np.meshgrid(mu, phi, indexing="ij")
    n = np.stack((np.sqrt(1-m*m)*np.cos(p),
                  np.sqrt(1-m*m)*np.sin(p), m), axis=-1).reshape(-1, 3)
    weights = np.repeat(w/64, 32)
    return n, weights

def test_isotropic_moments_by_spherical_quadrature():
    n, w = sphere()
    delta = np.eye(3)
    assert_allclose(np.einsum("a,ai,aj->ij", w,n,n), delta/3, atol=1e-14)
    expected = (np.einsum("ij,kl->ijkl",delta,delta)
              + np.einsum("ik,jl->ijkl",delta,delta)
              + np.einsum("il,jk->ijkl",delta,delta))/15
    assert_allclose(np.einsum("a,ai,aj,ak,al->ijkl",w,n,n,n,n),
                    expected, atol=1e-14)

def test_energy_projection_and_deviatoric_stress():
    n,w = sphere()
    # Gradient G_ij = d_i q_j, arbitrary compressive/shear field.
    G = np.array([[.7,.2,-.4],[-.1,.3,.8],[.5,-.6,-.2]])
    v, tau = 3., .07
    # Even kinetic correction -tau*(C theta_t + 3 n_i n_j d_i q_j),
    # using energy conservation C theta_t = -tr G.
    dg = -tau*(-np.trace(G)+3*np.einsum("ai,ij,aj->a",n,G,n))
    assert_allclose(np.dot(w,dg), 0, atol=1e-14)
    assert_allclose(v*np.einsum("a,ai,a->i",w,n,dg),0,atol=1e-14)
    stress = v*v*np.einsum("a,ai,aj,a->ij",w,n,n,dg)
    expected = -v*v*tau/5*(G+G.T-2*np.eye(3)*np.trace(G)/3)
    assert_allclose(stress, expected, atol=1e-14)
    assert_allclose(np.trace(stress),0,atol=1e-14)
    # Omitting the energy time derivative produces the historical trace error.
    wrong = -3*tau*np.einsum("ai,ij,aj->a",n,G,n)
    assert_allclose(np.dot(w,wrong),-tau*np.trace(G),atol=1e-14)

def test_longitudinal_coefficient_from_angular_integral():
    n,w=sphere()
    # q_z has unit z-gradient; the conserving correction is tau*(1-3 n_z^2).
    pi_zz = np.dot(w,n[:,2]**2*(1-3*n[:,2]**2))
    assert_allclose(pi_zz,-4/15,atol=1e-14)
    wrong = np.dot(w,-3*n[:,2]**4)
    assert_allclose(wrong,-3/5,atol=1e-14)

def test_1d_closure_conventions_are_same_when_longitudinal_length_is_fixed():
    a,b,tau_r=2e-5,12000.,1e-9
    ell_h=80e-9
    ell_ce=1.5*ell_h
    p=2j*np.pi*np.geomspace(1e4,2e8,40)
    z_h=relaxation_impedance(p,b,tau_r,3*ell_h**2/a)
    z_ce=relaxation_impedance(p,b,tau_r,(4/3)*ell_ce**2/a)
    assert_allclose(z_h,z_ce,rtol=1e-14)
