"""Validation of the axisymmetric FDTR twin with Fourier, Cattaneo and GK layers.

(a) Fourier limit equals src/fdtr.py.
(b) k = 0 (large-spot) limit equals the 1D lambda_eff model of forward_model.py
    and a 1D quadrupole chain built from the repository quadrupoles.
(c) Cattaneo half-space phase -pi/4 + arctan(omega tau)/2 (Camacho et al. 2025).
(d) Independent ODE shooting (solve_ivp) of the Hankel-space layer system.
(e) alpha drops out when ell = 0; GK tends to Cattaneo as ell -> 0.
(f) Interface flux continuity, layer energy balance, dissipation identity.
(g) Convergence of the Hankel quadrature.
Also: Hennessy-Myers (2021) half-space closed form, thin-film Poiseuille limit
of the in-plane conductivity, 200-digit transfer matrices where growing
exponentials reach 1e100, and the inadmissible normal-gradient condition.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.special import erfcx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import fdtr  # noqa: E402
import forward_model as fm  # noqa: E402
import quadrupoles as q  # noqa: E402
import twin_fdtr as tw  # noqa: E402

AU = tw.Material(150.0, 2.49e6)
SAPPHIRE = tw.Material(35.0, 3.03e6)
FREQ = np.geomspace(1e4, 2e8, 9)


def gk_aln(tau_r=2e-11, ell=50e-9, alpha=1 / 3):
    return tw.Material(60.0, 2.41e6, tau_r=tau_r, ell=ell, alpha=alpha)


def entropic_like(material, weight):
    """Kinetic temperature with a fraction of the entropic coefficients."""
    w = material.ell ** 2 / material.lam_z
    return tw.GKBoundary(temperature="kinetic", c_n=weight * (1 + material.alpha) * w,
                         c_t=weight * material.alpha * w)


# --------------------------------------------------------------------------
# (a) Fourier limit
# --------------------------------------------------------------------------

def _fourier_pair():
    films = [fdtr.Film(fdtr.Medium(150.0, 2.49e6), 80e-9, 1e-8),
             fdtr.Film(fdtr.Medium(60.0, 2.41e6, 120.0), 500e-9, 2e-8)]
    stack = tw.Stack((tw.Layer(AU, 80e-9, 1e-8),
                      tw.Layer(tw.Material(60.0, 2.41e6, 120.0), 500e-9, 2e-8)), SAPPHIRE)
    return films, fdtr.Medium(35.0, 3.03e6), stack


@pytest.mark.parametrize("k", [0.0, 1e4, 3e5, 2e6, 1e7])
def test_fourier_kernel_matches_fdtr(k):
    films, sub, stack = _fourier_pair()
    assert_allclose(tw.kernel(FREQ, k, stack), fdtr.radial_impedance(FREQ, k, films, sub),
                    rtol=1e-12)


def test_fourier_gaussian_matches_fdtr_and_adaptive():
    films, sub, stack = _fourier_pair()
    ref = fdtr.gaussian_response(FREQ, films, sub, 10e-6, 8e-6, epsrel=1e-13)
    assert_allclose(tw.gaussian_response(FREQ, stack, 10e-6, 8e-6), ref, rtol=1e-10)
    assert_allclose(tw.gaussian_response_adaptive(FREQ, stack, 10e-6, 8e-6), ref, rtol=1e-10)


def test_bare_fourier_halfspace_closed_form():
    lam, C, wp, wr = 60.0, 2.41e6, 10e-6, 8e-6
    f = np.array([1e-2, 1e2, 1e4, 1e6, 2e8])
    eta = (wp * wp + wr * wr) / 8
    exact = erfcx(np.sqrt(2j * np.pi * f * C / lam) * np.sqrt(eta)) / (4 * lam * np.sqrt(np.pi * eta))
    actual = tw.gaussian_response(f, tw.Stack((), tw.Material(lam, C)), wp, wr)
    assert_allclose(actual, exact, rtol=1e-11)


# --------------------------------------------------------------------------
# (b) k = 0 and large-spot limits
# --------------------------------------------------------------------------

@pytest.mark.parametrize("alpha", [1 / 3, 2.0])
@pytest.mark.parametrize("slip", [0.0, 1.0, np.inf])
@pytest.mark.parametrize("gk_substrate", [False, True])
def test_gk_k0_equals_forward_model(alpha, slip, gk_substrate):
    film = gk_aln(tau_r=3e-11, ell=60e-9, alpha=alpha)
    sub = (tw.Material(35.0, 3.03e6, tau_r=1e-11, ell=20e-9, alpha=alpha)
           if gk_substrate else SAPPHIRE)
    bc = tw.GKBoundary(slip=slip)
    stack = tw.Stack((tw.Layer(film, 500e-9, 2e-8, bc),), sub, bc)
    sample = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9,
                       sub_lam=35.0, sub_rho_c=3.03e6, contact_resistance=2e-8,
                       relaxation_time=3e-11, nonlocal_time=film.tau_l,
                       sub_relaxation_time=sub.tau_r, sub_nonlocal_time=sub.tau_l)
    f = np.geomspace(1e3, 1e10, 15)
    assert_allclose(tw.kernel(f, 0.0, stack), fm.response(2j * np.pi * f, sample), rtol=1e-12)


@pytest.mark.parametrize("convention", ["local", "entropic", "kinetic"])
def test_k0_stack_matches_quadrupole_chain(convention):
    """At k = 0 the S modes decouple and theta = (1 + c_n C p) T on GK faces."""
    film = gk_aln(tau_r=2e-11, ell=80e-9)
    if convention == "kinetic":
        bc = entropic_like(film, 0.1)
        c_n = bc.c_n
    else:
        bc = tw.GKBoundary(slip=1.0, temperature=convention)
        c_n = 0.0 if convention == "local" else (1 + film.alpha) * film.ell ** 2 / film.lam_z
    stack = tw.au_aln_sapphire(film, boundary=bc)
    f = np.geomspace(1e4, 2e9, 12)
    p = 2j * np.pi * f
    c = 1 + c_n * film.rho_c * p
    dmat = np.zeros(p.shape + (2, 2), complex)
    dmat[..., 0, 0], dmat[..., 1, 1] = c, 1.0
    dinv = dmat.copy()
    dinv[..., 0, 0] = 1 / c
    b = np.sqrt(60.0 * 2.41e6)
    xi1 = 500e-9 / np.sqrt(60.0 / 2.41e6)
    chain = (q.homogeneous_wall(p, 150.0, 2.49e6, 80e-9) @ fm.interface_resistance(p, 1e-8)
             @ dmat @ fm.relaxation_wall(p, b, xi1, film.tau_r, film.tau_l) @ dinv
             @ fm.interface_resistance(p, 2e-8))
    z1d = q.front_face_temperature(chain, q.semi_infinite_impedance(p, np.sqrt(35.0 * 3.03e6)))
    assert_allclose(tw.kernel(f, 0.0, stack), z1d, rtol=1e-11)


def test_large_spot_converges_to_forward_model_at_second_order():
    """Gaussian average of the 3D twin -> 1D forward_model response of the same observable,
    2 Z_1D/(pi (wp^2 + wr^2)), with error O(1/w^2) (measured order 1.997-2.000)."""
    film = gk_aln(2e-11, 50e-9)
    stack = tw.Stack((tw.Layer(film, 500e-9, 2e-8),), SAPPHIRE)
    sample = fm.Sample(film_lam=60.0, film_rho_c=2.41e6, thickness=500e-9, sub_lam=35.0,
                       sub_rho_c=3.03e6, contact_resistance=2e-8, relaxation_time=2e-11,
                       nonlocal_time=film.tau_l)
    f = np.array([1e6, 1e7, 1e8])
    radii = np.array([30e-6, 100e-6, 300e-6, 1e-3])
    err = []
    for w in radii:
        with pytest.warns(tw.NonPassiveObservableWarning):          # bare GK surface, T(0)
            h3 = tw.gaussian_response(f, stack, w, w)
        h1 = fm.response(2j * np.pi * f, sample) * 2 / (np.pi * 2 * w * w)
        err.append(np.max(np.abs(h3 / h1 - 1)))
    err = np.array(err)
    order = np.log(err[:-1] / err[1:]) / np.log(radii[1:] / radii[:-1])
    assert np.all(np.abs(order - 2) < 0.01), order
    assert err[-1] < 1e-5


# --------------------------------------------------------------------------
# (c) Cattaneo half-space phase
# --------------------------------------------------------------------------

def test_cattaneo_halfspace_phase_closed_form():
    tau = 1e-10
    omega = np.logspace(-3, 3, 13) / tau
    stack = tw.Stack((), tw.Material(60.0, 2.41e6, tau_r=tau))
    phase = np.angle(tw.kernel(omega / (2 * np.pi), 0.0, stack))
    assert_allclose(phase, -np.pi / 4 + 0.5 * np.arctan(omega * tau), atol=1e-12)
    f = np.array([1e6, 1e7, 1e8, 1e9])
    big = np.angle(tw.gaussian_response(f, stack, 1e-3, 1e-3))
    assert_allclose(big, -np.pi / 4 + 0.5 * np.arctan(2 * np.pi * f * tau), atol=2e-5)


# --------------------------------------------------------------------------
# (d) Independent ODE shooting
# --------------------------------------------------------------------------

def _rhs(m, k, p):
    g = 1 + m.tau_r * p
    C = m.rho_c
    if m.ell > 0:
        lam, l2, al = m.lam_z, m.ell ** 2, m.alpha
        lp = lam + (1 + al) * l2 * C * p

        def rhs(z, y):
            T, Qz, Qr, Dr = y
            return [-((g + l2 * k * k) * Qz + l2 * k * Dr) / lp,
                    -C * p * T - k * Qr,
                    Dr,
                    ((g + l2 * k * k) * Qr - (lam + al * l2 * C * p) * k * T) / l2]
    else:
        def rhs(z, y):
            T, Qz = y
            return [-g * Qz / m.lam_z, -(C * p + m.lam_r * k * k / g) * T]
    return rhs


def _theta(m, bc, y, k, p):
    if m.ell == 0:
        return y[0]
    c_n, c_t = bc.coefficients(m)
    return y[0] + c_n * (m.rho_c * p * y[0] + k * y[2]) - c_t * k * y[2]


def _extra(m, bc, y, face):
    a, b = bc.slip_weights()
    sign = 1.0 if face == "top" else -1.0
    return a * y[2] - sign * b * m.ell * y[3]


ADMITTANCE = 1e7     # W/(m^2 K): temperature scale of the shooting unknowns


def _unscale(m, yhat):
    """Scaled unknowns (Y T, Q_z, Q_r, ell D_r) -> physical state."""
    y = np.array(yhat, complex)
    y[0] /= ADMITTANCE
    if y.size == 4:
        y[3] /= m.ell
    return y


def _normalise(m, y):
    s = np.array(y, complex)
    s[0] *= ADMITTANCE
    if s.size == 4:
        s[3] *= m.ell
    return y / np.linalg.norm(s)


def _shoot(stack, f, k):
    """Integrate the physical ODEs upward; impose interface conditions on subspaces.

    Unknowns are scaled so that null vectors have comparable components; without
    this the SVD fixes temperatures only to absolute precision and loses digits.
    """
    p = 2j * np.pi * f
    sub = stack.substrate
    g = 1 + sub.tau_r * p
    gam = np.sqrt((sub.rho_c * p * g + sub.lam_r * k * k) / sub.lam_z)
    states = [np.array([1.0, sub.lam_z * gam / g], complex)]   # Fourier substrate
    lower, lower_bc = sub, stack.substrate_boundary
    for layer in reversed(stack.layers):
        m, bc, R = layer.material, layer.boundary, layer.resistance_below
        states = [_normalise(lower, y) for y in states]
        n_u, n_l = (4 if m.ell > 0 else 2), len(states)
        rows = []
        # unknowns: scaled upper bottom state (n_u), coefficients of lower states (n_l)
        def lin(fun):
            row = np.zeros(n_u + n_l, complex)
            for i in range(n_u):
                e = np.zeros(n_u, complex)
                e[i] = 1
                row[i] = fun(_unscale(m, e), None)
            for j in range(n_l):
                row[n_u + j] = fun(None, states[j])
            return row
        rows.append(lin(lambda u, l: u[1] if u is not None else -l[1]))
        rows.append(lin(lambda u, l: (_theta(m, bc, u, k, p) - R * u[1]) if u is not None
                        else -_theta(lower, lower_bc, l, k, p)))
        if m.ell > 0:
            rows.append(lin(lambda u, l: _extra(m, bc, u, "bottom") if u is not None else 0.0))
        if lower.ell > 0:
            rows.append(lin(lambda u, l: 0.0 if u is not None else _extra(lower, lower_bc, l, "top")))
        basis = null_space(np.array(rows))
        new = []
        # absolute tolerance on the physical scale of each component (scaled norm <= 1);
        # a vanishing atol makes roundoff-level components stall the step control
        atol = np.abs(_unscale(m, np.full(n_u, 1e-15)))
        for col in basis.T:
            y0 = _unscale(m, col[:n_u])
            sol = solve_ivp(_rhs(m, k, p), (layer.thickness, 0.0), y0, method="DOP853",
                            rtol=1e-12, atol=atol)
            assert sol.success
            new.append(sol.y[:, -1])
        states, lower, lower_bc = new, m, bc
    assert len(states) == 1 or stack.layers[0].material.ell > 0
    if stack.layers[0].material.ell > 0:
        raise NotImplementedError("bare GK top surface is checked against Hennessy-Myers")
    y = states[0]
    return y[0] / y[1]


@pytest.mark.parametrize("bc", [tw.GKBoundary(slip=0.0), tw.GKBoundary(slip=1.0),
                                tw.GKBoundary(slip=np.inf),
                                tw.GKBoundary(slip=1.0, temperature="entropic"),
                                tw.GKBoundary(slip=0.3, temperature="kinetic",
                                              c_n=4e-17, c_t=-2e-17)],
                         ids=["noslip", "diffuse", "freeslip", "entropic", "kinetic"])
@pytest.mark.parametrize("f,k", [(1e5, 3e5), (2e7, 1e6), (2e8, 0.0), (1e6, 2e6)])
def test_kernel_matches_ode_shooting(bc, f, k):
    # Single shooting amplifies integration errors by about exp(2 gamma_S d):
    # ell = 200 nm keeps gamma_S d ~ 2.5 (measured agreement 3e-9 to 2.4e-8;
    # 6e-8 to 1.9e-7 at ell = 120 nm). The arbitrary-precision test below covers
    # gamma_S d ~ 250.
    film = gk_aln(tau_r=5e-11, ell=200e-9)
    stack = tw.au_aln_sapphire(film, boundary=bc)
    assert_allclose(tw.kernel(f, k, stack), _shoot(stack, f, k), rtol=1e-7)


# --------------------------------------------------------------------------
# External closed form: GK half-space with slip (Hennessy and Myers 2021)
# --------------------------------------------------------------------------

def _hennessy_myers(k, f, m, C):
    """Eqs. (11)-(13) of the accepted manuscript, converted to exp(+i omega t).

    Per unit flux; general alpha through tau_l = (1+alpha) ell^2/a (the paper uses 2).
    """
    p = 2j * np.pi * f
    g = 1 + m.tau_r * p
    lp = m.lam_z + (1 + m.alpha) * m.ell ** 2 * m.rho_c * p
    gp = np.sqrt(k * k + m.rho_c * p * g / lp)
    gs = np.sqrt(k * k + g / m.ell ** 2)
    if np.isinf(C):
        phi = gs / (k * k * gp - gs * gp * gs) * gs
    else:
        phi = gs * (1 + C * m.ell * gs) / (k * k * (1 + C * m.ell * gp) - gs * gp * (1 + C * m.ell * gs))
    return -(1 / m.lam_z) * g / (1 + m.tau_l * p) * phi


@pytest.mark.parametrize("alpha", [1 / 3, 2.0])
@pytest.mark.parametrize("C", [0.0, 1.0, 5.0, np.inf])
def test_gk_halfspace_matches_hennessy_myers(alpha, C):
    # Silicon values of their Table 1 at 311 K, used only as a numerical benchmark.
    m = tw.Material(150.0, 1.692e6, tau_r=42e-12, ell=185e-9, alpha=alpha)
    stack = tw.Stack((), m, tw.GKBoundary(slip=C))
    f = np.geomspace(1e5, 1e10, 11)[:, None]
    k = np.array([0.0, 1e5, 1e6, 5e6, 3e7])[None, :]
    assert_allclose(tw.kernel(f, k, stack), _hennessy_myers(k, f, m, C), rtol=1e-12)


# --------------------------------------------------------------------------
# (e) alpha and the ell -> 0 limit
# --------------------------------------------------------------------------

@pytest.mark.parametrize("temperature", ["local", "entropic"])
@pytest.mark.parametrize("slip,min_3d_effect", [(0.0, 1e-2), (1.0, 5e-3), (np.inf, 1e-4)])
def test_alpha_enters_the_plane_wave_limit_only_through_tau_l(temperature, slip, min_3d_effect):
    """At fixed L^2 = (1+alpha) ell^2 (same tau_l) the k = 0 kernel cannot see alpha,
    for both temperatures (theta = (1 + tau_l p) T at k = 0); at k > 0 it can, through
    the tangential boundary layers of thickness ell (measured 2e-4 to 3e-2)."""
    ell1 = 90e-9
    ell2 = ell1 * np.sqrt((1 + 1 / 3) / 3)
    bc = tw.GKBoundary(slip=slip, temperature=temperature)
    a = tw.au_aln_sapphire(tw.Material(60.0, 2.41e6, tau_r=2e-11, ell=ell1, alpha=1 / 3), boundary=bc)
    b = tw.au_aln_sapphire(tw.Material(60.0, 2.41e6, tau_r=2e-11, ell=ell2, alpha=2.0), boundary=bc)
    assert a.layers[1].material.tau_l == pytest.approx(b.layers[1].material.tau_l, rel=1e-14)
    f = np.geomspace(1e5, 1e9, 9)
    assert_allclose(tw.kernel(f, 0.0, a), tw.kernel(f, 0.0, b), rtol=1e-13)
    diff = np.max(np.abs(tw.kernel(f, 3e5, a) / tw.kernel(f, 3e5, b) - 1))
    assert diff > min_3d_effect


@pytest.mark.parametrize("slip,order_cattaneo", [(0.0, 1.0), (1.0, 1.0), (np.inf, 2.0)])
def test_gk_tends_to_cattaneo_as_ell_vanishes(slip, order_cattaneo):
    """alpha enters through alpha ell^2: its effect is O(ell^2). Finite slip leaves
    an O(ell) boundary layer (in-plane Poiseuille reduction), free slip does not."""
    f, k = np.array([1e5, 1e7, 2e8]), 2e5
    cat = tw.kernel(f, k, tw.au_aln_sapphire(tw.Material(60.0, 2.41e6, tau_r=2e-11)))
    ells = np.array([20e-9, 10e-9, 5e-9, 2.5e-9])
    d_alpha, d_cat = [], []
    for ell in ells:
        z = {al: tw.kernel(f, k, tw.au_aln_sapphire(gk_aln(2e-11, ell, al),
                                                    boundary=tw.GKBoundary(slip=slip)))
             for al in (1 / 3, 2.0)}
        d_alpha.append(np.max(np.abs(z[2.0] / z[1 / 3] - 1)))
        d_cat.append(np.max(np.abs(z[1 / 3] / cat - 1)))
    o_alpha = np.log2(d_alpha[-2] / d_alpha[-1])
    o_cat = np.log2(d_cat[-2] / d_cat[-1])
    assert abs(o_alpha - 2) < 0.05, o_alpha
    assert abs(o_cat - order_cattaneo) < 0.05, o_cat


# --------------------------------------------------------------------------
# (f) Interface continuity, energy balance and dissipation
# --------------------------------------------------------------------------

BCS = [tw.GKBoundary(slip=0.0), tw.GKBoundary(slip=1.0), tw.GKBoundary(slip=np.inf),
       tw.GKBoundary(slip=1.0, temperature="entropic"),
       tw.GKBoundary(slip=0.0, temperature="entropic")]


@pytest.mark.parametrize("bc", BCS)
@pytest.mark.parametrize("f,k", [(1e5, 3e5), (2e7, 1e6), (2e8, 0.0), (1e6, 1e3)])
def test_continuity_energy_and_dissipation_identity(bc, f, k):
    stack = tw.au_aln_sapphire(gk_aln(2e-11, 80e-9), boundary=bc)
    d = tw.dissipation_balance(f, k, stack)
    assert max(i["flux_jump"] for i in d["interfaces"]) < 1e-12
    assert max(d["energy_residuals"]) < 1e-12
    assert d["residual"] < 1e-12
    # the depth quadrature is converged: doubling the points per panel changes nothing
    d32 = tw.dissipation_balance(f, k, stack, n=32)
    assert max(abs(a / b - 1) for a, b in zip(d["bulk"], d32["bulk"])) < 1e-12


@pytest.mark.parametrize("slip", [0.0, 1.0, np.inf])
def test_entropic_convention_has_nonnegative_dissipation(slip):
    stack = tw.au_aln_sapphire(gk_aln(5e-11, 150e-9),
                               boundary=tw.GKBoundary(slip=slip, temperature="entropic"))
    for f, k in ((1e4, 1e5), (3e6, 2e6), (2e8, 0.0), (2e8, 5e6)):
        d = tw.dissipation_balance(f, k, stack)
        for i in d["interfaces"]:
            assert abs(i["nonequilibrium"]) < 1e-12 * d["re_impedance"]
            assert min(i["tangential_above"], i["tangential_below"]) >= -1e-12 * d["re_impedance"]
    z = tw.kernel(np.geomspace(1e3, 1e10, 30)[:, None], np.geomspace(1e2, 1e8, 30)[None, :], stack)
    assert np.all(z.real > 0)
    assert tw.passivity_margin(np.geomspace(1e3, 1e10, 30)[:, None],
                               np.geomspace(1e2, 1e8, 30)[None, :], stack) > 0


def test_local_convention_interface_term_has_no_sign():
    """Documented property: the local temperature does not make each interface dissipative."""
    stack = tw.au_aln_sapphire(gk_aln(2e-11, 80e-9))
    d = tw.dissipation_balance(2e8, 0.0, stack)
    assert d["interfaces"][0]["nonequilibrium"] < -0.01 * d["re_impedance"]
    assert d["residual"] < 1e-12


def test_normal_gradient_condition_is_inadmissible():
    film = gk_aln(2e-11, 80e-9)
    bad = tw.au_aln_sapphire(film, boundary=tw.GKBoundary(tangential="normal_gradient"))
    with pytest.raises(np.linalg.LinAlgError):
        tw.kernel(1e6, 0.0, bad)
    d = tw.dissipation_balance(2e7, 1e6, bad)
    tangential = sum(i["tangential_above"] + i["tangential_below"] for i in d["interfaces"])
    assert tangential < -0.5 * d["re_impedance"]          # boundary terms of the wrong sign
    assert d["residual"] < 1e-12
    good0 = tw.kernel(1e6, 0.0, tw.au_aln_sapphire(film))
    assert abs(tw.kernel(1e6, 1.0, bad) / good0 - 1) > 0.01   # k -> 0 limit is not the 1D model


# --------------------------------------------------------------------------
# Physics limits
# --------------------------------------------------------------------------

def test_fourier_resonance_depends_on_temperature_convention():
    """tau_R = tau_l makes a GK film under a transducer indistinguishable from
    Fourier at k = 0 with the local temperature; the entropic temperature puts
    (1 + tau_l p) factors on both film faces and removes the resonance
    (measured: 2.6 % in |Z|, 0.64 degree in phase, ell = 80 nm)."""
    film = gk_aln(tau_r=1.0, ell=80e-9)
    film = tw.Material(60.0, 2.41e6, tau_r=film.tau_l, ell=80e-9)      # tau_R = tau_l
    f = np.geomspace(1e5, 1e9, 9)
    zf = tw.kernel(f, 0.0, tw.au_aln_sapphire(tw.Material(60.0, 2.41e6)))
    local = tw.au_aln_sapphire(film)
    entropic = tw.au_aln_sapphire(film, boundary=tw.GKBoundary(temperature="entropic"))
    assert_allclose(tw.kernel(f, 0.0, local), zf, rtol=1e-12)
    ratio = tw.kernel(f, 0.0, entropic) / zf
    assert np.max(np.abs(ratio - 1)) > 1e-2
    assert np.max(np.abs(np.degrees(np.angle(ratio)))) > 0.5


@pytest.mark.parametrize("slip", [0.0, 1.0])
def test_thin_film_poiseuille_limit(slip):
    """At low frequency a GK film acts as an anisotropic film with the slip-flow
    in-plane conductivity lam [1 - (2 ell/d) tanh(d/2ell)/(1 + C tanh(d/2ell))]."""
    d, ell, lam = 500e-9, 50e-9, 60.0
    t = np.tanh(d / (2 * ell))
    lam_r = lam * (1 - (2 * ell / d) * t / (1 + slip * t))
    f = np.geomspace(1e4, 3e5, 5)
    gk = tw.gaussian_response(f, tw.au_aln_sapphire(gk_aln(2e-12, ell), aln_thickness=d,
                                                    boundary=tw.GKBoundary(slip=slip)), 10e-6, 8e-6)
    iso = tw.gaussian_response(f, tw.au_aln_sapphire(tw.Material(lam, 2.41e6, tau_r=2e-12),
                                                     aln_thickness=d), 10e-6, 8e-6)
    sur = tw.gaussian_response(f, tw.au_aln_sapphire(tw.Material(lam, 2.41e6, lam_r=lam_r, tau_r=2e-12),
                                                     aln_thickness=d), 10e-6, 8e-6)
    dev = np.degrees(np.abs(np.angle(gk / sur)))
    sig = np.degrees(np.abs(np.angle(gk / iso)))
    assert np.max(dev) < 2e-3
    assert np.max(dev) < 0.01 * np.max(sig)


# --------------------------------------------------------------------------
# (g) Hankel quadrature
# --------------------------------------------------------------------------

def test_trapezoid_rule_converges_and_matches_adaptive():
    stack = tw.au_aln_sapphire(gk_aln(5e-11, 150e-9),
                               boundary=tw.GKBoundary(slip=1.0, temperature="entropic"))
    ref = tw.gaussian_response(FREQ, stack, 10e-6, 8e-6, step=0.0125)
    errs = [np.max(np.abs(tw.gaussian_response(FREQ, stack, 10e-6, 8e-6, step=h) / ref - 1))
            for h in (0.4, 0.2, 0.1)]
    assert errs[0] > errs[1] > errs[2]
    assert errs[1] / errs[2] > 1e3                      # exponential, not algebraic
    default = tw.gaussian_response(FREQ, stack, 10e-6, 8e-6)
    assert np.max(np.abs(default / ref - 1)) < 1e-12
    low = tw.gaussian_response(FREQ, stack, 10e-6, 8e-6, x_low=-30.0)
    assert np.max(np.abs(default / low - 1)) < 1e-13
    adaptive = tw.gaussian_response_adaptive(FREQ, stack, 10e-6, 8e-6)
    assert np.max(np.abs(adaptive / ref - 1)) < 1e-10


# --------------------------------------------------------------------------
# Arbitrary precision where exp(gamma_S d) ~ 1e108
# --------------------------------------------------------------------------

def _mp_kernel(stack, f, k, dps=220):
    """Shooting with exact propagators expm(-A d) in 220-digit arithmetic."""
    import mpmath as mp
    mp.mp.dps = dps
    p = mp.mpc(0, 2 * mp.pi * f)
    kk = mp.mpf(k)

    def generator(m):
        g = 1 + mp.mpf(m.tau_r) * p
        C = mp.mpf(m.rho_c)
        if m.ell > 0:
            lam, l2, al = mp.mpf(m.lam_z), mp.mpf(m.ell) ** 2, mp.mpf(m.alpha)
            lp = lam + (1 + al) * l2 * C * p
            return mp.matrix([[0, -(g + l2 * kk ** 2) / lp, 0, -l2 * kk / lp],
                              [-C * p, 0, -kk, 0],
                              [0, 0, 0, 1],
                              [-(lam + al * l2 * C * p) * kk / l2, 0, (g + l2 * kk ** 2) / l2, 0]])
        return mp.matrix([[0, -g / mp.mpf(m.lam_z)], [-(C * p + mp.mpf(m.lam_r) * kk ** 2 / g), 0]])

    au, aln, sub = stack.layers[0], stack.layers[1], stack.substrate
    bc, m = aln.boundary, aln.material
    a_w, b_w = (mp.mpf(x) for x in bc.slip_weights())
    c_n, c_t = (mp.mpf(x) for x in bc.coefficients(m))
    Cp = mp.mpf(m.rho_c) * p
    ys = mp.sqrt((mp.mpf(sub.rho_c) * p + mp.mpf(sub.lam_r) * kk ** 2) / mp.mpf(sub.lam_z)) * mp.mpf(sub.lam_z)
    R2, R1 = mp.mpf(aln.resistance_below), mp.mpf(au.resistance_below)
    ell = mp.mpf(m.ell)
    # AlN bottom: admissible states parametrised by substrate amplitude and Q_r (or D_r)
    T1 = (1 + R2 * ys) / (1 + c_n * Cp)
    s1 = mp.matrix([T1, ys, 0, 0])
    if b_w != 0:
        s2 = mp.matrix([-(c_n - c_t) * kk / (1 + c_n * Cp), 0, 1, -a_w / (b_w * ell)])
    else:
        s2 = mp.matrix([0, 0, 0, 1])
    phi = mp.expm(-generator(m) * mp.mpf(aln.thickness))
    w1, w2 = phi * s1, phi * s2
    e1 = a_w * w1[2] - b_w * ell * w1[3]
    e2 = a_w * w2[2] - b_w * ell * w2[3]
    w = w1 - (e1 / e2) * w2
    theta = w[0] + c_n * (Cp * w[0] + kk * w[2]) - c_t * kk * w[2]
    yau = mp.matrix([theta + R1 * w[1], w[1]])
    top = mp.expm(-generator(au.material) * mp.mpf(au.thickness)) * yau
    return complex(top[0] / top[1])


@pytest.mark.parametrize("temperature", ["local", "entropic"])
def test_extreme_boundary_layer_against_arbitrary_precision(temperature):
    film = gk_aln(tau_r=2e-11, ell=2e-9)             # gamma_S d ~ 250
    stack = tw.au_aln_sapphire(film, boundary=tw.GKBoundary(slip=1.0, temperature=temperature))
    for f, k in ((2e8, 1e6), (1e5, 2e5)):
        assert_allclose(tw.kernel(f, k, stack), _mp_kernel(stack, f, k), rtol=1e-12)


# --------------------------------------------------------------------------
# Instrument, synthetic data and input validation
# --------------------------------------------------------------------------

def test_instrument_and_synthetic_data():
    stack = tw.au_aln_sapphire(gk_aln())
    f = np.geomspace(1e4, 2e8, 6)
    h = tw.gaussian_response(f, stack, 10e-6, 8e-6)
    obs = tw.instrument(f, h, log_gain=0.1, phase_offset_deg=2.0, delay=1e-9)
    assert_allclose(obs / h, np.exp(0.1 + 1j * (np.deg2rad(2.0) - 2 * np.pi * f * 1e-9)), rtol=1e-14)
    a = tw.synthetic_measurement(f, stack, 10e-6, 8e-6, 0.1, 2.0, 1e-9, seed=20261002)
    b = tw.synthetic_measurement(f, stack, 10e-6, 8e-6, 0.1, 2.0, 1e-9, seed=20261002)
    assert np.array_equal(a["observed"], b["observed"])
    assert_allclose(a["truth"], obs, rtol=1e-14)
    assert a["status"] == "SYNTHETIC ONLY"


@pytest.mark.parametrize("kwargs", [dict(lam_z=-1.0, rho_c=1.0), dict(lam_z=1.0, rho_c=0.0),
                                    dict(lam_z=1.0, rho_c=1.0, lam_r=2.0, ell=1e-9),
                                    dict(lam_z=1.0, rho_c=1.0, alpha=-0.5),
                                    dict(lam_z=1.0, rho_c=1.0, tau_r=np.nan)])
def test_nonphysical_material_rejected(kwargs):
    with pytest.raises(ValueError):
        tw.Material(**kwargs)


def test_invalid_boundary_and_frequency_rejected():
    with pytest.raises(ValueError):
        tw.GKBoundary(tangential="stick")
    with pytest.raises(ValueError):
        tw.GKBoundary(temperature="hot")
    with pytest.raises(ValueError):
        tw.GKBoundary(slip=-1.0)
    with pytest.raises(ValueError):
        tw.GKBoundary(c_n=1e-17)                       # coefficients need "kinetic"
    with pytest.raises(ValueError):
        tw.kernel([0.0], 0.0, tw.Stack((), SAPPHIRE))
    with pytest.raises(ValueError):
        tw.gaussian_response([1e6], tw.Stack((), SAPPHIRE), -1.0, 1e-6)


# --------------------------------------------------------------------------
# Configurations without a transducer or with a GK substrate, against an
# independent arbitrary-precision solver (eigen-decomposition modes, physical
# boundary functionals, no equilibration; naive and stable referencing agree
# to <= 6e-17). Campaign file experiments/mpref_regression_constants.py.
# --------------------------------------------------------------------------

def _gk_layer(tau, ell, alpha, d, r, slip, temp, c=None, lam=60.0, rho_c=2.41e6):
    m = tw.Material(lam, rho_c, tau_r=tau, ell=ell, alpha=alpha)
    if temp == "kinetic":
        w = ell ** 2 / lam
        bc = tw.GKBoundary(slip=slip, temperature="kinetic", c_n=c[0] * (1 + alpha) * w,
                           c_t=c[1] * alpha * w)
    else:
        bc = tw.GKBoundary(slip=slip, temperature=temp)
    return tw.Layer(m, d, r, bc)


def _gk_substrate(tau, ell, alpha, slip, temp):
    return (tw.Material(35.0, 3.03e6, tau_r=tau, ell=ell, alpha=alpha),
            tw.GKBoundary(slip=slip, temperature=temp))


def _reference_stack(name):
    au = tw.Layer(AU, 80e-9, 1e-8)
    specs = {
        "A1": ([_gk_layer(2e-11, 80e-9, 1 / 3, 500e-9, 2e-8, 1.0, "entropic")], None),
        "A2": ([_gk_layer(5e-10, 300e-9, 2.0, 200e-9, 1e-8, 0.0, "local")], None),
        "A3": ([_gk_layer(1e-11, 40e-9, 1 / 3, 1e-6, 5e-9, np.inf, "kinetic", (0.1, -0.2))], None),
        "C1": ([au, _gk_layer(2e-11, 80e-9, 1 / 3, 500e-9, 2e-8, 1.0, "entropic")],
               _gk_substrate(1e-11, 30e-9, 1 / 3, 1.0, "entropic")),
        "C2": ([au, _gk_layer(1e-10, 150e-9, 2.0, 2e-6, 2e-8, 0.0, "local")],
               _gk_substrate(1e-10, 100e-9, 2.0, np.inf, "local")),
        "D1": ([_gk_layer(3e-11, 100e-9, 1 / 3, 300e-9, 1e-8, 1.0, "entropic")],
               _gk_substrate(1e-11, 50e-9, 1 / 3, 1.0, "entropic")),
        "D2": ([_gk_layer(2e-10, 250e-9, 2.0, 1e-6, 0.0, np.inf, "local")],
               _gk_substrate(5e-11, 20e-9, 1 / 3, 0.0, "local")),
    }
    layers, sub = specs[name]
    if sub is None:
        return tw.Stack(tuple(layers), SAPPHIRE)
    return tw.Stack(tuple(layers), sub[0], sub[1])


REFERENCE_POINTS = [(1e5, 0.0), (1e5, 3e5), (2e7, 1e6), (2e8, 5e6), (1e9, 2e6)]
REFERENCE = {
    "A1": [complex(1.0111984686697931e-07, -8.896295068064875e-08),
           complex(8.979766073244263e-08, -2.0681835037561133e-08),
           complex(4.884965639935458e-09, -6.381176676635955e-09),
           complex(1.5995885432901783e-09, -1.7717308067210862e-09),
           complex(1.9287330491540555e-10, -6.635633082432542e-10)],
    "A2": [complex(9.491738479328286e-08, -8.713695827674422e-08),
           complex(9.646300678339567e-08, -2.608064523569842e-08),
           complex(6.147231413721437e-09, -9.691768568861905e-09),
           complex(2.351177491901836e-10, -1.727668991785549e-09),
           complex(6.004513040827772e-11, -3.4799485374912277e-10)],
    "A3": [complex(8.506412965459948e-08, -8.69649768777347e-08),
           complex(7.037880408043134e-08, -1.5983433129603456e-08),
           complex(5.533165462516569e-09, -4.6861664199598995e-09),
           complex(1.8498439199867326e-09, -1.2986500115512212e-09),
           complex(5.393616766845674e-10, -8.276985774973959e-10)],
    "C1": [complex(1.0921520288574914e-07, -8.946226580714873e-08),
           complex(9.070127395151189e-08, -1.7913626206158228e-08),
           complex(1.0166461817520941e-08, -6.935728794168709e-09),
           complex(2.006819045079791e-09, -1.244902101187782e-09),
           complex(2.604304131005964e-10, -8.086030744173444e-10)],
    "C2": [complex(9.588726153635419e-08, -9.071991634583848e-08),
           complex(7.367812936511044e-08, -1.3430628203611196e-08),
           complex(1.0032392632315086e-08, -6.566461901386824e-09),
           complex(1.963274854631115e-09, -1.1915067713414448e-09),
           complex(2.6830520609655836e-10, -8.081469259738396e-10)],
    "D1": [complex(9.406724507477666e-08, -8.733858656656772e-08),
           complex(8.957398957469209e-08, -2.2842244257952454e-08),
           complex(5.9732272850451465e-09, -7.603537140681453e-09),
           complex(1.3348900052344055e-09, -1.9294630842191908e-09),
           complex(1.317295775087311e-10, -5.607579878993926e-10)],
    "D2": [complex(8.134047811319391e-08, -8.604011876609039e-08),
           complex(6.863305894708222e-08, -1.6209147280544353e-08),
           complex(2.945859853544597e-09, -5.43066528744785e-09),
           complex(1.4032628505685796e-10, -8.085653064027124e-10),
           complex(8.643869090885035e-11, -1.628208546144895e-10)],
}


@pytest.mark.parametrize("name", sorted(REFERENCE))
def test_configurations_without_transducer_or_with_gk_substrate(name):
    stack = _reference_stack(name)
    got = [complex(tw.kernel(f, k, stack)) for f, k in REFERENCE_POINTS]
    assert_allclose(got, REFERENCE[name], rtol=1e-12)


# --------------------------------------------------------------------------
# Diagnostics: anisotropic Fourier layers, underflow, wave-like layers
# --------------------------------------------------------------------------

@pytest.mark.parametrize("f,k", [(1e5, 1e6), (2e7, 3e5), (2e8, 5e6)])
def test_anisotropic_fourier_layer_balances(f, k):
    """The radial flux of a Fourier layer is lam_r k T / g: the lateral term of the
    energy balance and the bulk dissipation |Q_r|^2/lam_r depend on it."""
    stack = tw.Stack((tw.Layer(AU, 80e-9, 1e-8),
                      tw.Layer(tw.Material(60.0, 2.41e6, lam_r=150.0), 500e-9, 2e-8)), SAPPHIRE)
    d = tw.dissipation_balance(f, k, stack)
    assert max(d["energy_residuals"]) < 1e-12
    assert d["residual"] < 1e-12


@pytest.mark.parametrize("d,k", [(4e-6, 5e8), (1e-6, 1e9)])
def test_dissipation_balance_with_underflowing_fields(d, k):
    stack = tw.au_aln_sapphire(gk_aln(2e-11, 100e-9), aln_thickness=d,
                               boundary=tw.GKBoundary(temperature="entropic"))
    out = tw.dissipation_balance(1e6, k, stack)
    assert out["residual"] < 1e-12
    assert np.all(np.isfinite(out["energy_residuals"])) and max(out["energy_residuals"]) < 1e-10
    assert max(i["flux_jump"] for i in out["interfaces"]) < 1e-10


@pytest.mark.parametrize("f", [7.13e8, 1e9])
def test_dissipation_balance_resolves_wave_like_layers(f):
    """omega tau_R ~ 50 in a 4.28 um film (~80 wavelengths): default settings must close."""
    film = tw.Material(60.0, 2.41e6, tau_r=8.59e-9, ell=4.21e-9)
    stack = tw.au_aln_sapphire(film, aln_thickness=4.28e-6,
                               boundary=tw.GKBoundary(temperature="entropic"))
    out = tw.dissipation_balance(f, 0.0, stack)
    assert out["residual"] < 1e-12
    assert max(out["energy_residuals"]) < 1e-12


# --------------------------------------------------------------------------
# Passivity and the observable of stacks without a transducer
# --------------------------------------------------------------------------

def _kinetic_counterexample(temperature="kinetic"):
    """Member found by an adversarial search: 50 nm Au transducer, Re Z/|Z| = -0.94."""
    ell, alpha, lam = 9.7e-7, 1.95, 60.0
    w = ell ** 2 / lam
    film = tw.Material(lam, 2.41e6, tau_r=1.8e-11, ell=ell, alpha=alpha)
    if temperature == "kinetic":
        bc = tw.GKBoundary(slip=7.0, temperature="kinetic", c_n=0.05 * (1 + alpha) * w, c_t=0.97 * w)
    else:
        bc = tw.GKBoundary(slip=7.0, temperature=temperature)
    return tw.Stack((tw.Layer(tw.Material(150.0, 2.49e6), 50e-9, 1.13e-9),
                     tw.Layer(film, 9.13e-6, 1.13e-9, bc)), SAPPHIRE)


def test_kinetic_temperature_can_be_non_passive_and_is_flagged():
    ks = np.geomspace(1e4, 1e8, 801)
    bad = _kinetic_counterexample()
    assert tw.passivity_margin(5.83e8, ks, bad) < -0.9
    with pytest.warns(tw.PassivityWarning):
        h = tw.gaussian_response([5.83e8], bad, 2e-6, 2e-6)
    assert_allclose(h, [5.1349 - 55.009j], rtol=1e-4)        # value of the adversarial search
    with pytest.raises(tw.PassivityError):
        tw.gaussian_response([5.83e8], bad, 2e-6, 2e-6, passivity="raise")
    good = _kinetic_counterexample("entropic")
    assert tw.passivity_margin(5.83e8, ks, good) > 0
    with warnings.catch_warnings():
        warnings.simplefilter("error", tw.PassivityWarning)
        tw.gaussian_response([5.83e8], good, 2e-6, 2e-6)


def test_gk_top_surface_temperature_is_not_the_passive_variable():
    """Bare GK film (entropic): Re T(0) < 0 while the bounded variable theta_e(0) > 0."""
    film = tw.Material(60.0, 2.41e6, tau_r=1e-13, ell=1e-6, alpha=2.0)
    stack = tw.Stack((tw.Layer(film, 5.01e-9, 1e-11,
                               tw.GKBoundary(slip=0.13, temperature="entropic")),), SAPPHIRE)
    assert stack.top_is_gk
    zT = complex(tw.kernel(1.38e8, 1.4e3, stack))
    zth = complex(tw.kernel(1.38e8, 1.4e3, stack, observable="theta_e"))
    assert zT.real < -0.5 * abs(zT)
    assert zth.real > 0.5 * abs(zth)
    d = tw.dissipation_balance(1.38e8, 1.4e3, stack)
    assert d["surface_input"] == pytest.approx(zth.real, rel=1e-10)
    assert d["residual"] < 1e-12
    with pytest.warns(tw.NonPassiveObservableWarning):
        tw.gaussian_response([1.38e8], stack, 2e-6, 2e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("error", tw.NonPassiveObservableWarning)
        tw.gaussian_response([1.38e8], stack, 2e-6, 2e-6, observable="theta_e")


# --------------------------------------------------------------------------
# Quadrature step chosen per frequency inside the Gaussian window
# --------------------------------------------------------------------------

def test_quadrature_step_follows_the_gaussian_window():
    """tau_R = 10 ns: the wave branch point narrows the strip to 0.008 rad but lies far
    outside the window of 10/8 um beams, so the default rule stays small; with 100 nm
    beams it lies inside and the narrow strip is kept."""
    stack = tw.au_aln_sapphire(tw.Material(60.0, 2.41e6, tau_r=1e-8))
    f = np.geomspace(1e4, 1e9, 60)
    eta = (10e-6 ** 2 + 8e-6 ** 2) / 8
    assert tw.quadrature_strip(f, stack).min() < 0.01
    assert tw.quadrature_strip(f, stack, eta=eta).min() > 0.5
    assert tw.trapezoid_nodes(f, stack, 10e-6, 8e-6)[0].size < 1000
    assert_allclose(tw.gaussian_response(f, stack, 10e-6, 8e-6),
                    tw.gaussian_response(f, stack, 10e-6, 8e-6, step=0.01), rtol=1e-12)
    half = tw.Stack((), tw.Material(60.0, 2.41e6, tau_r=1e-8))
    for w in (100e-9, 240e-9):
        assert_allclose(tw.gaussian_response([1e9], half, w, w),
                        tw.gaussian_response_adaptive([1e9], half, w, w, epsrel=1e-12), rtol=1e-11)
