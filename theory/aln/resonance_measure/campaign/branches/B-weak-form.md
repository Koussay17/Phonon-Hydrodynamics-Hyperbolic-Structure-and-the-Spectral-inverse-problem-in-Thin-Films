# B — Entropy weak form: a constructive finite scheme and a positive-weight obstruction

Independent first pass, 2026-09-27. Read only the campaign question and assumptions. **First pass complete, 2026-09-28. Finite conditional derivations and bounded checks; independent adversarial audit not yet incorporated. Not an AlN result, novelty claim, or general convergence theorem.**

## 1. Domain, represented quantities, and finite equation

Let the mode domain be finitely many copies of a periodic Brillouin torus, with continuous physical energy epsilon >= epsilon_0 > 0. Fix a scalar Bose kinetic approximation and one consistent reversible-event counting convention. A reaction r has q_p=q_a+q_b modulo the reciprocal lattice and Delta(r)=epsilon_p-epsilon_a-epsilon_b. Initially restrict to smooth regular zero sets, |grad Delta| >= g_0 > 0, and a nonnegative integrable interaction weight. Acoustic zeros, crossings and critical resonances are excluded here.

Choose a real finite-dimensional space V=span{phi_i}, containing the **physical energy function itself**: epsilon=sum_i e_i phi_i. Adding epsilon as an extra basis function is permitted; interpolation of its nodal values alone need not achieve this. Let Q_v be a positive volume quadrature with nodes z_l, weights m_l>0, and full-column-rank evaluation matrix Phi_li=phi_i(z_l). Choose positive reaction quadrature weights omega_r and evaluation points p_r,a_r,b_r; initially these are on the true physical resonance surface. Define

    xi_alpha(z)=sum_i alpha_i phi_i(z),
    n_alpha(z)=1/(exp(xi_alpha(z))-1),
    b_ri=phi_i(p_r)-phi_i(a_r)-phi_i(b_r).

Work on the open coefficient domain where xi_alpha is positive at all volume and reaction evaluation points. Occupations are reconstructed by interpolating xi=log((1+n)/n), **not by interpolating n**. At a reaction put

    A_r=n_p(1+n_a)(1+n_b),  B_r=(1+n_p)n_a n_b,
    Lambda_r=(A_r-B_r)/(log A_r-log B_r)>0,

with the continuous value Lambda=A=B when A=B. Then A-B=-Lambda b_r.alpha exactly. Define retained moments U_i=Q_v(phi_i n_alpha) and

    M_ij(alpha)=Q_v[n_alpha(1+n_alpha) phi_i phi_j],
    K(alpha)=sum_r omega_r Lambda_r b_r b_r^T,
    dU/dalpha=-M,   dot U=K alpha,   M dot alpha=-K alpha.       (1)

The sign agrees with parent loss: dot U_i=-sum_r omega_r(A_r-B_r)b_ri. M is positive definite and K is symmetric positive semidefinite. On the stated open domain both are smooth, giving a unique local ODE solution. No global domain-invariance theorem is asserted.

## 2. Exact finite identities and the obstruction

Use the quadrature entropy and energy

    S=Q_v[(1+n)log(1+n)-n log n],  E_Q=Q_v[epsilon n]=e^T U.

Since the entropy derivative is xi,

    dot S=alpha^T K alpha >= 0,
    dot E_Q=e^T K alpha,
    e^T K e=sum_r omega_r Lambda_r [b_r.e]^2.                 (2)

Consequently, for strictly positive active weights:

**Finite lemma.** K e=0 if and only if b_r.e=0 for every active reaction. Because epsilon belongs to V, b_r.e is the true physical Delta at that evaluation triple. Hence true resonance nodes give exact energy conservation for every state in the local solution domain. At any beta>0, alpha=beta e gives the physical Bose equilibrium at every evaluation point and is stationary.

**Obstruction.** Keeping conventional event evaluations on off-shell physical tuples and merely fitting their nonnegative weights cannot enforce this nullspace unless all off-shell weights vanish. At the physical Bose state,

    dot E_Q=beta sum_r omega_r Lambda_r Delta_r^2,

which is strictly positive whenever an off-shell tuple has positive weight. Opposite signs of Delta cannot cancel this drift. Thus even nonlinear equilibrium stationarity fails in this class; a symmetric positive linearized matrix alone does not establish it. This is a support obstruction, not a prohibition on all structure-preserving methods.

At a resonant equilibrium the linearized equation is

    M_* dot(delta alpha)=-K_* delta alpha,
    C=M_*^(-1/2) K_* M_*^(-1/2) >= 0.

The energy null vector is M_*^(1/2)e. Additional collision invariants, a sparse reaction graph, or an undersampled basis can produce extra null vectors; no uniqueness of thermal equilibrium is claimed.

With imperfect roots the same calculation gives

    |dot E_Q| <= [sum_r omega_r Lambda_r Delta_r^2]^(1/2)
                 [dot S]^(1/2).

If |Delta_r|<=tau, the first factor is at most tau sqrt(sum omega Lambda). At the Bose state the energy drift is quadratic in tau; the full equilibrium vector-field defect need only be first order.

## 3. Physical energy, surrogate energy, and an interpolation counterexample

Equation (2) conserves a quadrature-defined energy using the physical dispersion. It is not automatically the exact continuum integral integral epsilon n_alpha. Using exact volume integrals instead of Q_v gives a finite-dimensional variational equation conserving that continuum functional; using finite Q_v requires a separate volume-quadrature consistency check.

If one replaces epsilon by a finite interpolant epsilon_h and solves Delta_h=0, the same argument conserves **surrogate energy** and gives Bose equilibrium for epsilon_h. This does not prove exact conservation or detailed balance for epsilon. Similarly, roots of the true Delta do not help if the energy used by the test space evaluates as a different epsilon_h.

A minimal exact counterexample isolates occupation interpolation. Take beta=log 2; daughter energies are 1,1 and the parent energy 2 is the midpoint of interpolation nodes of energy 1 and 3. Physical resonance and interpolated energy are exact. Nodal Bose populations are 1 and 1/7. Linear occupation interpolation gives n_p=4/7 rather than 1/3, and the equilibrium decay bracket is

    n_p(1+2 n_a)-n_a^2=3(4/7)-1=5/7 != 0.

Interpolating xi instead gives xi_p=2 log 2 and the correct n_p=1/3, so the bracket vanishes. For the two nodal entropy coefficients, the event row is b=(-3/2,1/2), e=(1,3), and b.e=0. With unit volume weights, M_*=diag(2,8/49), Lambda_*=4/3; the two eigenvalues of C are 0 and 85/24. This is a finite interpolation example, not a physical two-mode reaction: fractional weak-form redistribution replaces individual nodal-event incidence.

## 4. Positive surface quadrature and an analytic normalization test

A constructive way to respect the obstruction is to generate off-grid points on Delta=0 and integrate the coarea measure

    a(r) delta(Delta(r)) dr
      = a(r)/|grad Delta(r)| dH^(dimension-1)(r).

Positive surface weights then give (1)-(2). Root location, surface Jacobians and coverage are independent accuracy requirements; setting Delta=0 does not determine the measure normalization.

For an exactly solvable periodic benchmark take q_a,q_b in [0,2pi), q_p=q_a+q_b modulo 2pi, daughter energies identically 1, and parent energy

    epsilon_p(q)=2+d+u sin q,   u>0,  |d|<u,  2+d-u>0.

Then Delta=d+u sin q_p. If t_1,t_2 are its two roots, an exact change of variables gives

    integral f delta(Delta) dq_a dq_b
      =1/sqrt(u^2-d^2) sum_{j=1,2}
           integral_0^(2pi) f(x,t_j-x) dx.                  (3)

For unnormalized torus Lebesgue measure the total mass is 4pi/sqrt(u^2-d^2). Equally spaced x_k=2pi k/N, each root receiving weight 2pi/[N sqrt(u^2-d^2)], give positive quadrature with exact total mass for every N. For d=0 and f=cos^2(q_a), the value is 2pi/u and this rule is exact for N>=3. For any fixed smooth periodic f, ordinary periodic quadrature gives convergence; no statement about evolving numerical solutions follows from this alone.

The benchmark has an empty resonance set for |d|>u and a divergent mass as |d| approaches u from below. At |d|=u the regular-surface hypotheses fail. A method that keeps a finite, smoothly varying limiting mass there would fail this test. With normalized torus volume the above expressions must be divided by (2pi)^2. Energies u,d carry energy units; delta and the total measure factor carry inverse-energy units.

## 5. Convex interpretation and limits of the candidate

The convex potential

    Psi(alpha)=-Q_v log(1-exp(-xi_alpha))

satisfies grad Psi=-U and Hess Psi=M. Locally its entropy dual has gradient alpha, so (1) is the finite entropy-gradient flow dot U=K grad_U S. This provides a precise variational reason to use entropy interpolation. It does not prove invariant geometric probability families in a Fock model; the scalar kinetic equation remains an approximation.

Positive moment-fitting quadrature is useful **after** restricting candidate nodes to the resonance surface. Without resonant support the nonnegative square constraint in (2) makes exact energy preservation infeasible. Projecting an arbitrary matrix to remove its energy direction can preserve linear PSD, but changes the event weak form and provides neither nonlinear Bose balance nor convergence to the physical resonance measure by itself.

For consistency, one still needs dense approximation spaces, convergent volume and surface quadratures, bounded positive reconstructed xi, and uniform control of the regular-surface Jacobian. The present candidate proves finite structural identities and supplies a falsifiable normalization benchmark. A general error estimate, global ODE-domain preservation, branch-crossing treatment, singular-resonance limit, and material calculation remain open.


## 6. Completed bounded checks — 2026-09-28

Reproducible artifacts: [B_weak_form_check.py](../scripts/B_weak_form_check.py) and [B_weak_form_check.json](../runs/B_weak_form_check.json). The script ran successfully with SymPy and 70-decimal-digit mpmath arithmetic.

- Exact symbolic residuals are zero for the dual first derivative, dual Hessian, and Bose forward/reverse ratio.
- The finite interpolation example has occupation-interpolation drift 5/7, energy residual exactly (0,0), and entropy-generator eigenvalues 0 and 85/24.
- Physical mismatches Delta=-0.1 and +0.1 both give positive thermal energy drift, respectively 0.0098043593 and 0.0087339718 in the chosen unit. The positive-square identity agrees to below 9e-73.
- For the coarea benchmark u=0.5,d=0.2, the analytic mass is 27.42206883389030149289273. Root residual is below 1e-71; mass and cos-squared quadrature discrepancies are below 3e-70 for N=3,4,8,16,32 points per root.
- For the non-polynomial integrand exp(cos(q_a)), relative quadrature errors at N=4,8,16,32 are 4.32e-3, 1.57e-7, 1.17e-18 and 1.41e-45, respectively. The reference is the constant Fourier coefficient's convergent positive power series.

These numerical errors describe this regular analytic benchmark and the selected arithmetic precision. They do not validate a singular limit, a time integrator, global positivity of the finite ODE, or a material collision operator. The support obstruction and finite structural identities rest on the displayed derivations; numerical agreement alone is not their proof. No peer report was consulted for this first pass.
