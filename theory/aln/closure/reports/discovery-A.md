# Discovery A: dressed variables and a passive infrared continuum

Independent second-generation proposal, 2026-09-17. Bounded DISCOVERY stage; no substantial implementation is authorized by this proposal. Read only the campaign question, assumptions, checkpoint 1, and branch A projection report, together with the research protocol. No discovery-B or red-team report was read. No novelty claim is made. Both proposals require an independent counterexample attack before implementation.

## 0. Scope and preliminary decision

The physical system is a small entropy-coordinate perturbation of a bulk periodic crystal at 300 K. Assume detailed balance, C=C* >= 0, bounded Hermitian velocity operators V_i, and exact energy conservation Ce=0. Write L(k)=C+i k_i V_i. Streaming conserves the entropy norm, while collisions decrease it. This reversible/irreversible split is sufficient for a dissipative Galerkin formulation; a canonical Hamiltonian or Poisson structure is not supplied by these assumptions.

Two directions are retained:

1. Dress the slow coordinates by collision minimization, carrying the induced entropy metric instead of assuming unit inertial weight.
2. Represent a gapless current-coupled infrared sector by a positive continuum of auxiliary relaxation fields, retaining its nonanalytic frequency and potentially spatial response.

Neither identifies real-AlN coefficients from diagonal rates. Both change what must be measured or computed. Direction 1 asks for a few collision cell solves and residuals. Direction 2 asks for a current-weighted low-rate spectral measure and its interaction with streaming.

## 1. Direction A1: collision-minimizing coordinates and their entropy metric

### Precise finite-dimensional construction

Let E have orthonormal columns containing all relevant exact invariants; Q=I-EE*. Define A=E*CE, K=QCE, D=QCQ restricted to Q, and assume D>0 for this construction. The lift

    T = E - D^{-1} K

is understood as a map into the full space, with its second term in Q. For prescribed raw coordinates E*y=a, the vector Ta minimizes y*Cy over the eliminated coordinates. Direct multiplication gives

    CT = E A_eff,                 A_eff = A-K*D^{-1}K >= 0,
    G = T*T = I+K*D^{-2}K >= I,  J_i = T*V_i T = J_i*.

Here T and G are dimensionless; A_eff has units 1/s; J_i has units m/s. The pulled-back equation

    G partial_t a + J_i partial_i a + A_eff a = T*s

has entropy (1/2) integral a*G a, with decay rate integral a*A_eff a on a periodic domain. Set b=G^{1/2}a for an orthonormal entropy variable. If an invariant is a column of E, its K column vanishes and the lift preserves it. Conservation must also be checked with the source and chosen boundary conditions.

The exact Schur operator from the projection report has, on a finite grid, the derivatives

    S(0,0)=A_eff,
    partial_z S(0,0)=G,
    partial_{k_i} S(0,0)=i J_i.

The last identity includes the term K*D^{-1}(QV_iQ)D^{-1}K, which is absent if one merely sets the eliminated resolvent to D^{-1} before expanding in k. It is essential when K is not small. These are algebraic derivative identities, not an accuracy theorem on a finite frequency/wave-vector band.

For C=[[1+epsilon,1],[1,1]], E=(1,0), the lift is T=(1,-1), G=2, and A_eff=epsilon. The reduced decay rate epsilon/2 agrees with the leading small eigenvalue; retaining unit G predicts epsilon. This example only explains the correction; it does not certify arbitrary moment choices.

### Hypothesis A1 and intended use

HYPOTHESIS: physically chosen energy/drift coordinates in AlN may admit a finite, mesh-converged dressing D^{-1}K even when the rest of D has no uniform gap. On a specified finite response band, a small projected dynamical residual could then certify a useful entropy-stable model without certifying an operator-wide gap.

The finite-grid derivative identities are derived. Existence, mesh convergence, and useful small residuals for AlN are unresolved. In a continuum, D^{-1}K must exist as a Hilbert-space vector, not merely as a formal or distributional inverse. A finite static Schur form alone does not ensure a finite G.

Let U=T G^{-1/2}, P_U=UU*, L_U=U*L(k)U, and R=(I-P_U)L(k)U. For data initially in range(U), contractivity gives the elementary sufficient estimate

    ||exp[-t L(k)]U - U exp[-t L_U]|| <= t ||R||.

For Re z=sigma>0, the corresponding lifted resolvent error is bounded by ||R||/sigma^2. These conservative bounds do not assume a gap; they lose usefulness for long times or sigma approaching zero. Observable-specific improvement is a research target, not established here. The model also omits order-k^2 transport unless enough flux fields are retained or a passive memory channel is added. Pure energy retention alone gives no diffusion in this first-order system.

### First attacks; stop conditions before implementation

1. Construct D with eigenvalues tending to zero and K satisfying finite K*D^{-1}K but divergent K*D^{-2}K. This would make static elimination finite and the dressed inertia ill-defined. PSD of the full C and inclusion of energy must both be enforced.
2. Add a velocity coupling from range(T) into a collision-slow direction invisible to K. This should break any certificate based only on collision dressing; test the full k-dependent R.
3. Use nonzero eliminated initial data. A valid closure must retain its forcing term or state a prepared-data restriction.
4. Compare first derivatives with the exact Schur complement for a 3- or 4-mode PSD collision matrix and a noncommuting Hermitian velocity. Attack signs and the additional K*D^{-1}W_iD^{-1}K term.
5. Reject A1 as a quantitative closure if G or the relevant residual does not converge under infrared mesh refinement, or if the desired band makes the residual bound vacuous. Do not rescue it by fitting G independently of its cell problem.

## 2. Direction A2: passive continuum memory and fractional corrections

### Minimal spectral hypothesis

For the energy invariant, put E=e/||e|| and b_i=QV_iE. At k=0 define the positive matrix measure dM_ij(lambda)=b_i* dE_D(lambda) b_j, with any extra zero modes projected out. Suppose, only as a testable hypothesis, that

    dM_ij(lambda)/d lambda = A_ij lambda^alpha + o(lambda^alpha),
    0 < alpha < 1,              A=A* >= 0,

near lambda=0, with sufficient domination for the asymptotic integrals. Then the diffusivity kernel H_ij(z)=integral dM_ij(lambda)/(lambda+z), equal to conductivity divided by volumetric heat capacity in this normalization, has units m^2/s, finite H(0), and a divergent first z derivative, and

    H(z) = H(0) - [pi/sin(pi alpha)] A z^alpha + o(|z|^alpha),
    M(t) = integral exp(-lambda t) dM(lambda)
         ~ Gamma(1+alpha) A t^(-1-alpha).

The power uses the causal branch Re z>0. These formulas describe a fractional correction to an otherwise finite DC response. They do not imply fractional leading diffusion, an isolated hydrodynamic pole, or an observed AlN exponent. In a diagonal 3D acoustic illustration, lambda proportional to omega^p gives alpha=3/p-1; p=2 gives alpha=1/2. Applying that relation to exact collective collision modes without verification is invalid.

### Physical realization and spatial obstacle

Retain a positive low-rate integral, rather than using the fractional expansion at all frequencies:

    H_tail(z)=integral_0^Lambda A lambda^alpha/(lambda+z) d lambda.

This is a continuum of damped auxiliary flux fields, with a positive spectral weight. It has nonnegative Hermitian real part for Re z>0. A naive truncation H(0)-c A z^alpha need not remain passive away from its asymptotic range.

To include spatial streaming one would need a stronger, separate infrared hypothesis: the current-coupled low-rate sector must admit a joint asymptotic description by rate lambda and bounded ray velocity v. If that sector is represented by a positive velocity-dependent matrix measure dA_ij(v), the candidate kernel is

    H_tail,ij(z,k)=integral dA_ij(v) integral_0^Lambda
                  lambda^alpha/[lambda+z+i k.v] d lambda.

Its leading nonanalytic correction is proportional to

    -integral dA_ij(v) (z+i k.v)^alpha.

This represents fractional derivatives along phonon rays, not merely a fractional time derivative at fixed position. Reciprocal pairs of v are needed for the appropriate real/parity structure. The hypothesis is not implied by a k=0 spectral measure because D and QV_iQ generally do not commute.

For a diffusion-scale measurement z of order |k|^2, streaming of order |k| can dominate the argument of the nonanalytic correction. Under the joint-spectrum hypothesis this suggests a directional correction of order |k|^(2+alpha) to the energy Schur operator, rather than the |k|^(2+2 alpha) inferred by inserting z~k^2 into a frequency-only kernel. This is a conditional asymptotic prediction for the causal response, not a dispersion-pole claim.

### Hypothesis A2 and first attacks

HYPOTHESIS: the converged current-weighted collision spectrum for AlN may have a gapless regularly varying tail in the finite-DC/infinite-slope range. If so, a positive continuum memory model may converge where a fitted GK memory time cannot. Its spatial extension requires independently testable rate/velocity information.

1. Fit nothing initially. Test whether integral_0^ell dM scales as ell^(1+alpha), whether the semigroup current correlation has the predicted algebraic tail, and whether causal finite-frequency response has the matching power and coefficient under independent infrared refinement.
2. Build two collision/velocity pairs with identical k=0 current spectral measure but different commutators [D,QV_iQ]; attack any asserted spatial kernel identification from frequency-only data.
3. Test marginal exponents alpha=0 and alpha=1, multiple angular exponents, non-power tails, and a small finite spectral cutoff. Logarithms and crossover scales can mimic fractional powers over a short band.
4. Test positivity after any finite approximation. A continuum represented by positive quadrature weights and nonnegative relaxation rates can preserve passivity; an unconstrained power-law fit can fail it.
5. Abandon the ray model if low-rate collisions mix velocity directions strongly enough that the joint asymptotic description fails. Retain only the k=0 passive spectral representation in that case.
6. Reject a continuum exponent inferred from one finite grid. At fixed minimum lambda>0 the response is analytic sufficiently close to z=k=0; the mesh limit and low-frequency limit do not generally commute.

## 3. Candidates considered but not promoted

- Full canonical Hamiltonian/GENERIC construction: the antisymmetric streaming and symmetric collision split is already useful, but a nonlinear Poisson bracket with a verified Jacobi identity, energy, and degeneracy conditions is absent. Naming a Hamiltonian framework does not supply missing collision matrix elements.
- Pure symmetry reduction: reciprocal parity and the uniaxial point group reduce tensor components and can separate thermal and momentum-flux sectors. They do not fix low eigenvalues, spectral weights, or the coupling between retained and discarded spaces. Use symmetry to reduce the proposed tests, not as a closure theorem.
- Topological phonon transport: no Berry-curvature, magnetic/time-reversal-breaking, coherent-polarization, or edge-state input is present. A population-only real dissipative model cannot establish such effects. This direction changes the physical model before the available evidence warrants it.
- Classical stochastic jump-process representation: a symmetric PSD collision matrix need not have nonpositive off-diagonal entries in the given coordinates, so it need not generate a positivity-preserving Markov chain. The positive scalar/matrix spectral measure above is legitimate without inventing a microscopic random-walk interpretation.
- Entropy maximization over collision matrices matching diagonal rates: it selects one representative of an underdetermined operator class. It supplies a prior, not material identification or an error certificate.
- Bare fractional Fourier law: a k=0 nonanalytic conductivity does not establish its spatial continuation; its all-frequency truncation may also violate passivity. Only the positive tail representation and its separately testable spatial hypothesis are retained.
- Standard grey GK fitting or another uniform-gap estimate: neither resolves non-invariant coordinate inertia, a divergent weighted moment, or absent collision information.

## 4. Status and next decision

DERIVED UNDER FINITE-GRID ASSUMPTIONS: the lift, entropy metric, Schur derivatives, dissipation identity, and elementary residual bounds in A1.

DERIVED CONDITIONAL ASYMPTOTICS: the A2 power law from a specified regularly varying current spectral density. The joint rate/velocity extension is a separate model assumption.

CONJECTURED FOR AlN: convergent dressed coordinates and a useful residual certificate; a current-coupled spectral tail in the indicated exponent range; asymptotic joint rate/velocity structure.

UNRESOLVED: actual collision action; physical realizability of attack matrices; continuum convergence; useful measurement band; boundary/interface layers. No significant implementation has been performed.

NEXT GATE: independent counterexample attack on A1's continuum metric and k-dependent residual, and A2's streaming extension. Only surviving narrowly stated claims should receive computational development. Literature search must precede any novelty language; established search families include harmonic extension/static condensation, Mori-Zwanzig, passive model reduction, Stieltjes relaxation spectra, distributed-order memory, fractional material derivatives, and kinetic anomalous transport.
