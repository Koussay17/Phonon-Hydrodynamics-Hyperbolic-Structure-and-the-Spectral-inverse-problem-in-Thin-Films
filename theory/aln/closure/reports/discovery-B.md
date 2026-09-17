# Discovery B: positive spectral measures and a common observable subspace

Independent second-generation discovery, 2026-09-17. Read only the campaign question, checkpoint 1, the research protocol, and first-pass B report. No discovery-A or red-team reports were consulted. An early draft was saved before the calculations. The branch owns only this file.

**Status: DERIVED UNDER STATED ASSUMPTIONS; AWAITING COUNTEREXAMPLE ATTACK AND PROOF AUDIT.** No novelty claim, new AlN numerical result, or extensive implementation is made.

## 1. Selected directions and alternatives

Two bounded directions merit independent attack:

1. **Response envelopes:** optimize over a positive observable spectral measure constrained by a few actual response measurements or collision solves. Report response intervals rather than coefficients whose inverse moments may be unstable.
2. **Common observable subspace:** generate a space jointly under collision and streaming from the requested sources and exact invariants. Global slow modes outside this space can be irrelevant, but modes invisible to uniform current can become relevant under streaming.

| Route considered | Decision |
|---|---|
| Temperature/momentum variable changes | Do not provide missing collision information; useful only after finding an invariant or approximately invariant space. |
| Dual/variational geometry | Selected: positive-measure constraints and linear majorants give rigorous response bounds. Lifting measures back to microscopic collision matrices is not assumed. |
| Spectral/transform methods | Selected: positive Laplace shifts avoid singular inverse moments and preserve the exact Stieltjes structure. |
| Hamiltonian methods | Streaming is skew-adjoint after Fourier transformation; collision is dissipative. Passive-system reduction is useful, whereas a purely Hamiltonian model drops essential physics. |
| Symmetry | Reduces spaces, but streaming connects parity sectors. The example below preserves inversion symmetry. |
| Dimensional/asymptotic limits | Selected for falsification: simultaneous small frequency, wave vector, and hidden rate expose a nonuniform expansion. |
| Conserved/monotone quantities | Energy conservation and entropy decay certify structure, not identifiability or a gap. |
| Probability | A positive spectral measure needs no Markov interpretation. A jump-process interpretation would impose extra off-diagonal sign conditions. |
| Topology | Graph connectivity can count invariants, but gives no quantitative gap. No useful topological response constraint was identified. |

## 2. Precise setting

Initially use a finite-dimensional real entropy-coordinate space. Let C=C^T>=0, with exact invariant columns H satisfying C H=0. Let b be a nonzero current source perpendicular to ker(C), and M=||b||^2. A physical current can be b=V e, with e the energy vector and V a velocity operator. For Re(z)>0 define

    F(z)=b^T(C+zI)^(-1)b.

Accretivity gives existence and uniqueness of this inverse. Write K=F(0)>0 and tau=-F'(0)/K. Finite-dimensional DC uses orthogonality to every zero mode, not just the declared energy invariant. Physical conductivity differs by the known factor 1/(k_B T^2).

All spatial statements concern bulk Fourier modes. No boundary model or realization by AlN force constants is claimed.

For a positive self-adjoint continuum operator, finite scalar DC means b belongs to Dom(C^(-1/2)) on its positive spectral part. A finite-energy strong DC solution requires the stronger b in Dom(C^(-1)), equivalent to finiteness of the numerator defining tau. These conditions must not be conflated.

## 3. Direction 1: spectral-measure response envelopes

### 3.1 Exact representation

Let E_C be the spectral resolution and define mu(S)=b^T E_C(S)b. Then

    F(z)=integral (lambda+z)^(-1) dmu(lambda),
    integral dmu=M,
    K=integral lambda^(-1) dmu,
    K tau=integral lambda^(-2) dmu.

The measure is positive. These are exact collision-operator statements, not a diagonal-rate approximation. For p>0,

    (-1)^j F^(j)(p)=j! integral (lambda+p)^(-j-1) dmu >=0.

Domination by j! M/p^(j+1) justifies differentiation. At p=0 each inverse moment requires its own convergence test. Total modal rates do not determine mu.

### 3.2 Bounds from one positive-shift solve

Assume an independently justified spectral ceiling 0<=C<=L I, L>0. Suppose an actual solve or measurement gives S=F(p0)>0 at a real p0>0. Define

    dnu=dmu/[(lambda+p0)S],
    a=integral lambda dnu=M/S-p0.

This is a probability measure with 0<=a<=L, and

    F(p)/S=integral h_p(lambda) dnu,
    h_p(lambda)=(lambda+p0)/(lambda+p).

For p>p0 the function h_p is concave. Jensen and its endpoint chord give

    D(p)<=F(p)<=U(p),

    U(p)=M S/[M+(p-p0)S],
    D(p)=S[(1-a/L)p0/p+(a/L)(L+p0)/(L+p)].

For 0<p<p0 it is convex, so U(p)<=F(p)<=D(p). At p=p0 all three equal S. The one-point and endpoint measures make the bounds sharp in the relaxed measure class. If current-visible zero eigenvalues are forbidden, endpoint mass at zero is approached by small positive rates rather than necessarily attained.

If an exact finite DC value K is known, set dnu=dmu/(K lambda) and a=M/K. For every p>0,

    M/(L+p)<=F(p)<=M K/(M+pK).

The upper endpoint is a single visible rate M/K. When K>M/L, the lower endpoint can be a limiting envelope: vanishing spectral mass near zero carries finite DC weight.

**Scope:** these inequalities are for real positive Laplace shifts. Substituting p=i omega does not establish harmonic-response ordering. The finite-grid choice L=tr(C)=sum r_mu is valid, but can grow with mesh size; it need not provide a continuum-uniform ceiling. Sharpness is only in the relaxed positive-measure class, not the class with all prescribed diagonals, symmetries, and microscopic constraints.

### 3.3 Dual bounds from several measurements

Positive-shift samples F(p_j) give integrals of continuous rational functions against nu. Its mass and mean are already constrained. For Re(z)>0, a directional component of complex response is linear:

    Re(exp(-i theta)F(z)/S)
      =integral Re[exp(-i theta)(lambda+p0)/(lambda+z)] dnu.

Any linear combination of constraint functions that majorizes the target function on [0,L] produces an upper bound after integration. A minorant gives a lower bound. This weak-duality certificate needs no unproved optimizer or strong-duality assertion. Measurement error intervals can remain inequalities.

Every originally admissible C supplies a feasible measure. The converse is not assumed. Discarding fixed-diagonal and microscopic lifting constraints therefore produces conservative outer envelopes, not identified collision operators. A wide envelope is a useful failure result.

### 3.4 Internal counterexample: accurate response does not certify memory

Fix M,K>0 and a=M/K. For 0<t<1 set

    lambda_low=a t^2,
    lambda_high=a(1+t+t^2),

with spectral masses M t^3 and M(1-t^3). Then

    F_t(z)=M[t^3/(z+a t^2)
              +(1-t^3)/(z+a(1+t+t^2))].

Every member has total mass M and exact DC K, but

    tau_t=(1/a)[1/t+(1-t)/(1+t+t^2)] -> infinity.

Let F_*(z)=M/(z+a), x=z/a, and h=1+t+t^2. Direct subtraction gives

    (F_t-F_*)/K
      =-t(1-t)(1+t)^2 x/
         [(x+t^2)(x+h)(x+1)].

For Re(x)>=0, |x|<=|x+t^2| and |x+h|>=h. Consequently

    |F_t/F_*-1|
      <=t(1-t)(1+t)^2/h<=t.

Thus the entire closed right-half-plane response converges uniformly in relative error while its boundary derivative at zero diverges. The nearest pole approaches zero, so this does not contradict interior analytic derivative convergence.

The measure can be realized with fixed b in a two-dimensional dissipative space by rotating eigenvectors so their squared overlaps with b are M t^3 and M(1-t^3), and appending an orthogonal energy invariant. **This does not preserve mode-basis diagonal rates.** It must not be cited for that stronger conclusion.

The implication to reject is: an accurate single-pole fit determines a physical memory time. Controlling tau requires control of low-rate weight in the inverse-square measure.

### 3.5 Bounded proposal and falsification gates

Before implementing an optimizer:

1. Independently audit the envelopes and two-atom identities.
2. On a small available positive matrix, compare exact response with one-shift bounds at three positive target shifts; do not infer harmonic ordering.
3. Check whether obtainable AlN data supply actual F(p0), rather than an RTA surrogate.
4. Reject practical usefulness if outer envelopes remain too wide in the requested window. Reject a memory interpretation while inverse-square weight is uncontrolled.

Useful next data are a few residual-controlled positive-shift collision solves or response measurements with errors. Dense C storage is unnecessary. The presently supplied diagonal rates do not provide these measurements.

## 4. Direction 2: common collision/streaming observable subspace

### 4.1 Exact finite-dimensional construction

Let V_j=V_j^T be streaming operators. Let B contain requested sources/outputs and every exact invariant the reduced state must retain. Define

    K_B=span{W(C,V_1,...,V_d)b:
             W a finite word, including the empty word;
             b a column of B}.

All generators are self-adjoint, so K_B and its orthogonal complement are invariant under every generator. This is the smallest common invariant subspace containing B. It is not asserted to be the smallest possible parametric input-output realization under every definition of reduction.

If U is an orthonormal basis of K_B, then for real k and Re(z)>0,

    B^*(zI+C+i sum k_j V_j)^(-1)B
      =(U^*B)^*(zI+U^*CU+i sum k_j U^*V_jU)^(-1)(U^*B).

The proof is direct block diagonality. No collision gap outside K_B is needed. For arbitrary initial data, retain its projection onto K_B; the complementary part cannot contribute to the specified outputs.

Orthogonal compression preserves positivity and skew-adjoint streaming. Unforced Fourier dynamics obey

    d(||y||^2/2)/dt=-y^*Cy<=0.

Included exact invariants retain their collision nullspace. These facts certify structure and exactness on the common subspace, not that this space is small or hydrodynamic. Generic operators can generate the full mode space.

### 4.2 Internal counterexample: collision-only visibility fails spatially

In an orthonormal basis (e,h,d), take gamma,epsilon,v,w>0 and

    C=diag(0,gamma,epsilon),
    V=[[0,v,0],[v,0,w],[0,w,0]],
    b=V e=v h.

Energy e is the unique collision invariant. Parity R=diag(1,-1,1) satisfies R C R=C and R V R=-V. Diagonalizing the real symmetric V supplies a physical basis with diagonal velocities; the components of e can all be chosen positive.

At k=0,

    F(z,0)=v^2/(z+gamma),

independent of epsilon. The collision-only cyclic space generated by b is span(h). Including energy gives span(e,h), but V couples h to the hidden even mode d.

Eliminating e and d gives

    F(z,k)=v^2/[z+gamma+k^2 v^2/z+k^2 w^2/(z+epsilon)].

Projection onto span(e,h) omits the last term. Set w=v, z=epsilon, and k v=sqrt(gamma epsilon). Then

    F_full/F_projected
      =(epsilon+2 gamma)/(epsilon+(5/2)gamma) ->4/5

as epsilon/gamma->0, although z/gamma and |k|v/gamma both approach zero. A uniform-current observable gap gamma therefore cannot certify this spatial reduction. The hidden rate remains relevant in this simultaneous hydrodynamic scaling.

This example preserves inversion symmetry and admissible streaming. It is not a fixed-diagonal pair or an AlN microscopic model.

### 4.3 A positive-shift residual certificate

For an approximate real orthonormal U containing B, put A(k)=C+i sum k_j V_j. At Re(z)>=sigma>0 compute

    X_U=U(zI+U^*A(k)U)^(-1)U^*B,
    R=B-(zI+A(k))X_U.

Accretivity implies ||(zI+A(k))^(-1)||<=1/sigma, hence

    ||B^*(zI+A(k))^(-1)B-B^*X_U||
      <=||B|| ||R||/sigma.

This a posteriori certificate needs operator actions and a computed residual, not a global gap. It is useful only when the residual is small enough for the selected sigma. Taking sigma to zero destroys the estimate. Pure harmonic response needs a separate resolvent or singular-value bound or a justified limiting argument.

Check (I-UU^*)CU and (I-UU^*)V_jU separately. Streaming residuals can activate hidden modes even when the collision residual vanishes exactly. Small residual norm alone does not certify a uniform low-frequency approximation.

### 4.4 Bounded proposal and falsification gates

1. Independently attack the three-mode example, parity, diagonalization of V, conservation, and simultaneous scaling.
2. On a small available model, generate only a few word-space layers from the required columns. Check rank and collision/streaming residuals separately.
3. Test a small positive-shift (z,k) window with exact response where affordable and the residual certificate.
4. Abandon efficient reduction if the space fills rapidly or its certificate remains loose at all useful shifts. Do not extrapolate to z=i omega or zero.

The modest hypothesis is that some requested response windows admit a small certifiable common observable space. Existence and dimension for real AlN remain unresolved.

## 5. Prior art and scope of search

These are established mathematical languages, not claimed inventions. A short primary-source cross-check identified:

- Zimmerling, Druskin and Simoncini, *Monotonicity, Bounds and Acceleration of Block Gauss and Gauss-Radau Quadrature for Computing B^T phi(A) B*, Journal of Scientific Computing 103, article 5 (2025), [DOI 10.1007/s10915-025-02799-z](https://doi.org/10.1007/s10915-025-02799-z). It treats positive-shift resolvents, matrix-valued Stieltjes measures, and monotone two-sided block quadrature bounds. Its positive-definiteness and endpoint assumptions require checking before transfer to singular or continuum settings.
- Pozza and Pranic, *The Gauss quadrature for general linear functionals, Lanczos algorithm, and minimal partial realization*, [arXiv:1903.11395](https://arxiv.org/abs/1903.11395), revised 2020. This documents connections among moments, quadrature, Krylov methods, and partial realization.
- Gugercin, Polyuga, Beattie and van der Schaft, *Structure-preserving tangential interpolation for model reduction of port-Hamiltonian systems*, Automatica 48 (2012), 1963-1974, [DOI 10.1016/j.automatica.2012.05.052](https://doi.org/10.1016/j.automatica.2012.05.052), [arXiv:1101.3485](https://arxiv.org/abs/1101.3485). This provides prior art for reduction retaining passivity.

This is a narrow cross-check, not an exhaustive novelty search or verification of every theorem in these papers. Further terms: generalized moment problem, Stieltjes interpolation, positive-real realization, common reducing subspace, noncommutative Krylov space, block Gauss-Radau bounds, observable-dependent reduction. No novelty language is justified or needed.

## 6. Ledger and handoff

| Status | Finding |
|---|---|
| KNOWN | Positive spectral representations; substantial moment/quadrature and passive-reduction prior art. |
| DERIVED UNDER ASSUMPTIONS | One-shift envelopes; fixed-M, fixed-DC two-atom family with uniform relative response convergence and divergent memory; exact common-subspace restriction; three-mode streaming failure; positive-shift residual certificate. |
| OBSERVED NUMERICALLY | None. No material calculation or numerical convergence claim is made. |
| FAILED | Uniform response accuracy as memory certification; collision-only visibility as spatial certification; entropy decay or connectivity as a quantitative gap estimate. |
| UNRESOLVED | Microscopic/fixed-diagonal lifting; useful AlN constraints; small common observable space; harmonic boundary certificates; continuum convergence. |

**Earliest fragile step, direction 1:** identifying a feasible spectral measure with a collision matrix satisfying all original constraints. Only the forward map and conservative outer bounds are used.

**Earliest fragile step, direction 2:** silently replacing collision invariance by joint collision/streaming invariance. The three-mode example defeats that step.

An optional symbolic identity check was attempted, but automatic approval review rejected execution because the usage limit had been reached. No check ran; no symbolic or numerical verification is claimed. The calculations remain hand derivations awaiting independent counterexample attack and proof audit. No extensive implementation was performed.
