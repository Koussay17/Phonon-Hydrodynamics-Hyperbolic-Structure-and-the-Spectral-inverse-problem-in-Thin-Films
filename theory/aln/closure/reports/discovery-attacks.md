# Adversarial review of four discovery directions

Date: 2026-09-17. Reviewer: counterexample hunter. Bounded independent attack before implementation. Read discovery-A.md and discovery-B.md after they were complete, plus my own earlier counterexample review. No other red-team report was read. No mathematical novelty is asserted.

## Decision summary

| Direction | Main attack result | Gate before implementation |
|---|---|---|
| A1: collision-minimizing dressed coordinates | Finite static Schur form can coexist with an infinite entropy metric and no minimizing lift in the Hilbert space. Small absolute residual can coexist with order-one low-frequency error. | Reject unless the lift exists, its metric converges, and the full collision/streaming residual gives a useful error on the stated band and initial-data class. |
| A2: passive infrared continuum | The stated PSD coefficient allows A=0, which invalidates an unconditional divergent-slope assertion. A nonzero scalar tail does not identify streaming; the proposed spatial exponent can vanish/change for degenerate ray support. | Require nonzero directional leading weight, independently test the spatial hypothesis, and distinguish a continuum tail from a finite cutoff. |
| B1: positive-measure response envelopes | Finite DC plus exact mass, a positive-shift sample and a fixed spectral ceiling still permits unbounded memory. Outer-envelope sharpness is not microscopic sharpness. | Reject any memory inference; compute actual envelope width for the requested observable and certify dual inequalities over the entire spectral interval. |
| B2: common collision/streaming space | Distinct physical velocities and a nonzero energy component in every mode force the exact common space to be the full mode space, even for a trivial collision spectrum. | Count velocity distinctions before word expansion; reject exact reduction when the Vandermonde obstruction applies. Approximate reduction still needs a useful band-specific certificate. |

The finite-grid algebraic identities in A1, the positive real-Laplace envelopes in B1, and exact reducing-space identity in B2 survived these attacks as stated. This is not a proof audit or evidence that the proposals are useful for AlN.

## A1. Finite static elimination does not guarantee a dressed state

### Explicit continuum counterexample

Let the Hilbert space be span(e) direct-sum R direct-sum L2(0,1). The first coordinate e is an exact energy invariant. On the other two sectors define

    D f(lambda)=lambda f(lambda),
    K a(lambda)=sqrt(lambda) a,
    A=1+delta,                  delta>0.

Set C equal to zero on e and to the symmetric block [[A,K*],[K,D]] elsewhere. C is bounded and positive: for any scalar a and f in L2,

    <(a,f),C(a,f)>
      =delta |a|^2 + integral_0^1 |sqrt(lambda)f(lambda)+a|^2 d lambda >=0.

Its only null vector is e. Retain E=(e, raw scalar coordinate). Then

    K*D^(-1)K=1,
    A_eff=delta,
    D^(-1)K=lambda^(-1/2),
    ||D^(-1)K||^2=integral_0^1 d lambda/lambda=infinity.

Thus the static infimum is finite but its minimizer is not in L2. The claimed continuum lift T and entropy metric G do not exist as bounded finite-energy objects. On an infrared cutoff eta, G=1+log(1/eta) diverges while the scalar static Schur coefficient is delta+eta -> delta.

This enforces both PSD and the included energy invariant. It directly defeats any implication from a converged static Schur coefficient to converged dressed inertia. A1 already anticipates this failure; the actual AlN hypothesis remains entirely conditional.

### Absolute residual is insufficient at a shrinking band

Take C=diag(0,epsilon), V=[[0,v],[v,0]], retain only energy e1, and use k v=epsilon, z=epsilon>0. Collision dressing is exact: K=0, T=e1 and G=1. The streaming residual has norm epsilon ->0. Nevertheless,

    e1^T(z+C+i k V)^(-1)e1
      =1/[z+k^2 v^2/(z+epsilon)]
      =2/(3 epsilon),

whereas the retained model gives 1/epsilon. Relative error against the exact answer is 1/2, independent of epsilon. The stated residual certificate becomes correspondingly non-small; it is not contradicted.

For nonzero eliminated initial data the raw retained derivative includes -K* w(0), so an unforced reduced state initialized only by E*y(0) can also fail immediately. Prepared data or its forcing term is essential.

### Failed attacks and stop criteria

The finite-dimensional identities CT=E A_eff, G=T*T, and the additional K*D^(-1)W D^(-1)K term in the first streaming derivative have no apparent sign defect. Entropy stability does not fix accuracy. Stop if the continuum lift fails, G depends without convergence on the infrared cutoff, hidden streaming modes produce a non-small band error, or initial-data restrictions do not match the application. Pure energy retention has no first-order diffusion; treating it as a diffusion closure would change the proposal.

## A2. Tail and spatial-continuation attacks

### A genuine degenerate-coefficient counterexample

The hypothesis permits A=A*>=0, including A=0. For the scalar measure

    dM(lambda)=lambda^2 d lambda,       0<lambda<1,

choose alpha=1/2 and A=0. Then dM/dlambda=A lambda^alpha+o(lambda^alpha), exactly as allowed, but

    H(0)=1/2,
    H'(0)=-1,

so the first slope is finite. A nonzero leading coefficient is essential for divergence. With the matrix measure diag(lambda^(1/2),lambda^2)d lambda and leading A=diag(1,0), the slope diverges in the first channel and remains finite in the second. A PSD matrix tail does not make every directional memory divergent.

### Same uniform spectrum, different spatial response

Consider energy e, a current mode h and another even mode d, with

    C=diag(0,a,b),   a,b>0,
    V=[[0,v,0],[v,0,w],[0,w,0]],   v>0.

This respects reciprocal parity diag(1,-1,1). The k=0 current measure is v^2 delta_a, independently of w and b. The eliminated diffusivity kernel is

    H_w(z,k)=v^2(z+b)/[(z+a)(z+b)+k^2 w^2].

For w=0 it is v^2/(z+a); for w!=0 a hidden collision mode changes the spatial response. Identical temporal spectral data do not identify a joint rate/velocity measure. Adding a common continuum tail to both examples preserves this lack of identification. The ray proposal therefore has a genuinely additional assumption, already acknowledged in A2.

### A spatial-exponent degeneracy

The proposed positive ray measure permits support at v=0, or at velocities perpendicular to the Fourier direction n. Then the tail depends on z only in that direction. At z=c |k|^2, c>0, its nonanalytic contribution to the energy Schur operator is of order

    |k|^2 z^alpha = c^alpha |k|^(2+2 alpha),

rather than |k|^(2+alpha). The latter scaling requires a nonzero relevant non-grazing angular coefficient. Reciprocal pairs alone do not guarantee it. This does not defeat the conditional ray model; it defeats an unqualified directional exponent extracted from it.

### Passivity and finite-cutoff traps

For scalar A=1, alpha=1/2, Lambda=1, the exact positive integral has H(0)=2. Its asymptotic truncation 2-pi sqrt(z) is negative on the positive real axis once z>4/pi^2. It cannot be used as a passive model on arbitrary bands.

A small positive cutoff makes the exact response analytic near zero and can mimic a fractional law on an intermediate band. A finite set of exponent fits, including semigroup fits on one finite time range, does not establish a continuum tail. Stop if directional leading weight is zero, a cutoff explains the data equally well within errors, independent refinement changes the coefficient/exponent, or noncommuting streaming fails the joint description. The scalar k=0 positive spectral representation can remain valid while the ray extension fails.

## B1. Response envelopes cannot control memory

### Exact fixed-constraint adversarial family

Start with mu_0=(delta_1+delta_2+delta_3)/3. It has

    M=1,   K=11/18,   F(1)=13/36,   support in [0,3].

Add an atom of mass t^3 at lambda=t^2, t>0, and change the three original masses by d_1(t),d_2(t),d_3(t), chosen to solve

    d_1+d_2+d_3 = -t^3,
    d_1+d_2/2+d_3/3 = -t,
    d_1/2+d_2/3+d_3/4 = -t^3/(1+t^2).

The fixed three-by-three matrix is invertible. Thus d_j=O(t), all original masses remain positive for sufficiently small t, and M, exact DC K, the exact positive-shift value F(1), and ceiling L=3 stay fixed. But

    integral lambda^(-2) dmu_t >= t^3/t^4=1/t ->infinity.

This is a relaxed positive-measure counterexample, not a fixed-diagonal microscopic construction. It shows why adding one accurate collision solve does not turn the envelope procedure into a memory estimator. The same mechanism extends to finitely many independent moment constraints by adjusting sufficiently many interior atoms; no such general extension is needed for the explicit three-constraint example here.

The one-shift envelope formulas themselves pass the curvature check: h_p'' has sign p0-p, and the endpoint-versus-Jensen ordering reverses at p=p0. The internal two-atom counterexample in B is consistent with that failure mode.

### Practical and numerical falsifiers

- If F(p0) is only an RTA value, the feasible measure need not contain the true collision response. The purported bounds then lose their asserted connection to AlN.
- A ceiling L=tr(C) that grows with mesh size can make the lower bound tend to zero. Finite-grid validity is not continuum usefulness.
- A dual majorant checked only at a spectral grid is not a certificate over [0,L]. A missed sign change between nodes can reverse the bound. Rational denominators must stay positive and the entire inequality must be controlled.
- The relaxed extremal measure need not lift to a fixed-diagonal, fixed-event physical operator. The fixed finite event-cone obstruction in my earlier report makes this more than a cosmetic issue for unbounded-memory sequences. Outer bounds remain valid if the original class maps into the relaxed class; their sharpness is restricted.
- Real positive-Laplace ordering gives no componentwise ordering of complex harmonic values.

Stop practical development if measured/solved constraints are unavailable or their certified uncertainty produces an interval too wide for the requested decision. Reject memory or viscosity identification from these envelopes alone. No counterexample was found to weak duality when its pointwise hypotheses hold.

## B2. The exact common space can be forced to be full

### Physical-velocity Vandermonde obstruction

Suppose some real linear combination V_n of the diagonal physical velocity matrices has n distinct diagonal entries v_1,...,v_n, and the retained energy vector has e_i!=0 for every mode. Then

    [e, V_n e, ..., V_n^(n-1)e]

has determinant

    (PROD_i e_i) PROD_(i<j)(v_j-v_i) !=0.

Every column lies in the common word space containing energy. Therefore that space is the entire mode space, irrespective of C. Distinct vector-valued velocities admit such a separating direction except on finitely many direction hyperplanes.

A four-mode example with

    e=(1,1,1,1)/2,
    C=I-ee^T,
    V=diag(-2,-1,1,2)

has a unique energy invariant, reciprocal parity under index reversal, and only one positive collision rate. Yet det[e,Ve,V^2e,V^3e]=9/2, so exact common-space reduction gives no compression at all. The issue occurs even without a complicated collision spectrum.

This does not falsify the exact common-subspace identity, which is correct by block invariance. It is a concrete failure of the hoped-for efficiency. Approximate reduction is a separate band-dependent claim.

### Further failure regimes and gate

The A1 two-mode example above shows a small absolute streaming residual with order-one response error as the target shift approaches zero. B's own three-mode hidden-even-mode example similarly defeats collision-only visibility. Neither contradicts the stated positive-shift residual bound: the factor 1/sigma becomes large exactly where needed.

Adding boundary reflection/transmission operators generally changes the common invariant space; a bulk word space does not automatically reduce a film problem. Arbitrary omitted initial data decouples only for an exact common reducing space. Approximate spaces need an initial-state error estimate.

Before expensive word expansion, count distinct physical velocities and check the included energy support. Stop exact compression when the Vandermonde condition applies. For approximate compression, stop if separate collision and streaming residuals fail to certify the actual (z,k) window, the basis fills the mode space, or harmonic/boundary extrapolation is needed without its own estimate.

## Evidence and limits

All decisive examples above are explicit hand derivations. At this report's initial save no discovery attack experiment is claimed to have run. Small exact checks, if subsequently executed, will be recorded with their actual artifacts and outcomes. No material data, continuum convergence result, or physical AlN counterexample has been established.
