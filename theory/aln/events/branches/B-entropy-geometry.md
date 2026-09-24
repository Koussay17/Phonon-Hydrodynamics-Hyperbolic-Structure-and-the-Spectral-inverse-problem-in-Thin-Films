# Branch B: entropy Hessian and stoichiometric geometry

Independent first pass, completed 2026-09-18 after a usage interruption. Only this campaign's 00-question.md was read. No new peer reports, repository prototypes, or implementations were inspected. This branch owns only this report and has not edited implementation files.

**Status: DERIVED UNDER EXPLICIT ASSUMPTIONS; AWAITING INDEPENDENT PROOF AUDIT AND NONLINEAR VALIDATION.** No absolute material scattering rate, real AlN event list, or novelty is claimed.

## 1. Finite model and event convention

Let i=1,...,m label finitely many modes with energies E_i=hbar omega_i>0 on an equal-weight full grid. Fix T>0 and beta=1/(k_B T). Occupations n_i are positive real mean populations. The prescribed equilibrium is

    N_i=1/(exp(beta E_i)-1),
    D_i=N_i(1+N_i)>0.

An unoriented reversible event is stored once, with a chosen orientation c -> a+b. The parent energy is exactly E_c=E_a+E_b. Its stoichiometric vector is

    nu=e_a+e_b-e_c.

Repeated a=b is allowed and means nu=2e_a-e_c. Strictly positive energies exclude c=a or c=b for an exactly resonant event.

Define its forward and reverse Bose factors by

    P=n_c(1+n_a)(1+n_b),
    Q=(1+n_c)n_a n_b,

using powers when a=b. Choose a positive event coefficient kappa, whose counting and normalization convention must be stated separately. The nonlinear collision model is

    J=kappa(P-Q),
    dn/dt=sum_events nu J.

Here kappa has units of inverse time if occupations and the common grid weight have been normalized accordingly. It is not assigned a material value.

This is an assumed finite kinetic collision law. It is not the exact unitary dynamics of an isolated finite collection of oscillators. Obtaining irreversible rates from a microscopic Hamiltonian additionally requires an appropriate kinetic/weak-coupling framework and consistent phase-space normalization.

## 2. Entropy establishes the linearization

Use dimensionless bosonic entropy

    S(n)=sum_i [(1+n_i)log(1+n_i)-n_i log n_i],

so physical entropy is k_B S times the common state-counting weight. Its gradient and Hessian are

    eta_i(n)=log((1+n_i)/n_i),
    Hess S=-diag(1/[n_i(1+n_i)]).

For each event,

    log(P/Q)=nu^T eta(n).

Consequently

    dS/dt=sum kappa(P-Q)log(P/Q)>=0.

Each summand is nonnegative because log is increasing. This is entropy production within the specified model, without a linearization or a sign assumption on off-diagonal matrix entries.

At Bose equilibrium eta(N)=beta E. Exact resonance gives nu^T E=0, hence P(N)=Q(N). Write their common positive value as

    A=N_c(1+N_a)(1+N_b)
      =(1+N_c)N_a N_b.

A small perturbation x=n-N satisfies

    delta J=kappa A delta log(P/Q)
           =-kappa A nu^T D^(-1)x.

Thus the population generator for one event is

    dx/dt=-kappa A nu nu^T D^(-1)x.

Let y=D^(-1/2)x. The entropy deficit relative to its tangent plane is (1/2)||y||^2 to quadratic order. In these entropy coordinates,

    dy/dt=-C y,
    C=sum_events gamma g g^T,
    gamma=kappa A>0,
    g=D^(-1/2)nu.

Every event is one positive rank-one update. The coefficient multiplying the outer product includes the equilibrium Bose factor A. Replacing gamma by a bare rate without declaring what that rate includes is a normalization change.

Checks:

- C=C^T>=0 and y^T C y=sum gamma(g^T y)^2.
- d(||y||^2/2)/dt=-y^T C y<=0.
- C has units of inverse time.
- The population generator is generally nonsymmetric; its entropy-weighted transform is symmetric.
- Positive off-diagonal entries of C are allowed. Phonon event matrices are not generically graph Laplacians.

## 3. Conservation, nullspace, and rank

Let S_st=[nu_1 ... nu_R] be the stoichiometric matrix, with all listed coefficients strictly positive. Then

    rank(C)=rank(S_st),
    ker(C)={D^(1/2)ell : S_st^T ell=0}.

Indeed C=G diag(gamma)G^T, where G=D^(-1/2)S_st and D is invertible. Therefore x^T C x vanishes exactly when G^T x=0.

A linear population quantity ell^T n is conserved if S_st^T ell=0. Its entropy-coordinate invariant is D^(1/2)ell. In particular,

    e_E=D^(1/2)E,
    C e_E=0.

This also matches the temperature tangent:

    D^(-1/2) dN/dT=e_E/(k_B T^2).

The event model generally has more invariants than energy. Connected-looking mode support does not establish rank m-1. The condition for a unique collision null direction is exactly

    ker(S_st^T)=span(E).

Population number is not conserved by a 1 -> 2 event, since sum_i nu_i=1. Momentum components are linear invariants only if their actual signed sums vanish for every event; a reciprocal-lattice jump is not zero as an ordinary real-valued momentum sum.

Existence of N(T) as a stationary point follows from the stated exact-resonance assumption. Uniqueness across all populations is not claimed. On a fixed stoichiometric compatibility class, a strictly positive detailed-balanced equilibrium is unique if it exists: S is strictly concave, and S_st^T eta=0 is its stationarity condition on that affine class. Small event sets ordinarily define many such classes and extra equilibrium parameters. No global nonlinear convergence claim is needed here.

For the finite exactly resonant model, the polynomial vector field is locally Lipschitz. Nonnegative populations are forward invariant, and positive conserved energy bounds every occupation. This yields a unique global nonnegative solution from nonnegative finite data; positive data remain positive at finite times because on the bounded energy set every loss term is bounded by a constant times the affected population.

## 4. Counting and repeated indices

### 4.1 Reaction reversal is already included

The stored reversible event contains both P and Q. If the reverse orientation is stored with -nu and Q-P and is again evaluated as a full reversible event, its contribution duplicates the first event. At linear order (-g)(-g)^T=g g^T, so this error doubles the generator.

An alternative directed-event representation is valid if it stores only nu kappa P and -nu kappa Q separately. These two conventions must not be mixed.

Likewise a and b are interchangeable in a distinct-mode event. An ordered list containing both (a,b;c) and (b,a;c) doubles an unordered-event implementation unless its coefficient convention compensates.

### 4.2 The unambiguous multiplicity

For a=b,

    nu=2e_a-e_c,
    P=n_c(1+n_a)^2,
    Q=(1+n_c)n_a^2,
    A=N_c(1+N_a)^2.

The two entropy components are

    g_a=2/sqrt(D_a),  g_c=-1/sqrt(D_c).

Thus

    C_aa=4 gamma/D_a,
    C_ac=-2 gamma/sqrt(D_a D_c),
    C_cc=gamma/D_c.

Using nu=e_a-e_c for this event fails energy conservation and gives the wrong rank-one direction. The factor 2 belongs in the population change, even after coefficient symmetries have been resolved.

### 4.3 The coefficient multiplicity requires a declared Hamiltonian convention

Entropy geometry does not by itself determine the absolute prefactor kappa.

For illustration only, suppose the microscopic monomial is g (a_a^dagger)^2 a_c plus its Hermitian conjugate. Fock-state factors are (r_a+1)(r_a+2)r_c and r_a(r_a-1)(r_c+1). Under independent geometric Bose distributions of mean n_a,

    <r_a(r_a-1)>=2n_a^2,
    <(r_a+1)(r_a+2)>=2(1+n_a)^2.

Thus this monomial convention gives a factor 2 multiplying the Bose bracket, relative to the common Golden Rule factor proportional to |g|^2.

If instead the Hamiltonian is written as (1/2)sum_ab g_ab a_a^dagger a_b^dagger a_c with g_ab=g_ba, the collected monomial coefficient is g_ab for distinct a,b and g_aa/2 for a=b. Its repeated-event coefficient is then half the distinct-event base coefficient when those vertex values are equal: 2|g_aa/2|^2=|g_aa|^2/2. These are consistent descriptions using different definitions of g.

This example explains why there is no universal extra factor without the exporter convention. It does not derive an irreversible finite-system rate. It also shows why Fock integers r(r-1) must not be substituted directly for real mean occupations n(n-1): the required mean is 2n^2 under the stated Bose factorization.

## 5. Quadrature coordinates

The question's equal-weight grid avoids relative quadrature factors. A common weight can be absorbed consistently into entropy and event-flux normalization.

For clarity, a conserving extension to positive unequal weights W=diag(w_i) can be defined by applying each event to weighted populations p=W n:

    dp/dt=nu J,  dn/dt=W^(-1)nu J.

The entropy is S_W=sum w_i s(n_i), the energy is E^T W n, and the entropy coordinates are

    y=W^(1/2)D^(-1/2)x.

Linearization then gives

    C=sum kappa A g_W g_W^T,
    g_W=W^(-1/2)D^(-1/2)nu,

with energy invariant

    e_E,W=W^(1/2)D^(1/2)E.

This is a specified common-event-flux discretization, not a claim that arbitrary quadrature exporters already use it. Their event weights must be matched to this convention.

Naively keeping dn/dt=nu J while changing only entropy/energy quadrature weights is generally wrong. For example, E=(1,1,2), nu=(1,1,-1), and w=(1,2,1) give E^T nu=0 but E^T W nu=1. The naive update changes the weighted physical energy.

## 6. Momentum reversal/parity pairs

Assume the full-grid partner map i -> bar(i) is represented by an orthogonal involutive permutation P, with paired energies and weights equal. Then P commutes with D. This is an explicit finite-grid symmetry assumption; the actual mode mapping and event coefficients must be supplied by a material exporter.

A momentum-reversed companion has bar(nu)=P nu. If it is a distinct event and has the same coefficient, its Bose factor is also A. The pair contributes

    C_pair=gamma[g g^T+(P g)(P g)^T].

Define g_+=(g+Pg)/2 and g_-=(g-Pg)/2. They are orthogonal, and

    C_pair=2gamma[g_+g_+^T+g_-g_-^T].

Thus the pair commutes with P. Its nonzero eigenvalues are 2gamma||g_+||^2 and 2gamma||g_-||^2, for the components that are present. Cross-parity terms cancel.

A self-partner event, including one unchanged merely by exchanging a and b, is one event, not two. It must not be doubled by a symmetry expansion. For a genuine distinct companion in this positively oriented 1 -> 2 convention, g and Pg cannot be negatives: that would imply P nu=-nu, impossible because permutation preserves sum_i nu_i=1. Therefore a genuinely distinct stoichiometric companion gives pair rank two.

For two disjoint triplets, g is orthogonal to Pg. The even and odd nonzero eigenvalues are equal to gamma||g||^2. With six modes the pair has rank two and nullity four. Those extra zero modes are exact invariants of the deliberately small event model, not evidence that its eigensolver failed.

Reciprocity of reaction forward/backward rates and equality of momentum-reversed companion rates are separate statements. Neither should be inferred merely from the word reciprocal.

## 7. A concrete audit target with exact equilibrium populations

Take distinct a,b with energies E_a=E_b=E0 and E_c=2E0, and beta E0=log 2. Then

    N=(1,1,1/3),
    D=(2,2,4/9),
    A=4/3,
    g=(1/sqrt(2),1/sqrt(2),-3/2).

For one unordered reversible event,

    C=kappa [[2/3,2/3,-sqrt(2)],
             [2/3,2/3,-sqrt(2)],
             [-sqrt(2),-sqrt(2),3]].

Its spectrum is {0,0,13kappa/3}. The energy entropy-vector is proportional to (sqrt(2),sqrt(2),4/3) and is annihilated exactly.

An explicitly counted disjoint momentum-reversed triplet with equal coefficient gives two copies of this block: two eigenvalues 13kappa/3 and four zeros. In a parity basis the two nonzero directions are one even and one odd.

This target fixes dimensionless algebra and relative normalization only. The common kappa is a declared test parameter, not a measured or calculated AlN rate.

## 8. Nonresonant broadening breaks the equilibrium premise

Let Delta=nu^T E be nonzero while the same nonlinear Bose bracket is retained. At the prescribed Bose distribution,

    log(P(N)/Q(N))=beta Delta,
    J(N)=kappa Q(N)[exp(beta Delta)-1].

Therefore N is not stationary for that event, and its energy production at N is

    d(E^T n)/dt=Delta J(N)
      =kappa Q(N) Delta[exp(beta Delta)-1]>0.

This sign is positive for either sign of Delta. Consequently a positive weighted collection of such detuned reversible brackets cannot cancel this thermal-state energy production by simply including opposite detunings. A momentum-reversed companion has the same detuning anyway.

A broadening weight multiplying kappa does not repair detailed balance or energy conservation. At nonzero equilibrium affinity the Jacobian also acquires the derivative of the state-dependent mobility, so the equilibrium rank-one Hessian formula cannot simply be asserted for the actual nonlinear Jacobian.

One may construct a PSD outer-product model for detuned nu, but it has g^T e_E=Delta and does not conserve the prescribed energy. A corrected conserving discretization, changed equilibrium, or an explicitly modeled external energy exchange is a different model and needs its own derivation.

**Earliest fragile step:** assuming a broadened near-resonant event is stationary at the original Bose equilibrium merely because forward and reverse vertex coefficients agree.

## 9. Bounded validation specification and limitations

No implementation was edited and no computational checks were run in this branch. Suitable independent tests are:

1. Check the complete nonlinear vector field at N before its Jacobian. Stationarity must precede an equilibrium linearization claim.
2. Compare the full-population Jacobian action with -D^(1/2) C D^(-1/2)x, including every lower- and upper-mode update.
3. For these three-boson brackets the cubic terms cancel, leaving a quadratic polynomial. A central directional finite difference therefore has no truncation error in exact arithmetic; step-size sweeps primarily expose roundoff, underflow, and counting discrepancies. Perturbations must keep all n_i positive.
4. Verify energy conservation, entropy production, symmetry in entropy coordinates, and rank from the stoichiometric matrix separately.
5. Check unordered distinct indices, repeated indices, a self-partner event, a disjoint parity pair, and deliberate duplicate insertion. Expected rank and normalization must be declared in advance.
6. Use the exact-population target in Section 7 and a detuned event in Section 8 as positive and negative controls.

Before importing real AlN events one still needs: full-grid mode identifiers and weights; positive mode energies; an explicit partner map; consistent eigenvector/degeneracy conventions; counted event tuples and reciprocal-lattice vectors; exact or controlled energy-conservation treatment; matrix-element definitions and units; all combinatorial factors; temperature-dependent Bose weights; and an exporter contract stating whether rates are directed, reversible, ordered, or already symmetry-expanded.

A finite rank-two validation certifies the chosen algebra and counting. It cannot certify a real-material collision spectrum, completeness of scattering channels, a normal/resistive split, absolute rates, or a hydrodynamic window.

## 10. Status ledger

- DERIVED UNDER ASSUMPTIONS: nonlinear entropy production; equilibrium entropy-Hessian generator; weighted-population normalization; rank/nullspace characterization; parity decomposition; exact small audit target; nonresonant thermal-state energy production.
- KNOWN STRUCTURAL INPUTS: Bose equilibrium and bosonic entropy; the assumed reversible kinetic factorization; explicitly defined event coefficients.
- NOT VERIFIED COMPUTATIONALLY: every displayed identity remains subject to independent algebraic and finite-difference checks.
- FAILED INFERENCES: reaction reversal counted as an additional reversible event; missing repeated-index stoichiometry; arbitrary quadrature rescaling; graph connectivity used as a nullspace test; positive broadening used as a substitute for detailed balance.
- UNRESOLVED FOR MATERIAL IMPORT: exporter normalization, event completeness, momentum-reversal mapping, energy quadrature, and physical scattering prefactors.

An early save was attempted before the interruption but automatic approval review rejected execution because the usage limit had been reached. This report is now being saved under the resumed assignment.
