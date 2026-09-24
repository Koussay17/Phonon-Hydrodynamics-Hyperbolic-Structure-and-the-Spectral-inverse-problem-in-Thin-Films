# A: Direct nonlinear Bose-event linearization

Independent bounded first pass, completed 2026-09-18 after the interrupted save. Read only this campaign's `00-question.md`; no prototype, peer report, or new branch conclusion was consulted. Method: differentiate the nonlinear Bose gain/loss expression before constructing the linear operator. No implementation files were changed. No absolute material scattering rate is claimed.

## 1. Domain and event convention

Take finitely many scalar populations n_i>0 on an equal-weight full phonon grid at 0<T<infinity. All retained energies epsilon_i are positive. Reciprocal modes satisfy epsilon_bar-i=epsilon_i, with a consistent involutive mode permutation J. Coherences and ambiguous degenerate-mode pairing are outside this scalar model.

First take three distinct mode indices a,b,c and an exactly resonant reaction

    a <-> b+c,
    epsilon_a=epsilon_b+epsilon_c,
    q_a-q_b-q_c=G.

The reciprocal vector G refers to the chosen wave-vector representatives. Normal means G=0; Umklapp means G!=0. Count the forward/reverse reaction pair once, with daughters unordered. Define a nonnegative event coefficient kappa by the following nonlinear equation; it absorbs any microscopic amplitude, integration, symmetry-counting, and grid factors not specified here.

With r having entries (1,-1,-1) on (a,b,c), define

    A(n)=n_a(1+n_b)(1+n_c),
    B(n)=(1+n_a)n_b n_c,
    R(n)=kappa[A(n)-B(n)],
    dot n=-r R(n).                                      (1)

Forward decay decreases the parent and increases each daughter. The reverse reaction is already B(n); it is not an additional event. Occupations are dimensionless and kappa has units of inverse time under this definition.

## 2. Detailed balance and direct derivative

Set beta=1/(k_B T), n_i^0=(exp(beta epsilon_i)-1)^(-1), w_i=n_i^0(1+n_i^0), W=diag(w_i), and s_i=sqrt(w_i). Since 1+n_i^0=exp(beta epsilon_i)n_i^0,

    A(n0)/B(n0)=exp[beta(epsilon_b+epsilon_c-epsilon_a)]=1.

Write F0=A(n0)=B(n0)>0. Directly differentiating the two products gives

    delta(A-B)=F0 delta log(A/B)
              =F0[delta n_a/w_a-delta n_b/w_b-delta n_c/w_c]
              =F0 r^T W^(-1) delta n.

Here 1/n-1/(1+n)=1/[n(1+n)]. There is one F0, not 2F0: both gain and loss have already been differentiated. Therefore, with Lambda=kappa F0,

    dot(delta n)=-L_pop delta n,
    L_pop=Lambda r r^T W^(-1).                          (2)

The population generator is generally not Euclidean symmetric. It is symmetric in the entropy metric W^(-1), because W^(-1)L_pop=Lambda(W^(-1)r)(W^(-1)r)^T.

Define entropy coordinates y=W^(-1/2)delta n. Then

    dot y=-C_event y,
    u=W^(-1/2)r=(1/s_a,-1/s_b,-1/s_c),
    C_event=Lambda u u^T >=0,
    d(y^T y/2)/dt=-Lambda(u^T y)^2.                      (3)

The parent-daughter matrix entries are negative; the daughter-daughter entry is positive. The diagonal is positive for an active event. A common grid quadrature weight multiplies the entropy; rescaling all y by its square root leaves C_event unchanged. Unequal grid weights require a separate consistent event-update convention.

The Bose identity s_i=1/[2 sinh(beta epsilon_i/2)] implies that f_i=sinh(beta epsilon_i/2)delta n_i=y_i/2. This constant coordinate factor does not change the generator and is distinct from event-counting or lifetime/linewidth factors.

For a thermal-factor/sign cross-check, put S_i=sinh(beta epsilon_i/2). Resonance gives F0=1/(8 S_a S_b S_c), hence

    C_ab=-kappa/(2 S_c),  C_ac=-kappa/(2 S_b),
    C_bc=+kappa/(2 S_a),  C_aa=kappa S_a/(2 S_b S_c).

These use the kappa defined in (1), without claiming a first-principles prefactor.

## 3. Conserved quantities and normal/Umklapp distinction

An additive quantity h obeys d(h^T n)/dt=-(h^T r)R. If h^T r=0, then h is a left null vector of L_pop, W h is a right null vector of L_pop, and W^(1/2)h is a null vector of C_event. In general,

    C_event W^(1/2)h=Lambda u(r^T h).                   (4)

For energy, define E_i=epsilon_i s_i. Exact resonance gives C_event E=0. The uniform-temperature tangent is

    delta n_i=epsilon_i w_i delta T/(k_B T^2),
    y_i=E_i delta T/(k_B T^2).

Thus the energy perturbation on a unit-weight grid is E^T y; a common quadrature weight multiplies this expression.

For crystal momentum define P_(alpha),i=hbar q_(i,alpha)s_i. Then

    C_event P_alpha=Lambda u hbar G_alpha,
    P_alpha^T C_event P_beta=Lambda hbar^2 G_alpha G_beta. (5)

A normal event annihilates all momentum vectors. An Umklapp event relaxes additive phonon momentum along G and still preserves any component perpendicular to G. A network relaxes momentum only in the span of its reciprocal vectors. Momentum conservation modulo the reciprocal lattice is not conservation of the additive phonon momentum moment; the lattice receives the compensating momentum. The displaced-Bose tangent is y=sum_alpha P_alpha delta v_alpha/(k_B T).

Phonon number is not conserved: r^T 1=-1, so a forward decay increases total population count by one.

## 4. Reciprocal partner versus reverse reaction

The reciprocal event is bar-a <-> bar-b+bar-c, with bar-r=Jr and reciprocal vector -G. The reverse reaction of the original event is b+c -> a at unchanged mode labels, already included in (1). Reversing mode wave vectors and reversing a reaction are different operations. Microscopic time reversal combines both on an oriented transition; the event here includes both orientations from the start.

Assume reciprocal interaction coefficients kappa_bar=kappa. Bose equilibrium alone does not imply this equality for independently assigned coefficients. Since J commutes with W,

    bar-u=Ju,
    C_pair=Lambda[u u^T+(Ju)(Ju)^T],
    J C_pair J=C_pair.                                 (6)

For six distinct modes ordered (a,b,c,bar-a,bar-b,bar-c), this is diag(Lambda v v^T,Lambda v v^T), v=(1/s_a,-1/s_b,-1/s_c). Both orthonormal parity blocks equal Lambda v v^T. This is the physical action on both even and odd populations. A single event pair is highly rank deficient; it does not certify a full collision spectrum or material relaxation gap.

The distinct reciprocal pair satisfies

    P_alpha^T C_pair P_beta=2 Lambda hbar^2 G_alpha G_beta.

The opposite Umklapp vectors do not cancel dissipation. Nonlinear momentum changes as -hbar G(R-R_bar). Equal reaction imbalances cancel this change, but an odd linear perturbation has R_bar=-R.

If the reciprocal tuple is the same unordered event, count it once. Self-reciprocal modes or paired daughters can cause this coincidence. For overlapping events use the unique-event sum C=sum_e Lambda_e u_e u_e^T; do not assume the six-distinct-mode block structure. No crystal inversion symmetry is required beyond the stated reciprocal mode/coefficient pairing.

## 5. Identical daughters: multiplicity is not an arbitrary rate factor

Identical daughters mean the complete mode indices agree, including branch, not just their wave vectors. Positive-frequency resonance excludes a parent equal to either daughter.

For a <-> 2b, define an effective coefficient kappa_abb through

    R=kappa_abb[n_a(1+n_b)^2-(1+n_a)n_b^2],
    r=(1,-2),  dot n_a=-R,  dot n_b=2R.

At epsilon_a=2epsilon_b, the same direct derivative gives

    Lambda_abb=kappa_abb n_a^0(1+n_b^0)^2,
    u=(1/s_a,-2/s_b),
    C_event=Lambda_abb u u^T,
    C_bb=4 Lambda_abb/w_b,
    C_ab=-2 Lambda_abb/(s_a s_b).                       (7)

Accumulate repeated indices in r before forming the outer product. Omitting the cross term between the two daughter occurrences misses the factor four in C_bb.

The effective kappa_abb is not asserted equal to the distinct-daughter kappa. Exact occupation-number ladder factors are (N_b+1)(N_b+2) for pair creation and N_b(N_b-1) for annihilation. Under the quasifree/geometric single-mode population closure,

    <(N_b+1)(N_b+2)>=2(1+n_b)^2,
    <N_b(N_b-1)>=2n_b^2.

The common factorial factor can be absorbed into kappa_abb. Hamiltonian normalization and ordered/unordered sums may introduce compensating factors; only the microscopic convention determines them. Without this population closure, repeated-mode expectations involve higher moments rather than n_b alone.

For distinct daughters, adding both (b,c) and (c,b) with the coefficient of one unordered event doubles the rate. Use unique unordered events or a compensating counting convention. Do not insert a universal extra 1/2 into (7) without identifying the convention it converts. Forward/reverse pairing, daughter permutation, and wave-vector reciprocity are three separate counting operations.

## 6. Independent entropy and finite-difference checks

For Bose entropy S_B/k_B=sum_i[(1+n_i)log(1+n_i)-n_i log n_i], direct substitution yields

    dot S_B/k_B=kappa(A-B)log(A/B)>=0.                  (8)

At exact resonance the relative entropy

    Phi(n|n0)=sum_i[n_i log(n_i/n_i^0)
                  -(1+n_i)log((1+n_i)/(1+n_i^0))]

is nonincreasing. Its equilibrium Hessian is W^(-1), independently identifying the quadratic norm in (3).

A useful finite-difference detail: the cubic products cancel, leaving

    R=kappa[n_a(1+n_b+n_c)-n_b n_c]

for distinct daughters and R=kappa_abb[n_a(1+2n_b)-n_b^2] for identical daughters. These are quadratic. Consequently a central directional finite difference of the nonlinear drift about equilibrium exactly equals its derivative in exact arithmetic, for any step keeping populations positive. There is no O(h^2) truncation term to observe for this test; tiny steps instead amplify subtraction error. Compare its result with -Lambda r r^T W^(-1)delta n. A forward difference has linear-in-step truncation error.

Suggested acceptance tests: detailed balance at n0; nonlinear central differences in generic/even/odd directions; W^(-1/2)L_pop W^(1/2)=C; exact energy residual; normal momentum residual and the nonzero Umklapp identity (5); symmetry and nonnegative quadratic forms; repeated-index factor four; reciprocal commutation; and deliberate duplicate-event tests that reveal factor errors. High precision should separate roundoff from formula/counting mistakes. No numerical implementation was run in this branch.

## 7. Status and import limits

DERIVED UNDER STATED ASSUMPTIONS: equations (1)-(8), signs, coordinate transformation, detailed balance, nonlinear entropy production, exact energy conservation, N/U momentum structure, and reciprocal parity action. Independent campaign audit remains required before final acceptance.

ESSENTIAL ASSUMPTIONS: exact resonance, positive frequencies and temperature, a scalar quasifree Bose population closure, common grid weight, nonnegative specified event coefficient, unique-event counting, and reciprocal coefficients whenever reciprocal symmetry is claimed.

UNRESOLVED FOR REAL AlN: microscopic amplitude convention, absolute inverse-time units, quadrature/delta weights, symmetry multiplicities, complete mode and reciprocal maps, broadened integration, and continuum convergence. A broadened nonresonant event generally has A(n0)!=B(n0) and r.epsilon!=0, so the equilibrium and energy-null derivation fails without further justification. No ad hoc repair is proposed here.

A finite-event test may choose an arbitrary positive kappa as an explicitly defined test scale. Such validation establishes the algebra and event accounting; it does not establish a material scattering rate.