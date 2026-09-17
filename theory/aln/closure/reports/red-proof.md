# Independent red-team proof audit

Completed: 2026-09-17. Role: adversarial proof audit. Sources read: frozen red-team packet; branches A, B, C and D; pi-independent.md. No other red-team report was read. No agent conclusions were used in reaching these verdicts.

## Executive verdict

**The central finite-dimensional conclusions survive this audit.** I found no counterexample or missing algebraic step in D's eight-mode theorem, B1 under its exact hypotheses, B's projected-surrogate pseudoinverse, A's finite-dimensional Schur estimates, or PI's real-Laplace Jensen inequalities for a nonzero source.

**One subsidiary continuum inference is incomplete as written:** C turns a value asymptotic for the rate into a pointwise asymptotic for its pushforward density without a precise derivative hypothesis. The explicit smooth, strictly monotone counterexample below falsifies that inference under the explicitly quantified assumptions. It does not falsify C's moment thresholds, the elementary square-root-response example, or its grazing-angle calculation.

**Microscopic AlN remains unestablished.** D and B are theorems about their declared matrix classes. Neither constructs two allowed AlN event networks with the actual frequencies, momentum selection rules, crystal representations, matrix elements and measured diagonals. The branches explicitly acknowledge this limitation. Applying their finite-matrix existence theorems to that smaller microscopic feasible set would be an UNJUSTIFIED implication.

## 1. D: eight-mode construction — VALID in the stated class

Exact proposition: for every 0<t<1, C(t) is symmetric PSD with kernel span(e), physical diagonal 1, fixed e and current b, the stated reversal parity, K(0)=1 and tau_mem=1+t^(-2); the three-current extension has fixed DC tensor I_3.

Audit trail:

1. **Walsh transformation.** For H_ij=(-1)^(popcount(i AND j))/sqrt(8), distinct character columns are orthogonal. Thus H^T H=I. Congruence by H preserves symmetry, rank, eigenvalues and quadratic forms.
2. **Positive dissipative restriction.** Apart from the zero coordinate, A has blocks B=[[t^2,t],[t,2]], F=[[1,-t],[-t,1]], and scalars 2-t^2,1,1. B has positive first principal minor and determinant t^2; F has eigenvalues 1-t,1+t. All are strictly positive exactly throughout the stated open interval. No endpoint is being used as a positive-definite member.
3. **Physical diagonal.** Every squared Walsh entry is 1/8, while tr(A)=8. The two cross contributions are t S_i1 S_i2/4 and -t S_i4 S_i7/4. Character multiplication gives S_i1 S_i2=S_i3=S_i4 S_i7. They cancel identically. This is an entrywise argument, not a trace-only inference.
4. **Reversal.** Translation i -> i XOR 3 acts on column j by (-1)^(popcount(3 AND j)). Columns 1 and 2 have the same odd parity; 4 and 7 have the same even parity. Consequently both cross blocks commute with reversal; e is even and b=h2 is odd. Fixed diagonal velocities b_i/e_i are well defined because every e_i is nonzero.
5. **DC and memory.** Direct multiplication gives C(-h1/t+h2)=h2, and that candidate is perpendicular to e. Uniqueness on e-perp identifies it as C^+b. Orthogonality then gives K0=1 and ||C^+b||^2=1+t^(-2). This proof does not depend on a small-eigenvalue numerical inversion.
6. **Frequency response.** For d=p^2+(2+t^2)p+t^2, direct multiplication verifies (C+pI)[-t h1+(p+t^2)h2]/d=b. Taking the current inner product gives the stated K_t(p). For Re p>0, positivity ensures the denominator is nonzero.
7. **Three currents.** C h5=h5 and C h6=h6; these vectors are orthogonal to h1,h2. Therefore J^T C^+J=I_3 and J^T(C^+)^2J=diag(1+t^(-2),1,1). All three are odd under the same reversal. Their three binary characters are independent, so the claimed eight velocity sign combinations also check.
8. **Unboundedness and singular limit.** Every positive t is admissible, and 1+t^(-2) has no finite upper bound on that set. At t=0, rank changes and the physical zero-frequency limit must be reconsidered. K_t(p) -> 1/(p+2) for fixed p>0 while K_t(0)=1 is a valid noncommutation of limits, not a contradiction.

D's four-mode graph checks also survive: the common eigenvectors give the displayed eigenvalues, the equal-DC pair has K0=3/2 and memories 3/4,5/6, and the hidden-gap example has exactly constant scalar response. The four-mode unequal-energy perturbation argument has the stated local-interior and nonzero-derivative qualifications.

**Scope severity: high if omitted in synthesis.** PSD, one energy invariant and reciprocal parity are not a proof of realizability by resonant three-phonon/isotope events. The construction neither proves nor disproves nonidentifiability within the exact actual-AlN microscopic feasible set. No novelty verdict is supplied by this audit.

## 2. B1: local nonidentifiability — VALID under its exact hypotheses

The quantified statement is stronger than a dimension count, but its second-variation proof establishes it.

- The hypothesis ker(C0)=span(e) makes the restriction to e-perp strictly positive in finite dimension. Hence sufficiently small signed perturbations in T preserve positivity and the exact kernel. This is the required existence of an actual relative neighborhood; positive listed rates alone would not supply it.
- For each four-coordinate set, A1 f=0 requires f2=f3 and f1=f4; A2 f=0 requires f2=f4 and f1=f3. Their common kernel consists of constants. Under X=E^(-1) A E^(-1), the corresponding statement is that f_i/e_i is constant. Overlapping four-coordinate sets establish the same result for n>4. Thus the common kernel of T is exactly span(e).
- For nonzero b perpendicular to e, f=[(C0+pI)|_(e-perp)]^(-1)b is nonzero and perpendicular to e. Some X in T has Xf nonzero. Symmetry and Xe=0 ensure Xf also belongs to e-perp.
- The second variation is 2(Xf)^T[(C0+pI)|_(e-perp)]^(-1)(Xf)>0. Therefore the response cannot be constant on any relative neighborhood. The first derivative need not be nonzero; the proof correctly avoids needing it.

An explicit fragile-case check: for e=(1,1,1,1), C0=(4/3)[I-ee^T/4], f=(3,-1,-1,-1) and the displayed A1, the first derivative vanishes but the second derivative is 48. Thus replacing the second-variation argument by a generic first-derivative assertion would fail.

**Essential restrictions:** n>=4; every e_i nonzero; nonzero current; an existing feasible matrix with exactly one invariant; and the unshrunk tangent class. B1 does not automatically extend after imposing extra microscopic support, additional invariants or crystal constraints. For n=3 the off-diagonal preserving tangent space can be zero. These are correctly acknowledged limitations.

The feasibility lemma also checks: write C as a Gram matrix, use C e=0 iff the weighted generating vectors sum to zero, and apply the closed-polygon length criterion. The lemma only proves PSD feasibility, not the rank hypothesis of B1.

## 3. B projected surrogate — VALID; identification remains UNJUSTIFIED

For D>0, full-column-rank H, Q the orthogonal projector onto Ran(H)-perp, and C_D=D-DH(H^T D H)^(-1)H^T D:

1. Minimizing (y-H alpha)^T D(y-H alpha) proves C_D>=0 and kernel exactly Ran(H).
2. C_D H=0 and C_D D^(-1)Q=Q.
3. G=Q D^(-1)Q is symmetric, kills Ran(H), and satisfies C_D G=G C_D=Q. These identities imply the Moore-Penrose equations, so G=C_D^+.
4. For b perpendicular to H, b^T G b=b^T D^(-1)b. This is the surrogate's exact DC response, despite its changed actual diagonal.
5. The diagonal subtraction is strictly positive when the corresponding row of H is nonzero. Inclusion of an everywhere-nonzero energy vector guarantees that condition.
6. The complete odd-current response equals diagonal RTA only with the explicitly stated energy-only projection and reciprocal parity of D,e,b. It does not follow from the pseudoinverse identity alone.

No inverse fit preserving a supplied diagonal, and no fidelity of this surrogate to the true collision operator, follows from these identities.

## 4. A: signs, gap, bounds and slow-coordinate warning — VALID in finite dimension

Assumptions used: real wave vector; C=C*>=0; self-adjoint streaming; an orthogonal retained decomposition; D>=gamma I>0; and Re z=sigma>=0 (strictly positive for the inverse bounds).

- **A1/A2:** Elimination gives S=z+A+iU-Y R_Q X, with Y=K*+iB*, not X*. At K=0, -Y R_Q X=+B*R_Q B. The time convolution has the same positive sign. The fast initial-condition term remains present.
- **Resolvent/semigroup:** Re<w,(D+iW)w> >= gamma||w||^2 gives ||R_Q|| <= 1/(gamma+sigma) and the exponential semigroup bound. This requires no commutation of D and W.
- **A3:** The ordered resolvent identity R_Q-D^(-1)=-R_Q(z+iW)D^(-1) is exact. Submultiplicativity gives the displayed error. No unjustified Neumann expansion is being used.
- **A4/A5:** Expanding -Y D^(-1)X gives -K*D^(-1)K-i[K*D^(-1)B+B*D^(-1)K]+B*D^(-1)B. The collision Schur complement is PSD because the full collision block is PSD and D>0. The first-gradient coefficient is Hermitian. This proves accretivity and the displayed absolute Schur error.
- **A6:** Put w=-R_Q X a. The full block operator maps (a,w) to (Sa,0), so Re<a,Sa> >= sigma(||a||^2+||w||^2). Thus ||S^(-1)||<=1/sigma. The same follows for S0 from its PSD/skew decomposition. The inverse identity then proves the 1/sigma^2 error bound.
- **Interpretation:** This bound degenerates as sigma approaches zero. It gives no uniform relative error near a slow pole or a small observable. The branch explicitly includes this restriction and the fast-source correction.
- **Slow coordinates:** At k=0, differentiating the inverse gives the mass coefficient I+K*D^(-2)K. For the stated pair with epsilon>0, S=epsilon+2z+O(z^2), so the small decay rate is epsilon/2+O(epsilon^2). Static elimination's order-one error survives an arbitrarily small frequency-to-gap ratio. The warning is correct.
- **Normalization:** The displaced-Bose energy and drift equations give c^2=G_n S_P^(-1)G_n*/(k_B T^2 C_V), as written, when the momentum columns are independent and their coupling is nonzero. The viscosity sign and units also agree with the Fourier symbol.

These are finite-dimensional conditional theorems. Uniform continuum coercivity, an actual normal/resistive decomposition and boundary-layer closure are additional problems; the estimates do not supply them.

## 5. C: earliest incomplete inference and surviving statements

### C-density claim — INCOMPLETE as written; explicit counterexample to the value-only implication

Exact disputed inference: from r(q)=a q^alpha[1+o(1)] and radial weight b q^(n-1)[1+o(1)] dq, infer a pointwise pushed-forward density rho(r)~B r^(n/alpha-1) by changing variables as though r=a q^alpha exactly.

Take n=3, alpha=2, a=b=1 and, on 0<q<1/8,

    r(q)=q^2[1+q sin(1/q)].

This is smooth on the punctured interval, positive and strictly increasing. It satisfies r(q)/q^2 -> 1. Indeed

    r'(q)=q[2-cos(1/q)+3q sin(1/q)] > q(1-3q)>0.

The actual rate density is rho(r(q))=q^2/r'(q). The asserted leading density is (1/2)sqrt(r). At q_m=1/(m pi), their ratio is exactly

    rho(r(q_m))/[(1/2)sqrt(r(q_m))] = 2/[2-(-1)^m],

which alternates between 2 and 2/3. There is no limit of 1.

The branch's phrase "smooth enough for the stated uniform asymptotics" is not a quantified derivative hypothesis. If it is intended to exclude this example, the required exclusion must be stated. A sufficient additional condition is radial monotonicity plus r'(q,Omega)=alpha a(Omega)q^(alpha-1)[1+o(1)] uniformly, with the stated angular domination. Alternatively one can prove the needed response asymptotics from the cumulative rate measure without claiming a pointwise density asymptotic. That is an additional argument, not a consequence of the displayed change of variables.

**Severity: moderate, local proof gap.** The example does not contradict the power-moment criterion or necessarily the final long-time/resolvent asymptotics. It breaks their displayed pointwise-density route. The noninteger and integer Stieltjes expansions check if the asserted rate-density asymptotic is separately granted with the necessary remainder control.

### Other essential C claims — VALID under the stated uniform infrared/angular hypotheses

- With positive coupling on a set of nonzero angular measure, rate coefficient bounded above and away from zero, uniform value asymptotics and no other singularities, the radial integral is comparable to integral q^(n-1-m alpha)dq. Thus m alpha<n is precisely the convergence threshold. Differentiation at zero requires the corresponding inverse moment; finite DC alone is insufficient.
- In the exact q^2-rate benchmark, M1=1-1/N and M2=N-1, so tau_mem=N for N>1. The stated square-root response follows by direct integration. These are decisive finite-DC/divergent-memory examples independent of the density gap above.
- For fixed specularity p_s<1, S(K,p_s)~[(1+p_s)/(2(1-p_s))]/K. Integrating over epsilon<<mu<1 produces [3(1+p_s)/(4(1-p_s))]epsilon log(1/epsilon). The grazing region supplies the lower-order matching contribution. The limit is correctly stated as nonuniform near p_s=1.
- For diffuse walls, integrating the exponential term gives the displayed E3/E5 expression. The recurrence j E_(j+1)=exp(-epsilon)-epsilon E_j, together with E1=-gamma-log(epsilon)+O(epsilon), yields (3/4)epsilon[log(1/epsilon)+1-gamma]+O(epsilon^2), including the constant and sign.
- A finite quadrature with an exactly grazing transport weight w0 has the stated w0 floor by termwise limits. The continuum conclusion needs the stated zero-measure grazing set and integrability. It is not a statement about a finite-measure family of two-dimensional modes.
- The listed actual-AlN numerical moment/film values were not independently recomputed in this proof audit. Their continuum convergence cannot follow from these formal examples or from one physical grid.

## 6. PI Jensen bounds — VALID for b nonzero and real s>=0

Write a_j>=0, A=sum a_j>0, t_j>0, K0=sum a_j t_j. The first probability weights are a_j/A. For s>0, f(t)=t/(1+s t) has f''=-2s/(1+s t)^3<0; Jensen gives K(s)<=K0/(1+s K0/A). The second weights a_j t_j/K0 sum to one. For g(t)=1/(1+s t), g''=2s^2/(1+s t)^3>0; Jensen gives K(s)>=K0/(1+s M2/K0). Strictness gives the stated single-current-supported-time equality condition. At s=0 equality is automatic.

The standalone note should explicitly say b!=0: otherwise both defined times are 0/0. This is a minor domain omission, not a failure for a nonzero heat-current direction. Neither complex-frequency ordering nor a finite lower-pole parameter with divergent M2 is established, and the note does not claim either.

## 7. Evidence status and remaining work

- The verdicts above rest on the explicit algebraic arguments in this report, not on agent agreement or numerical agreement.
- Exact symbolic checks are additional reproducibility evidence, not the proof. They subsequently completed successfully; the addendum records their scope and environment.
- No independent microscopic AlN realization, continuum-grid convergence study, literature novelty determination, or full audit of the numerical scripts was performed in this bounded proof review.
- Earliest material proof gap located: C's rate-density change-of-variables inference. Earliest invalid extension to avoid: existence/nonidentifiability in the broad finite matrix class does not imply the same existence statement in the actual AlN microscopic class.

**Strongest justified conclusion:** the supplied exact finite-matrix structural data do not determine the response or memory in the declared broad operator class. The stronger actual-AlN reconstruction question remains UNKNOWN. The checked projection/Jensen theorems are conditional mathematical results; a controlled real-AlN closure is not established by them.


## Reproducibility addendum — completed

Independent script: `experiments/red_proof_checks.py`; recorded results: `experiments/red_proof_checks.json`. No campaign implementation is imported. All assertions passed with Python 3.14.7 and SymPy 1.14.0.

The checks use exact arithmetic to verify D's physical-coordinate conservation, diagonal, parity, response, three-current inverse moments and characteristic polynomial; B1's four-coordinate common kernel and the zero-first/nonzero-second-derivative case; all Moore-Penrose conditions for a rational two-invariant surrogate; a nonzero-K/nonzero-B Schur sign example; the alternating density ratios 2 and 2/3 in the C counterexample; and a strict two-pole Jensen instance. They supplement the general hand arguments rather than establish theorems by sampling.

Run:

    py -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260916-210742-aln-spectral-closure\experiments\red_proof_checks.py

No check essential to the stated audit verdict remains unexecuted. Actual-AlN numerical recomputation and microscopic realizability remain outside this bounded proof-audit scope.
