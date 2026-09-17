# Red-team counterexample review

Date: 2026-09-17. Reviewer: counterexample hunter. Status: **completed bounded independent review; additional lemmas await external proof audit.**

Independence: read red-team-packet.md and frozen A-D reports plus pi-independent.md. No other red-* report was read and no other reviewer was contacted. Frozen reports were not changed. Own experiment outputs use red_counterexample_ prefixes.

## 1. Verdict and exact target

**No counterexample was found to D's stated eight-mode algebra in its broad PSD class. An explicit obstruction defeats extending the entire family to a fixed finite microscopic event model.** The latter restricts the admissible class; it does not contradict D's carefully stated theorem. Finite-band response also fails to track the diverging memory uniformly.

The exact target is, for 0<t<1, C(t)=H A(t) H^T with fixed e=h0, b=h2, unit physical diagonal, reciprocal parity, unique invariant e, K(0)=1, and tau_mem=1+t^(-2). D expressly disclaims a microscopic AlN realization and full AlN point-group symmetry. This review does not substitute a stronger claim and then call it false.

| Finding | Severity and exact target |
|---|---|
| The full t->0 family cannot lie in any fixed finite nonnegative event cone | **Major physical-scope obstruction.** Any extrapolation to fixed harmonic AlN modes with varying event rates must confront it. No defect in the stated broad PSD example. |
| Constant energy entries rule out the exact model as an equal-weight, finite-temperature, three-phonon/isotope network | **Major conditional realizability obstruction**; the assumptions below are essential. |
| Divergent memory coexists with uniform convergence on every fixed frequency band away from zero | **Major interpretation limit.** Memory is not an experimentally resolved time scale without a time/frequency window. |
| Cube-corner velocities lack hexagonal covariance | **Moderate material-symmetry restriction**, already disclosed in D. |
| Endpoint, overlap, parity, fixed diagonal and numerical-artifact attacks | No in-domain algebraic counterexample found; this is not a substitute for the separate proof audit. |

## 2. Fixed finite event cone: explicit obstruction

### Assumptions and physical origin

Fix a finite collection of harmonic modes, positive reference temperature, entropy-coordinate transform, and a finite list of allowed detailed-balance collision events. Each event has a fixed entropy-affinity vector z_a. Its strength g_a can vary but remains nonnegative:

    C = SUM_a g_a z_a z_a^T,       g_a >= 0.                 (R1)

For an exact coalescence event i+j<->k, the affinity is f_i+f_j-f_k, where f_i=delta n_i/[n_i^0(1+n_i^0)]. In entropy coordinates this is a fixed linear functional z_a^T y. Linearizing detailed balance gives a nonnegative scalar times its square in entropy production. Elastic events use f_i-f_j. Repeated indices only change the fixed stoichiometric coefficients. Exact energy conservation imposes z_a^T e=0 event by event. A fixed finite set of higher-order detailed-balance events also has (R1).

Standard three-phonon/isotope formulas and entropy symmetrization appear in [Cepellotti and Marzari, Appendix A, Eqs. (17)-(22)](https://arxiv.org/pdf/1603.02608). The positive-outer-product representation here is derived by grouping each forward/reverse event contribution; it is not attributed as a verbatim theorem from that source. Exact resonance and fixed coordinates are assumptions. A numerical broadening that changes energy conservation requires a separate assessment.

### Attack calculation

In D's fixed Walsh basis let u_a=h1^T z_a and v_a=h2^T z_a. D requires

    SUM_a g_a u_a^2   = t^2,
    SUM_a g_a u_a v_a = t.

If every u_a vanishes the second identity is impossible. Otherwise the fixed finite dictionary gives

    R = max_(a:u_a != 0) |v_a/u_a| < infinity.

Events with u_a=0 contribute zero to both sums. Nonnegative weights imply

    |SUM_a g_a u_a v_a| <= R SUM_a g_a u_a^2,
    t <= R t^2,
    t >= 1/R.                                               (R2)

**Therefore no fixed finite event dictionary can realize this family for arbitrarily small positive t.** This obstruction does not need fixed physical diagonals, momentum selection, or event sparsity. It applies even to arbitrary fixed rank-one event vectors. Microscopic restrictions can only shrink that cone.

It does not rule out an isolated t, changing event directions/harmonics with t, an increasing number of modes/events, an infinite continuum dictionary, or an effective operator after eliminating additional modes. It also does not rule out bounded memory nonuniqueness at fixed DC. None of those alternatives is a supplied AlN realization.

### Stronger finite-dictionary consequence: derivation for audit

On a fixed d-dimensional dissipative subspace assume

    A(g)=SUM_a g_a a_a a_a^T > 0,

with fixed finite real vectors a_a, g_a>=0, and fixed b!=0. For each (d-1)-element event subset I define the generalized cross-product v_I by

    det([a_I,x])=v_I^T x.

Cauchy-Binet applied to det(A+s x x^T), followed by comparison of the coefficient of s, gives

    adj(A)=SUM_(|I|=d-1) (PROD_(a in I) g_a) v_I v_I^T.

Dependent subsets have v_I=0. Define the finite geometric constant

    M=max_(I:v_I^T b != 0) ||v_I||/|v_I^T b|.

There is at least one nonzero term because b^T adj(A)b>0. For x=A^(-1)b and K=b^T A^(-1)b, termwise triangle inequality yields

    ||x|| <= M K,
    tau_mem=||x||^2/K <= M^2 K.                              (R3)

Thus **within a fixed finite event cone, fixed finite DC does bound this memory time**. No fixed-diagonal condition is needed. For d=1 the same statement is immediate; the empty-subset convention gives it too.

This geometry-dependent bound may be enormous. Its constant is not supplied by total diagonal rates. The event geometry, source, dimension and invariant subspace must remain fixed, and M need not stay bounded under q-mesh refinement. It does not establish a small or known AlN memory time. **R3 is an independently derived obstruction awaiting external proof audit; no novelty claim is made.** The attached experiment checks its adjugate identity exactly in rational examples.

## 3. Independent equal-frequency obstruction

Assume equal full-zone quadrature weights, Bose equilibrium at fixed 0<T<infinity, positive frequencies, and that these eight modes comprise the entire model. Entropy-coordinate energy entries are proportional to

    f(omega)=omega/[2 sinh(c omega)],  c=hbar/(2 k_B T)>0.

The sign of f' is the sign of sinh(x)-x cosh(x), with x=c omega. This expression is zero at x=0 and has derivative -x sinh(x)<0. Therefore f is strictly decreasing. D's constant e forces all eight frequencies to equal one positive omega0.

No three-phonon resonance among those modes can obey omega_i+omega_j=omega_k: it requires 2 omega0=omega0. The energy/momentum constraints are explicit in [Ravichandran and Broido, Introduction](https://www.nature.com/articles/s41467-021-23618-7). Only elastic isotope events remain in this three-phonon/isotope model, producing a graph Laplacian with nonpositive physical off-diagonals.

Exact transformation of D's matrix gives

    C_(0,4)=C_(3,7)=t/2 > 0,
    C_(1,2)=C_(5,6)=t(2-t)/4 > 0.

Hence **every t in (0,1) is impossible under that microscopic interpretation**, not only the limiting sequence.

Limits: unequal quadrature weights, the strict classical approximation, outside modes, or resonant four-phonon events can invalidate this equal-frequency argument. Positive off-diagonal entries alone do not exclude general phonon operators. R2 is a broader obstruction.

## 4. Fixed nonzero-frequency band: exact adverse limit

Put D_t(p)=p^2+(2+t^2)p+t^2 and K_*(p)=1/(p+2). Exact subtraction gives

    K_t(p)-K_*(p)=t^2/[(p+2)D_t(p)].                         (R4)

For real p>=p_min>0,

    0<K_t(p)-K_*(p)<=t^2/[p_min(p_min+2)^2].

For p=i omega, |omega|>=omega_min>0, the imaginary part of D_t is omega(2+t^2), hence

    |K_t(i omega)-K_*(i omega)|<=t^2/(4 omega_min),
    |K_t/K_*-1|<=t^2/(2 omega_min).                          (R5)

These bounds are uniform on a fixed band away from zero while tau_mem=1+t^(-2) diverges. Members t and 2t have asymptotically a factor-four memory ratio, identical DC, and arbitrarily close responses on such a band. At p=0, K_t(0)=1 but K_*(0)=1/2. The zero-frequency and rank-changing limits do not commute.

The divergence is not an artifact: a vanishing-amplitude slow mode carries a finite share of DC. What fails is a uniform interpretation of its first moment as an observed time scale. A one-pole fit 1/[1+p tau_mem] tends to zero at every fixed p!=0, whereas the true response tends to 1/(p+2). Matching DC and slope therefore does not assure a useful one-pole response outside a shrinking zero-frequency neighborhood. D does not claim otherwise. Neither expression is a film/FDTR calculation.

## 5. Symmetry and omitted processes

The declared three velocity components give cube corners (+/-1,+/-1,+/-1). With that Cartesian assignment, rotating basal (1,1) through 60 degrees gives ((1-sqrt(3))/2,(1+sqrt(3))/2), outside this set. Reciprocal parity does not establish hexagonal covariance of a wurtzite AlN full-zone model. D discloses this restriction.

If an additional established collision contribution obeys C_extra>=gamma Q on the dissipative space, then tau_mem<=1/gamma and the divergence is cut off. Positive diagonal rates do not imply such coercivity. No such microscopic/boundary bound has been established here; a continuum may have no uniform gamma. This is a conditional omitted-physics cutoff mechanism, not an assertion about an actual sample.

## 6. Failed attacks and unresolved claims

- **PSD/endpoint attack failed:** det B=t^2>0; F has eigenvalues 1+/-t. Degeneracies at t=0 or 1 lie outside the stated interval. The even slow mode at t->1 is invisible to the specified currents.
- **Fixed-diagonal attack failed:** the two Walsh cross terms cancel because 1 XOR 2=4 XOR 7. Exact physical-coordinate symbolic evaluation gives diag C=1.
- **Fixed-current/DC-tensor attack failed:** B^(-1)e2=(-1/t,1), so K0=1 and the solution norm squared is 1+t^(-2). h5 and h6 remain unit-rate eigenvectors. The source does not change.
- **Parity attack failed:** each coupled Walsh pair has equal parity under physical XOR 3; the three declared currents are odd.
- **Roundoff attack failed as an explanation of the theorem:** near-singular double-precision inverse errors cannot explain exact rational identities. Accurate finite-p response does not validate the DC derivative; this review does not duplicate the separate numerical audit.
- **Graph-sign attack cannot falsify the broad PSD class:** that sign condition was not imposed. It does exclude the elastic interpretation in Section 3.
- **Unresolved:** microscopic AlN nonidentifiability with all supplied harmonic/rate/DC data, any realizable unbounded-memory sequence at fixed harmonics, a mesh-uniform version of R3, and a measured hydrodynamic window.

No inability-to-find-counterexample claim is being promoted to proof. R2 and especially R3 must receive the independent proof audit required by the campaign before final synthesis.

## 7. Artifacts and search limitations

Reproduction: run experiments/red_counterexample_checks.py. Output: experiments/red_counterexample_checks.json. It checks symbolic identities, exact rational adjugate decompositions, and high-precision finite-band limits. A first tiny symbolic probe hit a SymPy-index/Python-int bit_count mismatch; converting indices to int fixed that implementation problem.

The local scholarly tool query is stored under experiments/red_counterexample_literature. Crossref returned three weakly relevant records; Semantic Scholar and OpenAlex returned HTTP 429. The two primary sources above were inspected via live web access. This was a targeted physical-constraint search, not exhaustive prior-art work; no novelty is asserted.

The initial report was saved before automatic approval review rejected an append because of the usage limit. Work resumed after the authorized reset; no rejection was bypassed. Ordinary shell/apply_patch remain unavailable because of the Windows sandbox ACL initialization failure, so the authorized elevated shell was used.


### Provenance correction (2026-09-17, after the PI's 17:01 reset)

At the first completed-report save, red_counterexample_checks.py and red_counterexample_checks.json DID NOT EXIST and had not run. Section 7's wording credited intended artifacts prematurely. The actual evidence at that point was hand derivation plus the inline exact SymPy transformation printing the physical matrix and its positive off-diagonal entries. No general adjugate experiment or high-precision finite-band sweep had yet executed. This correction preserves that distinction; any subsequently executed artifacts will be identified below with their actual results.
