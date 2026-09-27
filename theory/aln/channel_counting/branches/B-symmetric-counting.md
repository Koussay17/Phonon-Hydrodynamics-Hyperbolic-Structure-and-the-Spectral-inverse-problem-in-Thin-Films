# B — Symmetric tensors and reversible-event counting

Independent first pass, 2026-09-26. Read only question/assumptions. **Conditional finite identities; awaiting audit. No claim that the tensor normalization below equals the pinned pp convention. No material or novelty claim.**

## 1. Declare a reference Hamiltonian before counting

Assume positive-frequency mode labels i, a specified reciprocal involution bar(i) preserving frequency, and a paired scalar gauge in which

    X_i=a_i+a_bar(i)^dagger,
    H3=(1/3!) sum_ordered(i,j,k) T_ijk X_i X_j X_k,
    T symmetric, T_bar(i),bar(j),bar(k)=T_ijk^*.

The X_i commute. A different Hamiltonian convention, especially one absorbing 1/3! into its tensor, changes the conversion below.

For c -> a+b, the all-incoming labels are (c,bar(a),bar(b)); their sum is a reciprocal vector exactly when q_c-q_a-q_b is one. The conjugate tuple (bar(c),a,b) supplies the reverse monomial. This sign assignment follows the displayed X convention, not a universal interpretation of an exported slot.

Assume an explicitly chosen Golden Rule/secular scalar-population approximation, with independent geometric occupation statistics for factorial moments. A finite isolated cubic Hamiltonian alone does not justify irreversible dynamics.

## 2. Permutation orbits and normalized symmetric coordinates

For an ordered triple with label multiplicities m_i, its S3 orbit has

    s=3!/product_i m_i!

distinct tuples. A normalized symmetric tensor coordinate is sqrt(s) times a symmetric ordered entry; therefore its squared norm contributes s|T_ijk|^2 to the ordered Frobenius sum. Deduplicating tuples without these factors changes the norm.

Exact positive-energy resonance E_c=E_a+E_b prevents the parent label c from coinciding with either reversed daughter label. Consequently only two cases occur:

| Daughters | Ordered S3 orbit | Collected coefficient g of a_a^dagger a_b^dagger a_c |
|---|---:|---:|
| a differs from b | 6 | T_c,bar(a),bar(b) |
| a=b | 3 | T_c,bar(a),bar(a)/2 |

“Same degenerate block” is not “same label”: distinct orthogonal modes use the first row. The second row requires identical momentum, branch, and chosen mode label.

If the stored tensor instead satisfies H3=sum Phi_ijk X_i X_j X_k, with Phi=T/6, these coefficients become 6 Phi and 3 Phi. A tensor convention must therefore be established before applying any factor to pp.

## 3. Fock matrix elements, Bose averages, and event coefficients

Write H_channel=g M+g* M^dagger. Apart from |g|^2, the forward Fock transition factors are

    r_c(r_a+1)(r_b+1)                    if a differs from b,
    r_c(r_a+1)(r_a+2)                    if a=b,

with reverse factors (r_c+1)r_a r_b and (r_c+1)r_a(r_a-1). These are integer Fock occupations r, not mean occupations n.

For a geometric mode distribution,

    <r(r-1)>=2n^2,
    <(r+1)(r+2)>=2(1+n)^2.

Let Kappa denote the common Golden Rule spectral multiplier, formally (2pi/hbar)delta(Delta E). After a justified integration/coarse-graining it gives a finite coefficient. It must not be evaluated as delta(0) for a discrete resonant tuple.

The reversible mean flux is

    J=kappa [n_c(1+n_a)(1+n_b)-(1+n_c)n_a n_b],
    kappa=(1+delta_ab) Kappa |g|^2
          =Kappa |T_c,bar(a),bar(b)|^2/(1+delta_ab).

Thus the repeated event has half the coefficient associated with the same reference T, not the same coefficient. Its incidence still has a factor two:

    nu=e_a+e_b-e_c;  nu=2e_a-e_c when a=b.

At equilibrium, the additional common Bose factor A=N_c(1+N_a)(1+N_b) gives the linearized event weight lambda=kappa A. Hamiltonian coefficient, Fock matrix element, spectral multiplier, kappa, and lambda are different objects.

A useful exact identity, for any symmetric daughter-pair integrand F, is

    (1/2)sum_ordered(a,b) F_ab
       =sum_unordered{a,b} F_ab/(1+delta_ab).

The factor 1/2 is valid for that complete ordered sum. Applying it again after deduplication, symmetry reduction, or another prefactor convention is not justified.

## 4. Reciprocal partners and reaction reversal are different quotients

A reversible event E=(c;{a,b}) already includes decay and absorption through one flux difference. Storing its reverse as another reversible event doubles its action.

Momentum reversal gives bar(E)=(bar(c);{bar(a),bar(b)}). It is self-reciprocal precisely when c=bar(c) and the daughter multiset is invariant under bar: the daughters can be individually fixed or exchange as a reciprocal pair. Equality is at complete mode-label level, not merely equal frequency or equal q-block.

For a reciprocal-closed set of M unordered events with F fixed events, the number of reciprocal orbits is (M+F)/2. Expanding J representatives gives 2J-F events. Blanket doubling fails for self-partners.

Even for two distinct equal-rate partners, replacing both incidences by twice one representative is wrong for a full population action:

    gamma(g g^T+bar(g)bar(g)^T) differs from 2gamma g g^T

in general. Disjoint triplets give rank two versus rank one. Scalar orbit weights can replace a pair only for a suitably invariant scalar contraction, not its missing mode updates.

## 5. Missing symmetry or phase data cannot be repaired by counting

If permutation symmetry of T is unverified, the monomial coefficient is the coherent sum of its ordered entries divided by 6. Six entries all +1 and six entries with three +1 and three -1 have identical squared entries but give coefficients 1 and 0 for a distinct triple. Averaging pp cannot enforce amplitude symmetry.

Likewise the reciprocal-label map above assumes a paired scalar gauge. At degenerate blocks, reversal can require a unitary sewing matrix, so q reversal alone does not provide the mode or amplitude map. Missing phases/sewing information must be reported.

For the pinned export, establish the exact Hamiltonian definition and event coverage before choosing numerical factors. With ordinary frequency nu, E=h nu and delta(Delta E)=delta(Delta nu)/h. The gamma half-linewidth convention is a separate conversion, not an event coefficient.

**Earliest fragile inference:** using one universal multiplier on pp while leaving tensor factorials, ordered coverage, repeated labels, reciprocal expansion, or spectral units unspecified.
