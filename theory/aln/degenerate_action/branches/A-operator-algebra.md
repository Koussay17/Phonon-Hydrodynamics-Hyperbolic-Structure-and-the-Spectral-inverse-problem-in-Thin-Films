# A — Operator algebra of degenerate populations

Independent first pass, 2026-09-25. Read only 00-question.md and 02-assumptions.md. No peer conclusions or implementation inspected. **Finite derivations below await independent audit; no novelty or material claim.**

## Domain and coordinates

Let H=direct_sum H_alpha, with finite dimensions d_alpha and harmonic energy E_alpha I on each exact-degenerate block; E_alpha>0 and T>0. Inside each block retain a full Hermitian occupation matrix rho_alpha>=0, not necessarily trace one. The linear perturbation space is V=direct_sum Herm(H_alpha), over the reals. Interblock coherences and their kinetic elimination are outside this declared domain.

The basis group is G=product U(d_alpha), acting by X -> U^dagger X U. Occupations on a chosen diagonal are n_i=rho_ii; block total is Tr(rho_alpha), and average occupation is that total divided by d_alpha. At Bose equilibrium rho_alpha=N_alpha I, with D_alpha=N_alpha(1+N_alpha)>0. Entropy coordinates are Y_alpha=delta rho_alpha/sqrt(D_alpha). Taylor expansion of bosonic matrix entropy at the scalar equilibrium gives the quadratic deficit (1/2)sum Tr(Y_alpha^2).

## Lemma 1 — Diagonal populations are not general gauge coordinates

Let Delta discard off-diagonal entries in the selected basis. A unitary maps the entire diagonal algebra onto itself **iff it is monomial**: a permutation times diagonal phases. Indeed conjugation must permute its minimal rank-one projections; the converse is immediate.

For diagonal X, the new diagonal is M_U diag(X), where (M_U)_ij=|U_ji|^2. This is not a group representation. For the two-dimensional Hadamard U, U^2=I but M_U^2=M_U differs from I. Projecting away coherences between rotations is a physical information loss, not an invertible coordinate change.

More directly, rho_plus/minus=nI +/- c sigma_x, with 0<c<n, have identical old populations. Hadamard rotation gives populations (n+c,n-c) and (n-c,n+c). Thus old populations do not identify new populations for general states. Observables must transform as full matrices: Tr(A rho) is invariant when both A and rho are conjugated.

## Lemma 2 — The invariant linear observables form the center

A linear functional Tr(A X) is invariant under every independent block unitary iff

    A=direct_sum a_alpha I_alpha.

Proof: invariance is equivalent to commuting with every block unitary; commuting with diagonal phases removes off-diagonal entries, and commuting with permutations equates the remaining diagonal entries within each block.

Hence the invariant linear data are block traces. Nonlinear spectral invariants also exist; this lemma does not classify them. The canonical orthogonal projection onto the center is

    E(X)_alpha=Tr(X_alpha) I_alpha/d_alpha.

It equals unitary group averaging. A selected population algebra is a basis-dependent maximal commutative algebra; their intersection over all block bases is precisely this center.

## Lemma 3 — Invariant observables need not have autonomous dynamics

Assume an explicitly supplied real-linear, self-adjoint, positive semidefinite entropy generator C on V:

    dY/dt=-C Y.

Define u_alpha=I_alpha/sqrt(d_alpha), q_alpha=<u_alpha,Y>_HS and E=sum |u_alpha><u_alpha|. Physical block totals differ from q only by known factors sqrt(d_alpha D_alpha).

For **every** initial perturbation, q obeys a closed constant linear equation iff

    E C (I-E)=0,

equivalently [E,C]=0 by self-adjointness. Then its unique matrix is Cbar_alpha,beta=<u_alpha,C u_beta>. Necessity follows by comparing initial data differing by ker(E); sufficiency follows directly from the projected equation and uniqueness of finite linear ODEs. Without self-adjointness only the displayed one-sided condition is required.

The same criterion applies to classical population vectors with E replaced by block averaging. E C E always defines a compression; it is not automatically an exact closure. Merely choosing initially uniform block populations does not remove subsequent leakage.

Full G-equivariance of a fixed C is sufficient: it implies commutation with group averaging E. **Basis covariance of a family C[V] is weaker:** conjugating both the amplitudes and operator can be consistent even when a fixed C[V] is anisotropic and does not commute with G.

## Counterexample — Conserving positive collision action without block closure

Take three distinct blocks: parent A of dimension two, and B,C of dimension one, with E_A=E_B+E_C. Keep one reversible event A1 <-> B+C; A2 is a spectator. Its valid classical entropy generator is

    C=gamma g g^T,   gamma>0,
    g=(-1/sqrt(D_A),0,1/sqrt(D_B),1/sqrt(D_C)).

It is PSD and annihilates the energy vector (E_A sqrt(D_A),E_A sqrt(D_A),E_B sqrt(D_B),E_C sqrt(D_C)). Yet Y=(1,-1,0,0) has all block sums zero while the parent block sum of -C Y is -gamma/D_A, nonzero. Arbitrarily small multiples are physical perturbations of positive Bose occupations. Energy conservation and exact degeneracy therefore do not supply block closure.

## Squared amplitudes: one sufficiency and two limitations

For three distinct blocks, a tensor transforms by a tensor product of unitary or conjugate-unitary factors. Its Frobenius block sum W=sum_ijk |V_ijk|^2 is invariant. Individual squared entries are not: vectors v_plus/minus=(1,+/-1)/sqrt(2) have identical squared entries, but Hadamard rotation gives (1,0) versus (0,1). The missing relative phases prevent determination of rotated population rates.

There is a precise limited sufficiency. Suppose a counted population-event model has a common known coefficient kappa A_Bose throughout this exactly degenerate channel. Each event has entropy vector g_ijk, and

    C=kappa A_Bose sum_ijk |V_ijk|^2 g_ijk g_ijk^T.

Projection onto normalized block sums sends every g_ijk to the same vector h, whose three nonzero entries are

    h_A=-1/sqrt(d_A D_A),
    h_B=+1/sqrt(d_B D_B),
    h_C=+1/sqrt(d_C D_C).

Therefore

    Cbar=kappa A_Bose W h h^T.

W identifies this **compression**, given the stated event convention and common prefactor. It does not establish the commutator condition in Lemma 3. The counterexample uses V=(v,0) and already defeats that inference. Repeated blocks/identical legs require their own symmetric-tensor counting convention and are not covered by this formula.

## Consequence for the importer

Keep separate objects: recorded-basis populations and squared amplitudes; invariant block sums; and full matrix occupations/amplitudes. Tag E C E as a compression unless closure is checked, for example by the exact finite commutator. Rotating event rates requires phase information; a density-matrix collision law additionally requires an explicitly justified kinetic model, not amplitudes alone.

**Earliest unsupported inference:** replacing “block-summed interaction weight is invariant” with “block populations have a closed collision operator.” Exact degeneracy identifies the basis freedom; it does not impose interaction isotropy.
