# B — Degenerate amplitudes and kinetic closure

Independent bounded first pass, 2026-09-25; saved after the usage interruption. Read only campaign question/assumptions and primary literature. No peer reports, PI experiments, or cross-agent conclusions were consulted. No material generator or novelty is claimed.

**Conclusion:** summing squared amplitudes over complete degenerate blocks identifies an invariant coupling strength. It does not by itself establish autonomous block-population kinetics. The latter requires a specified collision approximation and a closure condition or separately justified fast equilibration.

## 1. Variables and information lost by squaring

Take distinct exact-degenerate blocks A,B,C, with positive energies epsilon_A=epsilon_B+epsilon_C. Physical polarization blocks in a homogeneous crystal also share q: mixing equal-frequency states at different q while retaining unchanged momentum labels is not a harmless polarization gauge change. A resonant cubic term is

    H_int=sum_(alpha,beta,gamma) T_(beta gamma,alpha)
          b_beta^dagger c_gamma^dagger a_alpha + h.c.,
    T: H_A -> H_B tensor H_C.

T includes the chosen counting convention; repeated daughters require bosonic symmetrization. A finite Hamiltonian does not alone specify an irreversible collision rate.

Use one-body occupation matrices N_A, with entries <a_alpha'^dagger a_alpha>. Total population is M_A=Tr N_A; mean occupation is n_A=M_A/d_A. Under compatible block unitaries, T transforms as (U_B tensor U_C)T U_A^dagger and N_A as U_A N_A U_A^dagger. Consequently

    S=||T||_HS^2=sum |T_(beta gamma,alpha)|^2

is invariant. Individual diagonal populations and squared entries are basis dependent. Second-order kinetic contractions also require

    G_A=T^dagger T,
    G_B=Tr_C(T T^dagger),  G_C=Tr_B(T T^dagger).

P=|T|^2 determines their diagonals in the exported basis, but not their phase-sensitive off-diagonal entries. Polarized inputs further require maps such as T^dagger(Y_B tensor Z_C)T.

## 2. Conditional scalar-block equations

Assume a covariant weak-coupling, Markov, quasifree kinetic approximation exists, with common resonant weight k for this channel. At isotropic inputs N_X=n_X I, golden-rule contractions take the form

    F=n_A(1+n_B)(1+n_C)-(1+n_A)n_B n_C,
    dot N_A=-k F G_A,  dot N_B=k F G_B,  dot N_C=k F G_C.

Their instantaneous total flux is k F S. Remaining isotropic is a stronger requirement: for one active isolated channel, it requires G_X=(S/d_X)I. A complete network needs the corresponding weighted condition. Equilibrium F=0 masks this failure; perturbed dynamics must be tested.

If the isotropic manifold is invariant, dot M_A=-kFS and dot M_B=dot M_C=kFS, with n_X=M_X/d_X. Dimensions cannot be discarded. More generally, let Pi(N)_X=(Tr N_X/d_X)I. For a specified linearized matrix generator L, autonomous traces for arbitrary initial matrices require

    Pi L(I-Pi)=0.

Prepared isotropic data remain isotropic if (I-Pi)L Pi=0. Neither follows from block-summed P. These are algebraic closure tests, not a derived AlN generator.

## 3. Explicit insufficiency example within a stated bath model

Let d_A=d_B=2, d_C=1 and

    T1=[[1,1],[1,1]],  T2=[[1,1],[1,-1]].

Both have identical P (all ones) and S=4. Yet G_A,1 has eigenvalues (4,0), whereas G_A,2=2I. The first has a dark parent combination; the second does not.

For a conditional kinetic realization, hold the daughters in independent thermal reservoirs with short correlation times and equal channel spectral weights. Assume a weak-coupling degenerate-bath master equation, with no relevant noncommuting Lamb shift. Its parent covariance deviation obeys

    dot(delta N_A)=-(gamma/2){G_A,delta N_A}, gamma>0.

The thermal gain term fixes N_A^0=n_A^0 I; gamma includes Bose and dimensional factors. From delta N_A(0)=xI near positive-temperature equilibrium,

    Tr delta N_A,1(t)=x[1+exp(-4 gamma t)],
    Tr delta N_A,2(t)=2x exp(-2 gamma t).

Initial traces and first derivatives agree, but later traces differ. This establishes insufficiency in the stated bath model. It is NOT a derivation of closed three-phonon crystal kinetics or a claim that finite isolated modes obey a Lindblad equation. Exact Hamiltonian evolution generates higher correlations, so one-body closure also requires justification. [Davies' primary weak-coupling work explicitly concerns open systems coupled to a heat bath](https://link.springer.com/article/10.1007/BF01608389).

## 4. Physical limits, conserved quantities, importer

**Secular limit.** Discarding coherences needs a coarse-graining interval longer than bath memory, shorter than collision time, and satisfying |omega_i-omega_j| Delta t >>1 for discarded pairs. Exact degeneracy supplies no phase averaging. Splittings comparable to collision rates require retaining the coupled sector; near-degeneracy and weak coupling are distinct limits. One-body phonon matrix transport supplies established context, not the missing collision tensor. [Simoncelli et al., sections III–V](https://arxiv.org/pdf/2112.06897).

**When populations suffice.** A resolved nondegenerate secular basis or justified dephasing can support ordinary populations. Scalar degenerate-block occupations additionally need isotropy preservation or fast physical mixing. Symmetry can enforce isotropy on an irreducible representation when the full event set and state have the required covariance. Accidental degeneracy, repeated irreducible representations, or symmetry-breaking forcing invalidate that shortcut.

**Alternative variables and conservation.** Hermitian block charges are conserved by the resonant interaction when

    (X_B tensor I+I tensor X_C)T=T X_A.

This includes energy, normal-event momentum, and isolated-channel combinations M_A+M_B and M_A+M_C; Umklapp does not conserve additive phonon momentum. Kernels and singular subspaces of T identify dark and bright combinations. A kernel common to all channels can reveal extra conserved populations. Retaining those matrix components is a useful alternative to forced scalar-block closure.

**Importer consequence.** Preserve complex amplitudes with projector/gauge metadata, or sufficient covariant contraction maps. Label P-only output as specifying a chosen-basis population model only under an explicit population approximation. Record block dimensions and total/mean conventions. Test Gram isotropy, generator closure, symmetry, and conservation; invariant norms or equal linewidths are insufficient.

STATUS: amplitude algebra, closure criteria, and the stated bath example are derived conditionally. The actual degenerate phonon collision superoperator and its material validity remain unresolved. No numerical experiment or implementation edit was performed.