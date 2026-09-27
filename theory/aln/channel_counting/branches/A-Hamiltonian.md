# A — Hamiltonian orientation and channel counting

Independent first pass, 2026-09-26; saved after interruption. Read only this campaign's question/assumptions; no branch reports, PI results, or cross-discussion. All coefficients below depend on the declared convention and kinetic approximation. No material-rate or source-sewing consistency is inferred.

## 1. Canonical convention and physical orientation

Let lambda=(q,s), omega_lambda>0, epsilon_lambda=hbar omega_lambda, and [a_lambda,a_mu^dagger]=delta_lambda,mu. Choose a reciprocal gauge with an involution bar-lambda, omega_bar-lambda=omega_lambda and polarization e_bar-lambda=e_lambda^*. Displacements then use X_lambda=a_lambda+a_bar-lambda^dagger. An arbitrary exported gauge may instead require phases or a degenerate-block sewing matrix; squared amplitudes do not supply that map.

Define the ordered Hamiltonian

    H3=sum_(lambda,mu,nu) V_(lambda mu nu) X_lambda X_mu X_nu,
    V=W/3!,
    V_permutation=V,  V_(bar-lambda bar-mu bar-nu)=V_(lambda mu nu)^*.

W includes oscillator factors, masses, and crystal normalization. V has energy units, includes 1/3! already, and is supported on q_lambda+q_mu+q_nu=G. If an export represents W instead, use V=W/6. Actual complex permutation and sewing identities must be checked.

For a -> b+c, epsilon_a=epsilon_b+epsilon_c, use the all-incoming tuple (a,bar-b,bar-c), with q_a-q_b-q_c=G. Its conjugate tuple (bar-a,b,c) supplies b+c -> a:

    O=a_b^dagger a_c^dagger a_a,
    H_channel=A_(a;bc) O + A_(a;bc)^* O^dagger.

The positive-frequency channels of stored tuple (0,1,2) can be oriented as

    omega0=omega1+omega2: 0 -> bar-1+bar-2,
    omega2=omega0+omega1: bar-2 -> 0+1,
    omega1=omega0+omega2: bar-1 -> 0+2.

The last two use the reciprocal all-incoming tuple and therefore the conjugate amplitude under this convention. Global q reversal gives the reciprocal physical event, not the inverse reaction at unchanged labels.

## 2. Permutations and Fock-state matrix elements

Positive-frequency resonance excludes a=b, a=c, a=bar-b, and a=bar-c. If b!=c, all six permutations of (a,bar-b,bar-c) are distinct. If b=c, only (a,d,d),(d,a,d),(d,d,a), d=bar-b, are distinct. Thus

    A_(a;bc)=6 V_(a,bar-b,bar-c)   for b!=c,
    A_(a;bb)=3 V_(a,bar-b,bar-b)   for b=c.             (1)

Without verified symmetry, replace each expression by the coherent sum over DISTINCT ordered tuples producing O. Squaring entries before summing loses interference information. Equal magnitudes alone do not justify (1).

Forward squared Fock matrix elements are

    |A|^2 N_a(N_b+1)(N_c+1)       for b!=c,
    |A|^2 N_a(N_b+1)(N_b+2)       for b=c.

Reverse factors at a given Fock state are (N_a+1)N_b N_c and (N_a+1)N_b(N_b-1). A forward transition and its destination's inverse have equal matrix elements. A finite counting check: one parent with vacuum daughters couples with amplitude 6V to normalized |1_b,1_c>, and 3 sqrt(2)V to normalized |2_b>. Their squared strengths are 36|V|^2 and 18|V|^2.

The spectral golden-rule expression is

    rate_(i->f)=(2*pi/hbar)|<f|H_channel|i>|^2
                delta(epsilon_f-epsilon_i).            (2)

This is a transition density. A finite resonant Hamiltonian alone evolves coherently; it supplies neither an irreversible rate nor a numerical value for delta(0).

## 3. Additional approximation for mean-population kinetics

Assume weak anharmonicity, well-defined quasiparticles, justified continuum/Markov coarse-graining, and a resolved secular basis or separately justified suppression of coherences. Factorize different-mode correlations. Repeated daughters additionally require quasifree/geometric single-mode marginals:

    <N_b(N_b-1)>=2 n_b^2,
    <(N_b+1)(N_b+2)>=2(1+n_b)^2.

Their common factor two is a statistical factorial moment, not another inverse reaction. Without this closure, higher moments remain dynamical.

For one UNIQUE event with unordered daughters, accumulate repeated indices in r=e_a-e_b-e_c and write

    F(n)=n_a(1+n_b)(1+n_c)-(1+n_a)n_b n_c,
    dot n=-r K_e F(n).

Equations (1)-(2), together with the stated statistical closure, give

    K_e=72*pi*|V_(a,bar-b,bar-c)|^2
        *delta(epsilon_a-epsilon_b-epsilon_c)
        /[hbar*(1+delta_bc)]                            (3)
       =72*pi*|V|^2 delta(omega_a-omega_b-omega_c)
        /[hbar^2*(1+delta_bc)].

Identical daughters therefore have half the distinct-daughter coefficient in this convention. Their incidence remains r=(1,-2), so dot n_b=2 K_e F. Do not remove that multiplicity or add another factorial.

At Bose equilibrium, exact resonance gives the common forward/reverse factor

    F0=n_a^0(1+n_b^0)(1+n_c^0)
      =(1+n_a^0)n_b^0 n_c^0.

F0, K_e, A, and the Fock-state matrix element are distinct objects. With W_B,ii=n_i^0(1+n_i^0), y=W_B^(-1/2)delta n obeys

    dot y=-C_e y,
    C_e=K_e F0 u u^T,  u=W_B^(-1/2)r.

C_e is PSD and annihilates epsilon_i sqrt(W_B,ii). For repeated daughters C_bb=4 K_e F0/W_B,bb. Normal events conserve additive phonon momentum; Umklapp transfers hbar G to the lattice.

## 4. Units, duplicate counting, and application gate

If pp=|V|^2 in eV^2 is established in this convention, set F_THz=10^12 s^(-1), omega=2*pi*F_THz*f, and use hbar in eV s. A consistent integration weight g_f replacing the delta in numerical cyclic-THz coordinates gives

    K_e [s^-1]=36*pp*g_f/[hbar_eVs^2*F_THz*(1+delta_bc)]. (4)

This distribution/quadrature conversion does not certify energy conservation for a finite tuple. A half-linewidth gamma in cyclic THz has the usual lifetime rate 2 Gamma_angular=4*pi*F_THz*gamma. That rate is not automatically every diagonal entry of the complete collision action; repeated legs can add diagonal contributions.

The inverse reaction is already included in F(n). A reciprocal event is separate only if it is a different unordered physical event. Do not double a self-reciprocal channel, even when one all-incoming multiset produces both O and O^dagger. Ordered daughter exchange entered (1); fixed-root orbit weights require unfolding and reconciliation before forming unique-event rates.

A broadened energy mismatch generally breaks Bose detailed balance and exact harmonic-energy conservation of that tuple. Near-degenerate block rotations may require matrix occupations; a frequency cluster does not justify scalar closure or supply a sewing map. Applying (3)-(4) requires complex permutation/sewing checks, the 1/3! convention, units, energy integration, and unique-event bookkeeping. Otherwise only the conditional Hamiltonian/counting result is established. No implementation or material calculation was performed; independent audit remains required.