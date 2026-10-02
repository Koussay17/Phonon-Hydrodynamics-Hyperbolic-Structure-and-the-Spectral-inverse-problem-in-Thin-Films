# B2 — Scope of the conserving finite candidate and genuinely different alternatives

Second-generation cross-examination, 2026-09-28. Read A, C, D and L alongside B. **Bounded scope reconciliation and a candidate theorem for adversarial review; no new universal convergence or material claim. No additional computations.**

## 1. Four claims that must remain separate

**Declared physical dispersion versus a surrogate.** “Physical” here means the continuum function epsilon specified in the model, not an experimentally certified AlN dispersion. Including epsilon itself in the finite test space and evaluating it consistently gives an exact represented invariant. Including only an interpolant epsilon_h guarantees its represented invariant. If epsilon_h and epsilon agree at every volume node, their quadrature energies actually coincide; off-grid physical energy and Bose occupations can still differ. Full function inclusion is a convenient sufficient condition; the finite algebra only requires a coefficient vector e reproducing the energy at every volume and reaction evaluation point.

If ||epsilon_h-epsilon||_infinity<=a, a surrogate-resonant triple has physical mismatch at most 3a. This is a defect bound, not an exact nullspace. There is also a useful controlled-energy statement: if the conserved quadrature surrogate energy E_h uses epsilon_h>=epsilon_min>0, then, along any nonnegative solution,

    |E_Q,true(t)-E_Q,true(0)|
      <= a[N_Q(t)+N_Q(0)] <= 2a E_h/epsilon_min.

Thus a surrogate can give a controlled approximation without being renamed the physical invariant.

**Volume quadrature versus the continuum functional.** A/B conserve E_Q=Q_v(epsilon n_h). The exact functional E[n_h]=integral epsilon n_h obeys

    E[n_h(t)]-E[n_h(0)]
       =(E-E_Q)[n_h(t)]-(E-E_Q)[n_h(0)].

Uniform quadrature error eta along the trajectory bounds this difference by 2eta. Such a uniform bound requires control of the evolving reconstruction; a correct quadrature for one smooth test function does not supply it. Exact volume integrals would conserve E[n_h] in the finite-dimensional ODE, but actual numerical integration reintroduces this error.

**Finite algebra versus the correct resonance measure.** C's omitted-coarea-factor test preserves all invariants and equilibrium while giving a 20% rate error. D's vanishing-drift examples can miss essentially the entire measure. Hence algebraic structure and measure consistency need separate certificates. Fixed-integrand quadrature convergence is weaker than operator, spectral, or nonlinear evolution convergence.

**Local ODE versus global realizability.** Positivity of reconstructed occupations and invertibility of the capacity matrix hold on the declared open entropy-coefficient domain. They give local existence and uniqueness. They do not by themselves prove that the trajectory stays there for all time. Conserved positive energy bounds occupations at positive-weight volume nodes, but does not automatically control off-grid entropy evaluations or exclude approach to a zero-population boundary where coefficient coordinates diverge. Entropy-variable enrichment can introduce extrapolating evaluations. Nonnegative combinations of volume evaluation rows would preserve positivity at the additional evaluation points, but that extra compatibility has not been established for an energy-enriched space.

## 2. Exact finite theorem candidate for red review

This is the proposed theorem's complete finite scope; geometric convergence assumptions are deliberately separate.

**Hypotheses.**

1. Fix finitely many real basis functions, positive volume weights m_l, and a full-column-rank volume evaluation matrix Phi. Fix finitely many reaction triples and nonnegative, state-independent weights omega_r. All counting and units are declared once.
2. Let epsilon>0 at every evaluation point. A coefficient vector e represents this same energy at all volume nodes and all three legs of every reaction. The active triples obey epsilon_p-epsilon_a-epsilon_b=0.
3. For xi_alpha=sum_i alpha_i phi_i, let Omega be the open set where every required volume/leg value of xi_alpha is positive. Set n=(exp(xi)-1)^(-1). Take alpha(0) in Omega.
4. Use the same reconstruction on all legs. Define b_r,i=phi_i(p)-phi_i(a)-phi_i(b), the usual Bose factors A_r=n_p(1+n_a)(1+n_b), B_r=(1+n_p)n_a n_b, and their positive logarithmic mean Lambda_r. Define

       U=Phi^T diag(m) n,  M=Phi^T diag[m n(1+n)] Phi,
       K=sum_r omega_r Lambda_r b_r b_r^T,
       M dot alpha=-K alpha.                              (T)

**Candidate conclusions.** There is a unique maximal local solution in Omega. During its lifetime E_Q=e^T U is constant and dot S_Q=alpha^T K alpha>=0 for the quadrature Bose entropy. Every alpha=beta e, beta>0, is stationary and reproduces physical Bose occupations at all declared evaluation points, with each forward and reverse event individually equal. Equality with the physical Bose field everywhere additionally requires global function reproduction. The equilibrium entropy-coordinate generator

    C=M_*^(-1/2) K_* M_*^(-1/2)

is symmetric PSD and annihilates M_*^(1/2)e. Its nullity is exactly the nullity of the active row matrix b. In particular, unique thermal stationarity within this ansatz requires the additional rank condition ker b=span{e}; it is not supplied by resonance quadrature. Under that condition a prescribed E_Q>0 selects exactly one beta, since E_Q(beta) decreases continuously from infinity to zero.

**Short audit route.** dU/dalpha=-M, M>0, and Lambda is smooth on positive arguments, including A=B; these give local well-posedness. The two identities follow by testing with e and alpha. K e=0 holds for every alpha, so differentiating K at equilibrium produces no extra linearization term. Positive active weights give ker K=ker b. No global-domain, uniform gap, time-integrator, or approximation conclusion follows. Existence of a nonempty, accurately normalized reaction quadrature is a separate input; the zero collision operator also satisfies the finite conservation and entropy identities.

The earliest unsupported extension would be promoting this local finite theorem to globally realizable dynamics or to convergence of the physical collision operator. Regular C^2 surfaces, nonsingular charts, complete root coverage and correctly normalized coarea quadrature are needed for a separate integration theorem; they are not needed to prove the finite identities in (T).

## 3. Precisely what the positive-support obstruction forbids

With the same energy representation, unmodified event rows, positive common forward/reverse weights, and the usual Bose factors,

    e^T K e=sum_r omega_r Lambda_r Delta_r^2,
    dot E_Q at alpha=beta e = beta e^T K e.

Therefore fixed positive reweighting of unchanged off-shell triples cannot produce a conserving, thermally stationary model. Opposite detunings cannot cancel. This agrees with C/D without identifying a symmetrized surrogate matrix as an equilibrium Jacobian of a drifting nonlinear system.

It does **not** forbid changing the retained variables, event rows, mobility, dynamics, or conserved energy. For example, in fixed dimensionless moment coordinates put P=I-ee^T/(e^T e) and let K_0(U) be any symmetric PSD mobility. On a domain with a differentiable entropy S(U),

    dot U=P K_0(U) P grad S(U)

conserves e^T U, increases S, and fixes any state with grad S=beta e. This is a concrete projected alternative even when K_0 was assembled from off-shell samples. It modifies the event affinities and generally destroys the interpretation as the original labelled three-phonon events. Thermal stationarity here is not a proof of their microscopic eventwise Bose detailed balance. Global positivity and physical consistency remain separate tasks.

At fixed dimension, if K_0 already converges to a conserving operator in norm and the projector/metric is compatible, projecting can retain that limit. A small energy defect alone is insufficient, as C/D demonstrate.

## 4. Alternatives worth retaining, and what they cannot settle

**Dual and constrained variational formulations — useful.** B's convex potential Psi=-Q_v log(1-exp(-xi)) satisfies grad Psi=-U and Hess Psi=M; this identifies the entropy coordinates and the realizability problem exactly. L locates prior entropy-Galerkin work, so this principle is not a novelty claim. A distinct strategy works directly with nonnegative nodal populations: maximize Bose entropy minus a positive quadratic step cost on the compact set of fixed positive energy. A maximizer exists, preserves positivity/energy and cannot decrease entropy because the old state is feasible. This is a viable alternative to unrestricted coefficient evolution. Recovering the desired collision mobility as the timestep vanishes, and retaining additional invariants, would need a separate derivation; arbitrary step costs invent kinetics.

**Geometric charts — directly useful.** A's joint (parent,daughter) resonance manifold can be regular even when a fixed-parent slice is singular, as Delta=p-k^2 shows. This changes the integration strategy without pretending that an infinite pointwise rate is finite. Positive partitions of unity, periodic chart gluing and complete root coverage preserve the geometric measure. Chart existence does not supply a certified numerical atlas.

**Probabilistic state space — a genuinely different candidate.** A finite Fock-configuration Markov model can have nonnegative probabilities, exact pathwise energy conservation and reversible equilibrium on each connected energy sector. It retains factorial moments and correlations that a scalar population closure discards. Reducing it by conditional expectation can preserve a reversible quadratic form and represented conserved quantities, but exact retained dynamics generally has memory unless a lumpability/invariance condition holds. This avoids forcing microscopic dynamics into the entropy ansatz, at the cost of a larger state space and a new convergence/kinetic-limit problem. It does not make generic finite off-shell mode events conserve their original energies.

**Transform or spectral-energy formulations — distinct but presently speculative.** Fourier/time representations of delta are useful integration and asymptotic tools; a positive finite-time sinc-squared kernel remains off shell and obeys the same obstruction if inserted into an unchanged Markov event law. An enlarged model retaining spectral energy x=hbar*omega>0 as an independent variable could enforce x_p-x_a-x_b=0 exactly while using a broadened spectral distribution around each bare branch. Its conserved spectral energy and equilibrium n_beta(x) differ from a bare-mode population model. Deriving the spectral weights, interactions and zero-width limit would be new work; no physical validity is asserted here.

**Hamiltonian, symmetry, topology and asymptotics — limited but specific roles.** A finite isolated Hamiltonian supplies exact total energy and coherent dynamics; it does not by itself supply the dissipative scalar kinetic equation. Retaining interaction energy/coherences may change the legitimate conservation target, but cannot justify silently replacing it by bare-mode energy. Crystal/time-reversal group averaging can preserve PSD and a compatible energy invariant and reduce work; it cannot cancel positive off-shell squares. Umklapp does not furnish a globally additive real momentum invariant to import uncritically. Topological component tracking can expose missed resonance sheets, but neither their coarea density nor weighted integrability follows from topology alone. Weak-coupling/long-time limits and critical-root rescalings identify which model and integration scales are appropriate; D's examples rule out treating any single h-sigma relation as universal.

## Red-review targets and completion status

Attack: signs in (T); logarithmic-mean regularity; consistency of all energy evaluations; the omitted derivative-of-K term; exact kernel/rank claims; off-grid positivity; and any hidden substitution of quadrature energy for continuum energy. For approximation claims additionally attack missing roots, coarea normalization, critical weighted integrability, and spectral pollution/artificial invariants.

**Second-generation B scope report complete.** The leading result is a local finite structural theorem candidate, alongside explicit alternative models outside its obstruction class. Global realizability and a physical operator-convergence theorem remain unresolved. No additional material calculation, numerical experiment, or claim of novelty was introduced.
