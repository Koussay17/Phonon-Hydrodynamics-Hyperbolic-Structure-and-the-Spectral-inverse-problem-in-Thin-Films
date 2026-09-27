# B2 — Repeated daughters: population Jacobian versus linewidth

Second-generation reconciliation, 2026-09-26. Read A-D first passes. **Conditional finite derivations; awaiting audit. No phono3py bug, absolute material rate, or novelty is asserted.** Use one exactly resonant channel p <-> a+a, E_p=2E_a>0, dimensionless occupations, and rates in one fixed inverse-time unit. Linearize at finite inverse temperature 0<beta<infinity, so A=(exp(beta E_a)-1)^(-1)>0; the vacuum boundary, where W is singular, is excluded. All prefactors below refer to the same correctly counted unordered event.

## 1. The geometric mean model and its entropy factorization

Write parent/daughter means p,a and assume independent geometric mode distributions whenever evaluating the mean flux. Let kappa>0 include the repeated-event factorial/counting coefficient already derived in A-C. Then

    F(p,a)=p(1+a)^2-(1+p)a^2=p(1+2a)-a^2,
    J=kappa F,       (dot p,dot a)=(-1,2)J.

At Bose equilibrium (P,A), P=A^2/(1+2A). Put w_p=P(1+P), w_a=A(1+A), and Q=P(1+A)^2=(1+P)A^2. Detailed balance gives

    Q/w_p=1+2A,       Q/w_a=A-P,
    delta F=Q[delta p/w_p-2 delta a/w_a].

With incidence nu=(-1,2)^T, W=diag(w_p,w_a), the complete conditional Jacobian is

    dot x=-L_geo x,
    L_geo=kappa Q nu nu^T W^(-1)
          =kappa [[1+2A, -2(A-P)],
                   [-2(1+2A), 4(A-P)]].

For y=W^(-1/2)x,

    C_geo=kappa Q u u^T,    u=W^(-1/2)nu.

This is symmetric PSD, rank one, and has exactly the energy null direction e=(E_p sqrt(w_p),E_a sqrt(w_a)): u^T e=-E_p+2E_a=0. Positivity and the nullspace are consequences of the **declared closure**, not proof that the closure describes exact mean dynamics.

## 2. Why branch D finds ratios (1,2)

D's pinned ordered-linewidth arithmetic, expressed using this same repeated-event kappa, gives inverse-lifetime contributions

    r_p=kappa(1+2A),       r_a=2 kappa(A-P).

The complete Jacobian diagonals are therefore (r_p,2r_a). This is not a discrepancy in A-C's repeated-channel prefactor.

The elementary source of the extra contribution is the chain rule. Temporarily distinguish the two daughter occurrences x,y:

    F(p,x,y)=p(1+x+y)-xy,
    d/da F(p,a,a)=partial_x F+partial_y F=2(p-a).

The population update itself has multiplicity two. Both occurrences of the same physical population change together, producing 4 kappa(A-P), whereas a tagged one-particle response is a different derivative. Coincident occurrence terms that would be off-diagonal for different labels contribute to the same population diagonal.

Changing the common event prefactor by another factor 1/2 cannot fix both comparisons: it spoils the parent match. Replacing only the daughter diagonal while retaining the derived cross terms breaks energy conservation and PSD. Thus the source linewidth is not, without extra hypotheses, the diagonal of this full population operator.

## 3. A Fock master equation exposes the closure assumption

Choose an explicit Markov model with c=kappa/2>0 and integer state (m,n). Its jumps and rates are

    (m,n)->(m-1,n+2): c m(n+1)(n+2),
    (m,n)->(m+1,n-2): c (m+1)n(n-1).

Boundary rates suppress invalid states. Each trajectory preserves ell=2m+n, so each communicating sector is finite and the chain is nonexplosive. This is an assumed kinetic model after coarse-graining, not a derivation from a finite isolated Hamiltonian.

Its exact net event flux is

    J_exact=c[4 <mn>+2 <m>-<n(n-1)>].

On an independent geometric product with means (p,a), this reduces to kappa F(p,a). That assumption is sufficient, not necessary: other distributions can accidentally reproduce the same moment combination. In general means alone do not identify the flux. For example, parent geometric mean 1/3 and daughter fixed at n=1 have the same means as thermal geometric means (1/3,1), but J_exact=2c instead of zero.

The geometric product family is **not invariant**. Let H2=<n(n-1)>-2<n>^2. It vanishes on that family. Applying the exact jump generator at a geometric product yields

    dot H2=4c(1+2a) F(p,a).

At a=1 and p=1/3+epsilon, this is 36c epsilon. Thus even an arbitrarily small nonthermal tangent perturbation generates a nongeometric factorial moment at first order. Also,

    d Cov(m,n)/dt=4c(p-a) F(p,a)

at a product geometric state. Higher moments and correlations are dynamical variables, not fixed multiplicity constants.

At exact Bose equilibrium, the product geometric probability pi is reversible: pi(m,n)=pi(m-1,n+2), and the forward rate equals the destination's reverse rate. Its conditioning on ell is uniform; arbitrary mixtures of stationary sectors are also stationary. The full model consequently retains more conserved information than mean energy alone.

There is nevertheless a rigorous positive interpretation of C_geo. Let G be the backward Markov generator Gf=sum_{x'} rate(x,x')[f(x')-f(x)]. In L2(pi), use normalized number observables f_p=(m-P)/sqrt(w_p), f_a=(n-A)/sqrt(w_a). They are orthonormal, belong to the domain of G, and their Dirichlet matrix <f_i,-G f_j>_pi is exactly kappa Q u u^T: each jump changes these observables by +/-u, and each equilibrium directed mean rate is kappa Q. Hence C_geo is the exact equilibrium **compression** onto number variables, which need not form an invariant subspace.

## 4. A common quantum Markov example reconciles the two response types

An optional, explicitly additional model is

    L(rho)=c[D[M]rho+D[M^dagger]rho],
    M=p_hat (a_hat^dagger)^2,
    D[J]rho=J rho J^dagger-{J^dagger J,rho}/2.

Its Fock-diagonal restriction is precisely the chain above. On the algebraic occupation core its adjoint satisfies

    L^dagger(p_hat)=-c(2 N_a+1)p_hat,
    L^dagger(a_hat)=c(2 N_p-N_a)a_hat.

Perturb thermal product states by an infinitesimal coherent displacement. Using <N_a a_hat>=2A<a_hat> to first order gives initial amplitude damping coefficients c(2A+1) for the parent and 2c(A-P) for the daughter. Twice these coefficients are exactly (r_p,r_a) above.

This supplies a consistent model in which number-population slopes and twice amplitude slopes have the ratio (1,2). It does **not** establish that phono3py derives from this Lindblad model, that either exact response is a single exponential, or that an initial damping coefficient alone determines a spectral pole.

## 5. Reconciliation and alternative formulation

A-C count the channel consistently with D. The apparent disagreement starts when a single-particle half-linewidth is identified with every diagonal of a full, geometrically closed population Jacobian. The finite model distinguishes these observables, and the geometric family itself is not dynamically invariant.

For a probabilistic alternative, retain the reversible Fock master equation on finite energy sectors, its probability-entropy Dirichlet form, and the necessary factorial/correlation observables; or retain the exact projection memory when reducing to means. This avoids promoting a geometric ansatz to an invariant manifold. Material use still needs a justified kinetic limit, source conventions, mode sewing, and a specified closure.

## Verification record — 2026-09-27

An independent SymPy calculation applied the backward generator directly to n(n-1) and mn and averaged its polynomial outputs using geometric raw moments. The flux identity, both moment-defect derivatives, and both equilibrium entropy-factor identities had exactly zero symbolic residual. A separate NumPy ladder-matrix check used parent/daughter cutoffs 6 and 10, comparing only matrix elements with m<=3,n<=5; the two adjoint-operator residuals were 2.31e-14 and 1.95e-14. The cutoff check concerns the displayed algebraic identities away from artificial boundaries; it supplies no material rate or long-time closure validation.
