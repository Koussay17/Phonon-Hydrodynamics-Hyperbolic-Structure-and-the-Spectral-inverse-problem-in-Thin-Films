# Branch B: variational inverse problem for a conserving collision operator

Independent first pass, 2026-09-16. Only the campaign question, assumptions, protocol, and note 17 were read; no other branch reports were consulted. Results below are **DERIVED UNDER ASSUMPTIONS; AWAITING PROOF AUDIT**, not a claim of novelty or of a realizable alternative AlN force-constant model.

## 1. Precise problem and information granted

Work in a finite-dimensional real entropy-coordinate mode space R^n, for a bulk periodic crystal in linear response at fixed positive temperature. Grant full-zone harmonic data, including weights and the reconstruction of inversion/crystal stars. The normalized or unnormalized energy vector e has nonzero entries. For direction i, the heat-current source is b_i = diag(v_i)e. Inversion symmetry normally gives e^T b_i = 0. No boundary problem is inferred from this bulk problem.

Grant the strongest plausible interpretation of the supplied total rates: r_mu is the actual diagonal C_{mu,mu} of the entropy-symmetric collision matrix. This convention must be checked against the data format; a quoted lifetime/out-scattering parameter is not automatically the diagonal after a conserving model projection.

The minimal admissible class is

    F(r,e) = {C = C^T >= 0 : diag(C) = r, C e = 0}.

Known crystal symmetries, detailed process support, or further exact invariants can be imposed as additional constraints. Statements about F(r,e) do not assert that every member comes from the actual AlN three-phonon kinematics. Symmetric positive semidefiniteness and energy conservation are necessary structural conditions, not a complete characterization of microscopic phonon collision matrices.

For a fixed source b perpendicular to e and real p > 0, define

    F_p(C) = b^T (C + p I)^(-1) b.

Physical conductivity is F_p/(k_B T^2). The physical DC value is lim_{p down to 0} F_p: it is b^T C^+ b if b is perpendicular to ker(C), and +infinity otherwise. Writing a pseudoinverse alone would incorrectly suppress a ballistic divergence. If ker(C)=span(e), the DC solution exists and is unique after fixing y perpendicular to e. All derivations below are finite-dimensional; no continuum convergence is established.

### Feasibility lemma

A matrix in F(r,e) exists if and only if the nonnegative lengths l_mu = |e_mu| sqrt(r_mu) obey max(l_mu) <= sum of the others. This is a PSD-feasibility statement and does not assert that the only null vector is e.

Derivation: write C as the Gram matrix of vectors a_mu with lengths sqrt(r_mu). Then C e = 0 is equivalent to sum e_mu a_mu = 0. The triangle inequality gives necessity. Conversely, vectors of lengths l_mu forming a closed planar polygon provide such a Gram matrix; the usual polygon existence condition is precisely the displayed inequality. Degenerate zero lengths are allowed. If r_mu=0, PSD forces the entire corresponding row/column to vanish, so additional null vectors must be treated explicitly.

## 2. Variational characterization and guaranteed bounds

Let Q be an n by m matrix with orthonormal columns spanning the orthogonal complement of the declared exact invariants. Put A = Q^T C Q and beta = Q^T b, assuming b is perpendicular to those invariants. When A > 0, or when p > 0,

    beta^T (A+pI)^(-1) beta
      = max_z [2 beta^T z - z^T(A+pI)z].

Completion of the square establishes the identity and the unique optimizer z=(A+pI)^(-1)beta. At p=0 and singular A, the same supremum is infinite unless beta lies in the range of A; otherwise the maximizers are nonunique by ker(A).

Consequences:

1. The response is convex in the collision matrix on an affine admissible class. This is the standard matrix-fractional-function/Schur-complement structure.
2. For p >= 0, minimizing the extended response over a specified convex admissible class is an SDP: minimize t subject to the physical affine/PSD constraints and

       [[Q^T C Q + pI, beta], [beta^T, t]] >= 0.

   The Schur condition handles the singular DC range condition correctly. The upper endpoint is a convex maximization problem, not generally an SDP minimization. For p>0 compactness of F(r,e) ensures both endpoints exist.
3. A trial z gives a lower bound, not an equality unless it solves the collision equation. A diagonal RTA value is not automatically the value of this variational functional evaluated at its putative RTA solution.

With T_r = sum r_mu > 0, PSD gives lambda_max(C) <= T_r, hence

    ||b||^2/(p+T_r) <= F_p(C) <= ||b||^2/p       (p>0).

The finite DC value therefore obeys F_0 >= ||b||^2/T_r. Taking a one-coordinate trial y=t e_mu (here e_mu denotes a coordinate unit vector) in the full-space variational formula also gives

    F_0 >= max_{mu:r_mu>0} b_mu^2/r_mu.

These bounds are generally weak. There is no finite universal DC upper bound from positive diagonal entries alone (Section 4). An independently established dissipative gap gamma gives F_0 <= ||b||^2/gamma.

If F_0 is finite and nonzero, define the exact memory time

    tau_mem = ||C^+b||^2 / (b^T C^+b).

Cauchy--Schwarz gives tau_mem >= F_0/||b||^2, with equality exactly when b lies in one positive-eigenvalue eigenspace. A known gap gamma gives tau_mem <= 1/gamma. These are operator-level statements; replacing C by the modal diagonal is an extra model assumption.

Known harmonic data do fix heat capacity and the high-frequency coefficient ||b||^2/(k_B T^2). The next resolvent coefficient, -b^T C b/p^2, is generally not fixed by the diagonal. Likewise, sum r_mu fixes tr(C), but not its positive eigenvalues or eigenvectors.

## 3. General local nonidentifiability proposition

**Proposition B1.** Assume n>=4; all entries of e are nonzero; and there is C_0 in F(r,e) with ker(C_0)=span(e). For every nonzero b perpendicular to e and every fixed real p>=0, the response F_p is nonconstant in every relative neighborhood of C_0 in F(r,e). For p=0 all sufficiently small perturbations in the construction retain finite DC response. Thus diagonal entries, detailed balance/PSD, exact energy conservation, and even the assertion of a finite DC conductivity do not identify a nonzero heat-current response in this class.

**Derivation.** Define the tangent space

    T = {X=X^T : diag(X)=0, X e=0}.

Every sufficiently small signed perturbation C_0+tX remains PSD with exactly the same nullspace, since its restriction to e-perp is a small perturbation of a positive definite matrix. Its rates and invariant remain exact.

The common kernel of all X in T is span(e). To see this, put E=diag(e). For any four indices embed either of the symmetric matrices

    A1 = [[0,1,-1,0], [1,0,0,-1], [-1,0,0,1], [0,-1,1,0]],
    A2 = [[0,1,0,-1], [1,0,-1,0], [0,-1,0,1], [-1,0,1,0]].

Both have zero diagonal and zero row sum. Thus X=E^(-1) A_j E^(-1) lies in T. The common kernel of A1 and A2 on those four coordinates consists of constant vectors. Applying this to all four-index subsets shows that Xf=0 for every X in T implies f_mu/e_mu is constant, hence f is proportional to e. The reverse inclusion is immediate.

Now let f be the solution in e-perp of (C_0+pI)f=b. It is nonzero and cannot be proportional to e. Choose X in T with Xf nonzero. Differentiating the finite-dimensional inverse on e-perp gives

    dF_p(C_0+tX)/dt at 0 = -f^T X f,
    d^2F_p(C_0+tX)/dt^2 at 0
        = 2 (Xf)^T [(C_0+pI)|_(e-perp)]^(-1) (Xf) > 0.

So the response cannot be locally constant. No dimension-counting heuristic is used in the decisive step.

**Limits of B1.** Extra symmetry, sparsity, microscopic event constraints, or more invariants shrink T. For any such proposed class, the correct test is whether its admissible tangent directions satisfy Xf=0; the displayed proof does not silently grant directions excluded by new constraints. With only energy conservation, n=2 or n=3 can have a unique matrix once the diagonal is given, so n>=4 is material. The existence of a C_0 with precisely one null vector is an explicit hypothesis, not inferred from positive rates.

## 4. Four-mode constructive counterexamples with additional physical structure

Take equal mode energies/weights, e proportional to (1,1,1,1). Let a,b,c >= 0 and r=a+b+c:

    C(a,b,c) = [[r,-a,-b,-c],
                [-a,r,-c,-b],
                [-b,-c,r,-a],
                [-c,-b,-a,r]].

This is the generator of a reversible elastic-scattering graph. It is symmetric, conserves e, has the prescribed diagonal r, and has nonpositive off-diagonal entries, a constraint stronger than required of a general phonon collision operator. If a,b,c>0 the graph is connected, so the only invariant is e. It commutes with inversion R exchanging 1<->3 and 2<->4.

The orthonormal eigenvectors and eigenvalues are

    u0=(1,1,1,1)/2,       lambda0=0,
    ux=(1,1,-1,-1)/2,     lambdax=2(b+c),
    uy=(1,-1,1,-1)/2,     lambday=2(a+c),
    uz=(1,-1,-1,1)/2,     lambdaz=2(a+b).

The current vectors ux and uz are odd under R, whereas e and uy are even. Velocities proportional to the components of ux therefore supply a parity-consistent current b_source=ux, while all modal harmonic data remain fixed.

### B2: fixed rates, arbitrarily long relaxation and unbounded DC response

For fixed r>0 set a=r-epsilon and b=c=epsilon/2, with 0<epsilon<r. Then

    F_p = 1/(p+2 epsilon),
    F_0 = 1/(2 epsilon),
    tau_mem = 1/(2 epsilon).

The other two positive eigenvalues are 2r-epsilon. All four diagonal rates remain exactly r. The dissipative gap tends to zero even though every modal rate is bounded away from zero, each member has a unique equilibrium, and time-reversal symmetry holds. This identifies an actual slow collective eigenvector rather than a long modal lifetime.

The RTA value for this normalized source is 1/r. The exact response is larger for epsilon<r/2 and smaller for epsilon>r/2. Therefore the given structural assumptions do not make RTA a universal lower bound. Any stronger lower-bound claim needs extra hypotheses about the collision kernel; a generic appeal to a variational principle is insufficient.

### B3: even an exact DC conductivity does not identify memory

Use units r=1 and fix the odd source b_source=ux+uz=(1,0,-1,0). Two admissible connected graph operators are:

| (a,b,c) | (lambda_x,lambda_y,lambda_z) | F_0 | tau_mem |
|---|---|---|---|
| (2/5,1/5,2/5) | (6/5,8/5,6/5) | 5/3 | 5/6 |
| (1/2,1/4,1/4) | (1,3/2,3/2) | 5/3 | 13/15 |

Both have the same harmonic data, diagonal rates, energy invariant, inversion symmetry, and exact DC value. Their memory times differ. These rational identities were independently checked by exact Fraction arithmetic (no floating-point eigensolver). This example establishes nonuniqueness, not an unbounded memory at fixed DC.

**Physical scope.** Equal-energy reversible elastic redistribution is compatible with entropy production and energy conservation. These four-mode examples are not a claimed realization of the actual AlN force constants, momentum-selection rules, all space-group irreducible representations, or a measured sample. Their role is to falsify reconstruction based only on the specified information and structural axioms.

## 5. Measured diagonal versus a conserving relaxation parameter

For D=diag(d_mu)>0 and a full-column-rank invariant matrix H, the projected relaxation model is

    C_D = D - D H (H^T D H)^(-1) H^T D.

Its variational meaning is exact:

    y^T C_D y = min_alpha (y-H alpha)^T D(y-H alpha).

This gives PSD and ker(C_D)=Ran(H), and describes a weighted projection onto local equilibrium. It is related to a shorted positive operator/Schur complement.

Crucially,

    (C_D)_(mu,mu)
      = d_mu - d_mu^2 h_mu^T (H^T D H)^(-1) h_mu.

If H includes an energy vector with nonzero entries, this is strictly smaller than d_mu. Therefore entering measured diagonal rates as d_mu does not produce a collision matrix with those measured diagonals. The parameters d_mu belong to the surrogate before projection. Fitting them so that diag(C_D)=r is a separate nonlinear inverse problem, with existence/uniqueness not established here.

A useful qualification prevents an overcorrection of note 17: if Q_H is the orthogonal projector onto Ran(H)-perp, then

    C_D^+ = Q_H D^(-1) Q_H.

Indeed C_D D^(-1) b=b for b perpendicular to H, and the displayed expression chooses the solution perpendicular to H. Thus the static response of this projected model is exactly b^T D^(-1)b for b perpendicular to H, even though its diagonal is not D. This equality validates an observable of that surrogate, not identification of the true operator. With energy-only H, even D under inversion, and odd b, the projection term vanishes on the entire odd subspace, so the uniform current's complete Laplace response also equals RTA there. For more invariants, or absent that parity structure, this dynamic statement need not hold.

## 6. Viscosity and the missing collision information

Grant an orthonormal slow basis U and its orthogonal complement Q, and additionally assume the slow space is invariant under C (so U^T C Q=0). Let A=Q^T C Q have a verified gap gamma>0 and let V_j be the diagonal streaming velocity operators. The leading spatial Schur-complement coefficient has entries

    N_(ia,jb) = (Q^T V_i U_a)^T A^(-1)(Q^T V_j U_b).

For a Fourier direction k the contracted quadratic form is nonnegative. It is another inverse collision quadratic form, now driven by projected stress-like vectors rather than the heat current. Frequencies, velocities, and diagonal rates do not supply A^(-1), its eigenvectors, or its gap. If quasi-momentum is used to define U, coherent full-zone q-vectors are also required.

The formula itself is conditional: it arises from eliminating fast modes in p+C+i k dot V and expanding their resolvent about A. A sufficient finite-dimensional small parameter is

    (|p| + ||Q^T(k dot V)Q||)/gamma < 1,

with slow collision rates also small compared with gamma for a hydrodynamic interpretation. If collision slow/fast blocks do not vanish, the general Schur complement contains extra collision cross terms. A projected RTA rate split is not evidence that these hypotheses hold for the true operator.

Consequently there is no justified route from the supplied total diagonal rates alone to a numerical GK viscosity, a normal/resistive split, or a verified hydrodynamic window. Even exact DC conductivities add only a few scalar inverse constraints, as B3 illustrates; they do not supply the needed operator blocks.

## 7. Earliest fragile step, failed routes, and useful next data

**Earliest fragile inference:** treating total rates as independent relaxation eigenvalues, or assuming that a conserving projection fed those numbers retains the measured diagonal. Conservation/PSD do not fill the missing off-diagonal collision information.

Failed or limited routes:

- Counting unknown off-diagonal entries is suggestive but does not establish variation of a chosen observable. B1 replaces it with a tangent-space derivative argument.
- Positive modal rates do not lower-bound a collective collision gap. B2 is an exact failure.
- A generic variational principle does not certify the diagonal RTA conductivity as a lower bound. The trial functional includes the full C, and B2 crosses the RTA value.
- Fitting one DC value does not determine the first temporal moment. B3 is an exact failure.
- A fully conserving surrogate is still an extra constitutive assumption; PSD and invariants certify its consistency, not fidelity to AlN.
- The finite-dimensional examples do not settle which operators a fixed AlN interaction Hamiltonian can realize. That stronger question requires microscopic input.

Most useful next data: matrix-vector access to the symmetrized collision operator; exact invariants and symmetry action; event-resolved normal/resistive information if a momentum-based interpretation is intended; low eigenpairs or certified coercivity on the proposed fast subspace; and collision solves for heat-current and stress sources. Full dense storage is unnecessary if matrix-vector products and residual-controlled solves are available. Grid convergence of the relevant inverse moments remains a separate task.

## 8. Known frameworks, terminology, and limited literature cross-check

These calculations are not presented as novel. Useful equivalent search terms are: matrix fractional function; semidefinite inverse problems; PSD matrix completion with prescribed diagonal and kernel; correlation elliptope and realizable subspaces; Dirichlet/Thomson variational principle; effective resistance and graph bottlenecks; shorted positive operator; weighted least-squares projection; Schur complement/Feshbach reduction; Stieltjes resolvent and moment problem; observable controllability/identifiability.

Two primary checks were made during this independent pass:

- Boyd and Vandenberghe, *Convex Optimization*, Example 3.4, explicitly treats x^T Y^(-1)x as a convex matrix fractional function and gives the block-PSD epigraph used here: https://www.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf .
- Chaput, *Direct Solution to the Linearized Phonon Boltzmann Equation*, PRL 110, 265506 (2013), establishes the relevant pre-existing frequency-dependent spectral conductivity framework: https://doi.org/10.1103/PhysRevLett.110.265506 ; preprint https://arxiv.org/abs/1303.4062 .

The specific local nonidentifiability formulation B1 and the explicit examples still require the campaign proof audit. Literature coverage here is deliberately narrow; the designated literature branch should trace the broader matrix-completion, network, and phonon variational connections before any novelty language.

## 9. Reproducible exact check (standard Python only)

```python
from fractions import Fraction as F

def collision(a, b, c):
    r = a+b+c
    return [[r,-a,-b,-c], [-a,r,-c,-b],
            [-b,-c,r,-a], [-c,-b,-a,r]]

def matvec(A, x):
    return [sum(a*y for a,y in zip(row,x)) for row in A]

U = [[F(t,2) for t in signs] for signs in
     [[1,1,1,1], [1,1,-1,-1], [1,-1,1,-1], [1,-1,-1,1]]]

for a,b,c in [(F(2,5),F(1,5),F(2,5)),
              (F(1,2),F(1,4),F(1,4))]:
    A = collision(a,b,c)
    lam = [0,2*(b+c),2*(a+c),2*(a+b)]
    assert all(matvec(A,u)==[l*x for x in u] for l,u in zip(lam,U))
    K = 1/lam[1]+1/lam[3]
    tau = (1/lam[1]**2+1/lam[3]**2)/K
    print(K, tau)
# 5/3 5/6
# 5/3 13/15
```
