# Branch D: constructive inverse collision experiments

Date: 2026-09-16. Independent first pass; no other branch reports were read.

## Status and scope

**DERIVED UNDER ASSUMPTIONS, with exact symbolic and independent numerical checks.**
The elementary proofs below await the campaign's independent proof audit. No
novelty claim is made. The examples are finite-dimensional collision models,
not reconstructions of the actual AlN collision operator. The eight-mode
example uses the broad symmetric PSD, energy-conserving operator class; its
realizability by the allowed three-phonon/isotope channels of AlN has not been
established. The graph examples additionally have nonpositive off-diagonal
entries and admit an abstract reversible, equal-energy transition-network
interpretation.

**Main findings:**

1. Fixed diagonal rates and fixed current do not determine DC transport,
   memory, or the slow spectrum, even in connected weighted graph models.
2. Fixed diagonal rates **and fixed DC conductivity** do not determine memory:
   a four-mode graph pair gives exact memories 3/4 and 5/6 at the same DC value
   3/2.
3. In the broad operator class, an eight-mode family has identical unit diagonal,
   the same energy and current vectors, the same DC conductivity 1, and memory
   **1 + t^(-2)** for 0 < t < 1. Thus these data supply no finite uniform upper
   bound on memory. The same construction can preserve a complete three-current
   DC tensor while one direction's memory diverges.
4. Even the **entire uniform-frequency response of one current** need not
   identify slow modes orthogonal to that current.
5. Diagonal-rate positivity gives no lower bound on the nonzero collision gap.
   Diagonal-rate moments alone cannot certify controlled elimination of fast
   modes.

These statements distinguish what the supplied observables determine from what
additional physical restrictions might determine. They do not assert that two
actual AlN force-constant sets matching all Rao data have been constructed.

## Definitions and regime

Use entropy coordinates in a finite linearized kinetic model

    dy/ds = -C y,

with C = C^T >= 0, a prescribed normalized energy vector e, C e = 0, and a
prescribed current b perpendicular to e. For connected examples the nullspace
is exactly span(e). Define, on e-perpendicular,

    K(p) = b^T (C + p I)^(-1) b,       Re p > 0,
    K0   = b^T C^+ b,
    M2   = b^T (C^+)^2 b = ||C^+ b||^2,
    tau_mem = M2 / K0.

The physical conductivity contains the common prefactor 1/(k_B T^2), omitted
here. All numbers use the common diagonal collision rate as the frequency
unit. Multiplication C -> nu C multiplies K0 by nu^(-1), memory by nu^(-1),
and all eigenvalues by nu. No finite timestep, spatial grid, or boundary
condition is used; these are uniform bulk response experiments.

For a normalized orthogonal eigenbasis q_j, positive eigenvalues lambda_j,
and weights w_j = |q_j^T b|^2,

    K(p) = sum_j w_j / (lambda_j + p),
    tau_mem = [sum_j w_j/lambda_j^2] / [sum_j w_j/lambda_j].

Positivity ensures causal decaying modes. The Cauchy-Schwarz bound

    tau_mem >= K0 / ||b||^2

is valid, but it is a lower bound. None of the examples violates it.

## Construction 1: a fixed-degree four-mode graph

Let a,beta,c >= 0, a+beta+c=1, and

    L(a,beta,c) =
      [ 1    -a    -beta  -c   ]
      [ -a    1    -c     -beta]
      [-beta -c     1     -a   ]
      [ -c   -beta -a      1   ].

Each row sums to zero and every diagonal is exactly 1. The usual graph
quadratic form is

    x^T L x = sum_{i<j} w_ij (x_i-x_j)^2 >= 0.

Use the orthonormal vectors

    e  = (1, 1, 1, 1)/2,
    u1 = (1, 1,-1,-1)/2,
    u2 = (1,-1, 1,-1)/2,
    u3 = (1,-1,-1, 1)/2.

Direct multiplication gives eigenvalues

    0,  2(1-a),  2(1-beta),  2(1-c).

Consequently the constraints can be checked without an eigensolver. All these
operators also commute with the reversal permutation (1,2,3,4)->(4,3,2,1);
e is even and u1,u2 are odd.

### 1A. Same rates and current, arbitrarily different DC and memory

Choose a=0, beta=1-delta/2, c=delta/2 with 0<delta<=1. This graph is connected,
its positive spectrum is

    {2, delta, 2-delta},

and its diagonal is always (1,1,1,1). For the fixed current b=u2,

    K(p)=1/(delta+p),  K0=1/delta,  tau_mem=1/delta.

Thus the smallest collision eigenvalue tends to zero while no diagonal rate
changes. The same mode energies, heat-capacity weights and velocities can be
held fixed by fixing e and the diagonal velocity operator V with V e=b.

### 1B. Same rates and entire scalar response, different hidden gap

Keep the same L(delta), but prescribe b=u1. Then

    K(p)=1/(2+p),  K0=1/2,  tau_mem=1/2,

for every delta, although the collision gap is delta. Therefore a scalar
uniform-current response does not identify the whole collision spectrum.
An invisible slow mode need not obstruct a restricted closure if it is exactly
decoupled from the relevant dynamics; the example establishes missing spectral
information, not the necessity of retaining every such mode in every geometry.

### 1C. Same rates and DC, different memory

Prescribe the fixed current b=u1+u2=(1,0,0,-1). Compare

    C_A = L(1/3,1/3,1/3),
    C_B = L(1/2,0,1/2).

Both are connected graph Laplacians with diagonal 1 and the same energy/current.
Their responses are exactly

    K_A(p)=2/(4/3+p),
    K_B(p)=1/(1+p)+1/(2+p).

Hence

| Observable | C_A | C_B |
|---|---:|---:|
| Positive collision spectrum | {4/3,4/3,4/3} | {1,1,2} |
| K0 | 3/2 | 3/2 |
| M2 | 9/8 | 5/4 |
| tau_mem | 3/4 | 5/6 |
| K(1) | 6/7 | 5/6 |

The memory difference is 11.111...% relative to C_A. Their RTA diagonal data
are identical. Their finite-frequency responses differ despite equal DC.

## Construction 2: fixed DC and unbounded memory in eight modes

Number the physical modes and Walsh basis vectors by 0,...,7. Define

    H_ij = (-1)^(popcount(i AND j)) / sqrt(8).

Then H^T H=I, and h_j denotes its j-th column. Fix

    e=h0,   b=h2.

For 0<t<1 define the symmetric matrix A(t) in this basis by

    diagonal A = (0, t^2, 2, 2-t^2, 1, 1, 1, 1),
    A_12=A_21=t,
    A_47=A_74=-t,

with every other off-diagonal entry zero. Set

    C(t)=H A(t) H^T.

### Constraint proof

**Energy conservation:** A(t) has identically zero row and column 0, so
C(t)e=0.

**Positive semidefiniteness and unique invariant:** on e-perpendicular, the
nontrivial blocks are

    B(t) = [t^2  t; t  2],
    F(t) = [1   -t;-t  1],

and the scalars 2-t^2,1,1. The determinant of B is t^2>0 and its trace is
2+t^2>0; F has eigenvalues 1-t and 1+t. The scalars are positive for 0<t<1.
Thus C(t)>0 on e-perpendicular, with exactly one invariant.

**Fixed diagonal:** every diagonal contribution from A's diagonal is
trace(A)/8=1. The two off-diagonal pairs contribute

    (t/4) S_i1 S_i2 - (t/4) S_i4 S_i7 = 0,

where S=sqrt(8)H, because 1 XOR 2 = 4 XOR 7 = 3. Therefore every physical
diagonal entry of C(t) equals 1, for every t.

**Time reversal:** the physical permutation i->i XOR 3 commutes with C(t).
It fixes e and changes the sign of b. The coupled Walsh pairs have equal
parity under this permutation. Thus the example is not created by breaking
the elementary even-energy/odd-current symmetry.

**Fixed velocities:** since e_i=1/sqrt(8), choosing V_ii=b_i/e_i=+/-1 gives
the same V and b=V e for all t. No change in current weights or velocities
creates the effect.

### Exact response and memory

Because the source b=h2 is supported only in B(t),

    C(t)^+ b = -h1/t + h2.

Therefore, exactly,

    K0=1,
    M2=1+t^(-2),
    tau_mem=1+t^(-2).

The full response is

    K_t(p) = (p+t^2) / [p^2+(2+t^2)p+t^2].

This proves unbounded memory at fixed dimension, fixed diagonals, fixed energy,
fixed current and fixed DC within the stated broad operator class.

The two eigenvalues of B are

    lambda_+/- = [2+t^2 +/- sqrt(4+t^4)]/2.

For small t,

    lambda_- = t^2/2 - t^4/8 + O(t^8),
    w_- = |q_-^T b|^2 = t^2/4 + O(t^4).

The slow mode's vanishing current overlap is compensated by its vanishing
decay rate: w_-/lambda_- tends to 1/2, but w_-/lambda_-^2 diverges as t^(-2).
This is the mechanism for constant DC and divergent memory.

**Singular limit:** for fixed p>0, K_t(p)->1/(p+2) as t->0. For every t>0,
K_t(0)=1. At t=0 an extra invariant appears and the pseudoinverse DC is 1/2.
The zero-frequency limit and the rank-changing collision limit do not commute.
No assertion of uniform Taylor expansion in p near this limit is made.

### Extension: preserve a complete three-current DC tensor

Take the three fixed current vectors

    J=(h2,h5,h6).

They are mutually orthonormal, all odd under the same reversal permutation,
and generated by fixed diagonal velocities V_a e=J_a. Their entries describe
the eight sign combinations of three unit velocity components. Directly,

    J^T C(t)^+ J = I_3,
    J^T (C(t)^+)^2 J = diag(1+t^(-2),1,1).

Thus knowing all entries of the uniform DC tensor need not fix dynamic memory
in this operator class. This construction does not impose the full AlN point
group; its unit DC tensor is not a claim about AlN anisotropy.

### Physical realizability limitation

C(t) has positive off-diagonal entries for t>0, so it is not a conventional
nonnegative-weight graph Laplacian in the physical-mode basis. Symmetric PSD
entropy-coordinate phonon operators need not satisfy the graph sign condition,
but PSD, energy conservation and time reversal alone do not establish
realizability by any particular set of energy/momentum-resonant phonon events.
No frequencies, matrix elements, momentum selection rules, or AlN force
constants realizing this family were constructed. The unbounded-memory result
must be stated at the broad operator level.

## Why the ambiguity is not restricted to equal rates

Let all e_i be nonzero and write E=diag(e_i). Any symmetric four-by-four
matrix F with zero diagonal and F 1=0 generates

    Delta=E^(-1) F E^(-1),
    Delta e=0,  diag(Delta)=0.

An explicit two-parameter family is

    F = [ 0   x    y   -x-y ]
        [ x   0   -x-y  y  ]
        [ y  -x-y  0    x  ]
        [-x-y y    x    0  ].

It may be supported on any four nonzero-energy components of a larger system.
If a starting C0 is strictly positive on e-perpendicular, then C0+s Delta
has the same prescribed diagonal and energy invariant and remains positive on
that space for sufficiently small |s|. A sufficient condition is

    |s| ||Delta||_2 < lambda_min(C0 restricted to e-perpendicular).

For x0=C0^+ b,

    d K0(C0+s Delta)/ds at s=0 = -x0^T Delta x0.

This derivative need not vanish. Thus the ambiguity is an affine structural
freedom, not just a trick using equal rates. However this argument does not
prove that every given rate/current data set has multiple admissible
microscopically realizable operators, or that every direction changes every
observable. A known interior feasible operator and nonzero derivative are
required for that local claim.

For n>=3 and positive e, the linear constraints diag(Delta)=0 and Delta e=0
leave dimension n(n-3)/2. To see independence, use the invertible congruence
above: the row-sum map on off-diagonal symmetric matrices is the unsigned
incidence map of the complete graph. A dependence would require z_i+z_j=0
for every i!=j, whose only solution for n>=3 is z=0. There are n(n-1)/2
off-diagonal variables and n independent constraints.

## Computation and quantitative validation

Artifacts owned by this branch:

- `experiments/constructive_collision.py`
- `experiments/constructive_collision_results.json`

Reproduction from a PowerShell prompt:

    py C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260916-210742-aln-spectral-closure\experiments\constructive_collision.py

Environment used: Python 3.14.7, NumPy 2.5.3, SymPy 1.14.0, mpmath 1.3.0.
The JSON records actual versions and all diagnostics. The calculations use:

1. Hand block/eigenvector derivation as written above.
2. Exact symbolic validation of the physical-matrix diagonal, energy nullspace,
   parity, response formula and rational instances.
3. Dense physical-coordinate solves of (C+ee^T)x=b. Adding ee^T only lifts
   the known energy nullspace; it does not change the desired dissipative solve.
4. Independent eigendecomposition on e-perpendicular, calculating the spectral
   moments from current overlaps.
5. Physical-coordinate solves at 80 decimal digits for t=10^-2,10^-6,10^-10,
   10^-20.

No discretization error is present. Conditioning, invariant residuals, exact
moment discrepancies and the small-gap limit are the numerical concerns.

### Results

- Four-mode equal-DC pair: both DC values are 1.5 in double precision; the
  memories are 0.75 and 0.8333333333333335. Energy residuals are <=1.2e-16.
  Symbolic rational arithmetic obtains the exact values in the table above.
- Eight-mode family, t>=10^-3: dense DC absolute discrepancy <=1.84e-11 and
  memory relative discrepancy <=5.50e-11. Physical diagonal errors <=1.2e-16.
- Across t=0.75 down to 10^-8, the resolvent at p=0.3 agrees with the exact
  formula to <=2.23e-16 in absolute value. A well-resolved finite-frequency
  response therefore does not validate the DC slope in the near-singular regime.
- At t=10^-8, double-precision dense inversion reports DC 0.718356... instead
  of 1, and a memory 73.45% below the exact value. The reduced eigensolver
  reports DC 1.697852... and a memory 238.04% above the exact value. Both are
  untrustworthy: the exact slow eigenvalue is approximately 5e-17.
- At t=10^-20, the independent 80-digit physical-coordinate solve gives DC 1
  and memory 1+10^40 with DC absolute discrepancy 1.09e-42, memory relative
  discrepancy 3.27e-42 and relative solve residual 6.24e-62.

The divergence is exact, not an inference from unstable double-precision
eigenvalues. The unstable calculations are retained as a useful failure case.

## Failed or limited approaches

1. **Three modes with fixed diagonals and one nonzero energy vector:** there
   are no free symmetric off-diagonal directions preserving the energy
   invariant. A four-mode or larger construction is necessary for these
   constraints. This is a dimension statement, not a claim about all models.
2. **One fixed current eigenvector:** changing its eigenvalue changes DC and
   memory together, so it cannot demonstrate equal-DC/different-memory by
   itself. Two current-visible eigenmodes are used for the graph pair.
3. **Fixed-eigenvector, fixed-current weights:** finite K0 bounds every
   current-visible rate away from zero by lambda_j>=w_j/K0. Hence this form
   cannot create an arbitrarily slow mode with nonzero fixed overlap while
   keeping DC fixed. The eight-mode family instead makes the slow eigenvector's
   current overlap vanish along with its rate.
4. **Naive double precision near t=0:** small residuals and accurate response
   at p=0.3 coexist with very inaccurate DC moments. Precision must resolve
   the small eigenvalue, and exact formulas/independent high precision are
   needed here.
5. **Symbolic structural equality without simplification:** the first script
   run stopped at an unsimplified zero row-sum expression. Applying exact
   simplification made the identity explicit; this was a CAS comparison
   issue, not a failed conservation law.
6. **Physical realization:** no permitted phonon-event network for actual AlN
   was constructed. This remains an essential limit on extrapolating the
   strongest example from an operator theorem to a material-specific claim.

## Consequences for the AlN closure question

The supplied modewise diagonal rates fix the diagonal RTA surrogate. Even
adding a measured or iterative DC tensor does not, in the broad conserving
class, identify the off-diagonal collision couplings or time-response slope.
In particular, diagonal RTA second moments are properties of that surrogate,
not determined exact memory parameters of the full collision operator.

For a controlled hydrodynamic closure one needs information on current and
streaming couplings to the relevant collision eigenmodes, their rates, and an
error estimate for the eliminated subspace. Neither a positive minimum
diagonal rate nor agreement of DC conductivities supplies a collision gap.

The next discriminating computation should use an actual collision operator
or reproducible operator-vector action. It should compare DC and first
frequency derivative, inspect current overlaps with the slow eigenspace, and
test streaming into discarded modes. If only rates and a DC tensor are
available, obtain at least independent low-frequency response data or the
relevant projected resolvent; these constrain memory but still do not
automatically identify every hidden mode.

An independent proof auditor should check the fixed-diagonal cancellation,
the physical-coordinate constraints, the singular limiting argument, and the
scope distinction. A numerical analyst should inspect the stored near-singular
failure regime rather than infer reliability from the well-conditioned cases.
