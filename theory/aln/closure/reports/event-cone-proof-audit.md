# Proof audit of fixed-event-cone lemmas R2 and R3

Completed: 2026-09-17. This is controlled cross-examination after the four independent red-team reviews. Only the new event-cone R2/R3 claims and red-peer.md's independently derived event-library inequality are audited here. No novelty claim is made.

## Verdict

**R3 is VALID under its stated fixed finite-dimensional, fixed-dictionary, nonnegative-weight and positive-definite assumptions.** The adjugate identity and the resulting norm/memory bound have a complete exact proof below. Zero coefficients, dependent event subsets and vanishing current overlaps cause no gap.

**R2's substantive obstruction is VALID.** Its displayed reciprocal 1/R needs the trivial R=0 case separated: then no positive t is feasible. The independent red-peer inequality is the same valid termwise estimate.

Consequently, with a fixed finite event dictionary, fixed Euclidean/entropy inner product, fixed dissipative space and fixed nonzero current, **a uniform finite bound on DC response implies a uniform finite bound on the first response-moment time**, even if the positive collision gap tends to zero. A fixed diagonal is unnecessary. This strengthens the single-path obstruction; it does not provide a mesh-independent or numerically useful actual-AlN bound.

## 1. Exact statement and domains

Fix d>=1, finitely many vectors a_1,...,a_m in a fixed real d-dimensional inner-product space, and b!=0. Use an orthonormal basis. For every nonnegative weight vector g for which

    A(g)=sum_(a=1)^m g_a a_a a_a^T

is strictly positive definite, define

    x=A^(-1)b,  K=b^T A^(-1)b>0,  tau=||x||^2/K.

For each (d-1)-element subset I, order its columns by increasing event index and define v_I by

    det([a_I,y])=v_I^T y for every y.

When d=1, I is the empty subset, its weight product is 1 and v_I=1 in the chosen orientation. Set

    c_I=product_(a in I) g_a,
    M=max_(I:v_I^T b != 0) ||v_I||/|v_I^T b|.

The maximum is over the whole fixed dictionary, not a weight-dependent selection. Its index set is nonempty whenever any positive-definite A is feasible, as proved below. It is finite because there are finitely many finite nonzero ratios.

## 2. Audited adjugate proof

For s>=0 let Z_s be the d by (m+1) matrix with columns sqrt(g_a)a_a and sqrt(s)y. Then

    Z_s Z_s^T=A+s y y^T.

Cauchy-Binet gives its determinant as the sum of squares of all d-column minors. Terms not containing the last column sum to det(A). Terms containing it have the form

    s c_I det([a_I,y])^2 = s c_I (v_I^T y)^2.

There is no missing factorial: every subset occurs once, and the ordering sign disappears upon squaring. Thus

    det(A+s y y^T)=det(A)+s sum_I c_I(v_I^T y)^2.

Independently, determinant multilinearity in the columns gives

    det(A+s y y^T)=det(A)+s y^T adj(A)y.

Terms with two replaced columns vanish because those columns are proportional to y. Comparing the coefficient of s gives equality of these two quadratic forms for every y. Both matrices are real symmetric, so polarization proves

    adj(A)=sum_I c_I v_I v_I^T.                         (E1)

This polynomial identity does not itself require A to be invertible. Invertibility is required in the next step, not silently assumed during Cauchy-Binet.

## 3. Norm bound, positivity and zero terms

Because A>0, det(A)>0 and adj(A)=det(A)A^(-1)>0. Hence

    D_b := b^T adj(A)b
         = sum_I c_I (v_I^T b)^2
         = det(A) K > 0.                              (E2)

In particular, at least one subset has c_I>0 and v_I^T b!=0. This proves that M is defined. For each such subset put

    alpha_I=c_I(v_I^T b)^2/D_b.

The alpha_I are nonnegative and sum to one. Dividing (E1)b by D_b gives the exact convex-combination identity

    x/K = sum_(I:v_I^T b != 0) alpha_I v_I/(v_I^T b).   (E3)

Therefore

    ||x||/K <= sum_I alpha_I ||v_I||/|v_I^T b| <= M,
    ||A^(-1)b|| <= M b^T A^(-1)b,
    tau <= M^2 K.                                    (E4)

This is a complete proof. Equivalently one may use the branch's termwise triangle inequality; (E3) makes the normalization and nonnegative weights explicit.

Zero and rank cases:

- A subset of rank below d-1 has v_I=0 and contributes nothing.
- If any g_a in I vanishes, c_I=0 and that subset contributes nothing.
- If v_I!=0 but v_I^T b=0, the subset contributes zero both to adj(A)b and to D_b. Excluding it from M is legitimate; no vector contribution is lost.
- The active event vectors must span the full fixed space, since A>0. Some individual weights may nevertheless be zero.
- For d=1, M=1/|b| and (E4) is equality. No exceptional dimension invalidates the formula.
- There is no assumed lower bound on det(A) or the collision gap. They may approach zero along a family; the normalized convex combination still gives a bound independent of g.
- At a singular endpoint A^(-1) is undefined. A pseudoinverse version needs a restriction to a fixed range and the appropriate current range condition; that extension is not used here.

## 4. R2 and the independent event-library comparator

For fixed event vectors z_a and fixed vectors u,b, write U_a=u^T z_a, V_a=b^T z_a. If some U_a!=0, define

    R=max_(a:U_a!=0) |V_a/U_a|.

Term by term, nonnegative g_a gives

    |u^T C b| <= sum_a g_a |U_a V_a|
                <= R sum_a g_a U_a^2 = R u^T C u.

If every U_a=0, both quadratic expressions vanish. If R=0, the cross term vanishes even if the diagonal term does not. D's required positive cross term t is impossible in either case. If R>0, inserting u=h1, b=h2 and the D block gives

    t <= R t^2, hence t>=1/R for positive t.

Therefore no one fixed finite dictionary can represent D for a sequence t->0. This is exactly the independently obtained inequality in red-peer.md. Agreement between the reports is not the proof; the displayed termwise estimate is.

R2 alone rejects that path. R3 additionally rejects **every** positive-definite path in the same fixed finite cone with bounded K and divergent tau. Any earlier statement that the R2 argument alone does not establish a general memory bound remains correct; it is superseded by the additional R3 proof for this precise class.

## 5. Sharpness and failure of a dictionary-independent bound

The fixed-dictionary condition is essential, even at dimension two with just two events. Let epsilon>0,

    a_1(epsilon)=(epsilon,1)^T,  a_2=(0,1)^T,
    g_1=g_2=1,  b=(0,1)^T.

Then

    A=[[epsilon^2,epsilon],[epsilon,2]],
    A^(-1)b=(-1/epsilon,1)^T,
    K=1,  tau=1+epsilon^(-2).

For the one-element subsets, the relevant cross vector is v_1=(-1,epsilon)^T; v_2=(-1,0)^T has zero overlap with b and drops out. Thus

    M=sqrt(1+epsilon^(-2)),  tau=M^2 K.

The bound is sharp here, but M diverges because an event direction changes. This is also an explicit warning against treating M as a function of dimension or diagonal rates alone. The example does not violate R3, whose dictionary is fixed.

Under q-mesh refinement the dimension, event directions, entropy weights, invariant restriction and projected source may all change. Finiteness of each individual M gives no uniform bound on that sequence. An infinite event dictionary can likewise have unbounded ratios. No continuum memory estimate follows without an additional uniform geometric argument.

## 6. Scope and evidence status

- **PROVED under the exact hypotheses:** (E1)-(E4), the fixed-cone bounded-DC/bounded-memory consequence, and R2 with its R=0 case treated.
- **NOT PROVED:** a practical numerical M for AlN; a constant uniform under physical mesh refinement; finite continuum memory; uniqueness of memory or DC from diagonal data; any novelty claim.
- The proof assumes the event outer-product representation. Its microscopic validity for a particular discretization, broadening procedure or reduced effective operator must be checked separately.
- The theorem remains valid under extra constraints that merely restrict nonnegative g within the same fixed cone. It does not authorize changing the event geometry while calling the cone fixed.

Exact symbolic verification is supplementary. Execution details and zero/degenerate-subset checks are recorded in the addendum below.
