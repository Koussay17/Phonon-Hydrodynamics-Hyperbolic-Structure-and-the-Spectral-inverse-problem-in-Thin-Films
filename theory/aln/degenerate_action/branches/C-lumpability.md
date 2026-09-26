# C: reversible-event lumpability and weighted block projection

Independent first pass complete, 2026-09-25. Read only this campaign's question and assumptions. No other branch/PI result or discussion was read. Scope: a declared finite scalar-population event model, not an unspecified quantum generator.

## 1. Exact closure criterion

Let x be the population perturbation on an equal-weight mode grid:
\[
 \dot x=-Lx,\quad L=MW^{-1},\quad
 W_{ii}=n_i^0(1+n_i^0)>0,\quad
 M=\sum_e\lambda_e\nu_e\nu_e^T .
\]
The signed event incidence is nu_e; lambda_e includes its equilibrium Bose factor. Exact resonance and the stipulated reversible mean-population closure are assumptions.

For full-row-rank B, define observed block variables N=Bx. Autonomous linear dynamics for every initial x exist **if and only if**
\[
 BL=GB
\]
for some G. Necessity follows by differentiating at t=0; sufficiency follows by substitution. Equivalently, ker B must be L-invariant. This is the algebraic strong-lumpability criterion. The reaction generator need not have Markov transition signs, so Markov-chain theorems requiring those signs cannot be imported automatically.

Set
\[
 \Sigma=BWB^T,\quad H=WB^T\Sigma^{-1},\quad P=HB.
\]
Then BH=I and the only candidate reduced generator is
\[
 G=BLH=BMB^T\Sigma^{-1}.
\]
**Compression is not closure.** The exact condition is BL(I-P)=0. Since L is self-adjoint in the W^{-1} metric, this is equivalent to (I-P)LH=0. In entropy coordinates, both say the off-diagonal block of the symmetric matrix C=W^{-1/2}MW^{-1/2} vanishes between the retained space E=W^{1/2}B^T Sigma^{-1/2} and its orthogonal complement.

Thus, in this reversible class, closure for arbitrary microscopic states and preservation of the reconstructed manifold are equivalent. That equivalence need not hold for a general nonreversible linear generator.

## 2. Means are not totals

For disjoint block indicators and an exactly degenerate block A of size m_A, its common Bose variance is w_A. Then
\[
 \Sigma_{AA}=m_Aw_A,\qquad (HN)_i=N_A/m_A\quad(i\in A).
\]
Equal per-mode occupations correspond to dividing the block total by m_A. Mean perturbations a=D_m^{-1}N evolve, when closure holds, with generator D_m^{-1}GD_m, not G unchanged. The block susceptibility is m_A n_A^0(1+n_A^0), not N_A^0(1+N_A^0).

Even if x(0)=HN(0), compression always matches the initial derivative but can fail subsequently:
\[
 BL^2H-G^2=BL(I-P)LH.
\]
Equal mean occupations imposed initially therefore do not establish a preserved manifold.

## 3. Four-mode counterexample

Let p1,p2 be one degenerate parent block and a,b singleton daughters:
\[
 \epsilon=(\log6,\log6,\log2,\log3),\quad
 n^0=(1/5,1/5,1,1/2),
\]
\[
 W=\operatorname{diag}(6/25,6/25,2,3/4).
\]
The distinct-index events p1 <-> a+b and p2 <-> a+b have incidences
nu_1=(-1,0,1,1), nu_2=(0,-1,1,1).
Their common equilibrium factor is F=3/5. Set lambda_i=F Gamma_i.

Compare Gamma=(2,2) and Gamma=(3,1). Both have positive coefficients, identical total strength 4, and the same compressed block matrix. With B summing the parents and retaining each daughter, u=(-1,1,1),
\[
 BMB^T=(12/5)uu^T,\qquad
 G=\begin{pmatrix}
 5&-6/5&-16/5\\
 -5&6/5&16/5\\
 -5&6/5&16/5
 \end{pmatrix}.
\]

The equal-rate model satisfies BL=GB. The unequal-rate model does not: for h=(1,-1,0,0),
\[
 Bh=0,\qquad -BLh=5u=(-5,5,5).
\]
Hence x=0 and x=eta h have the same initial block totals but different block derivatives. For 0<eta<1/5 all perturbed occupations remain positive. Energy conservation holds in both models.

Equal initial parent occupations do not repair closure. In the unequal-rate model,
\[
 (I-P)LH=h(5/4,-3/10,-4/5),
\]
so an unresolved imbalance is generated immediately.

## 4. Relation to invariant sums

The rate allocations above can be assigned to squared amplitude vectors (sqrt2,sqrt2) and (sqrt3,1). A rotation by minus 15 degrees inside the degenerate parent block relates the vectors. Their squared norm is 4. Thus a basis-invariant block sum of squared couplings does not establish lumpability of the corresponding scalar population model.

More generally, Frobenius sums over a tensor block are invariant under unitary changes within its mode subspaces. If Bose factors and channel/counting prefactors are common within a block tuple, those sums can determine BMB^T and the compressed G. They do not force BL(I-P)=0.

This does not establish basis-dependent predictions for a correct quantum density-matrix generator. Diagonal populations in different degenerate bases are not related by rotating their diagonal entries alone; coherences matter. No full quantum generator has been specified here.

Also distinguish basis invariance from conservation. A conserved observable b^T x satisfies b^T L=0 and is trivially closed. Basis-invariant block totals or tensor norms are not automatically conserved. Conservation of total energy coexists with failure of closure for individual block populations in this example.

## 5. Status and importer implication

DERIVED UNDER ASSUMPTIONS: the intertwining criterion, weighted reconstruction, means-versus-totals normalization, and the obstruction from the explicit event network. Independent audit remains required.

Bounded computational check completed: experiments/C_lumpability_check.py and experiments/C_lumpability_results.json. All 13 checks passed (Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0). Exact rational algebra checks the intertwining defect, reconstructed-manifold leakage, second derivative, conservation, and a positive nonlinear population example. The amplitude rotation is also checked symbolically.

Direct matrix exponentials, starting from equal parent perturbations x(0)=(0.01,0.01,0,0), give:

| Model | Maximum block-propagator discrepancy from exp(-Gt), t=0.001,0.01,0.1,0.5 |
|---|---:|
| Equal rates | 1.66e-16 |
| Unequal rates | 0.0422073 |

For the unequal-rate model the selected initial-state block error is 3.04% at t=0.1 and 4.43% at t=0.5. The failure is not inferred solely from numerical propagation: the exact nonzero intertwining defect already establishes it. No timestep or spatial discretization is used.

Reproduce: py -B D:\ResearchLab\orchestration\campaigns\20260925-113533-aln-degenerate-collision-action\experiments\C_lumpability_check.py

No material calculation or novelty claim is made.

Importer implication: retain block multiplicities and mode-resolved incidence/weights long enough to test BL(I-P), or its entropy-coordinate counterpart. If only invariant block sums survive, label the result a projected initial-derivative model unless a separate symmetry/lumpability condition is established. Actual degenerate quantum transport may require density-matrix data and a specified generator.


