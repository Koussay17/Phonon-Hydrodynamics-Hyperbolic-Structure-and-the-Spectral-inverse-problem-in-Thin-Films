# Independent red-team proof audit — completed 2026-09-26

Scope: frozen candidate-for-review.md; A-operator-algebra.md; C-lumpability.md; D-phase-information.md; second-memory.md. No other red report or second-memory-check report was read. No material calculation was run. The algebra and exact checks were completed on 25 September; a usage-limit rejection interrupted only this final save.

## Verdict

**The central finite-algebra claims are VALID under their stated assumptions.** No incorrect memory sign, missing block multiplicity, or counterexample to the declared closure criterion was found. Independent exact SymPy 1.14.0 checks passed for C's network, reconstruction, curvature, memory blocks and Schur resolvent, and D's Hadamard examples. These support the explicit arguments below; they are not formal verification.

**Minor quantifier qualification:** D proves non-injectivity of squared-amplitude data. Its examples alone do not establish an almost-everywhere genericity theorem for every tensor format. A universal reading would be false: a 1x1x1 tensor's modulus identifies its local-unitary class. The frozen word "generic" must not be promoted to that stronger quantified theorem.

Candidate item 3's AlN amplitude/transport residuals were not recomputed in this bounded proof audit. They remain separate numerical evidence. No material ambiguity theorem or unique density-matrix generator follows from the finite examples; the candidate expressly disclaims these extensions.

## 1. D: phase information — VALID

K+ and K- have identical entrywise squares but Gram spectra {0,4} and {2,2}. Independent left/right unitaries conjugate the parent Gram matrix and preserve its spectrum. The tensors therefore are inequivalent under the stated group. Their Hadamard-rotated squares are exactly [[2,0],[2,0]] and [[2,0],[0,2]], proving the separate obstruction to reconstructing rotated squares.

The continuous family's Gram eigenvalues are 2 +/- |1+exp(i theta)|, matching the displayed formula. In the 2x2x2 examples each flattened row has squared norm 4 and their inner product is 4 or 2, giving {0,8} or {2,6}. Numerical permutation symmetry is correctly distinguished from physical equality of all leg energies.

D writes dot(x)=Lx with negative dissipative L, whereas C writes dot(x)=-Lx. This is a convention difference, not a sign error. The declared coherent covariance model is internally positive: E N(0) E+n0(I-E^2)>=0 for N(0)>=0 and G>=0. Its traces have slopes -2 and second derivatives 8 versus 4 as stated. That conclusion belongs to this chosen model.

## 2. C: weighted compression and closure — VALID

For full-row-rank B and W>0, Sigma=BWB^T>0 and BH=I. Autonomous block dynamics for every initial x imply BL=GB by differentiating at zero; conversely this identity gives closure. The unique candidate is G=BLH, and the exact residual condition is BL(I-HB)=0.

P=HB is orthogonal in the W^-1 metric. L=MW^-1 is self-adjoint in the same metric when M=M^T. Taking the metric adjoint makes BL(I-P)=0 equivalent to (I-P)LH=0. Reversibility is essential to that equivalence and is explicitly assumed.

For disjoint blocks with common variance w_A, Sigma_AA=m_A w_A and H distributes a total as N_A/m_A. The stated similarity transform for means is correct; no Bose variance of a summed occupation is substituted.

The exact network is resonant because log(6)=log(2)+log(3), with equilibrium flux factor 3/5. Independent rational calculations reproduce

    G=[[5,-6/5,-16/5],[-5,6/5,16/5],[-5,6/5,16/5]],
    -BL(1,-1,0,0)^T=(-5,5,5)^T,
    (I-P)LH=(1,-1,0,0)^T(5/4,-3/10,-4/5).

Equal rates satisfy BL=GB; unequal rates do not. A sufficiently small hidden perturbation preserves positive occupations, so the counterexample does not rely on inadmissible states. BL^2H-G^2=BL(I-P)LH is exact. Equal initial parent occupations consequently do not establish closure.

## 3. A: observables and symmetry — VALID on the declared domain

A unitary normalizing the diagonal algebra permutes its minimal rank-one projections, hence is monomial. The Hadamard example correctly shows that squared-modulus population maps are not a group representation.

On V=direct_sum Herm(H_alpha), invariant linear functionals are represented by block scalars; group averaging gives Tr(X_alpha)I/d_alpha. This classifies functionals on V, not arbitrary full-space observable matrices: off-block components are invisible on the declared domain.

For self-adjoint C, E C(I-E)=0 is equivalent to [E,C]=0. The spectator-parent example has zero block sums but the stated nonzero derivative; positivity and conservation do not imply closure. Common Bose/counting prefactors make projected event vectors identical, allowing the Frobenius sum to determine compression. Repeated blocks are explicitly excluded. Covariance of amplitudes and operator together is correctly distinguished from symmetry of a fixed generator.

## 4. Memory and zero modes — VALID

Solving h'=-K*a-Dh and substituting into a'=-Aa-Kh gives the **positive** convolution and **negative** hidden-initial-data term. For h(0)=0, a''(0)=(A^2+KK*)a(0). A selected state's vanishing second-order defect does not establish closure for every state; that requires K=0.

The Laplace Schur term is **minus** K(zI+D)^-1K*. Re z>0 ensures the required inverses. If Dv=0, positivity implies C(0,v)=0 and hence Kv=0. Hidden zero modes therefore contribute neither memory nor initial forcing. Completing the square with D+ gives A-KD+K*>=0, without implying invertibility.

For C's network, the exact audit gives A=(12/5)vv^T, K=-sqrt(3)v and D=5. Equal rates retain A,D but give K=0. Direct symbolic inversion verifies the Schur identity and KK* curvature. The real-p response ordering follows from inverse monotonicity of positive matrices. No gap or equilibrium initialization was silently required for exact elimination.

**Disposition:** accept the finite obstruction, weighted closure criterion and exact memory identities in their declared domains. Reject extensions from invariant compression to autonomous dynamics, covariance to fixed-generator symmetry, or these examples to a unique material quantum operator. Apart from the genericity qualification and unaudited material residuals, no unresolved central finite-algebra proof gap was identified. Only this report was written during finalization.
