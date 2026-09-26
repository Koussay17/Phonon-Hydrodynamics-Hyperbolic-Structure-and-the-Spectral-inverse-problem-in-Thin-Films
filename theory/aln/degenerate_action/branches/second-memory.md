# Second generation — Exact memory of retained block variables

2026-09-25. Cross-examination of A, C-lumpability, and D-phase-information. **Finite derivation under stated assumptions; awaiting audit. No material or novelty claim.**

Let dot(Y)=-C Y on a finite real or complex Hilbert space, with C=C*>=0. This includes a declared entropy-coordinate population generator or a separately specified Hermitian-matrix generator. Let U,Q be orthonormal bases of retained and hidden orthogonal subspaces, and write

    a=U*Y, h=Q*Y,
    A=U*CU, K=U*CQ, D=Q*CQ.

Thus a'=-Aa-Kh and h'=-K*a-Dh. Variation of constants gives the **exact** equation

    a'(t)=-A a(t)+integral_0^t M(t-s)a(s)ds-K exp(-Dt)h(0),
    M(t)=K exp(-Dt)K*.

The memory has a plus sign; the hidden-initial-data term has a minus sign. No equilibrium initialization, gap, or Markov approximation was used. Finite linear ODEs give unique global solutions. The positive feedback term is compatible with full entropy contraction because this is exact elimination of a contractive system.

For h(0)=0, compression a_c'= -A a_c matches the initial derivative, but

    a''(0)=(A^2+KK*)a(0),
    a(t)-a_c(t)=(t^2/2)KK*a(0)+O(t^3).

Hence initially equal block occupations do not suffice; this second-order defect vanishes for a selected state if K*a(0)=0. Closure for every initial state is equivalent to K=0.

## Resolvent and zero modes

For Re(z)>0 all Laplace transforms and both full/hidden resolvents exist, and

    [zI+A-K(zI+D)^(-1)K*] a_hat
       =a(0)-K(zI+D)^(-1)h(0).

No substitution at a pole is justified. Although D may be singular, positivity of C implies K ker(D)=0: a vector (0,v) with Dv=0 has zero C-quadratic form and therefore lies in ker(C). Hidden zero modes consequently contribute neither memory nor the initial-data term. In finite dimension the static Schur complement exists as

    A-K D^+ K* >= 0,

but may be singular because of full-system conserved modes. This does not assert a DC inverse or a continuum limit. If a conserved vector e satisfies Ce=0 and lies in the retained space, then A e_ret=0 and K* e_ret=0; its retained charge is preserved exactly.

## Concrete network from branch C

Use its unequal coefficients Gamma=(3,1), Bose variances W=(6/25,6/25,2,3/4), and event weights (9/5,3/5). Retain the normalized parent sum and the two daughters in entropy coordinates; hide the normalized parent difference. Direct substitution gives

    v=(-5/(2 sqrt(3)), 1/sqrt(2), 2/sqrt(3))^T,
    A=(12/5)vv^T, K=-sqrt(3)v, D=5,
    M(t)=3 exp(-5t)vv^T.

For Gamma=(2,2), A and D are unchanged but K=0. Thus identical compressed coefficients have respectively nonzero exponential memory and no memory. Uniform parent initialization sets h(0)=0, removing only the initial-data term; it leaves the convolution intact. These are algebraic identities, not new numerical results.

## Alternative formulations worth retaining

- **Conditional expectation/geometry:** the center projection identifies invariant block observables; it need not reduce C. Enlarging the retained space to span{C^j Ran(U)} gives an exact invariant realization, possibly the full space.
- **Variational projection:** useful for certified response bounds. For real p>0, Schur ordering gives U*(pI+C)^(-1)U >= (pI+A)^(-1). This bounds retained forcing/response, not arbitrary hidden initial data.
- **Spectral/dual:** M(t)=integral exp(-lambda t)dSigma(lambda), with dSigma=K dE_D K*, is a positive matrix measure. It identifies the missing information as coupling-weighted hidden spectral data, rather than another block norm.
- **Conserved charges:** retain them explicitly to prevent artificial relaxation. They do not imply closure of every block total.

Branch D's phase ambiguity can change the specified coherent generator and therefore this hidden spectral data. The memory formula requires that generator first; neither squared amplitudes nor invariant block sums determine it automatically.
