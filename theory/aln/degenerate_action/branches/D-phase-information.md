# D — phase information and explicit finite tensors

**Status: independent first pass complete, 2026-09-25.** Read only the campaign question and assumptions; no peer branches, PI experiments or AlN arrays. These are synthetic obstruction examples, not material counterexamples. Exact finite-dimensional algebra supplies the reference; NumPy checks it independently.

## 1. Identical squared entries can conceal inequivalent tensors

Let K map a two-dimensional parent space to a two-dimensional daughter space, with a singleton third leg: V[a,b,0] = K[b,a]. Parent frequencies are 2 and daughter frequencies are 1 and 1, in units hbar = k_B T = 1. The legs occupy distinct momentum or species sectors; equal daughter frequencies do not merge their mode spaces. Thus all frequencies and temperature are positive and all reactions are resonant. The interaction can include its Hermitian conjugate; no AlN force-constant realizability is asserted.

Take

\[
K_+=\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad
K_-=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
\]

Both have identical P = |K|², total squared norm 4, and diagonal parent Gram entries (2,2). Nevertheless G = K†K has spectra (0,4) and (2,2). Under independent basis changes K' = U_b† K U_a, G' = U_a† G U_a; hence these tensors cannot be related by such changes. Projectors, frequencies, P and even basis-specific diagonal rates therefore do not identify the interaction's unitary equivalence class.

This is a continuous ambiguity: replacing the bottom-right entry by exp(i theta) gives the same P and Gram eigenvalues 2 ± 2|cos(theta/2)|. Four checks from theta = 0 to pi match this exact formula within 4.44e-16. The inverse problem is non-injective, rather than an error repairable by higher precision. A zero Gram eigenvalue is a dark coupling direction, not a zero phonon frequency.

The obstruction also survives three nonsingleton legs and permutation-symmetric entries. For 2 x 2 x 2 tensors, compare all entries +1 with all entries +1 except V[1,1,1] = -1. Both have P identically 1 and squared norm 8. Every one-leg Gram spectrum is respectively (0,8) or (2,6). This is algebraic permutation symmetry, not a claim that three resonant phonon legs all have the same frequency.

## 2. Rotated squared amplitudes require interference information

For the parent Hadamard matrix H = [[1,1],[1,-1]] / sqrt(2),

\[
|K_+H|^2=\begin{pmatrix}2&0\\2&0\end{pmatrix},\qquad
|K_-H|^2=\begin{pmatrix}2&0\\0&2\end{pmatrix}.
\]

The same initial P and the same specified rotation produce different outputs. Squaring the rotation coefficients and multiplying P omits cross terms; its Frobenius error here is 2 for either case. All complete block sums remain 4. Gram covariance residuals are below 9e-17. Thus agreement of block norms is necessary under these unitary changes but is not sufficient to identify tensor equivalence. The examples use signs only: restricting amplitudes to be real does not remove the obstruction.

## 3. Block sums determine a projection, not necessarily closure

**Specified population model.** For independent incoherent Bose reactions a_i <-> b_j + c, define mode order (a0,a1,b0,b1,c), stoichiometric vectors s_ij = -e_ai + e_bj + e_c, equilibrium n_s = 1/(exp(omega_s)-1), and S = diag(n_s(1+n_s)). Linearizing the resonant forward-minus-reverse Bose flux gives

\[
\dot{\delta n}=L\delta n,\qquad
L=-\sum_{ij}P_{ji}s_{ij}s_{ij}^{T}S^{-1}.
\]

The common positive equilibrium flux factor has been absorbed into time. This assumes the diagonal, incoherent population approximation in the chosen basis; it is not a full density-matrix kinetic model.

Let B sum each of the three blocks and R distribute a block total equally among its modes, so BR = I. Autonomous closure for **every** perturbation exists precisely when

\[
BL(I-RB)=0.
\]

Indeed, closure BL = MB implies M = BLR and the displayed condition; conversely that condition gives BL = (BLR)B. Complete block sums determine the projected action BLR in this example because the frequencies and equilibrium susceptibilities are constant within each block. They do not determine the missing condition.

For the two Hadamard-rotated tensors, BLR is identical in the numerical calculation. Yet:

| Test | K+ rotated | K- rotated |
|---|---:|---:|
| norm(BL(I-RB), 2) | 27.06388005147698 | 0 |
| BL h, h = (1,-1,0,0,0) | 22.09756552866904 (-1,1,1) | (0,0,0) |

Here Bh = 0. The nonzero result proves that identical instantaneous block totals can have different derivatives. Exactly, the coefficient is 4/[n(2)(1+n(2))]. Scale h arbitrarily small to remain in the positive-population linear regime. The unrotated P is identical for both tensors and passes this particular classical closure test. Dropping coherences in different bases defines different approximations; this is not a violation of physical basis invariance.

Both constructed Onsager matrices are sums of nonnegative rank-one forms. Energy conservation follows from omega·s = -2+1+1 = 0; the largest computed conservation residual is 3.67e-15. Tiny negative computed entropy-operator eigenvalues, no smaller than -3.71e-15, are roundoff around exact nonnegativity.

## 4. A specified coherent observable distinguishes the phases

One additional model makes the information loss tangible: let a parent covariance relax as dN/dt = -{G,N-n0 I}/2, with n0 = n(2) > 0. This is an explicitly chosen thermal covariance model, not a derivation of the full phonon generator. It preserves positivity: N(t) = E N(0) E + n0(I-E²), E = exp(-Gt/2).

For delta N(0) = epsilon I/2, normalized total deviations are (1+exp(-4t))/2 for K+ and exp(-2t) for K-. Initial slopes are both -2; second derivatives are 8 and 4. At t = 1 the totals are 0.5091578194443668 and 0.1353352832366127. The phase information determines a dark direction and subsequent block relaxation despite identical original P, norms and initial scalar rates. Whether this model describes a particular material remains a separate question.

## Information limit and importer consequence

P determines the specified incoherent population operator **in its recorded basis**, once multiplicities and statistical factors are supplied. Complete block norms determine certain projected or fully contracted observables. Neither determines general rotated P, Gram spectra, coherent dynamics or autonomous block closure. A single Gram matrix also need not suffice for every three-leg collision action.

An importer should preserve the eigenbasis with P, distinguish a projected block action from a closed one, and test closure for its stated kinetic model. Arbitrary basis rotations require complex amplitudes or explicitly sufficient cross products, with their conventions. Numerical classification: algebraic obstruction with roundoff-limited checks; no spatial, temporal or material convergence question was tested. Matrix exponentials used spectral evaluation, with no time-stepping error.

## Reproduction

Owned artifacts at campaign root: `D-phase-information.py` and `D-phase-information.json` (26,870 bytes). Python 3.14.7, NumPy 2.5.3; environment and script SHA-256 are saved in JSON. No input files or absolute material paths are required. Example PowerShell command:

```powershell
& 'D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe' 'D:\ResearchLab\orchestration\campaigns\20260925-113533-aln-degenerate-collision-action\D-phase-information.py' --output 'D:\ResearchLab\orchestration\campaigns\20260925-113533-aln-degenerate-collision-action\D-phase-information.json'
```

The output argument must be absolute. The finite constructions establish an obstruction within the stated assumptions; they do not establish AlN degeneracy, dephasing, coupling phases or transport behavior.

