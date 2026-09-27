# A2 — Arbitrary-gauge sewing of an oriented cubic channel

2026-09-27. Second-generation derivation after reading A–D and the corrected PI permutation/sewing experiment. **Conditional Hamiltonian identities, not a material rate or population-closure result.** This report changes no other artifact.

## 1. Canonical field and the creation term

Assume an ordinary harmonic crystal with real force constants, positive frequencies, and a Fourier coordinate frame in which D(-q)=D(q)*. Include any reciprocal-lattice/atomic-position phase transport before identifying a stored representative with -q. Let E_q have orthonormal polarization columns, Omega_q=diag(omega_q), and F_q=sqrt(hbar/2) Omega_q^(-1/2). Use either complete spaces or complete matched eigenspaces. Define

    E_-q = E_q* S_q,       S_q = E_q^T E_-q.

Then S_q is unitary, S_-q=S_q^T, and harmonic consistency requires

    Omega_q S_q = S_q Omega_-q.                         (1)

For mass-weighted displacement, with common volume factors suppressed,

    u(q) = E_q F_q a_q + E_-q* F_-q a_-q^dagger
         = E_q F_q X_q,
    X_q = a_q + S_q* a_-q^dagger.                       (2)

Here a^dagger denotes a column of creation operators; * means entrywise conjugation. The second equality uses (1). Thus the creation coefficient in component j is **conj(S_q[j,b])**, not S_q[j,b]. Also X_q^dagger=S_q X_-q (column convention); X components commute, using S_-q=S_q^T.

Reference normalization: [Togo, Chaput and Tanaka, Eqs. (5), (9), (10)](https://arxiv.org/html/1501.00691) give the paired-gauge displacement and an ordered cubic coefficient already containing 1/3!. Equation (2) is the explicit sewing extension derived here.

## 2. Forward/inverse formulas and their proof

Use H3=sum_ordered V_ijk X_0i X_1j X_2k, with V in energy units already containing 1/3!. Assume verified permutation symmetry; otherwise perform the actual coherent permutation sum. The incoming wavevectors obey q0+q1+q2=G. Set S_l=S_ql and let Vminus be the tensor at (-q0,-q1,-q2), in the corresponding saved bases. Real force constants and consistent oscillator normalization imply

    Vminus[a,b,c] = sum_ijk conj(V[i,j,k])
                              S0[i,a] S1[j,b] S2[k,c]. (3)

This also follows directly from X_q^dagger=S_q X_-q and Hermiticity.

For distinct physical daughter modes beta=(-q1,b), gamma=(-q2,c), the coefficient of a_beta^dagger a_gamma^dagger a_(q0,i) is

    g_forward[i,b,c] = 6 sum_jk V[i,j,k]
                                  conj(S1[j,b]) conj(S2[k,c]).   (4)

The coefficient of its adjoint a_(q0,i)^dagger a_beta a_gamma is obtained from Vminus. Its parent-slot creation coefficient is conj(S_-q0[a,i])=conj(S0[i,a]), hence

    g_inverse[i,b,c] = 6 sum_a Vminus[a,b,c] conj(S0[i,a]).       (5)

Substitute (3) into (5): sum_a S0[l,a]conj(S0[i,a])=delta_li. Therefore

    g_inverse[i,b,c] = conj(g_forward[i,b,c]).                   (6)

**Both proposed index/conjugation formulas are correct under these assumptions.** This is an off-shell algebraic identity; resonance is a separate requirement. Wavevector reversal supplies the tensor used in (5), but selecting its parent creation term produces the inverse reaction. Reversing all physical labels instead produces the reciprocal event, which is a different operation.

For identical complete daughter labels beta=gamma, replace 6 by **3 in both (4) and (5)**. More generally, writing W for the contraction in (4) without 6,

    collected Hamiltonian coefficient = 6 W/(1+delta_beta,gamma),
    one-parent -> normalized daughter-vacuum-pair amplitude
                                      = 6 W/sqrt(1+delta_beta,gamma).

Equality of frequency or block membership is insufficient for this Kronecker delta. If both daughter legs use the same q, reuse one physical basis/sewing map; their symmetric tensor indices are not independently rotatable copies of a mode frame. Thermal factorial moments, spectral densities, and linewidth conversions remain separate from these Hamiltonian factors.

## 3. What is—and is not—an allowed gauge change

For basis changes E_q'=E_q U_q and E_-q'=E_-q U_-q,

    S_q' = U_q^T S_q U_-q.

Within exact degenerate eigenspaces (or consistent eigenvalue relabelings), (1) remains valid with diagonal frequencies. The oriented map transforms with U_q0 on the annihilated index and conjugated U_-q1, U_-q2 on the created indices. This is **passive basis covariance**, not evidence that the interaction possesses a physical symmetry or that diagonal populations stay diagonal. Rotating different q points also changes the momentum representation and cannot retain the original q labels.

A numerical frequency cluster is not an exact eigenspace. With unrestricted overlap matrices, the raw polarization contraction R obeys Rminus=conj(R) contracted with S0,S1,S2. Since V is proportional to R/sqrt(omega0 omega1 omega2), the general source-contraction audit must instead use

    Vminus[a,b,c] = sum_ijk conj(V[i,j,k]) S0[i,a] S1[j,b] S2[k,c]
                   *sqrt(omega0_i omega1_j omega2_k /
                         (omegaMinus0_a omegaMinus1_b omegaMinus2_c)). (7)

Unnormalize before transforming, then renormalize. Equation (7) reduces to (3) when (1) holds; arbitrary overlap unitarity alone is insufficient. Applying (4)–(5) to (7) while ignoring frequency incompatibility need not yield (6).

One can passively mix unequal frequencies consistently by also transporting the full, now nondiagonal harmonic frequency matrix and oscillator factors. The operator description remains valid, but scalar eigenmode frequencies/resonances and a diagonal-population kinetic interpretation no longer follow. Keeping the old diagonal frequencies after such a rotation is the inconsistent shortcut.

## 4. Bounded checks and an operator alternative

An in-memory NumPy check used seed 9272026, complex-normal random V of shape (2,3,2), and independent Haar-style unitary S_l from complex QR (positive-diagonal phase correction). Vminus was constructed by (3). Relative Frobenius residuals:

| Check | Residual |
|---|---:|
| (5) versus conjugate of (4) | 3.45e-16 |
| Independent gauge rotations on each paired space | 4.09e-16 |
| Deliberately omit the conjugation of S0 in (5) | 1.78 |
| Use (7), unrestricted S, and unequal frequencies 1,...,d while retaining (4)–(5) | 0.366 |

A scalar phase check gives absolute error 6.28e-16. For a symmetric 3x3 daughter tensor, the normalized unordered-pair strength sum_(b<=c) |6 W_bc/sqrt(1+delta_bc)|^2 equals 18||W||_F^2 and is invariant under a common daughter rotation to 4.27e-16 relatively. These are finite algebra checks, not material or kinetic validation; the proof is (1)–(6).

A useful alternative importer object is the **channel operator** T=P_two H3 P_one, mapping the one-parent space into H_beta tensor H_gamma, or Sym^2(H_beta) for coincident daughter spaces. Its entries are normalized Fock amplitudes, its reverse is automatically T^dagger, and its transformation is T'=(U_beta^dagger tensor U_gamma^dagger) T U_parent, with the symmetric restriction where needed. This avoids treating basis-dependent squared entries as autonomous events. Exact resonant projection satisfies

    (Omega_beta tensor I + I tensor Omega_gamma) T = T Omega_parent.

It expresses harmonic-energy conservation; signed wavevectors express crystal momentum modulo G. It does not supply a Markov generator, justify secular elimination inside a degeneracy, or establish closure of block traces. Such steps require additional kinetic assumptions and phase-sensitive contractions.

**Application gate.** Check coordinate-frame alignment, complete matched subspaces, unitarity, (1), tensor permutations and (3), plus physical-label counting. The PI artifact supports a corrected raw contraction/sewing identity for one tuple; its reported minimum admitted detuning is 0.00640228208505 THz, so it supplies no exactly resonant event in the stated thresholded first-leg-parent subset. Its near-roundoff agreement cannot replace these assumptions or certify a material rate.