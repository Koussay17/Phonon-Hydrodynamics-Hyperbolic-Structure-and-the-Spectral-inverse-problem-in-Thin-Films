# Primary conventions: one three-phonon channel

Independent literature first pass, 2026-09-26. Read only campaign question/assumptions and primary papers/code. No novelty claim or material calculation.

## 1. Verified primary convention

**SOURCE SAYS.** Togo, Chaput and Tanaka, *Distributions of phonon lifetimes in Brillouin zones*, Phys. Rev. B **91**, 094306 (2015), [DOI](https://doi.org/10.1103/PhysRevB.91.094306); inspected [arXiv:1501.00691v3](https://arxiv.org/html/1501.00691v3), whose title is singular. Equations (4), (5), (9), (10) give

\[
H_3=\sum_{abc}\Phi_{abc}A_aA_bA_c,\qquad A_a=a_a+a^\dagger_{-a}.
\]

The sum is ordered. The normal-mode tensor already contains $1/(3!\sqrt N)$ and three factors $\sqrt{\hbar/(2m\omega)}$. Its three polarization vectors are unconjugated; Fourier phases have positive signs. Momentum selection is $q_a+q_b+q_c=G$. Equation (5) uses $e^{iq\cdot r}W(q)$ and defines $-a=(-q,j)$. Equation (11), in angular frequency, is

\[
\Gamma_k(\omega)=\frac{18\pi}{\hbar^2}\sum_{ij}|\Phi_{-k,i,j}|^2
\left\{(n_i+n_j+1)\delta(\omega-\omega_i-\omega_j)
+(n_i-n_j)[\delta(\omega+\omega_i-\omega_j)-\delta(\omega-\omega_i+\omega_j)]\right\}.
\]

Here occupations are equilibrium Bose values. Equations (13), (18) distinguish lifetime $1/(2\Gamma_k)$ from its approximate use as the single-mode transport relaxation time.

**LOCAL PRIMARY CODE.** Installed phono3py 4.5.0, pinned commit `21fa8f3817fbcc603254656f525bb5aec113afb6`: [interaction.py L177-L190](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/interaction.py#L177) includes `1/36`, `1/8`, and `1/prod(mesh_numbers)` in conversion to eV^2. Its Python path, L1092-L1096, forms `pp = abs(fc3_normal)^2 * unit_conversion`. [reciprocal_to_normal.py L140-L174](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/reciprocal_to_normal.py#L140) contracts three unconjugated eigenvectors and divides by square-root frequencies/masses. [real_to_reciprocal.py L214-L239](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/real_to_reciprocal.py#L214) uses positive complex Fourier phases, including the basis-position prephase.

[imag_self_energy.py L200-L207, L775-L795](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L200) applies triplet weights and the Bose/delta combinations, then
$C_\Gamma=18\pi/[\hbar_{\mathrm{eV\,s}}(2\pi10^{12})]^2$.
Frequencies and numerical gamma use ordinary THz; the physical inverse lifetime is $4\pi10^{12}\gamma\ \mathrm{s}^{-1}$. Do not insert another mesh divisor or read pp as an unsymmetrized Fock transition probability. Native backend phase equivalence was not independently audited here.

## 2. Conditional channel translation, derived from the Hamiltonian

**ASSUMPTIONS.** Positive frequencies, exact resonance $\omega_k=\omega_i+\omega_j$, scalar modes, a consistently conjugate-sewn basis, and an actually permutation-symmetric cubic tensor. Decay $k\to i+j$ selects $\Phi_{k,-i,-j}$; absorption selects $\Phi_{-k,i,j}$. They are conjugates under these assumptions. Arbitrary numerical eigenvectors require verified phase maps, or unitary sewing within degenerate subspaces, with reciprocal-cell basis phases. The papers do not supply those maps for this AlN export.

Let $\phi=\Phi_{k,-i,-j}$ and $H_e=g_ea_i^\dagger a_j^\dagger a_k+\mathrm{h.c.}$; integer occupations are $m$, means are $n$.

| Daughters | Hamiltonian coefficient | Decay Fock matrix element | Absorption Fock matrix element |
|---|---|---|---|
| $i\ne j$ | $g_e=6\phi$ | $g_e\sqrt{m_k(m_i+1)(m_j+1)}$ | $g_e^*\sqrt{(m_k+1)m_im_j}$ |
| $i=j$ | $g_e=3\phi$ | $g_e\sqrt{m_k(m_i+1)(m_i+2)}$ | $g_e^*\sqrt{(m_k+1)m_i(m_i-1)}$ |

These count six versus three distinct ordered triples. Without permutation symmetry, replace $6\phi$ or $3\phi$ by the coherent sum of the distinct ordered amplitudes. Squared magnitudes alone cannot reconstruct that sum.

The spectral golden-rule expression is $2\pi|M|^2\delta(E_f-E_i)/\hbar=2\pi|M|^2\delta(\Delta\omega)/\hbar^2$. It is not a finite rate obtained by evaluating $\delta(0)$ for an isolated finite resonant Hamiltonian. A weak-coupling, Markov/Pauli kinetic approximation and a continuum or justified coarse-graining prescription are additional requirements.

**ADDITIONAL CLOSURE.** For a product of geometric occupation distributions, $\langle m(m-1)\rangle=2n^2$ and $\langle(m+1)(m+2)\rangle=2(1+n)^2$. These identities follow by summing the geometric distribution; they are not valid for an arbitrary state with the same mean. With this closure, one unordered reversible channel has

\[
\alpha_e=\frac{72\pi}{\hbar^2(1+\delta_{ij})}|\phi|^2\delta(\omega_k-\omega_i-\omega_j),
\quad J_e=\alpha_e[n_in_j(1+n_k)-(1+n_i)(1+n_j)n_k],
\quad \dot n=-s_eJ_e,\quad s_e=e_i+e_j-e_k.
\]

Thus $s_e=2e_i-e_k$ for repeated daughters. The equilibrium event weight is $\alpha_e\bar n_i\bar n_j(1+\bar n_k)$, distinct from $\alpha_e$, $g_e$, and pp. This translation is a conditional derivation, not an explicit repeated-mode theorem in either cited paper.

**LIMIT OF THE LINEWIDTH CHECK.** This coefficient reproduces the parent-mode decay contribution to $2\Gamma_k$: two ordered daughter terms when distinct, one when repeated. That check does not establish equality between the full population Jacobian and a linewidth-based diagonal on repeated modes. In particular differentiating a repeated-daughter flux also differentiates its second occurrence; an independent collision-generator audit is needed.

## 3. Independent collision-counting source and limits

**SOURCE SAYS.** Fugallo, Lazzeri, Paulatto and Mauri, *Ab initio variational approach for evaluating lattice thermal conductivity*, Phys. Rev. B **88**, 045430 (2013), [DOI](https://doi.org/10.1103/PhysRevB.88.045430); inspected [arXiv:1212.0470v2](https://arxiv.org/html/1212.0470v2). Equation (4) has a factor $1/2$ on the ordered decay sum and none on absorption. Equations (5)-(8) define signed vertices and occupation factors in their own normal-coordinate convention. After (13), detailed balance equates equilibrium decay and absorption weights using exact energy conservation. Section IV immediately after (30) explicitly notes that Gaussian broadening satisfies detailed balance only approximately. Its $V^{(3)}$ normalization has not been equated to Togo's $\Phi$ here.

**INTERPRETATION.** Daughter exchange is one quotient; forward/reverse flux is one reversible event. Equality of amplitudes under momentum reversal does not identify channels acting on different population indices. Neither source licenses discarding a time-reversed channel solely because its pp agrees. Broadening a finite off-resonant tuple does not preserve the exact Bose identity automatically.

## Scope, barriers, and next check

Queries: `Distributions phonon lifetimes Brillouin zones Togo Chaput Tanaka`; `Fugallo Lazzeri Paulatto Mauri variational lattice thermal conductivity 2013`. Local `research-search.cmd` queried Crossref and OpenAlex successfully (five results each/query); Semantic Scholar returned HTTP 429. DOI metadata and both complete primary arXiv HTML versions were inspected, plus the installed source files above. Logs: `C:\Users\Koussay\ResearchLab\literature\searches\20260926-channel-counting\` (`togo`, `collision-counting`). No exhaustive historical or novelty search; no older reference was attributed without inspection.

Unresolved for an actual import: coherent permutation consistency, negative-q sewing including reciprocal representatives, chosen repeated-mode occupation closure, and a legitimate resonant spectral measure. The minimal next step is a finite Fock-space one-channel enumeration and permutation/sewing check using complex amplitudes, followed by a separately stated kinetic approximation. A pp-only export or agreement with one linewidth cannot close these gaps.
