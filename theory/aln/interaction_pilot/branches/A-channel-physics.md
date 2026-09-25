# A — Channel physics for the pinned interaction pilot

Independent bounded first pass, 2026-09-24. Read only `00-question.md` and primary sources; no peer or PI results, installed environment, or prototype were inspected. Pin: phono3py 4.5.0, `21fa8f3817fbcc603254656f525bb5aec113afb6`. Assume positive-frequency scalar modes, reciprocal pairing, and permutation-consistent cubic interactions. This is a conditional derivation, not validation of an AlN export.

## 1. Interaction and channel orientation

With default conversions, `pp=|Phi_012|^2` is in eV^2. Its conversion already includes `1/(36*8*N_grid)` and inverse mode-frequency factors; it contains neither Bose factors nor energy-integration weights. Do not divide by N_grid again. Phi contains 1/3! in the ordered cubic Hamiltonian. [Pinned conversion](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/interaction.py#L177-L190); [Hamiltonian, equations 9–10](https://arxiv.org/pdf/1501.00691).

Stored triples satisfy q0+q1+q2=G. Write bar-i for the reciprocal physical mode and f_i for cyclic frequency in numerical THz. One consistent orientation, keeping physical mode 0 in every event, is:

| Weight | Resonance | Physical decay orientation |
|---|---|---|
| g0 | f0=f1+f2 | 0 -> bar-1 + bar-2 |
| g1 | f2=f0+f1 | bar-2 -> 0 + 1 |
| g2 | f1=f0+f2 | bar-1 -> 0 + 2 |

Each arrow includes its inverse through the Bose gain/loss difference. The global reciprocal event is separate unless it is the same unordered event. Absorption describes the fixed mode's participation, not N/U classification. Compute the physical reciprocal vector after orientation. [Channel definitions](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/triplets.py#L246-L252).

## 2. Bose law, rate conversion, and counting

For one unique event a <-> b+c,

    R=K_e[n_a(1+n_b)(1+n_c)-(1+n_a)n_b n_c],
    dot n=-r R,  r=e_a-e_b-e_c.

At exact resonance F0=n_a^0(1+n_b^0)(1+n_c^0)=(1+n_a^0)n_b^0 n_c^0. With W_ii=n_i^0(1+n_i^0), direct differentiation gives

    C_e=K_e F0 u u^T,  u=W^(-1/2)r,
    y=W^(-1/2)delta n.

Energy is annihilated; parent–daughter entries are negative and the daughter–daughter entry positive. Forward and reverse reactions are already included once.

The pinned self-energy formula is

    Gamma0=U_Gamma sum_(t,j1,j2) w_t pp_t
           [(n1+n2+1)g0+(n1-n2)(g1-g2)].               (1)

Gamma is a half-linewidth in numerical cyclic THz. With F=10^12 s^(-1), its usual lifetime rate is `4*pi*F*Gamma`. Numerically `U_Gamma=18*pi/[hbar_eVs^2*(2*pi*F)^2]`; record the actual runtime value. [Bose sum](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L767-L795); [conversion](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L200-L207); [lifetime arithmetic](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/utils.py#L74-L90).

For a fully unfolded unique unordered event with a consistent integration weight g_e in (1)'s convention,

    K_e=8*pi*F*U_Gamma*pp*g_e/(1+delta_bc).              (2)

Derivation: distinct daughters give amplitude 6 Phi from the six ordered Hamiltonian permutations, hence golden-rule weight `(2*pi/hbar)|6 Phi|^2 delta(epsilon_a-epsilon_b-epsilon_c)`. Equivalently, the parent linewidth sum counts (b,c) and (c,b). Identical daughters give amplitude 3 Phi; geometric Bose averages of pair creation/annihilation supply a common factor 2, leaving half the distinct-daughter coefficient. Accumulate repeated indices before the outer product: r=(1,-2), giving C_bb=4 K_e F0/W_bb. Gamma's lifetime diagonal need not be the complete collision diagonal for coincident indices.

Equation (2) is a conditional distribution/quadrature conversion, not a numerical prescription for delta(0). Verify counting, permutation symmetry, and runtime normalization before material use. Non-default conversion overrides invalidate the default pp interpretation.

## 3. Export and exactness requirements

- Triplet weights reproduce the ordered fixed-root sum, with sum_t w_t=N_grid. They aggregate swaps and crystal orbits. Multiplying one representative event by w_t does not reconstruct the full-grid populations changed. Unfold or explicitly represent the orbit, then deduplicate by parent plus unordered daughters. Do not add the same event again through absorption rows.
- Preserve mode/band indices, BZ/regular-grid and reciprocal maps, orbit maps/weights, runtime constants, backend/version, cutoffs and `g_zero` coverage, conversion overrides, and permutation/reciprocity checks. Skipped pp entries are not automatically physical zeros. Selected roots give an interaction sample, not a complete event network.
- Preserve individual channel weights or their complete integration recipe. `gamma_detail` combines Bose factors and channels; it is neither an event coefficient nor a collision operator.
- Gaussian weights include energy-mismatched tuples. Then A(n0)/B(n0)=exp[-beta(epsilon_a-epsilon_b-epsilon_c)] differs from one and energy changes by minus mismatch times R. Tetrahedron weights integrate a resonant surface; assigning them to nonresonant vertex energies does not make each tuple an exact event. Imposing positive outer products afterwards defines a modified model rather than the derivative of that nonlinear Bose law.

**Decision:** pp with verified mapping/counting can support subsequent reconstruction. First reproduce (1) and audit a unique exact-resonance event, retaining mismatch and integration provenance. A tiny broadened pilot or linewidth agreement does not certify a conserving material operator. The Rust implementation and actual pilot arrays were not checked. No repository files were edited.