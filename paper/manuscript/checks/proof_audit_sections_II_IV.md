# Proof audit — manuscript Sections II–IV (draft of 11–12 November 2026)

**Date:** 12 November 2026. **Object:** `paper/manuscript/sections/kinetics.tex`, `identifiability.tex`,
`frontface.tex` (compiled `main.pdf`, statement numbers as printed), checked against `SOURCES.md`,
`checks/verify_manuscript_algebra.py` and the audited sources `theory/aln/memory_gap/13-final-report.md`,
`theory/aln/memory_gap/review/proof-audit.md`, `.../review/counterexample-review.md`, `notes/14`, `notes/18`
and `paper/PRIOR_ART.md`. The `.tex` files were not modified.

**Priority.** Material that is not in the audited sources was audited first: (1) the dyadic-shell proof of the mesh
laws in Theorem 3(i); (2) Proposition 1; (3) Theorem 1(iv) and Remark 1; (4) the variational proof of Lemma 2;
(5) the matrix-determinant-lemma derivation of adj(A) in Theorem 2; (6) Eq. (6); (7) Proposition 5 with F > 0 and
its converse. Then the audited statements were compared with their sources in hypotheses, conclusions, status
labels and novelty wording.

**Method.** Each step was re-derived by hand. Fragile steps were checked by computations that differ from the
manuscript's check script: a different network, a different basis, a different solution route, exact arithmetic,
certified root isolation, and lattice sums compared with sharp constants derived here. The scripts are in
`checks/audit/` (§14). There are 88 checks and all pass. Four checks first failed because of defects in the scripts
themselves: a float produced by the integer division `0/F`, a first example search that was too narrow, and two
convergence tolerances tighter than the predicted pre-asymptotic rates. All four were corrected and rerun; none
involved the manuscript.

**Classification.** VALID · VALID UNDER EXTRA ASSUMPTION (the assumption is stated) · INCOMPLETE (the conclusion is
true but the argument as written does not reach it) · UNJUSTIFIED · FALSE (counterexample given) · UNKNOWN.

---

## 0. Verdict at a glance

| Item | Claimed status | Verdict | Defect | Repair |
|---|---|---|---|---|
| Thm 3(i), integral criterion | PROVED | **VALID** | — | — |
| Thm 3(i), mesh laws (dyadic shells) | PROVED | **VALID UNDER EXTRA ASSUMPTION** (gap) | The Riemann sums converge to 𝓜_m only if J²r⁻ᵐ is Riemann integrable away from Γ. (H1)–(H3) do not give this, so the claims "converges for 2α < d" and "the Riemann sums converge to 𝓜_m" fail without it (§1.10). The continuity of J₀, r₀ is sufficient for the shell sums, and some regularity is necessary, but continuity is stronger than needed and fails at acoustic axes (§1.12) | Add (H5): Φ_m = Σ_ac J₀²r₀⁻ᵐ is Riemann integrable on S^{d−1}, and J²r⁻ᵐ is Riemann integrable on BZ ∖ B_δ for every δ > 0. Write the lattice lemma of §1.13. The argument then gives sharp constants, confirmed numerically |
| Thm 3(ii) | DERIVED UNDER (H4) | **VALID UNDER** (H4) and uniform (H1) | "(H1)–(H3) uniformly in N₁" is not defined. The status note names only (H4) | Define the uniformity (§1.14). Name both assumptions in the label |
| Thm 3(iii) | DERIVED UNDER (H4) | **VALID** (conditional) | — | — |
| Prop. 1 (grey Callaway) | DERIVED UNDER ASSUMPTIONS | **VALID**; one step **INCOMPLETE** | "both O(n⁻¹)" is asserted without argument, and the limit is not stated | One sentence (§2, step 2.7). Verified: n·max(dev) → 22.8 |
| Thm 1(iv) (involution i ↦ i ⊕ 3) | PROVED | **VALID** | — | — |
| Remark 1 (slow pole) | — | **VALID** | — | Optional: η²/4 + η⁴/16; N_slow = η⁻² + 3/4 + o(1) |
| Lemma 2, inequality (variational proof) | PROVED | **VALID** | — | — |
| Lemma 2, "the constant is sharp" | PROVED | **INCOMPLETE** | Shown only for ς_max = 3 and 2, in a supplement that is still a skeleton. No example exists for ς_max = 4 in any source | Exact 1 ↔ 3 example and equality criterion derived here (§4.6) |
| Thm 2, adj(A) via the determinant lemma | PROVED | **VALID**; notational gap | ψ_I is defined only for I ∈ 𝔅, but Eq. (24) sums over all (n_H − 1)-subsets | Define ψ_I for every I with \|I\| = n_H − 1 |
| Extension to non-admissible g | — | **VALID** | — | — |
| Eq. (6) (drive at k = 0) | derived | **VALID** as the k → 0 limit; **INCOMPLETE** as worded | A uniform but time-dependent gradient with a uniform deviation violates the energy balance unless an external source acts | State the external-source setting, or the k → 0 limit of the H-projected plane-wave problem (§6) |
| Prop. 5 (F > 0, converse) | PROVED | **VALID**; definitional gap in the converse | The ratio g_γ/g^ref_γ is undefined where g^ref_γ = 0 (allowed for admissible g^ref) | Take the min over supp g^ref, and require τ > 0 |
| Remark 5, "physical currents have full support" | — | **FALSE** (inherited from the sources) | j^(i) vanishes where v_{ν,i} = 0 | Use time-reversal parity (§10.1). The conclusion survives |
| Text after Prop. 4, "By Proposition 3, large memory requires…" | — | **UNJUSTIFIED** as attributed | Prop. 3 gives sufficiency, not necessity | Derive the necessity from Eq. (22) (§10.2) |
| Remark 5, "transfers only if I* is a union of orbits" | — | **UNJUSTIFIED** ("only if") | The necessity is not shown | "The proof transfers if …" |
| Statements vs. sources | — | faithful (§8) | Lemma 2 at ς_max = 4; the mesh laws were numerical in the sources | as above |
| Novelty wording vs. `PRIOR_ART.md` | — | **faithful** (§9) | — | — |
| Sec. IV (Prop. 6, Cor. 2, Eqs. 30–40) | DERIVED UNDER (F1)–(F6) | **VALID** (independent checks) | — | — |

**Earliest failure in reading order** (Secs. II–IV): Eq. (6), a presentational gap (§6). **Earliest FALSE statement:**
Remark 5, "physical currents have full support" (§10.1). **Earliest incomplete proof:** the sharpness claim of
Lemma 2 (§4.6). **Weakest new step:** the convergence part of the mesh laws in Theorem 3(i) (§1.10, §1.12). None of
these defects propagates into another proof. Every repair is short.

---

## 1. Theorem 3(i): infrared criterion and mesh laws (new material: dyadic shells, continuity, additive error)

Setting: (H1) r = r₀(k̂)k^α[1 + o(1)] uniformly, with r₋ ≤ r₀ ≤ r₊. (H2) J → J₀ uniformly, J bounded, and
0 < A_m := Σ_ac ∫ J₀²r₀⁻ᵐ dΩ < ∞. (H3) r is bounded below away from Γ and on the optical branches. The mesh is
k = 𝖡n/N₁, and the mesh moments are 𝓜_m(N₁) = 𝒱⁻¹ Σ_ν J²r⁻ᵐ.

| # | Step | Class. | Reason |
|---|---|---|---|
| 1.1 | The integrand is bounded outside a ball k < k₀ | VALID | By (H3), r⁻ᵐ ≤ c(k₀)⁻ᵐ for m ≥ 0, and J is bounded |
| 1.2 | J²r⁻ᵐ = [J₀²r₀⁻ᵐ + ε(**k**)]k^{−mα} with sup_k̂ \|ε\| → 0 | VALID | J = J₀ + δ_J with δ_J → 0 uniformly, and \|J₀\| ≤ sup\|J\|. Then r⁻ᵐ = r₀⁻ᵐk^{−mα}(1 + o(1)) with r₀⁻ᵐ ≤ r₋⁻ᵐ, so ε = (2J₀δ_J + δ_J²)r₀⁻ᵐ(1 + o(1)) + J₀²r₀⁻ᵐ·o(1). **The additive form is necessary:** J₀ = k_BT v_s k̂_i vanishes on k̂_i = 0, where a relative error J²/J₀² − 1 is unbounded. The lattice test (§1.15) uses J² − J₀² = O(\|k\|), which does not vanish on that plane |
| 1.3 | The ball contributes between A/2 and 3A/2 times ∫k^{d−1−mα}dk | VALID | \|ε̄(k)\| ≤ n_ac\|S^{d−1}\| sup\|ε\|. Choose k₀ so that this is ≤ A/2 |
| 1.4 | 𝓜_m < ∞ ⟺ mα < d, with ln(1/k_min) at equality | VALID | Monotone convergence as k_min → 0 |
| 1.5 | "the memory moment [is finite] iff 2α < d" | VALID for α < d; wording INCOMPLETE for α ≥ d | For α ≥ d both 𝓜₂ and 𝓜₁ are infinite, so 𝓜₂/𝓜₁ is ∞/∞. Define τ^RTA(k_min) = 𝓜₂(k_min)/𝓜₁(k_min). It diverges like k_min^{−α} for α > d and like k_min^{−d}/ln(1/k_min) at α = d. Repair: one clause |
| 1.6 | Small-k mesh sum = N₁^{mα−d}Σ[φ(n) + o(1)\|𝖡n\|^{−mα}] (constant factor 1/V_cell) | VALID | The "o(1)" is ε(𝖡n/N₁). It is uniformly small only after k₀ has been chosen small, as in 1.3, and not uniformly in the shell index |
| 1.7 | Shell sums 2^{l(d−mα)}[c* + o(1)] with c* > 0 | VALID UNDER the added regularity | Substituting n = 2^l u gives a Riemann sum over the grid 2^{−l}ℤ^d of the homogeneous φ on the Jordan-measurable shell 1 ≤ \|𝖡u\| < 2. It converges for Riemann-integrable φ, for example J₀²r₀⁻ᵐ continuous with r₀ > 0. The o(1) is as l → ∞. For small l only the uniform bound Σ_shell φ ≤ C2^{l(d−mα)} holds, and that bound is what the summation needs. **Some regularity is necessary:** take J₀² = 1 on the irrational directions and 0 on the lattice directions 𝖡n/\|𝖡n\|, with r₀ = 1. Then A = \|S^{d−1}\| > 0, so (H1)–(H3) hold, but every mesh point has J₀ = 0 and the mesh laws fail |
| 1.8 | "the o(1) term is handled as in the integral" | VALID, INCOMPLETE as written | Needs Σ_shell \|𝖡n\|^{−mα} ≤ C′2^{l(d−mα)} for all l, with k₀ chosen so that sup\|ε\| ≤ c*/(2C′) (mα ≤ d). For mα > d, use dominated convergence instead: the sum is dominated by small \|n\|, where ε(𝖡n/N₁) → 0 |
| 1.9 | Σ_{l ≤ log₂N₁} gives orders 1, ln N₁, N₁^{mα−d} | VALID | Stolz–Cesàro for mα ≤ d. For mα > d the series converges to S_φ = Σ_{n≠0}φ(n) > 0, which is positive because the large-l shell sums are. The argument gives **sharp constants** (§1.13), confirmed in §1.15 |
| 1.10 | "for mα < d … the Riemann sums converge to 𝓜_m" | **UNJUSTIFIED from (H1)–(H3)**; VALID UNDER Riemann integrability of J²r⁻ᵐ on BZ ∖ B_δ | The shells near the origin are uniformly small, which is valid. The convergence of the rest needs Riemann integrability, and (H1)–(H3) assume only boundedness away from Γ. Counterexample: replace J by J·(1 − 𝟙_{mesh-rational points}) for \|k\| > k₀. This changes nothing a.e., so the hypotheses hold and 𝓜_m is unchanged, but the large-k part of every mesh sum vanishes. With a residue-class selection the mesh sums can even oscillate. The same gap affects the preamble's "\|j\|² → 𝓜₀, K_RTA → 𝓜₁" and the "\|j\|² → 𝓜₀" in the proof of (ii). Part (ii) needs only \|j\|² ≤ n_b sup J²/V_cell |
| 1.11 | Ratio 𝓜₂/𝓜₁ gives the four laws | VALID given 1.9 and 1.10 | The growth laws need only Θ-bounds on 𝓜₁(N₁), which hold without 1.10's assumption. Only "converges for 2α < d" needs 1.10 |
| 1.12 | Scope of "J₀ and r₀ continuous" | sufficient but too strong | In anisotropic crystals the frequency-sorted acoustic branches touch at acoustic axes (conical degeneracies) and may cross along curves. There the group velocity of a sorted branch depends on the direction of approach, so J₀ and r₀ are discontinuous on a closed null set of directions. Even the branch sum Σ_± v_{±,i}² is direction-dependent at a conical point. The argument needs only Riemann integrability of Φ_m = Σ_ac J₀²r₀⁻ᵐ on S^{d−1} (bounded, and continuous outside a closed null set), which these crystals satisfy |

**1.13 Repaired statement (proposed wording of the lattice lemma).** Let Φ ≥ 0 be Riemann integrable on S^{d−1} with
A = ∫Φ dΩ > 0, let s ≥ 0 and φ(n) = Φ(𝖡n/|𝖡n|)|𝖡n|^{−s}, and let |ε_N(n)| ≤ δ(|𝖡n|/N) with δ(t) → 0 as t → 0. Put
S_N = Σ_{0<|𝖡n|<k₀N}[φ(n) + ε_N(n)|𝖡n|^{−s}]. Then:
* for s < d, N^{s−d}S_N = |det 𝖡|⁻¹A k₀^{d−s}/(d − s)·(1 + O(sup_{t<k₀}δ)) + o(1);
* for s = d, S_N = |det 𝖡|⁻¹A ln N + o(ln N);
* for s > d, S_N → S_φ := Σ_{n≠0}φ(n) ∈ (0, ∞).

With |det 𝖡|V_cell = (2π)^d this gives:
* mα < d: 𝓜_m(N₁) → 𝓜_m, under Riemann integrability away from Γ;
* mα = d: 𝓜_m(N₁) = (2π)^{−d}A_m ln N₁ + o(ln N₁);
* mα > d: 𝓜_m(N₁) = N₁^{mα−d}S_φ/V_cell·(1 + o(1)).

For τ^RTA:
* 2α < d: → 𝓜₂/𝓜₁;
* 2α = d: (2π)^{−d}A₂ ln N₁/𝓜₁;
* d/2 < α < d: N₁^{2α−d}S_{φ,2}/(V_cell𝓜₁);
* α = d: (2π)^d S_{φ,2}N₁^d/(V_cell A₁ ln N₁).

**1.14 Part (ii): definitional gap.** "(H1)–(H3) uniformly in N₁" must mean that there is δ(k) → 0, independent of N₁,
with |r_{N₁}(k_ν)/(r₀(k̂_ν)k_ν^α) − 1| ≤ δ(k_ν) for all N₁ and all acoustic modes with k_ν < k₀. Only the upper side
r_{N₁} ≤ r₊k^α(1 + δ) enters the lower bound on N. The remaining steps are VALID:
* (H4) gives N ≥ c²Σ_{k_ν<k₀} j_ν²r_ν⁻²;
* Prop. 2 gives τ ≥ N^{1/2}/|j|;
* the one-sided rates are N₁^{α−d/2} and (ln N₁)^{1/2}.

The status note "(namely (H4))" should also name the uniform (H1) for diag 𝒞_{N₁}. By the paper's own definition,
that is also an unproved property of a specific operator family.

**1.15 Numerical confirmation** (`audit/audit_ir_lattice.py`). The test uses d = 2, 𝖡 = 2πI, odd N₁ ≤ 3201,
J² = (k̂_x² + 0.3|k|)e^{−k²/8} (additive error, with J₀ vanishing on a plane) and
r = (1 + 2.4k̂_x²k̂_y²)k^α(1 + 0.5k²). The sharp constants of §1.13 are reproduced:

| case | prediction | result at N₁ = 3201 |
|---|---|---|
| α = 0.6 (2α < d) | τ → 𝓜₂/𝓜₁ = 0.6291354 | rel. error 2.6·10⁻³, observed rate N₁^{−0.792} (predicted d − 2α = 0.8) |
| α = 1 (2α = d) | dτ/d ln N₁ = 0.3885388 | local slope / prediction = 0.99657, rising monotonically |
| α = 1.5 (d/2 < α < d) | 𝓜₂/N₁^{2α−d} → (2π)^{−2α}S₂ | within 2.5·10⁻³. 𝓜₁(N₁) → 𝓜₁ at the observed rate N₁^{−0.500} (predicted N₁^{−(d−α)}). Extrapolated τ ratio 0.99785 |
| α = 2 (= d) | d𝓜₁/d ln N₁ = (2π)^{−2}A₁, 𝓜₂/N₁² → (2π)^{−4}S₂ | 1.00087 and 1.00158. The combined τ ln N₁/N₁^d is still 1.057 times its limit (correction O(1/ln N₁)) |

NUMERICAL OBSERVATION relevant to Sec. V. Mesh laws approach their constants slowly: like N₁^{−(d−α)} through 𝓜₁
itself, and like 1/ln N₁ at α = d. Exponents fitted on coarse meshes are therefore biased. This is consistent with the
slow drift of the fitted exponent reported in the source audit (§10.1, α = 3.5: 3.37 → 3.44).

---

## 2. Proposition 1: grey conserving Callaway closure (derived in the manuscript)

| # | Step | Class. | Reason |
|---|---|---|---|
| 2.1 | Eq. (12): the normal term relaxes to the orthogonal projection onto span{e, P}; the resistive term relaxes to Π_e y | VALID (declared model) | The energy and momentum densities are eᵀy and P_iᵀy. The element of the span that matches them is the orthogonal projection, whatever the Gram matrix. The resistive term is the energy-conserving version (the compatibility conditions of `notes/14` §2) |
| 2.2 | j^(i) = v²P_i | VALID | ε_ν v_{ν,i} = ħv²k_{ν,i} |
| 2.3 | 𝒞_Cal j = j/τ_R | VALID | j ∈ range Π_N ∩ H. Not stated but true: 𝒞_Cal\|_H has rates 1/τ_R (multiplicity 3) and 1/τ_c (the rest), so ker 𝒞_Cal = span{e} and (K3) holds. Prop. 1 uses Def. 1, which requires (K3). One line would close this |
| 2.4 | λ(p) = λ/(1 + pτ_R), λ = Λ_iτ_R | VALID | K(p) = \|j\|²/(p + 1/τ_R) |
| 2.5 | Λ_i = ρc v²/3 | VALID UNDER cubic symmetry | Σ c_ν k̂_i² = Σ c_ν/3 (exact on O_h meshes) |
| 2.6 | τ_mem = τ_stat = τ_R | VALID | |
| 2.7 | r_ν = τ_N⁻¹ + τ_R⁻¹ − O(n⁻¹) | **INCOMPLETE** | The proof asserts "both O(n⁻¹)" without argument and without naming the limit. Repair: under refinement at fixed T and a fixed zone, use c_ν ≤ k_B and the boundedness of w_ν²k_{ν,i}² (which → 0 as k → 0), while Σc_ν = 𝒱ρc and Σ w²k_i² grow like n. The deviation is ≥ 0 |
| 2.8 | Lifetimes ≈ τ_c; (τ_c, λ) fixes (τ_N, τ_R) | VALID | Given ρc and v |
| 2.9 | Mutual orthogonality of e, P_x, P_y, P_z | VALID | It follows from inversion and the reflections, and is needed only for the projector diagonals in 2.7 |

Independent check (`audit/audit_callaway.py`): full O_h meshes N = 5 to 33 (n ≤ 35,936), Bose weights, matrix-free.
* j = v²P to 3·10⁻¹⁶, 𝒞j = j/τ_R to 10⁻¹⁴, and \|j_x\|²/\|e\|² = 1/3 to 12 digits.
* n·max(dev) = 20.7, 22.1, 22.5, 22.6, 22.7, 22.8, which is bounded.
* max\|1/r − τ_c\|/τ_c falls from 4.5·10⁻² to 1.7·10⁻⁴.

The manuscript's own check uses a 12-mode star mesh and does not test the factor 1/3.

---

## 3. Theorem 1(iv), the involution i ↦ i ⊕ 3, and Remark 1

| # | Step | Class. | Reason |
|---|---|---|---|
| 3.1 | 𝖯h_k = (−1)^{⟨3,k⟩}h_k | VALID | popcount((i⊕3)∧k) ≡ popcount(i∧k) + popcount(3∧k) (mod 2) |
| 3.2 | [𝖯, 𝒞_η] = 0 ⟺ 𝖳 couples only indices of equal sign | VALID | 𝖯 = 𝖧 diag(s_k)𝖧ᵀ |
| 3.3 | Signs: k = 1, 2 → −1; k = 4, 7 → +1; k = 0 → +1; k = 2 → −1 | VALID | Hence 𝖯e = e and 𝖯j = −j. The involution is fixed-point-free because 3 ≠ 0 |
| 3.4 | Remark 1: slow root −η²/2 + O(η⁴), residue η²/4 + O(η⁴) | VALID | Exactly: −η²/2 + η⁴/8 + O(η⁶) and η²/4 + η⁴/16 + O(η⁶). The residues sum to \|j\|² = 1 |
| 3.5 | "about half of K₀ and all of the divergence of N" | VALID | The K₀ share is 1/2 + O(η²). N_slow = η⁻² + 3/4 + o(1), against N = 1 + η⁻² |
| 3.6 | "within a relative η²/(2ω_min) of 1/(p + 2)" | VALID | This is (v). On 324 grid points the maximum of LHS/RHS is 0.99995 (sharp as ω ≫ η²) |

The other parts were re-verified independently (`audit/audit_thm1_walsh.py`, 29 checks):
* (i) from the exact characteristic polynomial of 𝒞_η;
* (iii) from the full 8 × 8 resolvent over ℚ(p, η), with no reduction to the block, and from a regularized solve;
* Lemma 1's time-domain form, by 60-digit quadrature;
* the spectral data, by 60-digit eigen-decompositions at η = 0.2, 0.05, 0.01;
* Example 1, by eigen-decomposition; the graph-Laplacian claim (both operators in the cone of the δ_μ − δ_ν) also holds;
* Remark 3, on 200 random rational cones.

**Theorem 1: VALID.**

---

## 4. Lemma 2: K₀ ≥ K_RTA/ς_max (second proof by the variational principle)

| # | Step | Class. | Reason |
|---|---|---|---|
| 4.1 | K₀ = max_y(2jᵀy − yᵀ𝒞y) for j ∈ range 𝒞 | VALID | The maximum over ℝⁿ is attained at 𝒞⁺j. The functional is invariant under y → y + te |
| 4.2 | ζζᵀ ⪯ s diag(ζ∘ζ), with s = \|supp ζ\| | VALID | Cauchy–Schwarz on the support, whatever the coefficients (−1, 2, …) |
| 4.3 | 𝒞 ⪯ ς_max diag 𝒞 | VALID | Only active events contribute, and s_γ ≤ ς_max. Exact char.-polynomial sign test on 48 admissible random networks (3-phonon, repeated daughters, 2 ↔ 2) |
| 4.4 | The test vector y = R⁻¹j/ς_max | VALID | R is invertible because r_ν > 0 under (K3). y need not lie in H |
| 4.5 | K₀ ≥ 2K_RTA/ς_max − K_RTA/ς_max | VALID | Independent of the source route (b_μ² ≤ r_μD_μ plus double counting) |
| 4.6 | "the constant is sharp" | **INCOMPLETE** | See below |
| 4.7 | "4 for four-phonon events" | imprecise | Four-phonon events with a repeated mode (2p → a + b, p + q → 2a) touch 3 modes. Write "up to 4" |

On 4.6. The proof cites ratios 1/3 (eleven modes) and 1/2 (a chain), deferred to Supplemental Sec. S2, which is
currently a skeleton. For ς_max = 4 no example exists in any source. The counterexample review's adversarial search
over mixed 3 + 4 networks reached K₀ς/K_RTA = 1.0000172, which is an approach, not an attainment, and pure 2 ↔ 2
networks are never admissible.

**Repair (derived and checked here, `audit/audit_lemma2.py`).**
* *Equality criterion.* K₀ = K_RTA/ς for some j ∈ H ∖ {0} ⟺ there is y ≠ 0 with 𝒞y = ς diag(𝒞)y. The current is
  then j = ς diag(𝒞)y. Forward direction: equality in 4.5 forces yᵀ(ς diag 𝒞 − 𝒞)y = 0, and since ς diag 𝒞 − 𝒞 ⪰ 0
  this gives (ς diag 𝒞 − 𝒞)y = 0. Reverse direction: substitute j = ς diag(𝒞)y.
* *Consequence.* Each event term of yᵀ(ς diag 𝒞 − 𝒞)y is ≥ 0, so each vanishes. An event with fewer than ς modes
  forces y = 0 on its support. A ς-mode event touching such a mode forces y = 0 on its own support. An admissible
  network is connected and covers every mode (otherwise energy restricted to a component is an extra invariant), so
  y = 0 everywhere. **Equality requires every active event to touch exactly ς_max modes.** For networks that mix
  three-phonon and 2 ↔ 2 events the inequality is strict but approached: min_j 4K₀/K_RTA = 1.2198, 1.0234, 1.00236,
  1.000024 for 3-phonon rates δ = 10⁻¹, …, 10⁻⁵.
* *Exact examples* (non-Bose rational weights, unequal rates, exact arithmetic):

| ς | modes | energies | events | K₀/K_RTA |
|---|---|---|---|---|
| 2 | 5 | 1, 2, 4, 8, 16 | chain 2a → a | 1/2 exactly |
| 3 | 8 | 6, 6, 6, 3, 3, 3, 3, 3 | 7 events 6 → 3 + 3 | 1/3 exactly |
| **4** | **9** | **6, 6, 6, 2, 2, 2, 2, 2, 2** | **8 events 6 → 2 + 2 + 2 (1 ↔ 3 four-phonon processes)** | **1/4 exactly** |

Each example is admissible, with stoichiometric rank n − 1.

**Minimal wording:** "The constant is attained exactly for ς_max = 2, 3, 4 (… Supplemental Sec. S2); equality requires
every active event to touch ς_max modes, so for networks mixing three- and four-phonon events the inequality is strict
but not improvable."

---

## 5. Theorem 2: adj(A) by the matrix determinant lemma; extension; Proposition 3

| # | Step | Class. | Reason |
|---|---|---|---|
| 5.1 | det(A + yzᵀ) = det A + zᵀadj(A)y | VALID | Holds for every square A, singular included (symbolic 3 × 3 check) |
| 5.2 | A + yzᵀ = [Z, y]diag(g, 1)[Z, z]ᵀ | VALID | |
| 5.3 | Cauchy–Binet split: subsets without the last column give det A; subsets with it give Σ g_I det[ζ_I, y]det[ζ_I, z] | VALID | The same column order in both factors, so no sign issue |
| 5.4 | adj(A) = Σ_{\|I\|=n_H−1} g_Iψ_Iψ_Iᵀ | VALID; **notational gap** | The theorem defines ψ_I only for I ∈ 𝔅, but the sum runs over all (n_H − 1)-subsets. Independent non-basis subsets contribute: 39 of 126 in the test, for a current in the span of two events. Define ψ_I for every I |
| 5.5 | ψ_Iᵀj ≠ 0 ⟺ I ∈ 𝔅, and then ψ_I/(ψ_Iᵀj) = y_I | VALID | Checked against a direct solve of Eq. (21) for every subset |
| 5.6 | Eqs. (22)–(23) | VALID | Exact for a generic j and for j in the span of two events |
| 5.7 | Independence from the orthonormal basis | VALID (implicit) | ψ_I changes by a common sign. adj is invariant under orthogonal changes of basis |
| 5.8 | Extension to non-admissible g with j ∈ range 𝒞(g) | VALID | Example with rank 𝒞(g) = 3 < 5. w_θ → 𝒞(g)⁺j, u_θ → −(Π_Q E Π_Q)⁻¹Π_Q E w and K₀(g_θ) → K₀(g), all with errors O(θ); w ⊥ u_∞ and \|w\| ≤ \|w + u_∞\| ≤ M K₀(g) exactly |
| 5.9 | Prop. 3: det A = Θ(ε), K₀ ~ c/ε | VALID; wording | "at least one such subset is a basis" should read "some I* ∪ {γ}, γ ∉ I*, is a basis", which is what gives the lower bound ε·det(Z_S)². That basis exists because I* is independent and Z spans H. Numerically det A/ε → 106.117, K₀ε → 9.824·10⁻⁵ and τ/(M²K₀) → 1 (1.000000 at ε = 10⁻¹⁴) |

Also checked (`audit/audit_thm2_cone.py`, 22 checks):
* Remark 5's attained case: j = ζ_β gives x/K₀ = y_{J*} exactly, for g supported on J* ∪ {β}.
* Prop. 4 on all 119 bases of a network with non-uniform W.
* LMI (27), including PSD at K* = K₀(g) for the non-admissible g.

All checks are in mode space, using det_{ℝⁿ}[ζ_I, y, e] = \|e\|·ψ_Iᵀy, so no orthonormal basis is needed.

---

## 6. Eq. (6): drive at zero wave vector (derived in the manuscript)

| # | Step | Class. | Reason |
|---|---|---|---|
| 6.1 | The streaming of the local equilibrium has entropy coordinates Σ_k j^(k)∂_kT/(k_BT²) | VALID | ∂_T f⁰ = w²ε/(k_BT²) and 𝖵_k e = j^(k) |
| 6.2 | ∂_t f⁰ ∝ e drops out on H, and its e-projection is the energy balance | VALID | |
| 6.3 | A uniform deviation removes the H-streaming term | VALID | |
| 6.4 | Self-consistency of "uniform ∇T(t)" with a uniform deviation | **INCOMPLETE** (presentational) | T = T₀ + r·G(t) with Ġ ≠ 0 and a uniform q violates ρc∂_tT + ∇·q = 0 unless an external heat source maintains T(r, t). Eq. (6) is exact as the k → 0 limit, at fixed gradient amplitude, of the H-projected plane-wave problem (p + 𝒞 + ikΠ_H𝖵_xΠ_H)y = −ik j_xθ/(k_BT²), whose e-component (the energy balance) decouples |
| 6.5 | Solution (7) and transform (8): sign and kernel | VALID | q(t) = −G∫₀ᵗΦ for a step gradient |

Check (`audit/audit_prop5_drive.py`), on a 12-mode time-reversal-symmetric system:
* −q/(ikθ) → jᵀ(p + 𝒞)⁻¹j with errors 2.48·10⁻⁶, 2.48·10⁻⁸, 2.48·10⁻¹⁰ for k = 10⁻¹, 10⁻², 10⁻³. The order is 2.000,
  because time reversal removes the O(k) term.
* The e-projection reproduces the energy balance with the source s_ext ≠ 0.

**Repair (one sentence):** "…switched on at t = 0 and maintained by an external heat source (equivalently, the limit
k → 0 at fixed gradient amplitude of the plane-wave problem, whose energy component decouples)."

---

## 7. Proposition 5: per-event prior, F > 0, converse

| # | Step | Class. | Reason |
|---|---|---|---|
| 7.1 | 𝒞(g) − F⁻¹𝒞(g^ref) ⪰ 0 | VALID for every F > 0 | It is a sum of PSD terms with coefficients ≥ 0 |
| 7.2 | Γ_min(𝒞(g)) ≥ Γ_min(𝒞(g^ref))/F > 0, hence admissibility | VALID | Courant–Fischer on the invariant subspace H |
| 7.3 | τ ≤ 1/Γ_min(𝒞(g)) | VALID | Lemma 1 |
| 7.4 | F > 0 instead of F ≥ 1 | VALID | F < 1 is exactly what the converse needs when every ratio g_γ/g^ref_γ exceeds 1 |
| 7.5 | Converse via F = max_γ g^ref_γ/g_γ | VALID with a **definitional repair** | Admissible g^ref may vanish on some events, where g_γ/g^ref_γ is undefined. Take the min (and the max) over supp g^ref. If g_γ = 0 for some γ in supp g^ref, the min is 0 and the converse is trivial. Require τ > 0 |

Check: a 7-mode, 12-event network with three zero reference rates. Γ_min(𝒞_ref) = 0.99153889740 was enclosed to
10⁻³⁰ by exact real-root isolation. The bound was certified in 30 cases with F ∈ {1/3, 1/2, 1, 3, 10}. The converse
was certified on 40 random rate vectors, 12 of them with ρ = 0.

---

## 8. Faithfulness of the audited statements to their sources

| Manuscript | Source | Hypotheses | Conclusion | Label |
|---|---|---|---|---|
| Ex. 1 | `notes/18` §3.1 | same (L(1/3, 1/3, 1/3), L(1/2, 0, 1/2) written explicitly) | same | PROVED ✓ |
| Thm 1 (i)–(iii), (v) | `notes/18` §3.2 (t → η) | same | same; the consequence τ ∈ (2, ∞) is derived | PROVED ✓ |
| Prop. 2 | 13 P4(a) + audit §10.2 | same | same, plus \|x\|/\|j\| | PROVED ✓ |
| Lemma 2 | 13 L3, audit §5 | same (r = diag 𝒞, distinct modes) | the sharpness claim is wider than its support (ς_max = 4) | PROVED for the inequality; sharpness see §4.6 |
| Def. 5, Thm 2 | 13 T1, audit §2 | same (Z spans H, g admissible, j ≠ 0) | same (active bases) | PROVED ✓ |
| Extension | audit §2.8 | same | same | ✓ |
| Prop. 3 | 13 T2, audit §3 | same | same (approached, not attained; rates diverge) | PROVED ✓ |
| Rem. 5 | 13 §2, audit §3.4 | same | same; "full support" inherited (FALSE, §10.1) | OPEN ✓ |
| Prop. 4 | 13 P3, audit §4 | same | same | PROVED ✓ |
| Cor. 1 | 13 P6, audit §8 | same | same | PROVED ✓ |
| LMI (27) | 13 P5, audit §7 (extended value) | same | same; "no convexity of τ" kept | PROVED ✓ |
| Prop. 5 | 13 P7, counterexample review §3.3 | **weakened**: F > 0 instead of F ≥ 1 (valid) | same, with a strict converse (valid) | PROVED ✓ (§7.5 repair) |
| Thm 3(i), criterion | 13 T7(i) | **generalized**: asymptotic (H1) instead of an exact power law | same | PROVED ✓ |
| Thm 3(i), mesh laws | 13 T7(i); audit §10.1 checked them **numerically only** | continuity of J₀, r₀ added | same four laws | PROVED only with §1.13 (new proof) |
| Thm 3(ii) | 13 T7(ii), audit §10.2 | uniform (H1)–(H3) and continuity added | same one-sided rates; two-sided rates disclaimed | DERIVED UNDER (H4) ✓ (§1.14) |
| Thm 3(iii) | 13 T7(iii), audit §10.4 | same | plus "M → ∞ if K₀ bounded", which is trivially valid | DERIVED UNDER (H4) ✓ |
| Ex. 2 | `notes/18` §6 | same | plus the expansion 1 − (π/2)√p + O(p), which is valid | PROVED ✓ |

No hypothesis was dropped. The two generalizations (asymptotic (H1); F > 0) are covered by the manuscript's own
proofs, which were re-derived here.

## 9. Status labels and novelty wording

* **Labels.** These are faithful, with three exceptions:
  * Theorem 3(i) is PROVED only once (H5) is added (§1.10, §1.12).
  * The sharpness half of Lemma 2 rests on a supplement that is not yet written, and it lacks ς_max = 4.
  * The parenthetical "(namely (H4))" of Theorem 3(ii) omits the uniform (H1).
* **Novelty.** Every novelty sentence in Secs. I–IV matches `paper/PRIOR_ART.md`:
  * C2a (Remark 2, Intro item 2): "not found … in the literature searched";
  * C2b (Remark 4, Intro item 3): inequality attributed to Ben-Tal–Teboulle, application "potentially novel within the
    literature searched";
  * C5 (after Prop. 6, Intro item 5): "elementary dimensional analysis … not found";
  * C1 and C3: background;
  * Theorem 3: "no novelty claimed";
  * Prop. 1: "a definition and an elementary check".

  No novelty is claimed for Lemma 2 or Prop. 5, which is appropriate (the source marks Lemma 2 "novelty not checked";
  𝒞 ⪯ ς_max diag 𝒞 is the hypergraph analogue of the bound L ⪯ 2·Deg for graph Laplacians). The search scope is
  stated (Intro, status paragraph). No overclaim was found.

## 10. Other findings

**10.1 Remark 5, "Physical currents (4) have full support": FALSE.** j^(i)_ν = 𝒱^{−1/2}w_νε_νv_{ν,i} vanishes wherever
v_{ν,i} = 0: on the plane k_i = 0 of the isotropic/cubic model (`audit_callaway.py`: 24, 80, …, 1088 zero entries =
N² − 1), on mirror planes perpendicular to i, and for the optical modes at Γ. The same wording is in
`memory_gap/13-final-report.md` §2 and `review/proof-audit.md` §3.4. **The conclusion survives.** Time reversal gives
j_ν̄ = −j_ν with ε_ν̄ = ε_ν, so supp j is a union of pairs {ν, ν̄} of equal-energy modes. A three-mode event has odd
support, and a repeated-daughter event joins modes of different energy, so neither can be parallel to j. *Repair:*
"Physical currents are odd under time reversal, so their support consists of pairs of equal-energy modes and they are
never parallel to a three-phonon event vector."

**10.2 Text after Prop. 4, "By Proposition 3, large memory requires rates concentrated on such a sub-network".**
Proposition 3 shows that concentration is *sufficient*. The necessity follows from Eq. (22): if τ_mem ≥ θ²M²K₀, then
Σϖ_I|y_I| ≥ θM, so a convex weight of at least (θ − θ′)/(1 − θ′) lies on bases with |y_I| ≥ θ′M. The weight is in
the products g_I(ψ_Iᵀj)², not necessarily in individual rates. *Repair:* "By Eq. (22), large τ_mem/K₀ requires convex
weight on bases with large |y_I|; Proposition 3 shows that concentrating the rates on one such basis attains the bound."

**10.3 Remark 5, "Proposition 3 transfers only if I* is a union of orbits".** "If" holds: g^ε is then G-invariant. "Only
if" is not shown: another maximizing basis may be a union of orbits, and the symmetric supremum has not been computed.
*Repair:* "the proof of Proposition 3 transfers when some maximizing I* is a union of orbits; otherwise the symmetric
supremum is not known."

**10.4 Sec. IV.** Independent checks (`audit/audit_frontface.py`, 6 checks):
* Z from a direct solution of the two-layer boundary-value problem equals Eq. (35) to 3·10⁻⁴¹ at five complex p.
* σ → −σ leaves Z unchanged.
* The gauge 𝒢_μ leaves Z unchanged in the direct solution.
* Cor. 2 gives the ratio 3/2.
* The half-space phase (40) holds on the decaying branch.
* The entropy production equals (q² + L²q_z²)/(λT₀²) plus a third-order remainder qθ_z(1/T₀² − 1/T²).

The no-extra-boundary-condition claim was also checked: three constants (two in the film, one decaying in the
substrate) meet three conditions. **Sec. IV: VALID**, within (F1)–(F6).

## 11. Earliest failures in reading order

1. Eq. (6): the energy-balance/external-source gap (presentational; §6.4).
2. Prop. 1, step 2.7: O(n⁻¹) without argument (§2).
3. Lemma 2: sharpness incomplete for ς_max = 4, and its supplement is not written (§4.6).
4. Theorem 2: the scope of the definition of ψ_I (§5.4).
5. Remark 5: "full support" FALSE (§10.1); "only if" unproved (§10.3).
6. Text after Prop. 4: misattributed necessity (§10.2).
7. Prop. 5: the domain of the ratio in the converse (§7.5).
8. Theorem 3(i): Riemann-integrability gap and the scope of the continuity hypothesis (§1.10, §1.12, the most
   substantive); "finite iff" for α ≥ d (§1.5); compressed shell bookkeeping (§1.8).
9. Theorem 3(ii): the undefined uniformity, and the label (§1.14).

## 12. Minimal repairs (consolidated)

* Thm 3: add (H5): "Φ_m = Σ_ac J₀²r₀⁻ᵐ is bounded and continuous outside a closed null set of directions, and
  J²r⁻ᵐ is Riemann integrable on BZ ∖ B_δ for every δ > 0." This replaces "J₀ and r₀ continuous" and covers acoustic
  axes. State §1.13 as a lemma (Supplemental Sec. S4), with the sharp constants. Define τ^RTA for α ≥ d through the
  cutoff. In (ii), define "uniformly in N₁", replace "|j|² → 𝓜₀" by "|j|² bounded", and name the uniform (H1) in the
  status note.
* Prop. 1: "(as N₁ → ∞ at fixed T: c_ν ≤ k_B and w_ν²k_{ν,i}² are bounded while Σ_νc_ν = 𝒱ρc and Σ_νw_ν²k_{ν,i}² grow
  like n)". Optionally add "𝒞_Cal|_H has the rates 1/τ_R (×3) and 1/τ_c, so (K3) holds."
* Lemma 2: replace the sharpness sentence by the wording of §4.6, and put the ς = 4 example in S2. Write "up to 4
  for four-phonon events".
* Theorem 2: "for every I ⊂ 𝓔 with |I| = n_H − 1 define ψ_I …".
* Eq. (6): add the external-source / k → 0 sentence of §6.
* Prop. 5: "min over γ with g^ref_γ > 0", and "τ > 0".
* Remark 5 and the text after Prop. 4: the wordings of §10.1–10.3.
* Prop. 3 proof: "some I* ∪ {γ} is a basis".

## 13. Candidates for machine checking

* **Lemma 2** (inequality 𝒞 ⪯ ς_max diag 𝒞, the variational bound, and the equality criterion of §4.6). It is short,
  central (it is the feasibility test of the SDP side) and finite-dimensional, and was already recommended by the
  source audit (§14 there). It is recommended for machine checking in Lean 4 + Mathlib (`Matrix.PosSemidef`,
  Cauchy–Schwarz on finite sums).
* The lattice lemma of §1.13 is the weakest new step but is costly to formalize (lattice Riemann sums). A written proof
  in the Supplement is the appropriate remedy.
* Theorem 1 and Eq. (24) are finite exact identities, checked here in exact arithmetic by routes independent of the
  manuscript's script. They are low priority for formalization.

## 14. Check files (`paper/manuscript/checks/audit/`, research Python, single thread, < 300 MB)

| file | content | result |
|---|---|---|
| `audit_thm1_walsh.py/.out.txt` | Ex. 1; Thm 1 (i)–(v) by char. polynomial, full resolvent over ℚ(p, η), regularized solve; involution; Remark 1 (series, 60-digit eigen-decomposition, time-domain quadrature); Remark 3 | 29/29 |
| `audit_callaway.py/.out.txt` | Prop. 1 on O_h meshes N = 5…33 with Bose weights; support of j_x | 8/8 |
| `audit_lemma2.py/.out.txt` | 𝒞 ⪯ ς diag 𝒞 exactly; variational bound; exact equality examples ς = 2, 3, 4; mixed networks strict but approached; 2 ↔ 2 not admissible | 8/8 |
| `audit_thm2_cone.py/.out.txt` | determinant lemma; Eq. (24) over all subsets; Eqs. (21)–(23); extension; Prop. 3 asymptotics; Remark 5; Prop. 4; LMI | 22/22 |
| `audit_prop5_drive.py/.out.txt` | Prop. 5 certified by root isolation (F < 1 included; converse); Eq. (6) k → 0 limit (O(k²)); Eq. (7) | 11/11 |
| `audit_ir_lattice.py/.out.txt/.json` | Thm 3(i) mesh laws against the sharp constants of §1.13 (d = 2, N₁ ≤ 3201, additive error, J₀ vanishing on a plane) | 4/4 |
| `audit_frontface.py/.out.txt` | Sec. IV by direct BVP solution; branch; gauge; Cor. 2; Eq. (40); entropy production | 6/6 |

Run each with `python -B <file> > <file>.out.txt` from `checks/audit/`.

## 15. Overall verdict

**The mathematical content of Sections II–IV survives. No conclusion stated in a theorem or proposition is false.**

* Theorems 1 and 2, Propositions 1–5, Lemma 2 (the inequality), Corollary 1, the extension paragraph, the LMI, both
  Examples and the Section IV results are valid proofs or derivations under their stated hypotheses. The new material
  was checked: the parity involution, the slow-pole remark, the variational proof of Lemma 2, the
  matrix-determinant-lemma derivation of adj(A), the grey Callaway operator, the drive (6) and Prop. 5 with F > 0.
* **Theorem 3(i)** is correct for the integral criterion. The dyadic-shell proof of the mesh laws is correct in
  structure and even yields sharp constants, which the lattice sums confirm. Its convergence statements need a
  regularity hypothesis that (H1)–(H3) do not supply, and the added continuity of J₀ and r₀ is stronger than needed
  and fails at acoustic axes. With (H5) and the lemma of §1.13 written out, the label PROVED is justified.
* **Lemma 2's sharpness** is supported in the paper and its sources for ς_max = 2 and 3 only. The case ς_max = 4 is
  supplied here by an exact 1 ↔ 3 network, together with the criterion that equality requires every active event to
  touch ς_max modes.
* One **false factual statement** (Remark 5, "full support") and two **misattributed or unproved logical claims**
  (§10.2, §10.3) occur outside the formal statements. Each has a one-sentence repair that keeps the conclusion.
* The status labels and novelty wording are faithful to the sources and to `PRIOR_ART.md`, subject to the three label
  qualifications of §9.
