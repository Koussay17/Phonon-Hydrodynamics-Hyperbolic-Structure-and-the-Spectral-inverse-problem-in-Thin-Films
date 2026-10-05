# Proof audit — memory-gap campaign (20261004-memory-gap)

**Date:** 4 October 2026. **Object:** every statement labelled PROVED, PROVED here, PROVED (standard),
PROVED (trivial) or DERIVED UNDER ASSUMPTIONS in `13-final-report.md`, with the definitions of
`02-assumptions.md`. I also audited the numerical claims that are stated alongside these proofs or used to
interpret them (`07-experiments.md`, `09-failed-approaches.md`, `results/`), the code that implements the
definitions (`scripts/memgap.py`, `bt_exact.py`, `debye_events.py`, `aln_events.py`, `aln_geometry.py`,
`symmetry.py`, `fastops.py`, `optim_multi.py`, `certify.py`, `certify_aln.py`) and the Cauchy–Binet section
of note 18 (§4.1, thesis repository, read only).

**Method.** I re-derived each step independently. For each step I asked what licenses it, and I tested the
fragile steps in exact rational arithmetic or with rigorous float bounds. The check scripts and their raw
outputs are in this directory (index in §15). The thesis repository and the campaign scripts were not
modified.

**Classification legend.** VALID · VALID UNDER EXTRA ASSUMPTION (the assumption is stated) · INCOMPLETE
(the conclusion is true but the argument given does not reach it) · UNJUSTIFIED (no argument, or the
argument does not support the statement) · FALSE (counterexample given).

---

## 0. Verdict at a glance

| Statement (13-final-report) | Claimed status | Audit |
|---|---|---|
| T1 cone bound τ ≤ M²K | PROVED | **VALID** (exact identity check). It extends to non-admissible g with b ∈ range C(g) (§2.8) |
| M = 16284.0067 for 1D N = 9 ("exact", E1) | numerical, exact | The original singularity test was heuristic (float cond cut). **Value certified here**: integer-rank + time-reversal certificate over all 9.66e6 subsets, M = 16284.006723562 (§2.9) |
| T2 sharpness without lifetimes | PROVED here | **VALID** (sup approached, not attained; rates on I* → ∞). The "new" label is not supported |
| "cone bound is sharp **only** without the lifetime constraint" (§1) | asserted | **FALSE** as a general statement (exact counterexample, §3.4). Open for physical b |
| P3 rational-kernel form of M | PROVED | **VALID** (exact check on all 219 bases of a test geometry) |
| L3 K ≥ K_RTA/3 | PROVED | **VALID**, including repeated-daughter events. **Constant 3 is sharp**: attained exactly (§5.3) |
| P4 lower bounds on τ | PROVED | **VALID** |
| P5 convexity / LMI / SDP | PROVED (standard) | Convexity **VALID**. The set identity is **FALSE as written**, **VALID** with an extended-value K (§7.3). The remark on τ is not a proof |
| P6 gap ≤ M²\|b\|² | PROVED | **VALID** |
| Dimension count 2508 / 43,503 | PROVED (trivial) | **VALID UNDER** full-rank condition, **verified here at the reference**. The growth law m ~ N^{2d−1} is **FALSE** for the exact 1D model |
| T7(i) RTA moment criterion | PROVED | **VALID**. The mesh law at α = d needs a log (N^d/log N) |
| T7(ii) τ → ∞ for 2α ≥ d (event operator) | DERIVED UNDER ASSUMPTIONS | Divergence **VALID UNDER H4** after a one-line repair. Mesh exponents **UNJUSTIFIED**. AlN sentence **UNJUSTIFIED** |
| T7(iii) M²K₀ diverges | derived | **VALID as a conditional**: it inherits H4 |
| "M²K₀ ≥ 1e8 τ_ref in every geometry computed" (§1, E2) | numerical | **FALSE** on the campaign's own data (1D N = 15: 7.88e6; N = 21: 3.15e7) |
| "τ_min within a factor ≤ 1.65 of τ_CS in every geometry tested" (§1) | numerical | **FALSE** on the stored scans (2D N = 5 with full symmetry: 1.684) |

**Earliest failure** in reading order: §1, first bullet of `13-final-report.md`. Its numeric range "0.001–65 %"
(repeated as "≤ 1.65" in the Paper-1 sentence) is exceeded on the campaign's own data (§11.0). The earliest
failure of a **mathematical** claim is §1, third bullet: "sharp only without the lifetime constraint" (§3.4). The
same bullet also carries the wrong "≥ 1e8" number (§11.1). **Earliest defective step inside a proof:** the LMI set
identity of P5 (§7.3). It is definitional and harmless. The first substantive gap is T7(ii) (§10.2).
**The proved core (T1, T2, P3, L3, P4, P6) is sound as stated.**

---

## 1. Conventions fixed for the audit

* H = e^⊥ ⊂ R^n with e = Dε. Event vectors a_α = D^{-1}s_α ∈ H. C(g) = Σ g_α a_α a_αᵀ with g ≥ 0.
  Admissible means ker C(g) = span(e).
* b ∈ H \ {0}. For admissible g, x = C⁺b is the unique solution of Cx = b in H. Also K = bᵀx, N = |x|²,
  τ = N/K = −K'(0)/K(0).
* **Extended value** (used where needed): K(g) := bᵀC(g)⁺b if b ∈ range C(g), and +∞ otherwise, for every g ≥ 0.
  This is sup_y {2bᵀy − yᵀC(g)y}, a supremum of affine functions of g. It is therefore convex and lower
  semicontinuous on all of R^m_{≥0}.
* s(α) := |supp a_α| is the number of **distinct** modes touched by event α: 3 for p → a + b with a ≠ b, and
  2 for the repeated daughter p → 2a (stoichiometry −e_p + 2e_a).
* Basic solutions: for |I| = d − 1 = n − 2, y_I is the solution of [a_Iᵀ; eᵀ; bᵀ] y = (0, …, 0, 1) when the
  matrix is nonsingular. M := max_I |y_I|. Equivalently, M = max |v_I| / |v_Iᵀb| over I with v_Iᵀb ≠ 0, where
  v_Iᵀy = det_H[a_I, y].

---

## 2. T1 — cone bound τ ≤ M²K (Ben-Tal–Teboulle; Cauchy–Binet proof of note 18 §4.1)

| # | Step | Classification | Reason |
|---|---|---|---|
| 1.1 | For admissible g and b ⟂ e: C⁺b = (C\|_H)^{-1}b | VALID | C maps R^n into H and C\|_H ≻ 0 iff ker C = span(e). The code's regularized solve (C + γêêᵀ)x = b returns the same x: applying êᵀ gives γ êᵀx = êᵀb − (Cê)ᵀx = 0, hence Cx = b with x ∈ H. With b ⟂ e only to 1e-16 in floats, the error is (êᵀb)/γ, which is negligible |
| 1.2 | adj_H(A) = Σ_I c_I v_I v_Iᵀ with c_I = Π_{α∈I} g_α | VALID | Cauchy–Binet for the (d−1)-th compound. Dependent subsets give v_I = 0 |
| 1.3 | bᵀadj(A)b = det(A)K > 0 | VALID | Needs b ≠ 0 and A ≻ 0 on H |
| 1.4 | x/K = Σ θ_I y_I with θ_I = c_I(v_Iᵀb)² / (det A · K), θ ≥ 0, Σθ = 1 | VALID | The sum runs over I with c_I > 0 and v_Iᵀb ≠ 0, so the convex hull is over **active** bases only, which is sharper |
| 1.5 | v_Iᵀb ≠ 0 ⇔ a_I independent and b ∉ span a_I ⇔ [a_Iᵀ; eᵀ; bᵀ] nonsingular, and then y_I = v_I/(v_Iᵀb) | VALID | Exact check: `A_basis_iff_VIb_nonzero` and `A_yI_equals_VI_over_VIb` are both true |
| 1.6 | \|x\| ≤ MK, τ ≤ M²K | VALID | Convexity of the norm |

**Hypotheses actually used:** g admissible, and b ∈ H \ {0}. The spanning of H by the full event set is needed
only for admissible g to exist. If the events span only V ⊊ H, the statement must be made on V, with b ∈ V and
(dim V − 1)-subsets. Every geometry used has kernel dimension 1, so this does not arise. Rank-deficient
subsets are correctly excluded: they contribute zero to both the numerator and the denominator. There is no
pseudo-inverse subtlety for admissible g. The bound holds for symmetric (orbit-tied) rates as well, as a
special case. d = 1 is covered by the empty subset.

**Exact verification** (`check_cone_bounds.py`, test geometry n = 7, m = 10 with repeated daughters,
219 bases, M² = 1142947629646/111419102025). The identity x/K = Σ_I c_I V_I(V_Iᵀb)/Σ_I c_I(V_Iᵀb)² holds
**exactly** for 8 rate vectors, including 3 with two zero rates (48 active bases). In all cases τ ≤ M²K.

**2.8 Extension (not claimed in the report, but used implicitly by the SDP side).** If g ≥ 0 is not admissible
but b ∈ range C(g), then |C(g)⁺b| ≤ M·K(g). *Proof:* g_θ = g + θ𝟙 is admissible for θ > 0. Write
E = C(𝟙) ≻ 0 on H, and let P and Q be the projectors onto range C(g) and onto ker C(g) ∩ H. The two blocks
of C(g_θ)x_θ = b give u_θ = −(QEQ)^{-1}QEP w_θ, and w_θ → C(g)⁺b =: w. Hence x_θ → w + u_∞ with u_∞ ⟂ w and
K(g_θ) → K(g). Then |w| ≤ |w + u_∞| ≤ M·K(g), by T1 and continuity. VALID.

**2.9 The value M = 16284.0067 for 1D N = 9 (E1).** *Gap in the original argument:* subsets are declared
singular when cond([a_Iᵀ; bᵀ]) ≥ 1e11. In float64, an exactly singular subset and a nonsingular subset with
|y_I| ≳ 1e15 are indistinguishable. The clean histogram gap (no log₁₀ cond in [8, 16)) is suggestive, but it
does not exclude such a subset. The rational-kernel recomputation re-evaluates only the maximizer I*, not the
maximality. *Independent certification* (`certify_M_1d9.py`, 608 s) classifies all 9,657,700 subsets
rigorously:

1. **Integer rank.** S_I has entries 0, ±1, 2. If rank S_I = 14, then det(S_IᵀS_I) is an integer ≥ 1, which
   forces σ₁₄ ≥ L_I := 1/Π_{i≤13}σ_i.
   * 6,082,937 subsets have computed σ₁₄ ≤ 1.2e-15 with σ₁₄/L_I ≤ 5.3e-12, so their rank is < 14.
   * The other 3,574,763 have σ₁₄ ≥ 0.062, so their rank is 14.
2. **Time-reversal certificate.** T (k → −k) fixes ε and D, and b is exactly odd. So bᵀDφ_I = 0 whenever
   Tφ_I = φ_I. 8,430 rank-14 subsets have |Tφ̂ − φ̂| ≤ 8.4e-15. By the Hadamard bound, a non-even integer
   kernel vector would give ≥ 1/H_max = 1.8e-7, so these subsets are exactly even and **exactly singular**.
   They are precisely the "rank 14 but bᵀDφ = 0" subsets that a condition-number cut cannot certify.
3. **Bases.** The remaining 3,566,333 subsets (exactly the campaign's count of bases) have |Tφ̂ − φ̂| ≥ 0.66 and
   relative |bᵀDφ_I| ≥ 3.4e-6. All their |y_I| are therefore finite and well determined, with
   log₁₀|y_I| ∈ [−2, 5).
4. **Maximum.** The 40 largest values lie on one hyperplane, with
   φ = (−932, 2001, 1334, 667, 167, 334, −399, 668, −1064, −798, −532, −266, 134, 268, 402, 536) ⟂ ε.
   Recomputed exactly (integer kernel, 60-digit Bose scales, energies exactly 2πk/N):
   **M = 16284.00672356239098…**. The campaign's 50-digit value 16284.00672356293 used float-rounded energies,
   a 3e-14 relative difference. The next hyperplane gives 6725.96.

**E1 is confirmed.** With this certificate, the label "exact" is justified up to the stated float error
bounds, which lie many orders of magnitude inside the margins. Without it, the label rested on a heuristic
gap in a condition-number histogram.

---

## 3. T2 — sharpness without the lifetime constraint

| # | Step | Classification | Reason |
|---|---|---|---|
| 2.1 | I* attaining M exists | VALID | Finite set |
| 2.2 | g^ε = 1 on I*, ε elsewhere, is admissible for ε > 0 | VALID UNDER the spanning hypothesis | All events active. Spanning is checked numerically for each geometry |
| 2.3 | c_J = ε^{\|J \ I*\|}, and every J ≠ I* has \|J \ I*\| ≥ 1 | VALID | Equal cardinalities |
| 2.4 | 1 − θ_{I*} = O(ε) and x/K → y_{I*} | VALID | v_{I*}ᵀb ≠ 0, and the y_J are finitely many and bounded |
| 2.5 | Rescaling to K = K₀ (K and τ of degree −1, \|x\|/K of degree 0) | VALID | C(λg) = λC(g) |
| 2.6 | sup = M²K₀ (upper bound from T1) | VALID | |

**3.1 Attained or approached?** Approached, and in general not attained. det A = O(ε) while det A · K → (v_{I*}ᵀb)²,
so K(g^ε) ~ c/ε. The rescaled rates λg^ε are therefore ~1/ε on I* (→ ∞) and O(1) on the complement. The
report's parenthesis "rates on the complement of I* must vanish" is true only relative to the I* rates. At
fixed K₀, the complement rates tend to the finite value c/K₀ (K·ε → 0.0193 in the example). The
sequence is unbounded, its diagonal diag C → ∞ on the modes touched by I*, and its formal limit (support I*)
is not admissible: d − 1 events cannot span H. Every member of the sequence is admissible, so the
**supremum** statement is correct. Exact illustration (test geometry): the deficit M² − τ/K equals 10.25,
9.93, 5.46, 0.885, 0.0941, 0.00946 for ε = 1e-1, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8. The convergence is O(ε) only
for ε ≲ 1e-6 (deficit/ε ≈ 9.5e5). In this 7-mode example, coming within 1 % of M² needs a rate ratio of
about 1e7 between I* and the other events.

**3.2 Fixing K but not diag C.** Essential. The approaching sequence drives the lifetimes of the modes in I* to
zero, so T2 says nothing about F(r, K₀).

**3.3 Novelty.** The argument is the standard weight-concentration argument for weighted least squares.
x/K = argmin{yᵀC(g)y : bᵀy = 1, y ∈ H} is a diagonally weighted least-squares solution, and the basic
solutions are the vertices of the Ben-Tal–Teboulle hull. Closely related: the geometry of the set of scaled
projections (Hanke & Neumann, Linear Algebra Appl. 190 (1993) 137, which proves it is a finite union of
polytopes) and the D-independent bounds for scaled pseudoinverses (Stewart, LAA 112 (1989) 189; O'Leary,
LAA 132 (1990) 115, doi:10.1016/0024-3795(90)90056-I, which shows the bounds are computable). I have not
established whether the fixed-right-hand-side corollary is stated verbatim there. **"New" is unsupported,
not refuted.** A literature check is required before claiming it.

**3.4 The remark "the cone bound is sharp only without the lifetime constraint" (13 §1).**
This is a negative claim (with lifetimes fixed, sup_F τ < M²K₀), and no proof is given. **It is FALSE in the
stated generality.**

*Exact counterexample* (`check_cone_bounds.py`, instance B). Take b = a_β, a current parallel to one event
vector (b ∈ H). Every subset containing β has b ∈ span a_I, so M is attained by some J* ∌ β. B = J* ∪ {β} is a
basis of H. For **any** g supported on B, adj(A)b = c_{J*}v_{J*}(v_{J*}ᵀb), so x/K = y_{J*} exactly and
τ = M²K. With r := diag C(g) and K₀ := K(g), F(r, K₀) contains g, so the lifetime-constrained supremum equals
M²K₀ and is **attained**. Numbers: r = (915/16, 7, 9/2, 29/9, 17/16, 2891/640, 153/80), K₀ = 4,
τ = M²K₀ = 44367253/1074614 exactly. The equality holds on the whole face (4 further rate vectors).

*Scope.* A physical current b = D(vε) has full support, so it is never parallel to a 3-mode event vector. For
such b the question is **open**. The generic indication is a strict and huge gap: on a 6-mode, 8-event instance
with P(r) a polygon (instance C), a dense scan of F(r, K₀) gives max τ = 0.503 against M²K₀ = 2.58e4
(ratio 2e-5). This is numerical and illustrative only. The obstruction to equality is visible in the proof:
an admissible limit point gives strict inequality unless b is parallel to an event vector. A non-admissible
limit point requires b ⟂ the extra invariant, and then the u-component stays bounded (§2.8). The simple
configurations I analysed again need b in the span of a single event vector. A general proof is missing.

*Minimal repair:* "Sharpness is proved without the lifetime constraint. With lifetimes fixed, sup_F τ can
equal M²K₀ in degenerate cases (b parallel to an event vector). For physical currents, whether it is strictly
smaller is open."

*Symmetric class.* g^ε is G-invariant only if I* is a union of orbits, so T2 does not transfer to F_G. For
symmetric rates (all AlN headline numbers), the relevant supremum is ≤ M²K₀, and equality is not shown.

---

## 4. P3 — rational-kernel form of M

| # | Step | Classification | Reason |
|---|---|---|---|
| 3.1 | a_I independent ⇔ S_I full column rank | VALID | D invertible |
| 3.2 | ε ∈ ker S_Iᵀ (conservation), dim ker S_Iᵀ = 2, so ker = span{ε, φ_I} | VALID | n − (n − 2) = 2 |
| 3.3 | y = Dφ with eᵀy = 0 and bᵀy = 1 gives y_I = P_H Dφ_I / (bᵀDφ_I) | VALID | Uses b ⟂ e. Independent of the representative φ_I + cε |
| 3.4 | φ_I rational | VALID UNDER integer (rational) stoichiometry | This is the 1D exact model only. For the d ≥ 2 Debye and AlN events, l_p = √((ε_a+ε_b)/ε_p) is irrational and φ_I is real. The formula itself still holds |

Exact check: the formula agrees with the direct solve for **all 219 bases** of the test geometry
(`A_P3_all_bases_exact = true`). The interpretation (M measures the weakest current coupling to an extra
invariant of a (d−1)-event sub-network) is correct.

---

## 5. L3 — K ≥ K_RTA/s (three-phonon: κ ≥ κ_RTA/3)

| # | Step | Classification | Reason |
|---|---|---|---|
| L.1 | b = Cx, i.e. b_μ = Σ_α g_α a_{α,μ}(a_αᵀx) | VALID | Needs b ∈ range C: admissible g with b ⟂ e, or the extended class |
| L.2 | b_μ² ≤ r_μ D_μ with D_μ = Σ_{α∋μ} g_α(a_αᵀx)² | VALID | Cauchy–Schwarz with weights g_α a_{α,μ}². The **coefficients** a_{α,μ} are absorbed into r_μ, so unequal entries such as (−1, 2) are harmless |
| L.3 | r_μ > 0 for all μ | VALID (admissible g) | If r_μ = 0, the unit vector e_μ ∈ ker C. But ker C = span(e) and e = Dε has full support (ε > 0, Γ zero modes excluded). So there is no D^{-1}, 0/0 or zero-mode issue |
| L.4 | Σ_μ D_μ = Σ_α s(α) g_α(a_αᵀx)² ≤ s_max · xᵀCx = s_max K | VALID | Double counting over the **distinct** modes of each event |
| L.5 | K ≥ K_RTA(r)/s_max | VALID | |

**5.1 Repeated labels.** The worry ("two identical daughters → incidence (1, −2), entries not bounded the same
way") does not apply to this proof. The matrix form is zzᵀ ⪯ |supp z|·diag(z²), which holds for every z by
Cauchy–Schwarz. For z = (−1, 2): 2diag(z²) − zzᵀ = [[1, 2], [2, 4]] ⪰ 0 (singular). A repeated-daughter event
has support 2 and satisfies the bound with constant 2 ≤ 3. Exact random tests (300 admissible networks with
repeated daughters, random current b ⟂ e): no violation, min K/K_RTA = 0.669 (`check_rta_bound.py`).

**5.2 Degenerate modes.** The inequality refers to the diagonal in the basis in which the events have support
≤ 3. A rotation inside a degenerate subspace changes both the supports and diag C. The campaign's operator
has a fixed band gauge, with orbit-averaged rates, so L3 applies to it as stated. A κ_RTA computed in a
different gauge or with different degeneracy averaging is a different number.

**5.3 The correct constant.** The sharp constant is **s_max = max over active events of the number of distinct
modes**: 3 for three-phonon events, 2 if every active event has a repeated daughter, 4 for four-phonon events.
It is independent of the stoichiometric coefficients. **It is attained exactly**, not just approached:

* *Equality condition* (derived in this audit): suppose there is ψ with s_{α,μ}ψ_μ = c_α for every μ ∈ supp α
  and every active α, and all active supports have size s. Then with y = Dψ, Cy = sRy. So b := sRy gives
  b ⟂ e, K = s·yᵀRy and K_RTA = s²·yᵀRy, i.e. K = K_RTA/s.
* For s = 3 this requires a "bipartite" sign structure: in each event both daughters carry ψ = c and the parent
  carries −c. Writing φ = σ ⊙ φ', the invariants of such a network are the solutions of the **zero-sum** system
  φ'_p + φ'_a + φ'_b = 0 on its triples. So an admissible example is exactly a 3-uniform hypergraph whose
  zero-sum space is spanned by a nowhere-zero integer w (energies |w|, signs sign w).
* *Exact example:* 11 modes with energies (3, 5, 3, 2, 1, 1, 2, 3, 4, 6, 2) and 12 distinct-daughter events.
  The network is admissible (rank 10 = n − 1), conserves energy and has b ⟂ e. It gives K = 28671/28,
  K_RTA = 86013/28, **K/K_RTA = 1/3 exactly**, and again 1/3 for 5 further rate vectors.
* A network made only of repeated-daughter events (a chain with energies 1, 2, 4, 8, 16) gives
  **K/K_RTA = 1/2 exactly**, so the constant 2 is sharp for that sub-class.

So "κ ≥ κ_RTA/3" is correct and **cannot be improved** in the class "three-phonon events, arbitrary current
b ⟂ e". With the physical current the attainable ratio is much larger (AlN: K₀/K_RTA = 1.257 (x),
1.432 (z), reproduced here). A sharper physical constant would need constraints on b, such as the time-reversal
parity and the velocity structure.

Novelty was not checked by the campaign. The inequality is the hypergraph analogue of the graph-Laplacian
bound L ⪯ 2·Deg. This short, central lemma is a good candidate for machine checking (§14).

---

## 6. P4 — lower bounds on the memory

| # | Step | Classification | Reason |
|---|---|---|---|
| 4a | K = bᵀx ≤ \|b\|\|x\|, so τ ≥ K/\|b\|² | VALID | Needs K > 0. It also holds with pseudo-inverses. Equality iff x ∥ b, i.e. b is an eigenvector of C\|_H |
| 4b.1 | N = \|x − Cv\|² + 2bᵀv − \|Cv\|² ≥ 2bᵀv − \|Cv\|² | VALID | Needs Cx = b |
| 4b.2 | \|Ce_μ\| ≤ Σ g\|a_α\|\|a_{α,μ}\| ≤ ρ_μ r_μ | VALID | ρ_μ is a max over all events containing μ, which is conservative |
| 4b.3 | Optimizing over t gives N ≥ b_μ²/(ρ_μ r_μ)² | VALID | |

Exact checks: P4a and P4b hold on 6 random rate vectors (τ/τ_CS = 1.05–1.16). Labelling: E5 calls 4b "Lemma
L4 in 13"; the report calls it P4(b). This is cosmetic.

---

## 7. P5 — convex structure, Schur complement, SDP

| # | Step | Classification | Reason |
|---|---|---|---|
| 5.1 | g ↦ K(g) convex | VALID | Matrix-fractional function composed with an affine map on {C\|_H ≻ 0}, which is convex. Or: a supremum of affine functions (extended value) |
| 5.2 | K nonincreasing | VALID | ∂K/∂g_α = −(a_αᵀx)² ≤ 0 |
| 5.3 | {g ∈ P(r): K(g) ≤ K₀} = {g ∈ P(r): [[C(g), b], [bᵀ, K₀]] ⪰ 0} | **FALSE as written**; VALID with the extended-value K | The generalized Schur lemma gives LMI ⇔ b ∈ range C(g) and bᵀC⁺b ≤ K₀. That admits **non-admissible** g, while 02-A defines K only for admissible g |
| 5.4 | K_lo(r) = min_{P(r)} K is an SDP | VALID with extended K | The minimum is attained (P(r) compact, K lsc). It is the infimum over admissible g, provided P(r) contains an admissible point: along g₀ → g₁ convexity gives K(g_θ) → K(g₀). The minimizer may itself be non-admissible: say "inf over admissible g" |
| 5.5 | F nonempty ⇒ K₀ ≥ max(K_lo, K_RTA/3) | VALID | The max is redundant: K_lo ≥ K_RTA/3 by L3 applied on the LMI set |
| 5.6 | "τ is not convex on F (indefinite term …)" | UNJUSTIFIED remark | F is a level set and not convex, so "convex on F" is undefined. The displayed term is not derived. Not load-bearing |

*Exact counterexample to 5.3 as written* (`check_lmi.py`). Two 4-mode energy-conserving sub-networks, bridged by
4 events, with b ⟂ e_A and b ⟂ e_B. Then g₀ (bridges off) has rank C = 6 < 7, so it is not admissible.
Yet K(g₀) = 406.29 is finite, and the LMI matrix at t = K(g₀) equals [I; x₀ᵀ]C[I, x₀], which is PSD
**exactly**. r := diag C(g₀) admits an admissible g₁ ∈ P(r). Along the segment, K(g_θ) − K(g₀) ≈ 0.335θ → 0
(θ = 1e-1 … 1e-5).
Such points are physically relevant: the AlN x-max witness approaches a hidden invariant with current
overlap 6.6e-19 (E9).

*Minimal repair:* define K with extended values on P(r) and replace "min" by "inf over admissible g". No
numerical result changes.

---

## 8. P6 — τ_max/τ_min ≤ M²|b|²

VALID: it combines T1 (on F, τ ≤ M²K₀) with P4a (τ ≥ K₀/|b|²), and the ratio is scale invariant. The
1D N = 9 value is reproduced: M²|b|² = 7.31e10, with |b|² = 275.64 and τ_CS = 63.64
(`check_bt_numbers.json`).

---

## 9. Dimension count (2508 symmetric, 43,503 unconstrained)

**Counting principle.** At a point of the **open** orthant where the constraint map g ↦ (Wg, K_cd(g)) has a
Jacobian of rank ρ, F is locally a manifold of dimension m_var − ρ (implicit function theorem). W does not
depend on g. So "m − n − n_K" is exact iff W has full row rank and the K-gradients are independent modulo
row(W). Otherwise the true dimension is **larger**. The report's "when the constraints are independent" names
the condition, but the condition was not checked. "PROVED (trivial)" is therefore premature, but the
principle is VALID UNDER the rank condition.

**Verified here at the reference rates** (`check_dimension_count.py`, AlN 5×5×3):

| map | shape | rank | σ_min/σ_max | gap after the W rows |
|---|---|---|---|---|
| W | 897 × 44406 | 897 | 1.96e-2 | — |
| [W; ∇K_ab] (6 components) | 903 × 44406 | 903 | 9.0e-6 | 1.3e3 |
| W_γ (orbits) | 117 × 2627 | 117 | 2.76e-2 | — |
| [W_γ; ∇K_xx, ∇K_zz] | 119 × 2627 | 119 | 5.9e-6 | 3.5e3 |

So dim F = 43,503 and dim F_G = 2,508 near the reference. Both numbers are confirmed. The reference has all
rates > 0, and its diagonal r is constant within each mode orbit to 2e-14, which validates the one-row-per-orbit
reduction. The K-gradients are only **weakly**
independent of the lifetime rows (relative singular value ~6e-6): near the reference, fixing the lifetimes
almost fixes K. This is consistent with L3, and it explains why the extremal witnesses must move far from the
reference.

**Growth law** "m grows like N^{2d−1}". It is correct for the tolerance events, where m/N³ ≈ 1.6–1.8 (2D) and
m/N⁵ ≈ 1.2–1.3 (3D) from the E3 table. It is **FALSE for the exact 1D model**: m = 26, 124, 294, 1032 for
N = 9, 21, 33, 63, i.e. m ≈ 0.26 N², not ~N, because collinear normal processes are resonant for all same-sign
pairs (09, item 7). The conclusion that the lifetimes fix a vanishing fraction of the rates survives.

---

## 10. T7 — infrared criterion

**10.1 Part (i), RTA moments: VALID**, with two precisions.

* "b² → b₀² > 0" is literally false. In entropy coordinates b → kT·v, so b₀(Ω)² = (kTc)²Ω_x², which vanishes on
  the plane Ω_x = 0. The criterion needs only ∫ b₀(Ω)² a(Ω)^{-m} dΩ ∈ (0, ∞), which holds when a is bounded
  above and below. Also needed: r bounded below away from q = 0, on every branch.
* The measure q^{d−1}dq dΩ is correct. The exponents are correct: M_m < ∞ ⇔ mα < d, with a log at equality.
* Mesh laws, checked with exact lattice sums in d = 3 up to N = 161 (`check_ir.py`):
  * α = 1: τ converges (0.486 → 0.5095).
  * α = 3/2: τ/log N = 0.188–0.191.
  * α = 2: local exponent 1.004 → 1.001.
  * α = 3.5: exponent 3.37 → 3.44 (→ 3.5 slowly).
  * **α = d = 3:** τ·log N/N³ = 0.00288 → 0.00280. So τ_RTA ~ N^d/log N, **not** N^α as stated.
  * E7's fitted exponents (1.29 for α = 2 at N ≤ 11) are pre-asymptotic, as the campaign says.

**10.2 Part (ii), event operator under H4 (|x_μ| ≥ c|b_μ|/r_μ for |q| < q₀, uniformly in the mesh).**

| # | Step | Classification | Reason |
|---|---|---|---|
| 7.1 | N ≥ c²·Σ_{small q} b²/r², so N → ∞ for 2α ≥ d | VALID UNDER H4 (uniform c) | |
| 7.2 | "τ = ∞" for d/2 ≤ α < d | **INCOMPLETE** | τ = N/K → ∞ needs K bounded above. L3 gives only a **lower** bound on K, and no upper bound is proved for the event operator |
| 7.3 | α ≥ d: "both diverge but τ still diverges, ~ N^α" | **INCOMPLETE / UNJUSTIFIED** | No argument is given. N^α is an RTA computation, and even for RTA it fails at α = d (log) |
| 7.4 | Mesh exponents τ ~ N^{2α−d}, log N for the event operator | **UNJUSTIFIED** (two-sided) | H4 is one-sided, so it gives only lower bounds |

*One-line repair of 7.2–7.3:* K = bᵀx ≤ |b||x| gives **τ ≥ |x|/|b| = √N/|b|**. Then τ → ∞ whenever N → ∞
and |b| stays bounded, with no condition on K, for all α with 2α ≥ d. Under H4 this gives the rigorous mesh
lower bound τ ≳ N^{α − d/2} for 2α > d (√log N at 2α = d). It gives τ ≳ N^{2α−d} if additionally K stays
bounded (α < d). Two-sided laws need
an upper RTA-likeness hypothesis as well. The inequality was checked on AlN: τ/(|x|/|b|) = 1.336 (reference),
1.0002 (x-min), 18.97 (x-max 6.9e3 ps).

*H4 is a property of an operator, not of the class.* On the AlN mesh the lowest x-coupled modes lie at
3.70–4.35 THz (8 modes), and there is nothing below 3.7 THz. There, ρ_μ = x_μ r_μ/b_μ is:

* 1.13–1.21 at the reference (RTA-like);
* **0.294** at the x-minimizer, with the low-mode share of |x|² falling from 0.53 to 0.061;
* 0.263 at the z-minimizer;
* 0.87–23.5 at the 6.9e3 ps maximizer.

Whether c stays bounded away from 0 under refinement, for the reference or for any member of F, is
untested.

**10.3 The AlN sentence** ("the physical AlN case (α ~ 2, d = 3) is in the divergent regime: τ_mem of AlN is not a
mesh-converged constant") is **UNJUSTIFIED as an unconditional statement**, for three reasons:

1. α ≈ 2 is fitted for the 1D Debye reference model (02-B). It is not measured for AlN in this campaign, and
   note 18 states that the AlN exponent is unresolved.
2. H4 is unverified for AlN in the infrared, because the 5×5×3 mesh has no infrared.
3. "34 ps (event operator, 5×5×3) vs 149/232 ps (RTA, 24³)" compares different operators on different meshes.

*Repair:* make it conditional: "if α ≥ 3/2 and H4 holds uniformly …".

**10.4 Part (iii), "M²K₀ diverges": VALID as a conditional.** τ_ref ≤ M²K₀ holds on every mesh (T1), so M²K₀
diverges whenever τ_ref does. That divergence is derived only under H4. The scaling remarks are correct:

* b → sb gives M → M/s, so M²K₀ and M²|b|² are b-scale invariant.
* M depends only on the hyperplanes spanned by the a_I, not on their normalization.
* M²K₀ scales like τ under g → λg.

A possible unconditional route, which I have **not** proved, would use the lifetime-free class: construct
rates that are RTA-like on the small-q modes, then use M²K₀ ≥ (N(g)/K(g)²)·K₀ together with K₀ ≥ K_RTA/3.

---

## 11. Numerical statements attached to the proofs

0. **"The minimum over the feasible set exceeds τ_CS by only 0.001–65 %" and "within a factor ≤ 1.65 in every
   geometry tested" (13 §1, bullet 1 and the Paper-1 sentence; §2 "lie within 1.0003–1.65"): FALSE on the
   campaign's own `results/scan_summary.json`.** τ_min(found)/τ_CS = **1.684 for 2D N = 5 with full cubic
   symmetry** and 1.652 for 1D N = 15 (time reversal). E5 quotes 1.10 for 2D N = 5, but that is the untied or
   time-reversal case; with full symmetry the value is 1.684. Likewise "1.012 (3D N = 3)" is the untied case;
   with full symmetry it is 1.206. The proved part (τ ≥ τ_CS) is unaffected. Correct summary: "within a factor
   ≤ 1.69 in every configuration tested". Over all 32 stored configurations the maximum is 1.6839.
1. **"M²K₀ ≥ 1e8 τ_ref in every geometry computed" (13 §1, E2): FALSE on the campaign's own lower bounds.**
   M²K₀/τ_ref ≥ 2.75e10 (1D N=9), **7.88e6 (1D N=15)**, **3.15e7 (1D N=21)**, 3.4e11 (1D N=27), 7.1e10
   (2D N=3), 1.8e13 (2D N=5), 1.3e13 (3D N=3) (`check_bt_numbers.json`). Correct statement: "≥ 7.8e6
   everywhere, ≥ 2.7e10 except 1D N = 15, 21". The vacuity conclusion survives.
2. **E2 "certified lower bounds".** The json stores neither the subsets nor the conditioning. A re-run with
   another seed (`check_pivot_bases.py`) shows the search is strongly seed-dependent: 1D N=9 → 6725.96 instead
   of 16284, and 3D N=3 → **1.77e6**, larger than the reported 2.05e5. The bases themselves are genuine:
   cond(B_I) = 2e5–3e9, and 50-digit re-solves agree to ≤ 1.2e-8. *Repair:* store I and cond(B_I) and
   re-verify in high precision. The 1D N=15/21 bounds (240, 412) are evidently weak.
3. **Lower bounds rounded upward.** "≥ 3.6e2" (x): 6906.75/19.1969 = 359.8. "≥ 6.2e4" (z): 61,848.
   "≥ 9.6e5 ps" (z): 959,690. 2D N=9 "3.2e3": 3188. A certified lower bound must be rounded **down**.
4. **"τ_max(lifetimes fixed) ≥ 627.5" (1D N=9).** This is the untied scan maximum (`scan_summary.json`,
   d1_N9_none_phys). It is not among the 12 witnesses re-verified independently in E4, and the E3 table
   lists the time-reversal-tied value 509.8. Minor.
5. **Inner bounds vs exact data.** A witness with residual 1e-12 is exactly feasible only for perturbed
   data. Transfer to the exact data needs an implicit-function step. I checked it for the robust AlN witnesses
   (`check_witness_ift.py`):
   * x-max 6906.75 ps: residual 1.6e-14, full-row-rank Jacobian (119 rows, σ_min = 0.0102), Newton correction
     ≤ 1.6e-12 in log-rates, first-order δτ/τ ≤ 6e-13;
   * x-min 19.197 ps: correction ≤ 2.6e-11, δτ/τ ≤ 4e-12.

   VALID UNDER this check, which should be added for the other headline witnesses.
5b. **"K_lo by SDP and by mirror descent agree to 1e-6" (P5/E6).** `results/K_lower_sdp.json` gives relative
   differences of 1.0e-7 (1D N=9), 2.3e-7 (2D N=3), 2.0e-6 (2D N=5), 3.3e-6 (3D N=3) and 4.1e-3 (1D N=15,
   SDP status "optimal_inaccurate"). The statement should read "to ≤ 3.3e-6 where the SDP converged". Minor.
6. **3.58e7 ps witness.** The arb certification builds its balls from the **float** stoichiometric coefficients
   l_p, l_d. The certified operator therefore conserves energy only to about 1e-16, while λ_min(H) = 2.5e-15
   (cond 6.4e13). The statement "attained for data within 1.1e-11" should also say "for an event geometry whose
   energy conservation holds to float precision". The campaign already relies only on the 6.9e3 ps witness.

---

## 12. Earliest failures, in order

1. **13 §1, bullet 1** (also lines 28 and 103): "0.001–65 %" / "≤ 1.65" is a numeric summary contradicted by
   the stored scans (1.684, §11.0). The proved inequality τ ≥ τ_CS that it accompanies is valid.
2. **13 §1, bullet 3**: "sharp only without the lifetime constraint" (FALSE in general, §3.4) and "M²K₀ ≥ 1e8 τ_ref
   in every geometry" (FALSE on the stored data, §11.1). Neither is used by any proof.
3. **13 §1, bullet 4**: certified lower bounds rounded upward (§11.3).
4. **P5 set identity** (§7.3): definitional. It is repaired by the extended-value K, and no number changes.
5. **T7(ii)** (§10.2): the step to τ = ∞ is incomplete, which the repair τ ≥ √N/|b| fixes. The mesh exponents
   for the event operator are unjustified, and so is the unconditional AlN sentence (§10.3).
6. Cosmetic or numerical items: §9 (1D growth law), §11.2, §11.4–11.6.

None of these failures propagates into T1, T2, P3, L3, P4 or P6.

## 13. Minimal repairs (consolidated)

* State T1 and T2 with: admissible g; b ∈ H \ {0}; events spanning H. Add the extension §2.8 if the SDP
  side uses non-admissible points.
* T2: "supremum, approached by unbounded rates on I*; not attained in general". Drop "new" pending a
  literature check (Hanke–Neumann 1993 and the scaled-projection literature).
* Replace "sharp only without the lifetime constraint" by the §3.4 wording.
* L3: state s = number of distinct modes per event. Add sharpness (equality example) and the basis/gauge
  remark (§5.2).
* P5: use the extended-value K and "inf over admissible g". Delete or prove the τ-non-convexity remark.
* Dimension count: "local dimension at a regular point; rank verified at the reference". Correct the 1D growth law.
* T7: add τ ≥ |x|/|b|. Make all event-operator mesh laws one-sided lower bounds under H4. Fix the α = d
  log. Make the AlN sentence conditional.
* Numbers: replace "0.001–65 % / ≤ 1.65" by "≤ 1.69 (max 1.684, 2D N = 5 full symmetry)". Correct the 1e8
  statement. Round lower bounds down. Store the pivot bases. Add the implicit-function transfer for the
  headline witnesses.

## 14. Candidates for machine checking

* **L3** (finite sums, Cauchy–Schwarz plus double counting; the statement in §5, step L.5). It is short,
  central and used as a feasibility test (P5).
* **P4a** and the equality case of L3 (§5.3) can go in the same file.
* T1 needs the Cauchy–Binet adjugate identity on a hyperplane. It is classical and verified here in exact
  arithmetic, so it is a lower priority for formalization.

## 15. Check files (all in this directory; run with the research Python, `-B`)

| file | content |
|---|---|
| `check_cone_bounds.py/.json/.log` | T1 identity, P3, T2 family, P4, P6 (exact); counterexample to "sharp only without lifetimes"; generic lifetime slice |
| `certify_M_1d9.py/.json/.log` | rigorous certification of M for 1D N = 9 (integer-rank bound, time-reversal certificate, exact top values) |
| `check_rta_bound.py/.json/.log` | L3 random exact tests; exact equality 1/3 and 1/2; repeated-daughter matrix inequality |
| `check_lmi.py/.json` | LMI set ≠ admissible sublevel set (exact) |
| `check_dimension_count.py/.json` | ranks of W and [W; ∇K] at the AlN reference |
| `check_ir.py/.json/.log` | RTA mesh laws (exact lattice sums), H4 diagnostic on AlN, τ ≥ \|x\|/\|b\| |
| `check_bt_numbers.py/.json` | vacuity numbers from the stored lower bounds; P6 value |
| `check_pivot_bases.py/.json/.log` | re-run and 50-digit verification of pivoting lower bounds |
| `check_witness_ift.py/.json/.log` | implicit-function transfer of AlN witnesses to the exact data |

## 16. Overall verdict

**The mathematical core survives.**

* T1, T2, P3, L3, P4 and P6 are valid proofs under the stated hypotheses. Each was re-derived independently
  and checked in exact arithmetic on test geometries.
* L3's constant 3 is correct, applies unchanged to repeated-daughter events, and is sharp (attained exactly).
* The dimension counts are confirmed at the reference.
* The 1D N = 9 constant M = 16284.006723562 is now certified (§2.9). The campaign's own singularity test was
  heuristic but happened to classify every subset correctly (3,566,333 bases).

**What does not survive as written:**

* The claim that the cone bound is sharp *only* without lifetimes (false in general; open for physical currents).
* The "≥ 1e8 in every geometry" vacuity number (false for two geometries).
* The "τ_min within 1.65 τ_CS everywhere" summary (1.684 on the stored data).
* The literal LMI set identity (definitional).
* The T7(ii) step to τ = ∞ (incomplete, repaired here).
* The event-operator mesh exponents and the unconditional AlN divergence sentence (unjustified).

None of these affects the campaign's headline conclusions, which rest on proved lower bounds (τ_CS) and on
explicit feasible witnesses, not on the sharpness of the cone bound. With the repairs of §13, the statement list
can be labelled: T1, T2, P3, L3, P4, P6 PROVED; P5 PROVED (extended K); dimension count VERIFIED AT THE
REFERENCE; T7(i) PROVED; T7(ii)–(iii) DERIVED UNDER H4.
