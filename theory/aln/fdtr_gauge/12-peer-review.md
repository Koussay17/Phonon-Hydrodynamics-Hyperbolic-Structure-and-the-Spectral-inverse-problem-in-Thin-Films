# 12 — Peer review: FDTR detectability of the film memory-time range (A) and rank and gauge of the 1D non-Fourier front-face response (B)

Hostile review at Physical Review B / Journal of Applied Physics standard, 12 November 2026. Object: the two results
intended for paper 1 (theory/methods), as recorded in this campaign and in the draft manuscript
(`C:\Users\Koussay\these\paper\manuscript\sections\frontface.tex`, `fdtr.tex`). The review assumes nothing that was not
re-derived, recomputed or traced to a source.

---

## 0. Material reviewed and independent checks

**Read.** `00-question.md`, `02-assumptions.md`, `04-branches.md`, `05-evidence-ledger.md`, `06-disagreements.md`,
`07-experiments.md` (A1–A7), `10-proof-status.md`, `11-open-questions.md`, `01-known-results.md`, `08-literature.md`,
`proofs/B-analytic.md` with its addendum, `branches/B-proof-audit.md`, `branches/A-numerical-audit.md`,
`branches/B-counterexamples.md` (§0–2), `derivations/AB-asymptotic.md`, `derivations/B-symbolic.md`; upstream
`theory/aln/fc3_memory/07-experiments.md` (E26, E30–E37), `notes/07_Note_partie_IV.tex`, `src/twin_fdtr.py` (docstring,
`GKBoundary`, `kernel`), `paper/PRIOR_ART.md`, `paper/CLAIM_MAP.md`, manuscript `frontface.tex` and `refs.bib`.
The D: copies of 05, 07 and the proof are identical to the repository copies.

**Independent checks** (own code, no import of campaign scripts except where a cross-check against the twin is the
point; scripts and raw output in `review/`; research Python, two threads).

| # | Script | What it tests | Result |
|---|---|---|---|
| R1 | `r1_b_checks.py` | Own 60-digit implementation of Z; ranks on the strata; certified twin; thin-slab ratio; five-parameter ranks; Laurent coefficients c₀–c₃ by a Cauchy integral; the same under an alternative interface-temperature convention | All B claims reproduced under the local convention (details §2). Under the entropic convention the resonance degeneracy disappears (§5, B-2) |
| R2 | `r2_twin_entropic.py` | Closed form with load (R + 1/S)/(1 + τ_ℓp) against `twin_fdtr.kernel(k = 0)` with `GKBoundary(temperature="entropic")` | Agreement ≤ 6.8e-16 for both conventions; resonant GK film (τ_R = τ_ℓ = 25 ps) vs Fourier film, 500 nm AlN on sapphire, 10 kHz–200 MHz: max \|Δφ\| = 0 (local), 0.089° (entropic) |
| R3 | `r3_six_param.py` | Known Fourier transducer, unknown top and bottom resistances R1, R2: rank of (ξ₁, b, τ_R, τ_ℓ, R1, R2) | Rank 6 at generic, R2 = 0, b = b_s and Cattaneo points (normalised s₆ = 1.8e-2–4.8e-2); 5 on Δ; 4 on Δ ∩ {R2 = 0} and on D₂ ∩ {R2 = 0} |
| R4 | `r4_film_kernel.py` | Is the film-regularised kernel a single pole in the band? RTA 64×64×35 mode data (fc3-memory E27), same regulator | ⟨τ³⟩⟨τ⟩/⟨τ²⟩² = 1.30 / 1.41 / 1.69 / 1.88 (30/100/300/500 nm). At 200 MHz Im κ(ω)/κ(0) of 1/(1 + iωτ_mem) matches the multi-pole kernel to 0.007 % / 0.04 % / 0.2 % / 0.4 %; the real-part correction is 1.3–1.9× the single-pole one but ≤ 6 % of \|Im κ/κ₀\| |
| R5 | `r5_info_gain.py` | Threshold-free comparison: uniform prior of half-width w, Gaussian likelihood σ | §5 A-1 and §6 A-m1 |
| R6 | `r6_regulator_norm.py` | Normalisation c of the heuristic rate c\|v_z\|/L | τ_mem,z(500 nm) = 36.3 / 29.7 / 25.6 / 22.0 / 17.7 ps for c = 1/2, 3/4, 1, 4/3, 2 (RTA 64×64×35) |
| R7 | `r7_entropic_scan.py` | 1200 random points, both conventions | Local: every Δ point rank-deficient (s₄ at the finite-difference floor), D₂ exact; entropic: min s₄ = 1.3e-3 on Δ, 1.9e-3 on D₂, 9.3e-3 interior — no rank drop found |
| R8 | `r8_a_1d_fisher.py` | Own 1D Au/R1/AlN/R2/sapphire Fisher analysis; (R1, R2) map at 500 nm | 500 nm: peak \|∂φ/∂ln τ_R\| = 0.0208° at 65 MHz; σ(ln τ_R) = 0.145 / 0.189 / 0.305 / 1.96 (N0–N3, 0.01°), consistent with the 1D-limit offsets of A4. R1 dominates (§5, A-1) |

---

## 1. Summary of the claims

**Result A.** In the illustrative Au(80 nm)/AlN/sapphire FDTR twin (note-16 values, R1 = 1e-8, R2 = 2e-8 m²K/W,
beams 10/8 µm, 80 points 10 kHz–200 MHz, 1 % amplitude noise), the range of the film-regularised Cattaneo memory time
left by the static data (fc3 per-orbit class, three meshes, ≈ ±2.8–5.2 % at 30–300 nm, inner bounds) is not
resolvable when instrument and interface nuisances are free (N ⊇ N2): σ(ln τ_R)/w ≥ 4.2 at 0.01° and ≥ 34 at 0.1°;
marginal with instrument-only nuisances at 500 nm and 0.01°.

**Result B.** For a 1D Cattaneo/GK film on a Fourier half-space (front-face flux and temperature, local interface
temperature, exact data): injectivity of θ = (ξ₁, b, τ_R, τ_ℓ) ↦ Z off Δ = {τ_R = τ_ℓ}; Jacobian rank 4 off Δ, 3 on Δ,
2 on Δ ∩ {b = b_s}; rank 5 with an interface resistance (4 and 3 on Δ for R > 0 and R = 0); two generic frequencies give
local rank 4 but not global uniqueness (interval-certified twin); the depth-scaling gauge is the only non-identifiable
direction off Δ, broken by a known length-carrying quantity.

---

## 2. Classification of statements

| Statement | Status after this review | Basis |
|---|---|---|
| B: impedance, single-valued quadrupole (Lemmas 0–3, with D1) | PROVED | proof + audit; R1 reproduces Z and c₀–c₃ to 1e-60 |
| B: Theorem 1 (injectivity, fibres), local convention | PROVED (under Assumption M) | two proofs + audit; the explicit low-frequency inversion (AB-A B5, B-S §8) is a third, elementary route |
| B: Theorem 2 (ranks 4/3/2, Σ = Δ) | PROVED under Assumption M; **CONVENTION-DEPENDENT** | R1, R7: reproduced for the local convention; no rank drop on Δ or D₂ under the entropic convention of the twin |
| B: Theorems 1′, 2′ (resistance) | PROVED under M (Proof B rests on Lemma T for R > 0) | audit; R1 consistent (rank 5 / 4 / 3) |
| B: Theorem N2, Proposition TW | PROVED (N2 with repairs D1–D2); TW certified | R1: \|Z(θ*) − Z(θ′)\| = 2.9e-21 and 3.5e-21 at ω = 1/2 and 4 (limited by the 20 printed digits); 2.1 % apart at ω = 1.7 |
| B: Conjecture N3 | CONJECTURED | multistart only |
| B: Theorem G (gauge off Δ; breaking by known e, λ, ρc, a, ℓ) | PROVED (elementary) | dimensional analysis |
| B: gauge breaking by lateral diffusion (isotropic) | CONJECTURED (D7 open); anisotropic counterexample PROVED | — |
| B: interior source breaks the gauge | PLAUSIBLE (structural remark only) | — |
| B: extension to unknown top resistance with known transducer | NUMERICALLY SUPPORTED (this review, R3), not proved | — |
| A: σ(ln τ_R) tables, linearity in τ_R | NUMERICALLY DEMONSTRATED (model-conditional) | two independent implementations to ≤ 7.2e-6; R8 reproduces the 1D scale |
| A: not resolvable, N ⊇ N2, 0.1° | NUMERICALLY DEMONSTRATED within the declared model; robust to every variation tried (≥ 10× over R1, R2 ∈ [2e-9, 2e-8] m²K/W) | A4, audit §7, R8 |
| A: not resolvable, N ⊇ N2, 0.01°, "≥ 4.2×" | DERIVED UNDER ASSUMPTIONS (note-16 resistances, extrapolated 500 nm width); **not robust** as a general statement (margin 1.1–1.8 for R1 ≤ 5e-9) | R8, audit F8 |
| A: marginal for N0/N1 at 500 nm, 0.01° | NUMERICALLY SUPPORTED | A4, audit F1 |
| A: τ_mem ↔ single-pole Cattaneo τ_R adequate in the band | NUMERICALLY SUPPORTED (RTA kernel, this review, R4) | R4 |
| A: Cattaneo vs Fourier film not resolved with realistic nuisances and 0.1° | NUMERICALLY DEMONSTRATED (model-conditional) | A7, audit §4 |
| A: any statement about a real sample | UNSUPPORTED (by design, A7 of 02-assumptions) | — |

---

## 3. Strongest parts

1. **Numerical reliability of A.** Two implementations that share only the physics (twin; independent Julia transfer
   matrices with adaptive quadrature, dual numbers and QR/BigFloat Fisher algebra) agree on all 270 finite-beam entries to
   ≤ 7.2e-6; the 1D limit, quadrature, finite-difference step and Fisher algebra are converged; noiseless profiles and
   200-seed maximum-likelihood fits confirm the local Fisher numbers; the response is shown to be linear in τ_R (not in
   ln τ_R), which justifies reading σ(ln τ_R) as a relative change. My own minimal 1D code (R8) reproduces the scale.
2. **Pre-declared nuisance sets** (N0–N5) and an explicit falsification criterion in `00-question.md`; the audit's
   narrowing of the claim (F1, F8) has been adopted in the ledger. This is good practice.
3. **Mechanism.** The exact identity ∂ ln Z/∂ ln τ_R = −[τ_Rp/(1 + τ_Rp)] s_λ and the leading-order law
   δφ ≈ w ωτ_R W(ω), W = −∂ ln|Z|/∂ ln λ_z, explain the result in one line and locate the information (65–200 MHz).
4. **B is correct as stated** within Assumption M. Every load-bearing identity was rederived by a second route
   (Riccati recursion, two CAS, ball arithmetic); I reproduced the rank table, the thin-slab ratio −10/21, the twin and
   the Laurent coefficients with an implementation of my own. The audit (D1–D9) is careful and the repairs are minor.
5. **Honesty of scope.** The proof states that it is a uniqueness argument, not an inversion method; the practical
   conditioning (s₄ ≈ 1e-6–1e-4, E-B8) and the finite-band near-equivalence of GK and Cattaneo (E-B5) are recorded.

---

## 4. Fatal issues

None. Nothing found invalidates either result within its stated assumptions. The major issues concern the
generality, the physical reading and the novelty of what is stated.

---

## 5. Major issues

### A-1 (MAJOR). The headline margins are set by an extrapolated width and by the illustrative interface resistances; the 0.01° statement is not robust

*What supports the claim.* "≥ 4.2× at 0.01°, ≥ 34× at 0.1°" are both attained at **500 nm**, against the upper end
(8 %) of the trend bracket 5.5–8 %: N2 gives 0.340/0.08 = 4.25 and 2.73/0.08 = 34. No admissible range was computed at
500 nm (A4 of 02-assumptions). Where the range is computed (30–300 nm), the minima for N ⊇ N2 are 0.662/0.052 ≈ 12.7
(0.01°, 300 nm) and 5.36/0.052 ≈ 103 (0.1°, 300 nm). The quoted sentence combines the 30–300 nm widths with the 500 nm
margin.

*Robustness.* Within the declared physics the margin depends on inputs that are uncertain by O(1):
- **Top resistance R1** is the dominant lever (R8, 1D, 500 nm, N2): σ(ln τ_R) at 0.01° is 0.09–0.11 for R1 = 2e-9,
  0.10–0.14 for R1 = 5e-9, 0.30–0.44 for R1 = 1e-8, 1.7–2.0 for R1 = 2e-8 (R2 between 2e-9 and 2e-8). With
  R1 ≤ 5e-9 m²K/W (G₁ ≥ 200 MW m⁻²K⁻¹) the 1σ ratio σ/w at 500 nm falls to 1.1–1.8 against 8 %. The audit's finite-beam
  value (R1/R2 = 3e-9/5e-9, σ × 1/3.5) agrees. The note-16 value R2 = 2e-8 also appears pessimistic: TDTR conductances of
  epitaxial nitride heterointerfaces and of epitaxial metal/sapphire interfaces are reported at several hundred
  MW m⁻²K⁻¹ (seen at abstract level in this review; the authors should cite primary AlN/sapphire and Au(Ti)/AlN values).
- **Regulator normalisation** (R6): τ_mem,z(500 nm) = 17.7–36.3 ps for c|v_z|/L with c = 2 … 1/2; σ ∝ 1/τ_R (A6), i.e.
  ×0.70–1.45, and w(L) moves too (not computed).
- **Inner bounds**: E39 is still running; the single-mesh 5×5×3 film range was [0.909, 1.098] at 100 nm, the calibrated
  classes up to ±8 %, so a factor ~2 widening is plausible.
- **Design**: c-optimal allocation of the same 80 points gains 1.3–2.9 (audit §7).

Compounding these, the 0.01° verdict at 500 nm ranges from "resolvable at 1σ" to "4× below". The 0.1° verdict survives
every variation computed here (σ/w ≥ 10 over the R map; ≥ 88 with λ_z, C free, N3).

*Information measure (R5).* With a uniform prior on the admissible range, FDTR removes 1.8 % of the prior variance
(0.013 bit) for N2, 0.01°, 500 nm on the declared stack, 0.03 % at 0.1°, but 20 % (0.16 bit) with all the favourable
changes of the audit at 0.01°.

*Required.* (i) Compute w(500 nm) with the same optimiser, or restrict the headline to 30–300 nm (≥ 12× at 0.01°,
≥ 100× at 0.1°, ≥ 6.7× with all favourable audit changes at 300 nm). (ii) State the 0.1° result as the claim and the
0.01° result as stack-dependent and marginal. (iii) Give σ/w (or the information gain) as a map over a declared
plausible box of (R1, R2, Au thickness, σ_φ), and quote its worst case. (iv) Correct `07-experiments.md` synthesis
("every nuisance set N0–N5", "margin of more than an order of magnitude at every thickness"), which the audit (F1, F8)
and R8 contradict.

### A-2 (MAJOR). The compared quantity is not the FDTR-visible one in the paper's own closure

*What supports the claim.* The twin is Cattaneo (τ_ℓ = 0) or GK with τ_ℓ held fixed; the width compared is that of
τ_mem = τ_R alone. But Result B and AB-A (Proposition B3′) show that at first order in ωτ the response depends on the
times only through **δ = τ_R − τ_ℓ**. In the conserving closure the paper advocates, τ_ℓ = 4τ_c/5 with
τ_c = (τ_N⁻¹ + τ_R⁻¹)⁻¹ (manuscript Eq. `noresonance`), so τ_ℓ co-varies with τ_R and τ_N, and
dδ/dτ_R = 1 − (4/5)τ_N²/(τ_N + τ_R)² ∈ (1/5, 1). Since τ_ℓ/τ_R = 4τ_N/[5(τ_N + τ_R)], τ_ℓ is comparable with τ_R
(0.4 τ_R for τ_N = τ_R) unless τ_N ≪ τ_R; it is not negligible a priori.

Two consequences. (a) The per-unit sensitivity used in A is the largest possible (dδ/dτ_R = 1), so the Cattaneo variant
is conservative for the detectability side — this should be stated, since it strengthens A. (b) The width side is not
established: static data fix total lifetimes but not the normal/resistive split, so the class range of δ(L) can differ
from that of τ_R(L), possibly by more than the ±5 %. The conclusion would change only if the range of δ exceeded that of
τ_R by about the margins of A-1 (≥ 4× at 0.01°, ≥ 34× at 0.1°), which is not excluded by anything in the record.

*Required.* Compute the class range of τ_ℓ(L) (or δ(L)) with the same optimiser, or restrict the claim explicitly to the
Cattaneo projection and justify τ_ℓ = 0 for these films. Replace the GK variants with τ_ℓ = 2τ_R, which violate the
closure bound τ_ℓ < 4τ_R/5, or label them as phenomenological.

### A-3 (MAJOR). Conditionality and significance: the informative statement is missing and the novel content is a scaling bound

*Conditionality.* "The range left undetermined by static transport data" is the range in the fc3 per-orbit ±10 %
class, as an inner bound, with the heuristic |v_z|/L regulator. The paper's own central contribution 1 states that
τ_mem is unbounded in the positive-semidefinite class. The claim must name the class; otherwise "FDTR cannot resolve
the static-data gap" reads as a statement about all closures, which it is not.

*Meaning given inner bounds and the heuristic regulator.* Both work against the claim only through the width and τ_R
scales, and the record quantifies how much they would have to change: at N2, 0.1°, the range would have to reach
τ_R ≥ 3.7 τ_R,0 (500 nm) or the base τ_R would have to be ≈ 30× larger. Within the class, the conclusion at 0.1° is
therefore meaningful; at 0.01° it is not robust (A-1). The regulator also replaces boundary scattering by a uniform
bulk rate and omits Knudsen-layer temperature slip at the interfaces, i.e. the same interface physics that controls B
(B-2). The single-pole reduction itself is adequate (R4): the reactive part of the multi-pole film kernel matches
1/(1 + iωτ_mem) to 0.4 % at 200 MHz for 500 nm; this should be reported, because it removes an obvious objection.

*Significance.* The verdict follows from δφ ≈ w ωτ_R W(ω): the film weight W is ≤ 0.05 and falls to zero near the band
top (sign change at ≈ 180 MHz, AB-A A3) while ωτ_R ≤ 0.035, so that max |∂φ/∂ ln τ_R| ≈ 0.02° (R8: 0.0208° at
65 MHz); a few-percent change is then ≲ 1.6e-3°, below any per-point noise considered.
The extensive Fisher tables confirm a back-of-envelope estimate. The physically more significant result is A7: with
λ_z and C free (N3, the realistic minimum, since published AlN film conductivities carry ≥ 20 % uncertainty) or at
0.1°, even the existence of a 25 ps Cattaneo time against a Fourier film is not resolved (Δχ² ≤ 0.46). This bears
directly on claims of FDTR-extracted relaxation times and should be positioned against Beardo et al. 2020 (Si: the
FDTR-visible deviation is nonlocal, τ = 42 ps negligible), the PRISM preprint (τ held fixed as unidentifiable) and
Camacho de la Rosa et al. 2025 (sensitivity without noise or nuisances). The complementary, informative direction is
what FDTR excludes, which is the statement relevant to the unbounded class: since the response is linear in τ_R (to
8 % up to 10 τ_R,0, A7), at 500 nm, N2, 0.1° the data exclude τ_R ≥ 1 + 2σ_F ≈ 6.5 τ_R,0 at 2σ (σ_F = 2.73); with
λ_z and C free (N3) only τ_R ≳ 35 τ_R,0 (Fisher value; the profile is steeper there); at 30–100 nm nothing useful.

*Required.* Name the class in the claim; lead with the scaling law and with A7; add the exclusion statement; state the
novelty as "first Cramér–Rao estimate with instrument, interface and transducer nuisances, within the literature
searched" (A-K8), not as a new physical effect.

### B-1 (MAJOR). Prior art and novelty: the result is partly anticipated and its content is routine

1. **Equivalence to DPL.** For this boundary-value problem the 1D GK film is a dual-phase-lag film with τ_q = τ_R,
   τ_T = τ_ℓ (B-K6). The reduction to Fourier at τ_q = τ_T, i.e. the degenerate set Δ, is classical in the DPL
   literature (Tzou, J. Heat Transfer 117, 8 (1995), and monograph, secondary attribution; the "diffusion mode" of
   Tang & Araki, Int. J. Heat Mass Transf. 42, 855 (1999)). The manuscript attributes the resonance to Both et al. 2016
   and Kovács 2018 only; `refs.bib` contains neither Tzou nor Tang–Araki.
2. **Ordóñez-Miranda & Alvarado-Gil, J. Heat Transfer 132, 061302 (2010), doi:10.1115/1.4000748.** The Crossref record
   confirms the configuration (finite DPL layer on a semi-infinite layer, modulated surface flux) and the content:
   in a frequency range with strong oscillations at large effusivity contrast, explicit formulas determine both time
   parameters and further thermal properties of the finite layer. That is a constructive identifiability result on part
   of the parameter space, i.e. partial prior art for Theorem 1. Its reference list points to further measurement
   methods for lag times (Ordóñez-Miranda & Alvarado-Gil, Int. J. Therm. Sci. 2009; Basirat Tabrizi et al., Int. Commun.
   Heat Mass Transf. 2009), not read either. Only the abstract has been read; this is not acceptable for a novelty
   statement.
3. **Method.** Identifiability of constant coefficients of 1D parabolic systems from boundary data by Laplace-domain
   analyticity is standard (Kitamura & Nakagiri, SIAM J. Control Optim. 15, 785 (1977), with the 1983 comment;
   Pierce 1979). The low-frequency route (AB-A B5, B-S §8, Proof C) inverts five explicit coefficients and gives
   injectivity and the rank off Δ in a few lines; on Δ the rank is a two-line computation. The analytic-continuation
   ("second sheet") and Lindemann–Weierstrass proofs are elegant but not needed, except for R > 0.
4. **Thin-slab confound.** The rank-3 stratum on Δ ∩ {R = 0} is the classical indistinguishability of film resistance
   and interface resistance at first order in thickness (Lee & Cahill 1997, thickness series), made exact for a Fourier
   half-space by p Z_s² = 1/b_s².
5. **Gauge.** Dimensional analysis; the manuscript says so. The same structure appears for telegraph lines (Jaulent
   1982) and finite-length diffusion impedances.

Within the literature searched, the complete classification (fibres, rank table, matched effusivity non-degenerate off
Δ, resistance strata, local-but-not-global identifiability with two frequencies) does not appear to have been
published, but it is incremental. It is a correct and useful section of a methods paper, not a headline result.

*Required.* Read O-M & A-G 2010 and companions in full and state precisely what is new; cite Tzou and Tang–Araki for
Δ, Kitamura–Nakagiri/Pierce for the method, Lee–Cahill for the thin-slab confound; present B as a proposition with the
elementary low-frequency proof in the main text and the other proofs in an appendix; use "structural identifiability"
throughout.

### B-2 (MAJOR). The exceptional set is a property of the interface-temperature convention, not of the GK law

*The assumption.* (F4) imposes continuity of the energy-density temperature at z = e and is justified in the draft by
"the one-dimensional problem needs none, because after elimination the flux in each layer satisfies a second-order
equation". That argument concerns the number of conditions (well-posedness), not their physics. The
kinetic temperature-jump condition of Beardo et al. 2020, cited in the same paragraph for 3D, has a nonzero 1D limit:
T_above − T_GK = Rq − γ⁻¹(β − χ_nn) ∂_z q, and ∂_z q = −ρc p θ. The project's own twin implements the entropic variant,
θ = T − (1 + α)(ℓ²/λ) ∂_z q, i.e. θ = (1 + τ_ℓ p) T at k = 0, documented there as the natural condition of the
dissipation identity, whereas the local convention is documented as having an interface term without sign. Since the
correction is O(ℓ²), the same order as the GK term retained in the bulk, keeping one and dropping the other is not a
consistent truncation.

*Evidence (R1, R2, R7).* Under the entropic convention the load seen by the film is (R + 1/S)/(1 + τ_ℓ p), which
depends on t along Δ. The twin's k = 0 kernel agrees with this closed form to 6.8e-16. A GK AlN film at
τ_R = τ_ℓ = 25 ps differs from the Fourier film by 0.089° in phase (500 nm on sapphire, no transducer), against 0
under the local convention. The four-parameter rank is 4 on Δ (s₄ = 0.030), on D₂ (0.032) and at the Fourier point
(5.5e-4), and no rank drop appears at 1200 sampled points (minimum s₄ = 1.3e-3 on Δ). The R = 0 thin-slab ratio becomes
frequency-dependent and the five-parameter rank is 5 on Δ ∩ {R = 0}. Hence "Σ = Δ exactly", the ranks 3 and 2, and the
rank-3 stratum with resistance hold only for c_n = 0.

*Required.* State Theorems 1, 2, 1′, 2′ as conditional on the local convention (c_n = 0) in the theorem, not only in
§11. Add the classification for c_n = κℓ²/λ, κ ≠ 0 (entropic κ = 1 + α; kinetic coefficients of Beardo et al.), and for
c_n free (one additional parameter). Discuss which convention follows from the conserving closure. Remove "Fourier
resonance" as a physical degeneracy from any statement about measurements.

### B-3 (MAJOR). The structural statements have little bearing on the closure and the experiment the paper is about

1. **The closure lies off Δ.** The draft itself (Eq. `noresonance`) shows τ_ℓ/τ_R < 4/5 for the conserving closure at
   all positive times, so the rank drop never occurs for it; Δ matters only for the historical closure, at
   τ_R/τ_N = 4/5 or 9/5, which the draft calls a formal marker.
2. **Rank 4 is practically rank 3 in FDTR.** τ_R and τ_ℓ separate only at relative order (ωτ)², with ωτ ≤ 0.035 in the
   band: s₄ ≈ 1e-6–1e-4 for 500 nm AlN (E-B8), whitened cosine −0.9994 between the two time columns (AB-A), GK and
   Cattaneo with the same δ indistinguishable to 1.2e-8 on a thin band (E-B5). Matched effusivity is non-degenerate off
   Δ, but its leading minor is ∝ (τ_R − τ_ℓ)².
3. **Two frequencies.** Theorem N2 and Proposition TW concern exact data. The certified twin lies on the stratum b = b_s;
   one further, uncertified twin was found off it (θ* = (0.8, 0.5, 0.2, 0.9), τ_ℓ > τ_R), and the multistart found none
   for the interior and AlN-like targets. Its practical bearing is unclear.

*Required.* Put the practical statement first (only ξ₁, b and δ are accessible at ωτ ≪ 1, with singular values), next to
the theorem; state that the paper's closure is always off Δ; give a twin in the physically relevant region or say that
none was found.

### B-4 (MAJOR). Scope gap between B and the configuration of A

Theorems 1–2′ cover a bare film, and Corollary 7.1 a fully known transducer and top contact. Every FDTR measurement
has an unknown transducer/film resistance R1, which A frees in N2 and which dominates the sensitivity (A-1). The
"theory" and the "numerics" are therefore about different models. R3 shows that, under the local convention, a known
Fourier transducer with unknown (R1, R2) gives rank 6 generically, 5 on Δ and 4 on Δ ∩ {R2 = 0}, so the structure of
Theorem 2′ appears to persist.

*Required.* Prove or certify this extension, or state its absence prominently; connect the practical conditioning of
B to the N2 numbers of A.

---

## 6. Minor issues

**A-m1 (comparison metric).** σ against the half-width is fair and errs towards declaring resolvability, so "not
resolvable" is conservative under it; but its reading is ambiguous (1σ vs 2σ, half-width vs full width, audit F2).
σ = w corresponds to a 25 % reduction of the prior variance (0.21 bit, R5). Report the information gain, or the
posterior-to-prior variance ratio, instead of a threshold.

**A-m2 (single-pole mapping).** Report R4 (RTA kernel): the mapping τ_mem → τ_R is adequate in the band; the second
moment exceeds the single-pole value by 1.3–1.9×, a correction ≤ 6 % of the first-order signal, largely absorbed by
λ_z. Repeat with the event operator.

**A-m3 (regulator).** Report τ_mem(L) for c ∈ [1/2, 2] (R6) and recompute at least one width w(L) with c ≠ 1.

**A-m4 (record consistency).** `07-experiments.md` A1 still misquotes the 300 nm event/RTA ratio (0.973 instead of
1.0097) and presents the RTA mesh step +0.22 % as a convergence figure (audit F4: 1–2 % is the honest uncertainty);
A2 presents the 2.2e-14 normalisation offset as a quadrature error (F3); `paper/PRIOR_ART.md` item 3 still reads
"numerical audit pending" and "margins 34×–4.5·10⁴×".

**A-m5 (noise model).** A constant per-point phase noise is optimistic at 150–200 MHz, where the amplitude is smallest
and the τ_R information sits; systematic phase errors (≥ 0.1°, A-K2) are not modelled. Both strengthen the negative
claim; say so, and treat 0.1° as the realistic case.

**A-m6.** The A6 GHz design sweep uses Fourier Au without electron–phonon coupling or optical penetration; remove it or
flag it as model-internal. Keep the finite-beam GK rows (slip walls, no independent reference, audit F7) out of the main
text.

**A-m7 (base values).** The 12×12×8 event values are probably 2–4 % low (audit §6), the 500 nm value is interpolated,
and the in-plane λ_r on uniform meshes carries the q_z = 0 artefact; state these with E-A4.

**A-m8 (reproducibility).** Inputs live on an external drive. Deposit the scripts, raw JSON, the RTA chunk files
(7.5 MB), the REF blocks and the twin commit (e4287dd) with checksums.

**B-m1.** Incorporate the audit repairs D1–D9 in the main text (Lemma 0′ before Lemma 1; absolute convergence before
Cauchy–Binet; uniform disc for termwise differentiation; corrected status wording D4; TW digits and regularity D5; D6
numbers; D7 status; D9 wording).

**B-m2.** b_s need not be known at function level (z₋₁ = 1/b_s); say so, or explain that it is kept for finite bands.

**B-m3.** "Its breaking mechanisms" holds only for known e, λ, ρc, a or ℓ; the relation a = v²τ_R/3 is a grey-model
identity with no single measurable v for AlN; lateral diffusion (isotropic) is a conjecture (D7) and interior sources a
structural remark. The draft already defers lateral diffusion; keep it that way.

**B-m4.** For R > 0 the proof rests on Lemma T (compact Riemann surfaces). An elementary route (real-axis exponential
versus algebraic growth, or completing the moment minors on Q = 0, where the record stops at Q = P₁ = P₂ = 0) would
remove the only non-elementary step.

**B-m5.** Homogeneous film only; under Fourier a graded film has a functional gauge (Krapez 2023). State the
limitation next to Theorem G.

---

## 7. Decisive tests

1. Read Ordóñez-Miranda & Alvarado-Gil 2010 and companions in full; list which parts of Theorem 1 and Proposition B5
   they contain.
2. Rank and fibres of the 1D map for θ = T − c_n ∂_z q at the film/substrate interface, with c_n = κℓ²/λ (κ = 0,
   1 + α, kinetic coefficients) and with c_n free. Prediction from R1/R7: Σ = ∅ at sampled points for κ ≠ 0.
3. Known transducer, unknown R1 and R2: certified rank 6 at generic points and the exceptional strata (R3 gives the
   prediction), with the practical singular values in the FDTR band.
4. Admissible range of τ_mem at 500 nm (or removal of 500 nm from the headline); E39 multistart ranges propagated to A.
5. Class range of τ_ℓ(L) or δ(L) with the same constraints; comparison with σ(δ).
6. σ/w and information gain over a declared plausible box (R1, R2 ∈ [2e-9, 3e-8] m²K/W, Au 20–100 nm, σ_φ 0.01–0.1°,
   optimal allocation of points); quote the worst case.
7. Twin evaluated with the full multi-pole film kernel of the extreme witnesses of E33/E37, comparing phases directly.
   This tests the observable, bypasses the τ_mem ↔ τ_R mapping and answers the question as posed.

---

## 8. Unresolved questions

- Which interface temperature does the conserving closure imply at a film/substrate or transducer/film interface, and
  is the local convention consistent at O(ℓ²)?
- Outer bounds of the fc3-class range (E39), and the range of the normal/resistive split under the static constraints.
- Identifiability of (ξ₁, b, τ_R, τ_ℓ, a) from the beam-averaged response of an isotropic GK film (Conjecture L, D7).
- Global identifiability with three frequencies (Conjecture N3).
- Whether the in-plane memory time (q_z = 0 artefact, τ_x = 111–157 ps) matters with beams smaller than 10 µm.

---

## 9. Recommendation

**Major revision.** Neither result is wrong. A is numerically sound and, at 0.1° phase noise and wherever the
admissible range is actually computed (30–300 nm), robust by one to two orders of magnitude. As stated, however, its
headline numbers rest on an extrapolated width and on illustrative resistances that the record itself identifies as the
largest lever; the 0.01° statement is not robust; the comparison ignores the τ_ℓ co-variation that the paper's own
closure implies; and the scientifically interesting content (a scaling bound, and the invisibility of the memory time
itself under realistic nuisances) is under-sold relative to the narrower negative statement. B is proved and audited
within Assumption M, but it is partly anticipated (DPL resonance, Ordóñez-Miranda & Alvarado-Gil 2010, standard
Laplace-domain identifiability), mathematically routine, of limited practical bearing in the FDTR band, and its
exceptional set depends on an interface-temperature convention that is not the one with a sign-definite dissipation
identity in the project's own twin. Each issue is fixable with computations that the existing code supports; until the
prior-art reading (B-1), the convention analysis (B-2) and the 500 nm width or headline restriction (A-1) are done, the
results are not ready for submission.

---

## 10. Objection list

| ID | Severity | Objection | Required action |
|---|---|---|---|
| A-1 | MAJOR | Headline margins (≥ 4.2×, ≥ 34×) set at 500 nm against an uncomputed width and the note-16 resistances; 0.01° margin falls to 1.1–1.8 for R1 ≤ 5e-9 | Compute w(500 nm) or restrict to 30–300 nm (≥ 12×, ≥ 100×); claim at 0.1°; sensitivity box; fix `07` synthesis |
| A-2 | MAJOR | FDTR sees δ = τ_R − τ_ℓ; closure ties τ_ℓ to τ_R, τ_N; width of δ not established | Class range of τ_ℓ or δ, or restrict to the Cattaneo projection; state that Cattaneo is the most sensitive case |
| A-3 | MAJOR | Class-conditional claim stated generally; significance lies in the scaling bound and in A7, not in the narrow negative | Name the class; lead with the scaling law and A7; add the exclusion statement; position against Beardo 2020, PRISM, Camacho 2025 |
| B-1 | MAJOR | Partly anticipated (Tzou, Tang–Araki, O-M & A-G 2010) and routine (Kitamura–Nakagiri, Pierce, Lee–Cahill) | Read O-M & A-G in full; cite; elementary proof in the main text; "structural" wording |
| B-2 | MAJOR | Σ = Δ and the R = 0 stratum depend on c_n = 0; the entropic and kinetic conventions remove them | Theorem conditional on c_n = 0; classification for c_n ≠ 0; discuss the closure's interface condition |
| B-3 | MAJOR | Paper's closure never reaches Δ; rank 4 practically rank 3 in the band; certified twin on a special stratum | Practical statement first; closure off Δ; physical twin or none |
| B-4 | MAJOR | Theorems exclude the unknown top resistance present in every FDTR measurement | Prove or certify rank 6 with known transducer and unknown R1, R2 (R3 supports it) |
| A-m1 | MINOR | Half-width threshold ambiguous | Information gain or variance ratio |
| A-m2 | MINOR | Single-pole adequacy not reported | Report R4; repeat with the event operator |
| A-m3 | MINOR | Regulator normalisation untested | τ_mem(L) for c ∈ [1/2, 2]; one w(L) |
| A-m4 | MINOR | Record inconsistencies (F3, F4, PRIOR_ART item 3) | Correct |
| A-m5 | MINOR | Noise model optimistic at the band top; systematics absent | State; 0.1° as realistic case |
| A-m6 | MINOR | GHz sweep outside model validity; finite-beam GK rows unverified | Remove or flag |
| A-m7 | MINOR | Base-value caveats (mesh bias, interpolation, q_z = 0) | State with E-A4 |
| A-m8 | MINOR | Data on an external drive | Deposit with checksums |
| B-m1 | MINOR | Audit repairs only in an addendum | Merge into the proof |
| B-m2 | MINOR | b_s known is unnecessary at function level | State |
| B-m3 | MINOR | "Breaking mechanisms" overstated beyond trivial ones | Restrict wording |
| B-m4 | MINOR | Lemma T load-bearing for R > 0 | Elementary route |
| B-m5 | MINOR | Homogeneous film only | State limitation |
