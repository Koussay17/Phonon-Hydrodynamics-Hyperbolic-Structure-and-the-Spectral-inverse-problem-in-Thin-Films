# B-revision — counterexample campaign against `proofs/B-revision.md`

Campaign `20261111-fdtr-gauge-memory`, branch CX-R, 12 November 2026. Object: the identifiability claims of
`proofs/B-revision.md` (answers to objections B-2 and B-4 of `12-peer-review.md`), including its post-audit addendum
(D-1 … D-11): rank tables (T3, C2, C4, C5, mixed faces), global injectivity (T2, C1), finitely many frequencies (T4,
Conjectures T4′, T4″) and three-frequency twins (TW3). Method: an independent model written from the definitions (no import
of the revision-branch scripts), 50–110-digit arithmetic, interval (ball) certification, blind searches.
Scripts and raw outputs: `derivations/B-revision-cx/` (§10). Units as in the claims: b_s = 1 and the time unit of the
frequencies; transducer (ξ_A, b_A) = (0.2, 3) unless stated.

---

## 0. Verdicts

| Claim (B-revision) | Verdict | Evidence (this branch) |
|---|---|---|
| Lemma R / Proposition O | SURVIVES | independent model = `twin_fdtr.kernel(k=0)` to 1.5e-15 (stack), 1.0e-15 (bare, T), 1.1e-15 (bare, θ_e); local, entropic, kinetic ε = 0.3, −0.5, 2 (§1) |
| T3, stack rank tables (local 6/5/4; tied ε ≠ 0: 6 except Π∩{R₂=0}: 5/4; free c_n 7/6/5/4) | SURVIVES | 517 random 50-digit points on 55 (model, stratum) cells: every rank as claimed, retained s ≥ 9.2e-6, deficient s ≤ 1.1e-50, claimed kernels annihilated to ≤ 1.8e-51 (§2); blind exact-degeneracy search, 8 tied conventions × 14 faces: 990 exact solutions, all on claimed strata with claimed kernels (§3) |
| R₁ never in a kernel | SURVIVES | parity proof (D_c(1/B̃) is odd in s) and every stratum above |
| C2, bare-film table (read-out T(0), tied and free c_n, R known/unknown) | SURVIVES | 132 + part of the 517 points (§2) |
| C4 (read-out θ_e) as first stated | **BROKEN at one stratum** (error D-1 of the audit) | at ε = 1 on Π∩{b=b_s, R=0} the kernel is span{e_ξ, (0,0,−1,1[,0])} (rank 2 of 4 / 3 of 5): confirmed independently, both vectors annihilated to 1e-54 / 1e-51; for ε ∈ {0.5, 0, −1, 2} only e_ξ (§2.2) |
| C4 corrected (addendum), C5, mixed face coefficients | SURVIVES | 132-point test with the corrected prediction; E₀ and E₋₁ rediscovered blindly; independent hand derivation of the kernel conditions for general (ε₀, ε) gives exactly the listed sets (§3.1). Blind search over mixed faces not run (§9) |
| T2 / C1, global injectivity (function level) | SURVIVES | proof steps re-derived (§5.1); not refutable by finite computation |
| T4 (i)–(iv) | SURVIVES | closed forms of G′ (P, Q) reproduced to 13 digits by an independent Cauchy-integral computation; σ¹⁵ law: ratio −0.82, −0.94, −0.98 at σ = 1e-2, 3e-3, 1e-3 (§4) |
| Conjecture T4′ (three frequencies on {G′ = 0}) | SURVIVES (supported) | 16 points exactly on {G′ = 0}: 40 random triples each, median normalised s_min(J₃) 1e-4 – 2.4e-3 (§4) |
| §3.3(4): two generic frequencies give rank 4 on Δ (entropic) | SURVIVES; PRACTICALLY ILL-CONDITIONED for small ε | median normalised s_min 7e-3 (ε = 1), 7e-4 (ε = 0.084), 1e-4 (ε = 0.01), 2e-52 (ε = 0) (§4) |
| T3/T4 in the FDTR band (Au/AlN/sapphire, note 16) | PRACTICALLY ILL-CONDITIONED | three band frequencies: best normalised s_min 3e-6 (local) / 4e-5 (entropic) over 300 triples; 80 band frequencies: two singular values ≤ 8e-5 (§7) |
| TW3 (three certified local twins) | SURVIVES | reproduced; re-certified 3/3 by an independent Krawczyk implementation (§5.2) |
| "None [three-frequency twin] was found under the entropic convention" (§0.3, §5.4, §8) | **BROKEN** | entropic twins for 501 of 6 000 random targets; 5 Krawczyk-certified (§5.3) |
| Kinetic ε = 0.084 (Beardo-type) | twins exist | certified three-frequency twins for the branch's **own** targets and triples (all three cases), also ε = −0.5, −1 (§5.2) |
| "None for the AlN-like target" (§5.4) | **BROKEN** for the note-16 stack | 11 certified three-frequency twins of the note-16 target (τ_ℓ = 10 ps): 1 entropic, 10 local, one at physically plausible values (τ_R′ = 55 ps, τ_ℓ′ = 19 ps) (§7) |
| Conjecture T4″ (four frequencies, generic targets) | SURVIVES as a generic statement | 24 random targets × 2 000 starts: no exact twin; exceptional pairs form locally a 4-manifold whose projection on the targets has rank 4 (codimension 2) (§6) |
| §7 "four are conjectured to suffice globally"; §0.3 "none with four frequencies" | **BROKEN as worded** | certified four-frequency twin pairs (5 local, 5 entropic) and five-frequency pairs (5 entropic); none found for N = 6 (moderate box) or N = 7 (§6) |
| Statements for ε < 0 (C4(d), E_ε, "every fixed ε ≠ 0", negative-kinetic faces) | **PHYSICALLY VOID** (mathematically consistent) | Lemma S: for ε < 0, τ_R > 0 the bare film has a real pole p₀ > 0 (growing interface mode); stacks: two such poles in 58–60/60 random cases; AlN note-16: growth rates 1.8e14–1.8e18 s⁻¹. Addendum D-2 (no poles on iℝ₊ for ε ≤ 1) is correct but does not exclude these right-half-plane poles (§8) |
| D-3 (local convention: Re Z_T > 0 at k = 0) | consistent | min Re Z/|Z| = 0.013 (bare) / 0.056 (stack), 200 random sets each (§8) |

---

## 1. Independent model and cross-checks (`cxmodel.py`, `cx_01_crosscheck.py`)

Film quadrupole from C p T + q′ = 0, (1 + τ_R p) q = −λT′ + (1 + α)ℓ²q″, i.e. q = −(λ/g)T′ with
g = (1 + τ_R p)/(1 + τ_ℓ p); A = ξ₁ s √g, B = b s/√g, p = s². Faces: ϑ = N T on the GK side, N = 1 + τ p; top factor
N₀ = 1 + τ_o p, bottom factor N = 1 + τ_n p; jump ϑ_above − ϑ_below = R q. Bare film Z = N₀ Mob(M; L/N), L = R + 1/(b_s s);
stack Z = Mob(T_A; R₁ + N₀ Mob(M; L₂/N)). Möbius maps evaluated as (1/B)(1 + rE)/(1 − rE), r = (u − 1)/(u + 1),
E = e^{−2A}. Exact parameter derivatives by forward-mode dual numbers on three back ends (numpy, mpmath, python-flint balls).

Checks: dual derivatives vs 80-digit central differences, worst relative error 8.4e-42; numpy vs mpmath 1.7e-13; model vs
the twin `src/twin_fdtr.py` (`kernel(f, 0, stack)`, SI units, Au 80 nm / AlN 500 nm / sapphire and bare AlN on sapphire;
local, entropic and `GKBoundary(temperature="kinetic", c_n = ετ_ℓ/C)` with ε = 0.3, −0.5, 2; α ∈ {1/3, 2};
(τ_R, τ_ℓ) ∈ {(25,10), (25,25), (0,20)} ps; (R₁, R₂) ∈ {(0,0), (1e-8, 2e-8)} m²K/W; 29 frequencies 10 kHz–100 GHz): worst
relative differences 1.5e-15 (stack), 1.0e-15 (bare, T(0)), 1.1e-15 (bare, `theta_e`).

---

## 2. Rank tables at 50 digits

### 2.1 Random points on every stratum (`cx_02_rank_tables.py`, `cx_02b_bare_free.py`)

Jacobian of log Z at 40 frequencies log-spaced in [1e-3, 1e3], 50 digits, columns normalised (exactly-zero columns, which
occur on r ≡ 0 strata, are counted as kernel directions instead of being normalised); random points drawn independently
of the branch's test points; derived coordinates (k_s, exotic-set values, 1 − 2ε) computed in the working precision.

`cx_02`: 517 points, 0 disagreements. Cells — stack local: off Δ (generic, R₂ = 0, b = b_s & R₂ = 0, τ_R = 0 & R₁ = 0,
τ_ℓ = 0 & R₂ = 0), Δ with R₂ > 0, Δ ∩ {R₂ = 0} (b ≠ b_s and D₂), D₂ with R₂ > 0, Π ∩ {R₂ = 0}; stack tied for
ε ∈ {1, 0.5, 0.084, −0.5, 2.5, −1, 1/3, 2, −3, 1e-3}: interior, Δ (R₂ > 0 and = 0), D₂ ∩ {R₂ = 0}, Π (R₂ > 0, b ≠ b_s and
b = b_s), Π ∩ {R₂ = 0} (b ≠ b_s: rank 5, kernel (k_s,0,0,0,0,1); b = b_s: rank 4, kernel e_ξ and (0,0,1−2ε,1,0,0)),
τ_ℓ = 0 & R₂ = 0 (also b = b_s), τ_R = 0 & R₂ = 0, τ_R = ετ_ℓ & R₂ = 0; stack free c_n: off Δ (τ_n of either sign),
Δ (τ_n ≠ 0, R₂ = 0, τ_n = 0), Δ ∩ {τ_n = 0, R₂ = 0} (b ≠ b_s, D₂), Π ∩ {τ_n = 0, R₂ = 0} (b ≠ b_s; b = b_s: rank 4 with the
three claimed vectors), τ_n = τ_R, τ_n = τ_ℓ (also τ_R = 0), τ_R = τ_n = 0 & b = b_s; stack with mixed (ε₀, ε) ∈
{(0.3,0.7), (−0.4,1.3), (2,−1), (0.5,0), (0,0.5), (0.5,1)}: interior, Δ, Δ ∩ {R₂ = 0}, Π ∩ {R₂ = 0} (b ≠ b_s, b = b_s),
τ_ℓ = 0 & R₂ = 0; exotic sets E₀ ((ε₀,ε) = (1,0)) and E_ε ((1, ε < 0)) in the stack (rank 5, claimed kernels); bare film,
read-out T(0), ε ∈ {0, 1, 0.084, −0.5, 2.5, 0.5}, R known (0, 0.37) or unknown, strata of Table 0.2. Smallest retained
singular value 9.2e-6 (stack, ε = 1e-3, Π ∩ {b = b_s}, R₂ > 0: the lifting is ∝ ε there); largest deficient 1.1e-50;
largest residual of a claimed kernel vector 1.8e-51.

`cx_02b`: 132 points, 0 disagreements — bare film, free c_n (off Δ, Δ with τ_n ≠ 0, Δ and D₂ with τ_n = 0 for R = 0 known,
R > 0 known, R unknown at R > 0 and R = 0) and read-out θ_e with ε ∈ {0, 0.5, 1, −1, 2.5} (off Δ, Δ off the exotic sets,
Δ ∩ {R = 0}, E₀, E_ε, Π ∩ {R = 0}, D₂ ∩ {R = 0} known, Π ∩ {b = b_s, R = 0} known with the corrected prediction of D-1).
Smallest retained 1.1e-4, largest kernel residual 1.3e-51.

### 2.2 The audit correction D-1 (`cx_02c_C4_correction.py`)

Read-out θ_e = (1 + τ_ℓp)T(0), bottom coefficient ετ_ℓ, Π ∩ {b = b_s, R = 0}: for ε = 1 the rank is 2 of 4 (R known) and
3 of 5 (R unknown); e_ξ and (0,0,−1,1[,0]) are annihilated to ≤ 2.5e-54 and ≤ 1.2e-51 (two random ξ₁ each). For
ε ∈ {0.5, 0, −1, 2} the rank is 3 / 4, e_ξ is in the kernel and (0,0,−1,1) is not (residual 1.3e-4 – 4.9e-3). On D₂ with
t > 0, ε = 1, R = 0 known: rank 4. Independent derivation: at Π with R = 0 the condition (O) reduces to (ε₀ − ε)c = 0, which
is vacuous for ε₀ = ε = 1, and the p¹ condition of (E) is vacuous at b = b_s, so c is free with c_{τ_R} = (1 − 2ε₀)c = −c.
The original statement of C4(a) is therefore false at this stratum and the corrected one holds (it is C5 at ε = 1). The
first test set (`cx_02`) did not contain this stratum, so the error was not caught here before the audit.

---

## 3. Degenerate sets: derivation and blind search

### 3.1 Independent derivation

From D_cW ≡ 0 (exponential versus algebraic growth in the sector, re-derived) the kernel conditions for general face
times τ_o = ε₀τ_ℓ, τ_n = ετ_ℓ are (i) c_b = 0 and a/(1+τ_Rp) − c/(1+τ_ℓp) + 2ε₀c/(1+ε₀τ_ℓp) ≡ 0 (a = c_{τ_R}, c = c_{τ_ℓ}),
(O) and (E) as in Lemma K. Coincidences of {τ_R, τ_ℓ, ε₀τ_ℓ}: generic (a = c = 0; residual system degenerate only on
Π ∩ {R = 0} and, for ε = 0, on Δ ∩ {R = 0}); Δ with ε₀ = 0 (a = c, removed by (E) unless ε = 0); Δ with ε₀ = 1 (a = −c:
exotic sets exactly E₀ for ε = 0 and E_ε for ε < 0, none for ε > 0); τ_ℓ = 0 with ε₀ = 1/2 (non-degenerate); Π
(a = (1 − 2ε₀)c; for R = 0 the thin-slab vector, plus (0,0,1−2ε,1) if ε₀ = ε and b = b_s — for the θ_e read-out this is the
D-1 case); r ≡ 0 points (c_{R₁} = 0 by parity in s). This reproduces every degenerate set of §4 and Table 0.3 of the claims
(with D-1) and nothing else.

### 3.2 Blind search (`cx_07b_kernel_solve.py`)

The exact equations G(θ)c = 0, |c| = 1 (G: relative Jacobian at 25 frequencies in [1e-3, 1e3], c ∈ ℝⁿ all directions, θ
on a face: interior, R₂ = 0, τ_R = 0, τ_ℓ = 0, Δ, Δ ∩ {R₂ = 0}, b = b_s, b = b_s ∩ {R₂ = 0}, τ_ℓ = 0 ∩ {R₂ = 0},
τ_R = 0 ∩ {R₂ = 0}, Π, Π ∩ {R₂ = 0}, D₂, D₂ ∩ {R₂ = 0}) are solved by batched Levenberg–Marquardt in (log θ_free, c), with
c orthogonal to the claimed generic kernel of the face (projector recomputed at every iterate); 300 walkers per face; box
ξ₁, b ∈ [0.2, 5], times ∈ [0.01, 5] (or 0), resistances ∈ [0.01, 3] (or 0). Near-degenerate box corners stall at residual
≳ 1e-7; exact degenerate points converge to ~1e-16. Each exact solution is classified against the claimed kernel of the
stratum on which it actually lies.

Validation: on the Δ face of the (ε₀, ε) = (1, 0) stack the solver finds 14 points with t/(4ξ₁²(1−b²)) =
R₂/(2ξ₁(1−b²)/b) = 1.00000000 (E₀) and kernel (1, 0, −t/ξ₁, t/ξ₁, 0, (b²−1)/b); for (1, −1) one point with b/√2 =
t/(2ξ₁²) = R₂/(√2ξ₁) = 1.00000000 (E₋₁). A first design that minimised log s_r directly was discarded after failing this
validation (§9).

Result (tied conventions ε ∈ {0 (local), 1, 0.084, 0.5, 2.5, −0.5, 1/3, −2}, 112 faces): 990 exact solutions, every one on a
claimed stratum with a claimed kernel (local: Δ, Δ ∩ {R₂ = 0}, D₂ ∩ {R₂ = 0}; tied: Π ∩ {R₂ = 0} including its b = b_s
substratum). One solution was flagged only because b − 1 = 3.7e-9 exceeded the classifier tolerance; its kernel is
(0,0,1−2ε,1,0,0) for ε = −2, the claimed b = b_s vector. No missing degenerate set.

---

## 4. Low-frequency minors and finitely many frequencies (`cx_08`, `cx_08b`, `cx_14`)

Laurent coefficients of W − R₁ − 1/(b_s s) by trapezoidal Cauchy integrals in s (160 nodes, two radii inside the nearest
singularity — for b > b_s a pole lies on the negative real s-axis near −atanh(b_s/b)/ξ₁ — agreeing to ~1e-22), with exact
parameter derivatives. G′₁₂₃₄₅ = det[∇c₁; …; ∇c₅] (columns ξ₁, b, τ_R, τ_ℓ, R₂) equals the record's 2(τ_ℓ−τ_R)ξ₁⁷P/(135b²)
(local) and −2ξ₁⁷Q/(2025b²) (entropic) to all 13 printed digits at 6 random points. Leading order:
det J₃(σ², 4σ², 9σ²)/(σ¹⁵ t(2,3) G′) = −0.8217, −0.9430, −0.9806 at σ = 1e-2, 3e-3, 1e-3.

On the hypersurface: 8 local points with P = 0 and 8 entropic points with Q = 0, off Δ (|G′| ≤ 4e-45 there): function-level
rank 6 (normalised s_min 2.7e-3 – 4.1e-2), and for 40 random triples in [10^−1.5, 10^1.5]³ the median normalised s_min of
J₃ is 1.1e-4 – 2.4e-3 at every point (minimum over triples 1e-10 – 1e-5: isolated bad triples). Conjecture T4′ is supported
on the set that T4(iv) leaves open.

Proposition LF: c₋₁ … c₂ (c₂ at R = 0) agree with the closed forms to ≤ 7.7e-48 at 8 random points with b_s = 1.3 and
independent face times τ_o, τ_n ∈ [−0.3, 0.8]; hence also δ_ε = τ_R − τ_ℓ + ε(1 + b²/b_s²)τ_ℓ.

Two frequencies on Δ, bare film, read-out T(0), R known, 30 samples per case: median normalised s_min of the 4 × 4
Jacobian 6.8e-3 / 9.1e-3 (ε = 1, R = 0 / 0.3), 6.7e-4 / 9.0e-4 (ε = 0.084), 9.1e-5 / 1.3e-4 (ε = 0.01), 2e-52 (ε = 0). The
lifting is roughly linear in ε.

---

## 5. Twins

### 5.1 Global statements (T2, C1)

Re-derived: the two cases of Lemma G; r ≡ 0 iff R₂ = 0, b = b_s, τ_n = 0, g = 1; the step "ρ ≠ 0 ⇒ p(1+τ_R′p)(1+τ_ℓ′p) is a
square in ℂ[p]" (impossible, simple factor p for τ ≥ 0), with the sign of D-6. Z(s²) depends on g only (M is invariant
under (A, B) → (−A, −B)) and is meromorphic off p = −1/τ_ℓ. No gap found.

### 5.2 TW3 and continuation in the face coefficient (`cx_03`, `cx_04`, `cx_09`)

The three printed local twins are roots of the independent model (60-digit Newton moves them by ≤ 3e-16, the float
rounding of the data; s_min/s_max of J₃ 0.9e-3 – 3.2e-3). Independent Krawczyk certification (python-flint, 400 bits,
exact decimal data, box radius 1e-40): 3/3 certified (K(X) − m within 2.5e-34 – 8.9e-34 of the radius); a control with the
centre moved by 1e-30 fails each time.

Pseudo-arclength continuation of each twin in ε (ε₀ = ε, target and triple fixed): no branch reaches ε = 1 (folds at
ε ≈ 0.16 – 0.48 or exit through τ_ℓ′ < 0), but certified twins exist on the way (all seven certified, K-deviation ≤ 2.2e-33):

| target, triple | ε | certified twin θ′ |
|---|---|---|
| (0.8, 2, 0.5, 0.1, 0.3, 0.2), (1/2, 4, 1.7) | 0.084 | (1.10376369045, 1.74298509244, 0.315156252333, 0.0117022360415, 0.263782475432, 0.355917217032) |
| same | −0.5 | (0.642061797341, 2.24989085416, 0.627532049223, 0.118294343716, 0.327644046687, 0.126118125529) |
| same | −1 | (0.836134542759, 1.96318897568, 0.478523538203, 0.097983880452, 0.296315559875, 0.221272633378) |
| (0.8, 2, 0.5, 0.1, 0.3, 0.2), (0.2, 1.1, 6) | 0.084 | (0.818048484592, 1.97926677060, 0.484553080777, 0.093408229977, 0.297690723022, 0.208569009512) |
| same | −0.5 | (0.396703596359, 3.26565808277, 1.20228918938, 0.161915506222, 0.386842240706, 0.083576263001) |
| same | −1 | (0.521201547175, 2.63483326872, 0.809998299979, 0.128522761379, 0.357689517157, 0.115019872330) |
| (0.8, 1, 0.5, 0.1, 0.3, 0.2), (0.2, 1.1, 6) | 0.084 | (0.551743515020, 1.19075140647, 0.821981273341, 0.208919933046, 0.366939227572, 0.322890777290) |

### 5.3 Entropic three-frequency twins (`cx_05`, `cx_06`, `cx_09`)

Target-fixed multistart (20 000 starts per case, log-box ±3, half of the starts with times in [1e-4, 10]) reproduces the
branch on its targets: local twins = the certified ones, entropic none. The joint search over random targets (one start per
target, N = 3, ω = (1/2, 4, 1.7), targets log-uniform within ±1.5 of (0.8, 2, 0.5, 0.1, 0.3, 0.2)) finds entropic twins for
**501 of 6 000** targets; the first 40, refined at 60 digits, are admissible, regular (s_min/s_max 2.5e-6 – 1.8e-2) and
separated (log-distance 0.12 – 4.65; mismatch at ω = 1.3 between 6e-7 and 1.2e-1). Five certified with exact decimal targets:

| target θ* | certified entropic twin θ′ | K-dev/r |
|---|---|---|
| (0.65317179, 3.8887279, 0.74853435, 0.12214185, 0.98266289, 0.36116796) | (1.85052490194, 1.12065802009, 1.57906214503, 0.072444174454, 0.044478046806, 0.468190039578) | 6.9e-35 |
| (3.3176622, 0.82406472, 1.5634862, 0.052208832, 0.28493317, 0.073053301) | (1.88081464320, 0.725256624275, 1.65349573648, 0.053386567846, 0.038040316298, 0.981297584778) | 1.5e-34 |
| (2.0534403, 1.1653533, 0.54494992, 0.12324068, 0.14450787, 0.23678845) | (2.34172298803, 1.17473949510, 0.511665311988, 0.131488833601, 0.164672308566, 0.509417082634) | 2.7e-34 |
| (3.5413392, 0.50406275, 1.798024, 0.28152833, 0.35388833, 0.16745881) | (1.68978887685, 0.543271060143, 1.22540502038, 0.407109810513, 0.823341232554, 2.39173684587) | 2.8e-34 |
| (1.6976177, 4.1601777, 1.1677385, 0.022621446, 0.27572095, 0.11675949) | (1.88738449283, 4.21276181249, 1.01136458487, 0.029119441136, 0.293539072913, 1.96793099713) | 7.9e-34 |

HEURISTIC (mechanism). The entropic ϑ-quadrupole has B_ϑ² = b²p/((1+τ_Rp)(1+τ_ℓp)), symmetric under τ_R ↔ τ_ℓ, while
A² = ξ₁²p(1+τ_Rp)/(1+τ_ℓp) is not. In the 40 refined pairs, median |b′/b − 1| = 0.018, median τ_Rτ_ℓ/(τ_R′τ_ℓ′) = 1.001,
median (τ_R+τ_ℓ)/(τ_R′+τ_ℓ′) = 1.033, and 12 of 40 twins reverse the order of τ_R and τ_ℓ: the twins nearly preserve the
symmetric data, which makes them easy to find once targets are not restricted to three hand-picked points.

---

## 6. Four, five, six, seven frequencies (`cx_06`, `cx_09b`, `cx_13`)

Joint pair search with 12 − 2N target coordinates fixed (random values), square system of 2N real equations, batched LM,
double-precision regularity filter (s_min/s_max > 1e-9), resonance-continuum pairs removed, 60-digit refinement.

| N | ω | local | entropic |
|---|---|---|---|
| 4 | (1/2, 4, 1.7, 0.2) | 5 regular admissible pairs off Δ, **5/5 certified** (8 × 8 pair system, K-dev 2.7e-34 – 6.5e-33) | 131 regular pairs; 40/40 refined admissible, residual ≤ 1.3e-53, s_min/s_max 8e-5 – 9e-3; **5/5 certified** (K-dev 6e-35 – 4e-33) |
| 5 | (…, 9) | 12/12 refined admissible, regular | 39/39 refined admissible, regular (s_min/s_max 3.5e-8 – 5.6e-3); **5/5 certified** (K-dev 8e-34 – 1.4e-31) |
| 6 | (…, 0.05) | none (moderate box [0.05, 10]) | none in the moderate box; in the wide box [4.5e-5, 148] only singular asymptotic artefacts |
| 7 | (…, 25) | — | no regular separated pair (residual 4e-16 only at singular asymptotic points) |

Examples (certified): local N = 4, (0.70571034, 0.57425376, 0.77016345, 1.09323666, 0.51951975, 0.11944515) ↔
(4.45228361, 0.86254549, 5.59953574, 3.04126006, 0.49060604, 0.28199605); entropic N = 4, (3.99926879, 0.71592362,
2.17167053, 0.06047942, 0.10704822, 0.91221724) ↔ (1.83934075, 1.16219434, 0.41199756, 0.30734718, 1.45265093, 2.28937002),
off-data mismatch 9.2e-3; entropic N = 5, (1.84650892, 0.60868786, 0.80192239, 0.04039468, 0.29556685, 0.87330234) ↔
(0.88970241, 1.06300276, 0.14225427, 0.29215297, 1.18713817, 1.91760404). Several pairs again nearly swap τ_R and τ_ℓ.

Generic targets (`cx_13`): 12 random entropic and 12 random local targets (log-uniform within ±1.5 of
(0.8, 2, 0.5, 0.1, 0.3, 0.2)), 2 000 starts each, four frequencies: no exact twin (best distinct relative residual 2.7e-8
and 1.5e-7 — near-twins only). Local structure of the exceptional set (10 entropic pairs, 1 local): the 8 × 12 pair
Jacobian has rank 8, the pair set is a 4-manifold, and its tangent space projects onto the target coordinates with rank 4.
So the targets that have four-frequency twins form, locally, a 4-dimensional set (codimension 2): T4″ survives as a generic
statement, while the unqualified wording of §7 and §0.3 of the claims is false. The pattern 12 − 2N (pairs of dimension 6,
4, 2, 0 for N = 3, 4, 5, 6; none generically for N ≥ 7) is borne out.

---

## 7. AlN note-16 targets (`cx_05 aln`, `cx_09`, `cx_16`, `cx_12`, `cx_15`)

Model units: t₀ = e²/a = 3.1946 ns, b_s = 1. Target θ* = (1, 2.07026, 0.0078257, 0.0031303, 1.82200, 3.64399) (τ_R = 25 ps,
τ_ℓ = 10 ps, R₁ = 1e-8, R₂ = 2e-8 m²K/W), Au transducer (ξ_A, b_A) = (0.18236, 1.87668). Triples: band1 = (1, 20, 150) MHz,
band2 = (0.1, 5, 60) MHz, T1 = (0.5, 4, 1.7) in model units (25, 85, 200 MHz).

**Three-frequency twins of the note-16 target** (target-fixed multistart, 20 000 starts per case; all 11 roots certified by
Krawczyk, K-dev 2.9e-35 – 2.2e-31, controls fail):

| triple | convention | certified twin θ′ (ξ₁, b, τ_R, τ_ℓ, R₁, R₂) | s_min/s_max | off-data mismatch |
|---|---|---|---|---|
| band1 | entropic | (1.86018926708, 1.16805931703, 0.29181896725, 0.06724668870, 1.41624313478, 3.05800398604) | 6.4e-3 | 2.3e-4 |
| band1 | local | (0.98425370467, 2.10216389580, 2.16689347649, 2.17936083906, 1.82767359510, 3.65176840326) | 1.6e-5 | 2.4e-6 |
| band1 | local | (0.98788430723, 2.09331135065, 19.0249850502, 19.2841892436, 1.82762207208, 3.64562245329) | 4.6e-6 | 1.8e-5 |
| band2 | local | 4 twins, e.g. (1.02187429772, 2.03391851765, 4.53192757148, 4.87767224537, 1.82764281723, 3.62710765676) | 9e-8 – 6e-5 | 8e-8 – 2e-5 |
| T1 | local | (1.02619311366, 2.02101345673, 0.01718617252, 0.00591067237, 1.81394633879, 3.63573087889) | 3.3e-5 | 5.8e-6 |
| T1 | local | 3 further, e.g. (0.97448925301, 2.12085330056, 0.04035473622, 0.04222999679, 1.82967653533, 3.65189463094) | 1.4e-5 – 6.3e-5 | 6e-6 – 3e-5 |
| band2, T1 | entropic | none found | — | — |

The local T1 twin is physically plausible: τ_R′ = 54.9 ps, τ_ℓ′ = 18.9 ps against 25 and 10 ps, with ξ₁, b, R₁, R₂ changed
by +2.6 %, −2.4 %, −0.44 %, −0.23 %; the two responses differ by 6e-6 (relative) away from the data. Most other local twins
lie near the resonance continuum (τ_R′ ≈ τ_ℓ′, large); the entropic twin is a different film (ξ₁′ = 1.86, b′ = 1.17). The
branch's "none for the AlN-like target" concerns a different target (R₁, R₂ about 36 times smaller, transducer (0.2, 3))
and does not carry over to the note-16 stack.

AlN neighbourhood (`cx_15`, ξ₁ = 1 and b = 2.07026 fixed, τ's within [0.25, 4] and R's within [0.3, 3] times the note-16
values, Au transducer, band1, entropic): 151 regular separated pairs; the refined ones are admissible and regular, with twin
times τ_R′ ≈ 0.57 – 0.98 ns, τ_ℓ′ ≈ 0.19 – 0.23 ns and off-data mismatch ≈ 1e-3. Homotopy (`cx_12`) of 20 entropic
random-target pairs to the note-16 target, transducer and band1: none arrives as a twin (folds or exits; one merges with the
target).

**Conditioning in the band** (`cx_16`, 50 digits, columns ∂ log Z/∂ log θ): best normalised s_min of J₃ over 300 random band
triples 3.2e-6 (local, τ_ℓ = 10 ps) and 4.4e-5 (entropic); medians 6e-10 and 1.4e-8. With 80 band frequencies the singular
values relative to s₁ are 1, 0.58, 0.24, 1.3e-2, 7.9e-5, 1.9e-6 (local) and 1, 0.58, 0.24, 1.3e-2, 7.6e-5, 2.9e-5 (entropic);
the weakest direction is (τ_R, τ_ℓ) ∝ (−0.38, −0.92) (local) and essentially τ_ℓ (entropic). On Δ (τ_ℓ = 25 ps) the local
s₆ is 2.6e-52 (exact resonance), the entropic one 7.2e-5. Hence T4 ("three generic frequencies suffice") and "rank 6 on Δ
under the entropic convention" are structurally true but PRACTICALLY ILL-CONDITIONED in the band: three frequencies would need
relative data accuracy ≲ 1e-6, and with 80 frequencies two directions lie below 1e-4.

---

## 8. Stability and passivity of the conventions (`cx_10`, `cx_10b`, `cx_10c`)

**LEMMA S (instability for negative interface coefficients).** Bare film, read-out T(0), τ_n = ετ_ℓ with ε < 0,
τ_ℓ > 0, τ_R > 0; any ξ₁, b > 0, R ≥ 0. Then Z has a pole at a real p₀ > 0.
*Proof.* For real s > 0, A, B, L > 0 and Z = (L + N tanh A/B)/D with D(s) = B L tanh A + N. D is continuous on (0, ∞),
D(0⁺) = 1 (BL → b/b_s, tanh A → 0); since 0 < tanh A < 1 and BL = bs(R + 1/(b_s s))/√g = O(s) (g → τ_R/τ_ℓ > 0),
D(s) ≤ BL + 1 − |ε|τ_ℓs² → −∞. Hence D(s₀) = 0 for some s₀ > 0, where the numerator equals L(1 − tanh²A) > 0. □
Thick-film limit (τ_R = τ_ℓ = τ, R = 0, ξ₁ → ∞): p₀ → (1 + b/b_s)/(|ε|τ), confirmed to 6 digits for ε ∈ {−0.05, −0.2,
−0.5, −1}, b ∈ {0.5, 2}, ξ₁ ∈ {1, 3, 10, 30}; verified directly on Z at 40 digits (|Z| ∝ 1/|p − p₀|). With τ_R = 0 and R > 0
no real pole was found (the argument does not apply).

PHYSICAL INTERPRETATION. A zero of D at real p₀ > 0 is a nontrivial solution of the homogeneous problem (q(0) = 0) with time
factor e^{p₀t}, decaying into the substrate like exp(−(p₀/a_s)^{1/2}z): a localised growing interface mode. In the time domain
|ε|τ_ℓ ∂_tT(e) = T(e) − T_s − Rq(e) drives the film-side face temperature away from T_s + Rq (anti-damped relaxation), whereas
for ε > 0 the same law relaxes it towards T_s + Rq. Mechanism in the energy identity: the bottom boundary factor
Re[(1+τ_ℓp̄)(R + 1/S̄)/(1+ετ_ℓp̄)] is negative on the real axis p > 1/(|ε|τ_ℓ) when ε < 0.

Census (argument principle on |arg s| < π/4 − 1e-3, |s| ∈ [1e-3, 3e3]; 60 random parameter sets per case): right-half-plane
poles in 59–60/60 bare films and 58–60/60 stacks (two real poles, one per film face) for each ε ∈ {−2, −1, −0.5, −0.2, −0.05};
none for ε ∈ {0, 0.084, 0.3, 0.5, 1, 1.5, 2.5, 5}. AlN note-16: growth rates p₀/t₀ of 2.4e11 (ε = −1) … 2.3e13 s⁻¹ (ε = −0.01)
for the bare film with R = 0, 1.8e14 … 1.8e18 s⁻¹ for the stack; the rate diverges like 1/|ε| as ε → 0⁻.

Consequence: for ε < 0 the frequency-domain response is the analytic continuation of the transfer function of an
exponentially unstable system; no steady state exists to be measured. Proposition C4(d) (E_ε), the negative-kinetic exotic
strata, the ε = −0.5 and −1 checks, and statements quantified over "every fixed ε ≠ 0" are correct for the formal transfer
function but have no measurement meaning for ε < 0. Addendum D-2 (no poles on iℝ₊ for every ε ≤ 1, bare film, and
0 ≤ ε ≤ 1, both-face stack) is correct — an independent derivation of the T-identity agrees, and no axis pole was found for
ε ≥ 0 — but it does not exclude poles in the open right half-plane, which exist for every ε < 0 (bare film, τ_R > 0). The
physically admissible range of the 1D model is ε ≥ 0; the sign of β − χ_nn therefore matters.

Passivity on the axis (200 random sets, ω ∈ [1e-4, 1e4]): stack min Re Z/|Z| > 0 for ε ∈ [0, 2.5] (0.037 – 0.113), −0.0024
for ε = 5, −1.0 for ε = 20; on the axis the phase of the bottom factor at R = 0, π/4 + atan(εωτ_ℓ) − atan(ωτ_ℓ), exceeds π/2
for some ω iff ε > (1+√2)² ≈ 5.83 (bottom-face term alone). Bare film with read-out T(0): non-passive for every ε > 0 tested
(−0.25 to −0.46; T(0) of a GK top layer is not the passive variable), positive for ε = 0 (0.013), consistent with D-3.

---

## 9. Attacks not performed or incomplete

1. **Interrupted runs (machine memory pressure).** Not run: the blind exact-degeneracy search over **mixed** face
   coefficients and over the **free-c_n** faces (`cx_07b mixed`, `cx_07b free`; mixed and free strata are covered only by the
   random-point tests of §2 and the derivation of §3.1); the note-16 targets with τ_ℓ = 0 and 25 ps (`cx_05 aln2`); the
   AlN-neighbourhood pair search at N = 4 (`cx_15 entropic 4`); the local homotopy and the band2 homotopy (`cx_12`); the
   practical misfit profiles (`cx_11`); random frequency triples for the branch's generic target (`cx_05b`). The random-target
   test of T4″ was run with 12 targets per convention instead of 30.
2. **No re-run of the branch's quantifier elimination.** The degenerate sets were attacked by the hand case analysis and by
   the blind solver, limited to the box of §3.2; sets lying only outside it would be missed, and codimension ≥ 3 subsets of a
   face are found with lower probability (absence not certified).
3. **Function-level injectivity** (T2, C1) cannot be refuted by finite computation.
4. **Twin searches** used a few frequency sets per N; absence of twins (branch targets in the entropic convention, N = 6, 7)
   is an observation about the searched boxes, not a proof.
5. **k > 0 (finite beams)** not examined; the instability of §8 is shown at k = 0 only.
6. **Certification** relies on python-flint balls (acb sqrt/exp) and on dual numbers propagated through balls for the
   interval Jacobian; the 50-digit ranks are floating computations with a ≥ 45-decade gap, not certified bounds.
7. The first design of the adversarial singular-value search (`cx_07_adversarial_smin.py`, direct minimisation of log s_r
   or s_r/s_{r−1}) is attracted by near-degenerate box corners and missed E₀ in validation; it is not used as evidence.

---

## 10. Files (`derivations/B-revision-cx/`)

| File | Content |
|---|---|
| `cxmodel.py`, `cxtwin.py`, `cxcont.py`, `cxkraw.py`, `throttle.py` | independent model (dual numbers, numpy/mpmath/flint); twin solvers; pseudo-arclength continuation; Krawczyk test; job throttle |
| `cx_01_crosscheck.py/.json` | derivative, back-end and twin cross-checks (§1) |
| `cx_02_rank_tables.py/.json/.log`, `cx_02b_bare_free.py/.json/.log`, `cx_02c_C4_correction.py/.json/.log` | 50-digit rank tables; D-1 check (§2) |
| `cx_07b_kernel_solve.py`, `cx_07b_validate.json`, `cx_07b_tied2.json` | blind degeneracy search (§3); `SUPERSEDED_cx_07b_tied.*` were written concurrently by two runs and are invalid |
| `cx_07_adversarial_smin.py` | discarded first design (§9) |
| `cx_08_three_freq.py/.json/.log`, `cx_08b_laurent.py/.json`, `cx_14_two_freq_entropic.py/.json` | low-frequency minors, T4′, LF, two frequencies (§4) |
| `cases_tw3.json`, `cx_03_twin_continuation.*`, `cx_04_eps_arclength.*` | TW3 reproduction and ε-continuation (§5.2) |
| `cx_05_multistart.py`, `cx_05_multistart_branch.json/.log`, `cx_05_multistart_aln.json` | target-fixed multistart (§5.3, §7) |
| `cx_06_joint_pairs.py`, `cx_06_joint_{entropic,local}_N{3,4,5,6,7}*.json/.log` | joint pair searches (§5.3, §6) |
| `cx_09_certify.py`, `cx_09_certified_{entropic,kinetic,aln}_N3.json`, `cx_09b_certify_pairs.py`, `cx_09b_certified_{local,entropic}_N4.json`, `cx_09b_certified_entropic_N5.json`, `cases_*.json` | Krawczyk certification of twins and pairs |
| `cx_10_stability.py/.json/.log`, `cx_10b_poles_verify.py/.json/.log`, `cx_10c_instability_lemma.py/.json/.log` | right-half-plane poles, Lemma S, passivity (§8) |
| `cx_12_homotopy_to_aln.py`, `cx_15_aln_region_pairs.py`, `cx_16_aln_conditioning.py/.json` (+ `queue_*.log`) | AlN homotopy, neighbourhood pairs, band conditioning (§7) |
| `cx_13_four_freq_random_targets.py`, `cx_13_four_freq_{entropic,local}.log`, `cx_13_four_freq_*.json` | T4″ random targets and exceptional-set structure (§6) |
| `run_queue.py`, `run_queue2.py`, `queue_*.log` | sequential job queues used for the runs, with their logs |
