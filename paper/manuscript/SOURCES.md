# Paper 1 manuscript — source map

Draft of 11–12 November 2026 (`main.tex`, `sections/*.tex`, `supplement.tex`, `refs.bib`), revised on
12 November 2026 to apply the repairs of the proof audit `checks/proof_audit_sections_II_IV.md` (cited below as
"audit §n"; its scripts are in `checks/audit/`).
Every equation, statement and number of the manuscript is listed below with the file it comes from.
Paths are relative to the repository root `C:\Users\Koussay\these` unless they start with `ResearchLab/`
(`C:\Users\Koussay\ResearchLab`). "Derived here" means that the manuscript writes out a derivation that is
not in the sources; each such item is checked in `paper/manuscript/checks/verify_manuscript_algebra.py`.

Printed numbers refer to the build of 12 November 2026 after the audit repairs (15 pages; numbering unchanged
by the repairs). Statement numbering in the paper: Theorems 1–3, Propositions 1–6, Lemmas 1–2, Corollaries
1–2, Definitions 1–5, Examples 1–2, Remarks 1–5, Eqs. (1)–(40), hypotheses (K1)–(K3), (H1)–(H5), (F1)–(F6).
Supplemental Material: Lemma S1, Eqs. (S1)–(S7), Tables S1–S2 (Secs. S2 and S4 written, Sec. 12 below).

## 0. Build and checks

| Item | Command / file | Result |
|---|---|---|
| Build | from `paper/manuscript/`: `pdflatex -output-directory=.build main.tex`; in `.build/`: `BIBINPUTS="..;." bibtex main`; then `pdflatex -output-directory=.build main.tex` twice more; copy `.build/main.pdf` to `main.pdf` (`.build/` is ignored by git) | 15 pages; 0 overfull boxes; 1 unavoidable class warning (Sec. 11) |
| Supplement | `pdflatex -output-directory=.build supplement.tex` (twice); copy `.build/supplement.pdf` to `supplement.pdf` | 4 pages (Secs. S2 and S4 written, the others headings only); 0 overfull boxes, no warnings |
| Proof audit | `checks/proof_audit_sections_II_IV.md`; scripts and outputs in `checks/audit/` (run each with `python -B <file> > <file>.out.txt` from `checks/audit/`) | 88 of 88 checks pass; every repair of its §12 is applied in this revision (rows below) |
| Record size | `checks/trim_records.py` keeps only the bibliographic fields of registry answers (applied by `verify_dois.py`) | records about 150 kB in total |
| Algebra | `checks/verify_manuscript_algebra.py` → `checks/verify_manuscript_algebra.out.txt` | 74 of 74 checks pass (SymPy, exact rationals; 50-digit eigensolvers; two numerical phase checks). Added with the audit repairs: Callaway spectrum (K3); Lemma 2 sharpness for ς_max = 2, 3, 4 and the mixed family; cut-off moment ratios for α ≥ d |
| Bibliography | `checks/verify_dois.py` reads every DOI of `refs.bib` → `checks/doi_records.json`, `checks/doi_records.txt`; `checks/compare_bib.py` → `checks/compare_bib.out.txt` | 62 of 62 DOIs resolve to the intended works; 0 discrepancies in year, volume and first page |
| Bibliography (history) | first query list → `checks/doi_records_initial.json/.txt`; `checks/search_titles.py` → `checks/title_search_records.json`; `checks/make_bib.py` drafted entries from these records | 62 DOIs queried first: 61 resolved, of which 3 were unrelated works (guesses, replaced), and 1 did not exist (Sec. 9) |
| Word counts | `checks/word_count.py` → `checks/word_count.out.txt` | see report |

Required LaTeX packages added to MiKTeX on 11 November 2026: `amscls`, `natbib`, `textcase`.

## 1. Front matter

| Manuscript | Source |
|---|---|
| Title (the central question) | `paper/RESEARCH_PLAN.md`, "Candidate central question"; `paper/PRIOR_ART.md`, "Central contribution" |
| Author name | `git log` of the repository (commit author). Affiliation and author list: PENDING (`paper/RESEARCH_PLAN.md`: authorship to be agreed) |
| Abstract, sentences 1–2 | `paper/PROGRAMME.md` (paper 1 scope); `paper/RESEARCH_PLAN.md` |
| Abstract, Cattaneo time = first memory moment | `paper/PRIOR_ART.md`, central contribution 1; `theory/aln/memory_gap/02-assumptions.md` §A ("Moments") |
| Abstract, lower bound and PSD unboundedness | Prop. 2 and Theorem 1 below |
| Abstract, cone bound, prior | Theorem 2, Prop. 5 below |
| Abstract, infrared | Theorem 3 below; "explicit infrared hypotheses" (plural) because Theorem 3(ii) rests on (H4) and on the uniform form of (H1) (audit §1.14) |
| Abstract, gauge | Prop. 6, Cor. 2 below |

## 2. Section I (Introduction)

| Statement | Source | Notes |
|---|---|---|
| TDTR/FDTR as standard probes | Refs. Cahill2004, Schmidt2009 (DOIs verified) | general background |
| Frequency-dependent apparent conductivities; geometry dependence | KohCahill2007, Regner2013, WilsonCahill2014 | titles verified in registry records |
| Inversion for mean-free-path / relaxation-time spectra | Minnich2012, Forghani2016 (`theory/aln/fc3_memory/13-final-report.md` §5); Chen2025 (`ResearchLab/orchestration/campaigns/20261002-paper1-prior-art/08-literature.md` item 23) | Chen2025 is an arXiv preprint (DataCite DOI) |
| Cattaneo/GK laws used to interpret thermoreflectance | Beardo2020, HennessyMyers2021 (`notes/15_Comparaison_travaux_anterieurs.tex` §§2–3); CamachoDeLaRosa2025 (`notes/camacho_2025.tex`; `notes/conventions.md` §14) | wording "interprets … through these parameters" (not "fits") because Hennessy–Myers analyze the model |
| AlN: about 321 W m⁻¹K⁻¹ at room temperature, TDTR, films of 18 and 22.5 µm; not transferable to 500 nm | Cheng2020; local copy `notes/cheng2020aa.pdf` (p. 3: "Our measured thermal conductivity (321 W m-1 K-1)"; sample thicknesses 18 and 22.5 µm); `notes/conventions.md` §11.1 | |
| First-principles workflows benchmarked | McGaughey2025; `CORRECTIONS.md` (2 Oct 2026, AlN datasets) | |
| Callaway descriptions of AlN examined | MaLiLuo2014; `notes/conventions.md` §11.2, §11.5 | |
| Three-step route | `paper/PROGRAMME.md` (paper 2 chain) | |
| Conserving closure, coefficient 1/3, τ_c form, attribution | `paper/PRIOR_ART.md` row C1; prior-art trail `08-literature.md` items 1–3; `notes/14_Derivation_fermeture_cinetique.tex` §1 | C1 is background, not claimed |
| Coefficient 2 attributed to approximate inverse of the normal collision matrix | prior-art trail item 1 (Guo & Wang 2022, "source says"), item 3 (Hardy & Albers, abstract only) | Sendra et al. not cited for this clause: the statement was read in the Sendra thesis, not in the 2022 article |
| Struchtrup (2025) | prior-art trail item 5 (abstract only) | cited as "see also" |
| Relaxons; parity of conductivity and viscosity | prior-art trail items 10, 21; `notes/18_Audit_fermeture_spectrale.tex` §11 refs 1–2 | |
| Exact lumping; projection memory | `paper/PRIOR_ART.md` row C3; prior-art trail items 26–27; `notes/21_Degenerescences_et_fermeture.tex` §3 | C3 is background, not claimed |
| Zeng & Dong (2026) | prior-art trail item 22 (abstract) | cited only as an existing memory-kernel framework |
| Diagonal-rate descriptions can fail | prior-art trail item 21 (Malviya & Ravichandran 2023, abstract) | |
| Krapez & Rigollet: non-identifiability whatever frequencies and noise | `notes/07_Note_partie_IV.tex` §1 ("inconditionnelle"); prior-art trail item 33 | |
| Contribution 1 (Cattaneo time = τ_mem) | `theory/aln/memory_gap/02-assumptions.md` §A; `paper/PRIOR_ART.md` central contribution 1 | stated as a definition |
| Contribution 2 strength ("elementary … Stieltjes …; explicit statement not found") | `paper/PRIOR_ART.md` row C2a; prior-art `13-final-report.md` §2.1 C2(a) | |
| Contribution 3 strength ("inequality known (Ben-Tal & Teboulle 1990); application potentially novel within the literature searched") | `paper/PRIOR_ART.md` row C2b; prior-art `13-final-report.md` §2.1 C2(b) | |
| Contribution 4 (no novelty claimed; Vermeersch precedent) | `theory/aln/fc3_memory/13-final-report.md` §5 ("No novelty is claimed"); prior-art `13-final-report.md` §2.1 C2 and item 25 | precedent concerns free-path moments, worded as "analogous"; "explicit hypotheses on the infrared rates and solution" names the uniform (H1) and (H4) of Theorem 3(ii) (audit §1.14) |
| Contribution 5 strength ("partially known; extension potentially novel") | `paper/PRIOR_ART.md` row C5; prior-art `13-final-report.md` §2.1 C5 | |
| Fourier resonance named by Both et al., analyzed by Kovács | `paper/PRIOR_ART.md` (footnote under the table); `CORRECTIONS.md` (3 Oct 2026); prior-art trail item 12 | |
| Status labels | `README.md` (status vocabulary) | "proved" = complete proof given |
| Literature scope: Crossref and OpenAlex (about fifty queries each), arXiv full text, publisher/repository pages, forward citations of Ben-Tal–Teboulle and Stewart, October 2026 | `paper/PRIOR_ART.md` first paragraph; prior-art `13-final-report.md` §2 ("Search scope"); trail §1 (50/52 queries) | |

## 3. Section II (Linearized kinetics)

| Item | Source | Status / notes |
|---|---|---|
| (K1) mesh, modes, Γ modes excluded | `theory/aln/memory_gap/02-assumptions.md` §A; `notes/17_AlN_spectral_temperature.tex` §2 (weights 1/(N_q V_0)) | |
| (K2) detailed balance, symmetric PSD | `notes/18_Audit_fermeture_spectrale.tex` §2; Ziman2001 (variational formulation) | |
| Eq. (1) entropy coordinates | `memory_gap/02-assumptions.md` §A (y = D⁻¹dn) with the volume normalization of `notes/17` §4 (y = √W δn/√(n(1+n))) | notation choice (Sec. 10) |
| Eq. (2) | `notes/18` §2 | |
| Eq. (3) second-order entropy, entropy production | derived here; `CORRECTIONS.md` 15–16 Sep ("entropy production is quadratic near equilibrium") | dimension checks in the check script |
| (K3) energy only invariant | `memory_gap/02-assumptions.md` §A ("No other invariant is assumed") | |
| Eq. (4) current; time-reversal parity, j ⊥ e | `notes/17` §4 (b_{μi}); `memory_gap/02-assumptions.md` §A ("b orthogonal to e on a time-reversal-symmetric full grid") | |
| Eq. (5) ρc = |e|²/(k_BT²), ballistic weight | derived here; consistent with `notes/17` §2 (C_μ) | checked (dimensions) |
| Eq. (6) drive at zero wave vector | derived here from the local-equilibrium streaming term | setting stated explicitly (audit §6): a uniform gradient maintained by an external heat source, equivalently the zero-wave-vector limit at fixed gradient amplitude of the plane-wave problem, whose energy component decouples. Checked in `checks/audit/audit_prop5_drive.py` (approach O(k²); the energy projection needs the source) |
| Eq. (7) memory law (Gurtin–Pipkin form) | derived here; GurtinPipkin1968 for the form | |
| Eq. (8) dynamic conductivity; 1/p warning | `notes/18` §2 (K_i(z), κ_i = K_i/(k_BT²), "un courant couplé à un invariant produit un terme en 1/z"); `notes/17` §4; `notes/conventions.md` §18 (exp(+iωt)) | |
| Eq. (9) Stieltjes/relaxon form | CepellottiMarzari2016; `fc3_memory/13-final-report.md` §5 ("τ_mem = N/K is a κ-weighted mean relaxon lifetime") | |
| Definition 1, Eq. (10) | `memory_gap/02-assumptions.md` §A ("Moments"); τ_stat = τ_CS of `memory_gap/13-final-report.md` §1–2 (P4) | τ_CS renamed τ_stat |
| Lemma 1, Eq. (11) | `memory_gap/02-assumptions.md` §A (N/K = −K′(0)/K(0)); `notes/18` §2; time-domain and relaxon forms and Γ bounds derived here | PROVED; checks "equivalent forms" |
| RTA forms of τ_mem and τ_stat | `notes/17` §3 ("Deux moyennes pour deux observables") | |
| Definition 2 (Cattaneo time = matched single pole) | `memory_gap/02-assumptions.md` §A ("relaxation time of the single (Cattaneo) pole … matches the response and its first derivative at z = 0"); `paper/PRIOR_ART.md` central contribution 1 | DEFINITION |
| Eq. (12), Proposition 1 (grey Callaway operator: j = v²P, Cj = j/τ_R, λ = ρc v²τ_R/3, τ_mem = τ_stat = τ_R, r_ν = τ_N⁻¹ + τ_R⁻¹ − O(1/n)) | derived here under the assumptions of `notes/14` §§1–2 (Callaway model, displaced equilibrium, compatibility conditions) | DERIVED UNDER ASSUMPTIONS (exact in the grey model); check-script block "Callaway" (12-mode star mesh, symbolic weights). Added with the audit repairs: the line on the spectrum of C_Cal on H (1/τ_R on span{P}, τ_N⁻¹ + τ_R⁻¹ on span{e,P}^⊥), which gives (K3) (audit §2.3; check "Callaway: spectrum"); the O(n⁻¹) argument with its limit, refinement at fixed T (audit §2.7; `checks/audit/audit_callaway.py`: n·max(dev) → 22.8 on O_h meshes up to n = 35,936) |
| Eq. (13) conserving GK closure, β = 1/3, ℓ² = v²τ_cτ_R/5; original GK β = 2 with τ_N | `notes/14` §1 and §5; prior-art trail item 1 (Guo & Wang 2022, Eq. 7); `notes/conventions.md` §7, §7.1 | background (C1) |
| Validity conditions |p|τ_N ≪ 1, vτ_N/L_v ≪ 1 | `notes/14` §1 | |
| "Lifetimes ≈ τ_c, not τ_R; (τ_c, λ) fixes (τ_N, τ_R)" | derived here from Prop. 1 | |
| Definition 3; r_ν > 0 under (K3) | `memory_gap/02-assumptions.md` §A ("Data"); `memory_gap/review/proof-audit.md` §5 step L.3 | |
| Remark that r = diag C is not a self-energy linewidth | `memory_gap/02-assumptions.md` §A; `memory_gap/13-final-report.md` §2 (L3) | |
| Definition 4, Eq. (14) | `memory_gap/02-assumptions.md` §A ("Operators", "admissible"); `notes/18` §4 | |
| Eq. (15) event form, linearization of the net rate | `memory_gap/02-assumptions.md` §A ("Events") and §D (rate formula); `notes/18` §11 ref 5 (Fugallo et al.); the linearization written here | |
| Declared energy-conserving stoichiometry for broadened events | `memory_gap/02-assumptions.md` §C | Supplemental Sec. S1 |
| Eq. (16) feasible set | `memory_gap/02-assumptions.md` §A ("Feasible set") | |
| Lumping criterion BL = GB | `notes/21` §3; `paper/PRIOR_ART.md` C3 | background |

## 4. Section III (Static data do not determine the memory time)

| Item | Source | Status / notes |
|---|---|---|
| Example 1, Eq. (17) | `notes/18` §3.1; `theory/aln/closure/10-proof-status.md`; `ResearchLab/orchestration/campaigns/20261002-paper1-prior-art/experiments/check_note18_examples.py` | PROVED; matrices written explicitly instead of L(a, β, c); check-script block "four-mode" |
| "Both operators lie in one fixed cone (vectors δ_μ − δ_ν)" | derived here | |
| Theorem 1, Eq. (18) (i)–(iii) | `notes/18` §3.2 (family A(t), t renamed η); `theory/aln/closure/CLAIM_LEDGER.md` C1 ("PROVED in class"); `closure/10-proof-status.md` | PROVED; check-script block "Thm 1" |
| Theorem 1 (iv), involution i ↦ i ⊕ 3 | `notes/18` §3.2 states that the construction respects the odd-current/even-energy parity; the explicit involution is derived here | checked |
| Theorem 1 (v) band bound | `notes/18` §3.2 (|K_t(iω)/(1/(iω+2)) − 1| ≤ t²/(2ω_min)) | same inequality, rewritten |
| Remark 1 (slow root −η²/2, residue η²/4, half of K₀, all of N) | derived here | checked ("Thm 1 remark") |
| Remark 2 (Stieltjes moments; spectral gap bound; elliptope; attribution) | `paper/PRIOR_ART.md` C2a; prior-art `13-final-report.md` §2.1 C2(a) and trail items 19–21 (Saunderson; Golden–Papanicolaou; parity; Malviya–Ravichandran) | |
| Remark 3, Eq. (19) (no fixed cone contains the family) | `notes/18` §4 (triangle inequality with ratio R; "t ≤ R t²"); `closure/errata.md` item 3 (case R = 0) | PROVED |
| Proposition 2, Eq. (20) | `memory_gap/13-final-report.md` §2 P4(a); `memory_gap/review/proof-audit.md` §6 (4a) and §10.2 (τ ≥ |x|/|b|) | PROVED |
| τ_stat = λ/Λ_i; RTA reduction τ_mem ≥ τ_stat | `notes/17` §3 (τ_stat, τ_mem ≥ τ_stat by Cauchy–Schwarz) | |
| PENDING on minima of τ_mem/τ_stat | `memory_gap/13-final-report.md` §1 (factor ≤ 1.69 in every configuration tested) | deferred to Sec. V |
| Lemma 2 (K₀ ≥ K_RTA/ς_max, sharp) | `memory_gap/13-final-report.md` §2 L3; `memory_gap/review/proof-audit.md` §5 (sufficient equality condition; 11-mode network gives 1/3, repeated-daughter chain gives 1/2) | PROVED; the variational-principle proof is an independent route written here; check "Lemma 2". Audit repairs (audit §4.6–4.7): "up to 4" for four-phonon events (events with a repeated mode touch 3 modes); equality criterion (equality iff C y = ς_max diag(r) y for some y ≠ 0; currents j = C y) with its consequence (every active event touches ς_max modes; strict for mixed networks), proof written here from audit §4.6; sharpness now PROVED for ς_max = 2, 3, 4 by the exact networks of Supplemental Sec. S2 (Sec. 12 below), which replace the eleven-mode example of the sources; "constant cannot be improved" for mixed networks from the family of Supplemental Sec. S2 |
| Self-energy linewidth counterexample | `memory_gap/13-final-report.md` §2 L3 (chain 1, 2, 4, 8, 16; ratio 0.3333 < 1/2); `memory_gap/review/counterexample-review.md` §6 | Supplemental Sec. S2 |
| Definition 5, Eq. (21) | `memory_gap/13-final-report.md` §2 T1; `memory_gap/review/proof-audit.md` §1 ("Basic solutions") | |
| Properties of M (rate independence, hyperplanes, scaling) | `memory_gap/review/proof-audit.md` §10.4 | |
| Theorem 2, Eqs. (22)–(24) | `memory_gap/13-final-report.md` §2 T1; `notes/18` §4.1; `memory_gap/review/proof-audit.md` §2 (steps 1.1–1.6 VALID; exact check `check_cone_bounds.py`) | PROVED; matrix-determinant-lemma derivation of adj(A) written here; checks "Cauchy-Binet". ψ_I is now defined for every I with \|I\| = n_H − 1, not only for bases, because Eq. (24) sums over all such subsets (audit §5.4; `checks/audit/audit_thm2_cone.py`: 39 of 126 non-basis subsets contribute) |
| Extension to non-admissible g with j ∈ range C(g) | `memory_gap/13-final-report.md` §2 T1 (last sentence); `memory_gap/review/proof-audit.md` §2.8 | |
| Remark 4 (weighted least squares; Ben-Tal–Teboulle Thm 2.1/Cor. 2.1; Berg; Stewart; O'Leary; Todd; Hanke–Neumann; Vavasis–Ye; Tunçel) | prior-art trail items 14–18; prior-art `experiments/check_bental_teboulle_mapping.py` (200 random cones); `memory_gap/review/proof-audit.md` §3.3 | inequality attributed; application "potentially novel within the literature searched" |
| Proposition 3, Eq. (25) | `memory_gap/13-final-report.md` §2 T2; `memory_gap/review/proof-audit.md` §3 (steps 2.1–2.6, §3.1 approached not attained, K ~ c/ε, §3.2) | PROVED; check "Prop. 3" (monotone approach to 1 for ε = 10⁻¹…10⁻¹¹). Proof wording repaired (audit §5.9): "some I* ∪ {γ}, γ ∉ I*, is a basis of H" (it exists because I* is independent and Z spans H), which gives det A ≥ ε det(Z_S)² |
| Remark 5 (lifetimes fixed: equality possible if j ∥ event vector; open for physical currents; symmetric rates) | `memory_gap/13-final-report.md` §2 ("Lifetime-constrained sharpness — counterexample and open question") and its erratum of 12 November 2026; `memory_gap/review/proof-audit.md` §3.4 (instance B: 7 modes, 10 events, r = (915/16, 7, 9/2, 29/9, 17/16, 2891/640, 153/80), K₀ = 4, τ = M²K₀ = 44367253/1074614) | OPEN for physical currents; Supplemental Sec. S3. Audit repairs: "physical currents have full support" (FALSE, audit §10.1: j vanishes where v_{ν,i} = 0; `checks/audit/audit_callaway.py` finds N² − 1 zero entries) replaced by the time-reversal statement (support = pairs {ν, ν̄} of equal-energy modes, never parallel to a three-phonon event vector, whose support has three modes or two modes of different energies); "transfers only if I* is a union of orbits" (unproved "only if", audit §10.3) replaced by the proved direction and "otherwise the supremum over symmetric rates is not known" |
| Proposition 4, Eq. (26) | `memory_gap/13-final-report.md` §2 P3; `memory_gap/review/proof-audit.md` §4 | PROVED; check "Prop. 4" (every basis of a test network) |
| Text after Prop. 4: large τ_mem/K₀ requires convex weight on bases with large \|y_I\| (for 0 < ε ≤ 1/2, τ_mem ≥ (1 − ε)²M²K₀ puts at least half of the weight on bases with \|y_I\| ≥ (1 − 2ε)M); the weight resides in g_I(ψ_Iᵀj)², not in individual rates; Prop. 3 shows that concentration attains the bound | derived here from Eq. (22) (audit §10.2, with θ = 1 − ε, θ′ = 1 − 2ε in its bound (θ − θ′)/(1 − θ′)) | replaces "By Proposition 3, large memory requires rates concentrated …", which attributed a necessity to a sufficiency result (audit §10.2, UNJUSTIFIED as attributed) |
| Corollary 1 | `memory_gap/13-final-report.md` §2 P6; `memory_gap/review/proof-audit.md` §8 | PROVED |
| Convexity, LMI Eq. (27), SDP | `memory_gap/13-final-report.md` §2 P5; `memory_gap/02-assumptions.md` §A (extended-value K); `memory_gap/review/proof-audit.md` §7 | PROVED with the extended value; "no convexity of τ_mem is claimed" from §2 P5 |
| Proposition 5, Eq. (28) and converse | `memory_gap/13-final-report.md` §2 P7; `memory_gap/review/counterexample-review.md` §3.3 | PROVED; check "Prop. 5"; F > 0 instead of F ≥ 1 (the proof needs only F > 0). Converse repaired (audit §7.5): the minimum of g_γ/g^ref_γ is taken over g^ref_γ > 0 (admissible g^ref may vanish on some events) and τ > 0 is required; the proof treats a zero minimum separately (`checks/audit/audit_prop5_drive.py`: certified, 12 cases with zero minimum) |
| "Information that confines τ_mem is microscopic" | `memory_gap/13-final-report.md` §5; `memory_gap/02-assumptions.md` §E | |
| Eq. (29) moments; Riemann-sum normalization | `memory_gap/13-final-report.md` §2 T7 (M_m = ∫ b² r⁻ᵐ); `notes/17` §2 (weights) | mesh matrix 𝖡 with \|det 𝖡\|V_cell = (2π)^d and the mesh moments 𝓜_m(N₁) defined in the preamble (audit §1.13 notation). The previous "\|j\|² → 𝓜₀, K_RTA → 𝓜₁" asserted a convergence that needs (H5) (audit §1.10); it now reads as the identities \|j\|² = 𝓜₀(N₁), K_RTA = 𝓜₁(N₁) on the meshes |
| (H1)–(H3) | `memory_gap/13-final-report.md` §2 T7(i) ("r = a(Ω) q^α with a bounded above and below and ∫ b0² a⁻ᵐ dΩ finite and positive"); `memory_gap/review/proof-audit.md` §10.1 (b₀² may vanish on a plane; r bounded below away from q = 0 on every branch); `theory/aln/closure/02-assumptions.md` ("uniform power-law value asymptotics … angular rate coefficient bounded away from zero") | notation: a(Ω) → r₀(k̂), b → J; the angular integral of (H2) is named A_m |
| (H5) (Φ_m = Σ_ac J₀²r₀⁻ᵐ bounded and continuous outside a closed null set of directions; J²r⁻ᵐ Riemann integrable on the part k ≥ δ of the zone for every δ > 0) | audit §1.10, §1.12 and §12 (wording of §12) | hypothesis added by the audit repair; it replaces the continuity of J₀ and r₀ of the previous draft, which is stronger than needed and fails at acoustic axes (sorted branches touching conically or crossing); some such regularity is necessary (counterexamples of audit §1.7 and §1.10, Supplemental Sec. S4) |
| J₀² = (k_BT v_s)² k̂_i² | `memory_gap/review/proof-audit.md` §10.1 | |
| (H4) | `memory_gap/13-final-report.md` §2 T7(ii) ("|x_μ| ≥ c|b_μ|/r_μ at small q, c uniform in the mesh") | |
| Theorem 3 (i) | `memory_gap/13-final-report.md` §2 T7(i) (criterion and mesh laws, log N at 2α = d, N^{2α−d}, N^d/log N at α = d); `memory_gap/review/proof-audit.md` §10.1 (exact lattice sums in d = 3 up to N = 161); `notes/18` §6 | PROVED under (H5) (audit §9, §15). Written here: the dyadic-shell argument; an additive (not multiplicative) error term, needed because J₀ vanishes on a plane (source audit §10.1). Audit repairs: (H5) replaces the continuity of J₀, r₀ (audit §1.10, §1.12); mesh laws stated with the sharp constants of audit §1.13 (𝓜_m(N₁) → 𝓜_m for mα < d; (2π)^{−d}A_m log N₁ + o(log N₁) for mα = d; N₁^{mα−d}S_m/V_cell·[1 + o(1)] for mα > d), proved by Lemma S1 (Supplemental Sec. S4, written from audit §1.6–1.9 and §1.13); the "∞/∞" wording for α ≥ d replaced by the ratio of the cut-off moments, ~k_min^{−α} (α > d) and ~k_min^{−d}/log(1/k_min) (α = d) (audit §1.5; check "cut-off moments"); constants confirmed by `checks/audit/audit_ir_lattice.py` (audit §1.15) |
| Theorem 3 (ii) | `memory_gap/13-final-report.md` §2 T7(ii); `memory_gap/review/proof-audit.md` §10.2 (repair τ ≥ √N/|b|; one-sided mesh exponents; two-sided UNJUSTIFIED) | DERIVED UNDER ASSUMPTIONS (H4) and the uniform form of (H1) (status note repaired, audit §1.14, §9). "(H1)–(H3) uniformly in N₁" replaced by its definition (r₀ independent of N₁, the o(1) uniform in k̂ and N₁; only the upper side is used); (H3) is not needed; the continuity of J₀, r₀ replaced by that of Φ₂ outside a closed null set; "\|j\|² → 𝓜₀" replaced by "\|j\|² ≤ (number of branches) sup J²/V_cell" (audit §1.10, §1.14) |
| Theorem 3 (iii) | `memory_gap/13-final-report.md` §2 T7(iii); `memory_gap/review/proof-audit.md` §10.4 | DERIVED UNDER ASSUMPTIONS (H4) and the uniform form of (H1) (it uses the hypotheses of (ii); label made consistent with (ii)). The clause "if K₀ stays bounded, M → ∞" is added here to make `paper/PRIOR_ART.md` ("so M diverges under refinement") precise |
| Slow approach to the mesh laws (𝓜₁(N₁) converges like N₁^{−(d−α)} for d/2 < α < d; relative correction ~1/log N₁ at α = d; exponents fitted on coarse meshes biased) | audit §1.15 (NUMERICAL OBSERVATION); `checks/audit/audit_ir_lattice.py`; Supplemental Sec. S4, Table S2 | stated as an observation of the lattice test, not as a theorem; consistent with the drift of the fitted exponent in `memory_gap/review/proof-audit.md` §10.1 (α = 3.5: 3.37 → 3.44) |
| Example 2 | `notes/18` §6 (K(z) = 1 − √z arctan(1/√z); cutoff 1/N gives τ = N) | PROVED; expansion 1 − (π/2)√p + O(p) written here |
| d = 3: α < 3 for κ, α < 3/2 for memory | `notes/18` §6 | |
| Herring (1954) | `fc3_memory/13-final-report.md` §5 | cited for the existence of branch- and symmetry-dependent power laws only |
| H4 holds with c = 1 in the RTA; is a property of an operator | `memory_gap/13-final-report.md` §2 T7(ii); `memory_gap/review/proof-audit.md` §10.2 | "(H1)–(H4)" in the material question and its PENDING marker (and in Sec. V) now reads "(H1)–(H5)" |
| Regularized memory time motivation | `fc3_memory/13-final-report.md` §4 items 1–2 | |

## 5. Section IV (Front-face measurement model)

| Item | Source | Status / notes |
|---|---|---|
| (F1)–(F6) | `theory/aln/fdtr_gauge/02-assumptions.md` B1–B3; `notes/07` §1; `notes/conventions.md` §7 | |
| Eq. (30) 1D GK, L² = (1+β)ℓ²; β = 2 historical, 1/3 conserving | `notes/conventions.md` §7.1 (α renamed β); `notes/14` §6 | |
| Eq. (31) τ_ℓ | `notes/conventions.md` §7.1; `notes/07` §2.3 | |
| Interface conditions in 3D models (slip, temperature jumps) | `notes/15` §2 (Hennessy–Myers: tangential slip), §3 (Beardo: jumps and slip) | |
| No extra boundary condition in 1D (second-order equation for q) | derived here | checked (GK 1D block) |
| Eq. (32) λ_eff | `notes/conventions.md` §7; `notes/07` §2.2; `notes/14` §6 | checked |
| Eq. (33) σe and λ_eff σ; times exchange roles | `notes/conventions.md` §7; `notes/07` §2.3 | checked |
| Eq. (34) quadrupole, flux positive toward z > 0, AD − BC = 1 | `notes/conventions.md` §§1–4 (Maillet, S = 1); `notes/14` §6 | checked (rear-to-front map, det = 1) |
| Entries even in σ (branch immaterial in the film) | derived here | |
| Z_s = 1/(b_s√p) | `notes/conventions.md` §5 | |
| Eq. (35) Z = (AZ_s + B)/(CZ_s + D), quotient form | `notes/conventions.md` §5 (h = 0); `src/forward_model.py`, function `response` (z_in = (z + t/coeff)/(1 + coeff·z·t)) | checked |
| Eq. (36) effusivity ratio ϱ(p) | derived here | |
| f_c = 1/(2πξ₁²) = a/(2πe²) | `notes/conventions.md` §12.1; `notes/07` §7 | |
| Proposition 6, Eqs. (37)–(38) | `notes/07` §3 (gauge group and Euler identity with the ℓ term); `notes/conventions.md` §8; `theory/aln/fdtr_gauge/00-question.md` B3; prior-art `13-final-report.md` §3 ("FDTR gauge for Cattaneo/GK: DERIVED UNDER ASSUMPTIONS") | DERIVED UNDER ASSUMPTIONS (F1)–(F6), exact in the model; checked (gauge and Euler identity) |
| Krapez–Rigollet (Fourier gauge, unconditional) | `notes/07` §1; prior-art trail item 33 | |
| "Extension not found for non-Fourier laws" | prior-art trail items 34–37 | |
| Gauge breaking by a known length, or v with a = v²τ_R/3 | `notes/07` §3 (conclusion); `notes/conventions.md` §8 | the lateral-length mechanism is left to the pending analysis (`fdtr_gauge/00-question.md` B3) |
| Corollary 2 (β, ℓ invisible separately; ℓ ratio 3/2) | `notes/14` §6 (ℓ_CE = (3/2)ℓ_hist) and §7 ("une mesure 1D dépend de L² seulement"); `paper/PRIOR_ART.md` central contribution 3 | DERIVED UNDER ASSUMPTIONS |
| Entropy density, extended entropy flux, production (q² + L²q_z²)/(λT₀²) ≥ 0 | `notes/14` §7; `notes/11_Note_partie_II.tex` §4; `notes/conventions.md` §10; `tests/test_admissibility.py` | JouCasasLebon2010 for the framework |
| Laplacian and grad-div coefficients independent in nonequilibrium thermodynamics | prior-art trail item 11 (Kovács, Phys. Rep., read by text search) | |
| Fourier resonance and its restriction to the boundary-value problem | `notes/07` §6; `notes/conventions.md` §9 | attributed (Both2016, Kovacs2018) |
| Eq. (39); historical resonance at 4/5 (τ_c kept) and 9/5 (τ_c → τ_N); "formal marker" | `notes/14` §7; `notes/conventions.md` §7.1 and §9; `paper/CLAIM_MAP.md` ("Audit admissions, 2 October 2026"); `CORRECTIONS.md` (2 Oct 2026) | checked (closure-dependent resonance block) |
| Classical Fourier coating response for τ_R = τ_ℓ = 0 | `notes/conventions.md` §§2, 5 | |
| Eq. (40) half-space phase; Cattaneo −45° → 0° | `notes/conventions.md` §7 (Cattaneo half-space phase); `notes/15` §2.1 (Hennessy–Myers 1D limit, φ = −π/4 + [arctan ωτ_R − arctan(3ωℓ²/a)]/2) | checked (8 + 9 numerical points) |
| Naive substitution p → p(1 + τ_R p): wrong flux coefficient by (1 + τ_R p); phase toward −90° | `notes/conventions.md` §7; `notes/07` §2.4; `CORRECTIONS.md` (2 Oct 2026, note 00) | checked |
| Finite speed (a/τ_R)^{1/2}; GK parabolic, unbounded speed | `notes/11` §3; `notes/conventions.md` §10 | |
| Structural identifiability (Raue2009) | `fc3_memory/13-final-report.md` §5 | |
| Rank theorem | PENDING; `theory/aln/fdtr_gauge/00-question.md` workstream B (claims B1–B3) | not stated |

## 6. Sections V–VII and back matter (placeholders)

| Placeholder | Where the pending content will come from |
|---|---|
| Sec. V (AlN) items (i)–(iii) | `theory/aln/fc3_memory/13-final-report.md` (S6–S8; film-regularized ranges E33/E37 listed in `fc3_memory/README.md` and `fdtr_gauge/00-question.md` A); `theory/aln/memory_gap/13-final-report.md` §1–§3 (Debye and AlN 5×5×3 witnesses, M values) |
| Sec. VI (FDTR) | `theory/aln/fdtr_gauge/00-question.md` workstream A; `fdtr_gauge/02-assumptions.md` A1–A7; `src/twin_fdtr.py` |
| Sec. VII, abstract tail, acknowledgments, data availability, affiliation | to be written |

Numbers inside the placeholders (500 nm; 10 kHz–200 MHz; 0.01°–0.1°; 1 %; 30, 100, 300 nm; nuisance sets) are copied from
`fdtr_gauge/00-question.md` (A) and `fdtr_gauge/02-assumptions.md` (A4–A6).

## 7. Every number in the established sections

| Number | Where | Source |
|---|---|---|
| 321 W m⁻¹K⁻¹; 18 and 22.5 µm; room temperature; TDTR | Sec. I | Cheng2020 (local PDF, p. 3); `notes/conventions.md` §11.1 |
| about fifty queries per database; October 2026 | Sec. I | prior-art `13-final-report.md` §2; `08-literature.md` §1 |
| 1/3, 2 (closure coefficient); 1/5 (ℓ²) | Secs. I, II, IV | `notes/14` §§1, 5, 6; Guo & Wang Eq. (7) via trail item 1 |
| ρc v²τ_R/3 | Prop. 1 | derived here; `notes/14` §1 (λ = Cv²τ_R/3) |
| O(n⁻¹) diagonal correction | Prop. 1 | derived here; argument (c_ν ≤ k_B, w_ν²k_{ν,i}² bounded, Σc_μ = 𝒱ρc and Σw_μ²k_{μ,i}² growing like n, refinement at fixed T) from audit §2.7 |
| rates τ_R⁻¹ (on span{P_x, P_y, P_z}) and τ_N⁻¹ + τ_R⁻¹ of C_Cal on H | Prop. 1 proof | audit §2.3; check "Callaway: spectrum" |
| 4/3, 2/(p + 4/3), 1/(p+1) + 1/(p+2), 3/2, 3/4, 5/6 | Example 1 | `notes/18` §3.1 |
| n = 8; (0, η², 2, 2−η², 1, 1, 1, 1); ±η entries; 1 + η⁻²; (2, ∞); η²/(2ω_min) | Theorem 1 | `notes/18` §3.2 |
| −η²/2, η²/4, "about half of K₀" | Remark 1 | derived here (checked) |
| ς_max = 3, 2, up to 4 | Lemma 2 | `memory_gap/13-final-report.md` §2 L3; "up to 4" from audit §4.7 |
| ratio 1/ς_max attained for ς_max = 2, 3, 4 | Lemma 2 and proof | Supplemental Sec. S2, Table S1 (audit §4.6); replaces "1/3 (eleven modes) and 1/2 (chain)" of `memory_gap/review/proof-audit.md` §5.3 |
| ε ≤ 1/2, (1 − ε)²M²K₀, at least half of the weight, \|y_I\| ≥ (1 − 2ε)M | after Prop. 4 | derived here from Eq. (22) (audit §10.2) |
| seven modes, ten events | Remark 5 | `memory_gap/review/proof-audit.md` §3.4 |
| m = 0, 1, 2; α < d; 2α < d; log N₁; N₁^{2α−d}; N₁^d/log N₁; N₁^{α−d/2}; (log N₁)^{1/2} | Theorem 3 | `memory_gap/13-final-report.md` §2 T7; `memory_gap/review/proof-audit.md` §10 |
| \|det 𝖡\|V_cell = (2π)^d; (2π)^{−d}A_m log N₁; N₁^{mα−d}S_m/V_cell; S_m ∈ (0, ∞) | Theorem 3(i) | audit §1.13 (sharp constants); Lemma S1, Supplemental Sec. S4 |
| k_min^{−α} (α > d); k_min^{−d}/log(1/k_min) (α = d) | Theorem 3(i) | audit §1.5; check "cut-off moments" |
| N₁^{−(d−α)}; 1/log N₁ | after Example 2 | audit §1.15 (lattice test, Supplemental Table S2) |
| factors 1/2 and 3/2 in the proof of Theorem 3 | proof | written here |
| α < 3, α < 3/2 (d = 3) | after Example 2 | `notes/18` §6 |
| 1 − (π/2)√p + O(p); K₀ = 1 − 1/N₁; τ_mem = N₁ | Example 2 | `notes/18` §6 (expansion written here) |
| 3/2 (ℓ ratio) | Corollary 2 | `notes/14` §6 |
| 4/5, 9/5, 4τ_c/5, 9τ_c/5 | Sec. IV E | `notes/14` §7; `notes/conventions.md` §§7.1, 9 |
| −π/4; −45°, 0°, −90° | Sec. IV F | `notes/conventions.md` §7; `notes/07` §2.4 |

## 8. Equation index (printed number → label → source row above)

(1) eq:entropycoords · (2) eq:kinetic · (3) eq:entropy · (4) eq:current · (5) eq:harmonic · (6) eq:drive ·
(7) eq:memorylaw · (8) eq:dynamicK · (9) eq:stieltjes · (10) eq:moments · (11) eq:forms · (12) eq:callaway ·
(13) eq:gk · (14) eq:cone · (15) eq:events · (16) eq:feasible · (17) eq:four · (18) eq:walsh ·
(19) eq:obstruction · (20) eq:lower · (21) eq:basic · (22) eq:hull · (23) eq:conebound · (24) eq:adj ·
(25) eq:sharp · (26) eq:subnetwork · (27) eq:lmi · (28) eq:prior · (29) eq:Mm · (30) eq:gk1d · (31) eq:taul ·
(32) eq:lameff · (33) eq:combinations · (34) eq:quadrupole · (35) eq:impedance · (36) eq:rho ·
(37) eq:gaugegroup · (38) eq:euler · (39) eq:noresonance · (40) eq:halfphase.

Supplemental Material: (S1) eq:Seq · (S2) eq:mixed · (S3) eq:SN · (S4) eq:Sa · (S5) eq:shellbound ·
(S6) eq:shellsum · (S7) event-operator lower bound (no label); Lemma S1 lem:lattice; Tables S1 tab:sharp,
S2 tab:ir.

## 9. Bibliography: verification record

All 62 entries of `refs.bib` (60 cited) carry a DOI whose registry record was retrieved on 12 November 2026
(`checks/doi_records.txt`; Crossref, DataCite for Chen2025). Titles, authors, journal, volume, pages and year
were copied from those records; journal names abbreviated; accents and capitalization protected by hand;
`checks/compare_bib.py` finds no discrepancy in year, volume or first page. Remarks:

- Three DOIs first guessed from memory resolved to *unrelated* works and were replaced by bibliographic search
  (`checks/search_titles.py`): Krapez 2019 (correct 10.1016/j.ijthermalsci.2018.10.018), Krapez 2023 comment
  (correct 10.1063/5.0151465; 10.1063/5.0155055 is the Mandelis response), Kovács Phys. Rep. (correct
  10.1016/j.physrep.2023.11.001).
- Maillet et al., *Thermal Quadrupoles* (Wiley, 2000): no registered DOI found; cited as PENDING in Sec. IV C.
- Cattaneo (1948): no DOI; not cited; the Cattaneo law is cited through Joseph & Preziosi (1989).
- JouCasasLebon2010: the Crossref record of the 4th edition lists two authors (Jou, Lebon); the chapter records
  of the same book and of earlier editions list Jou, Casas-Vázquez, Lebon, which is used.
- HennessyMyers2021: online 2020-12-16, print 2021 (year 2021 used); Crossref gives no volume or editors.
- Kovacs2024: issued 2024 (vol. 1048); `08-literature.md` item 11 cites it as 2023.
- Entries present but not cited in the current text (printed only when cited): Krapez2023, LebonDauby1990.
- Content claims rest on the reading recorded in the prior-art trail (`08-literature.md`), the notes and the
  campaign reports; several are abstract-level (Hardy & Albers, Struchtrup, Zeng & Dong, Malviya & Ravichandran,
  Vermeersch et al.). The manuscript cites them only for what the abstract states.

## 10. Notation choices where the sources disagree

| Quantity | Manuscript | Sources | Reason |
|---|---|---|---|
| Laplace variable | p (p = iω, e^{+iωt}) | z in `memory_gap`, `notes/18` | `notes/conventions.md` (binding); z is the depth coordinate |
| Mode index | ν = (k, s) | μ in `memory_gap`, `notes/17` | μ is the gauge parameter (task and `notes/07`) |
| Wave vector | **k** | q in `memory_gap`, `notes/18` | q is the heat flux (`notes/14`, conventions) |
| Bose occupation | f⁰_ν | N_μ in `memory_gap`, n⁰ in `notes/21` | N is the second memory moment |
| Bose weight matrix | 𝖶 = diag(w_ν) | D = diag(d_μ) in `memory_gap` | d is the spatial dimension |
| Entropy coordinates | y = δf/(𝒱^{1/2} w) | y = D⁻¹ dn (`memory_gap`); √W δn/√(n(1+n)) (`notes/17`) | volume factor makes sums Riemann sums, λ = K/(k_BT²) without 𝒱 |
| Energy invariant | **e** (bold) | e (`memory_gap`, `notes/18`) | e (italic) is the film thickness (conventions) |
| Current vector | **j** | b (`memory_gap`, `notes/18`) | b is the effusivity (conventions) |
| Collision operator | 𝒞 (calligraphic) | C(g), 𝒞 | C is a quadrupole entry |
| Event index, vectors | γ, ζ_γ | α, a_α (`memory_gap`); z_α (`notes/18` §4) | α is the infrared exponent; a is the diffusivity; z is depth |
| Event rates | g_γ | g_α | — |
| Infrared exponent | α | α (`memory_gap` T7, `notes/18` §6); n (`fc3_memory`) | task wording "2α ≥ d" |
| Closure coefficient | β (β = 2 historical, 1/3 conserving) | α (`notes/14` §6, conventions §7.1) | α is the infrared exponent |
| Spatial dimension | d | d (`memory_gap` T7); n (`notes/18` §6) | — |
| dim of dissipative space | n_𝓗 = n − 1 | d (`memory_gap` 02 §A, `notes/18` §4.1) | d is the spatial dimension |
| Dissipative space | 𝓗 | H | 𝖧 is the Walsh matrix |
| Lower bound K₀/|j|² | τ_stat ("static time") | τ_CS (`memory_gap`, `fc3_memory`); τ_stat (`notes/17`, RTA) | descriptive name; coincides with `notes/17` in the RTA |
| Theorem 1 parameter | η | t (`notes/18`) | t is time |
| Theorem 1 block matrix | 𝖳(η) | A(t) | A is a quadrupole entry |
| Cofactor vectors, rate products, weights | ψ_I, g_I, ϖ_I | v_I, c_I, θ_I (`notes/18`) | v velocities, c capacities, θ temperature |
| Hidden invariant | φ_I (bold) | φ_I | φ (italic) is a phase |
| Max distinct modes per event | ς_max | s_max | s_γ is the stoichiometric vector |
| Relaxon rates | Γ_k | λ (eigenvalues) | λ is the conductivity |
| Smallest eigenvalue on 𝓗 | Γ_min | λ_min (`memory_gap` P7) | idem |
| Angular rate coefficient | r₀(k̂) | a(Ω) (`memory_gap` T7, `notes/18` §6) | a is the diffusivity |
| RTA moments | 𝓜_m | M_m (`memory_gap`, `notes/18`) | M is the cone constant |
| Mesh moments | 𝓜_m(N₁) | 𝓜_m(N₁) (audit §1.13) | — |
| Angular function of (H5) | Φ_m = Σ_ac J₀²r₀⁻ᵐ, A_m = ∫Φ_m dΩ | Φ_m, A_m (audit §1.12–1.13) | distinguished from the autocorrelation Φ(t), Φ_{ik}(t) by the order index m and the argument k̂ |
| Reciprocal-lattice matrix | 𝖡 | 𝖡 (draft proof, audit) | the same letter is used for the block-observable matrix of the lumping remark in Sec. II (overload present in the previous draft) |
| Lattice sum | S_m (main text), S_Φ (Lemma S1) | S_φ (audit §1.13) | φ is reserved for the lattice function φ(n) |
| Error bound in Lemma S1 | ϑ(t) | δ(t) (audit §1.13) | δ is the radius in (H5) and in the proofs |
| Concentration tolerance after Prop. 4 | ε, with θ = 1 − ε, θ′ = 1 − 2ε | θ, θ′ (audit §10.2) | θ is the temperature rise and the regularization parameter of the extension paragraph |
| Mesh points per direction | N₁ | N (`memory_gap`) | N is the second moment |
| Memory kernel (time domain) | κ(t) = Φ(t)/(k_BT²) (Gurtin–Pipkin form) | "K(t)" in the brief | K(0) would be ambiguous (K(t=0) = |j|², K(p=0) = K₀) |
| Conductivity | λ (film and bulk); λ(p) dynamic | κ in kinetic notes; λ in conventions | conventions (binding) for the film; one symbol throughout |
| Prescribed data | r*, K* | r, K0 | distinguishes data from functions |
| Film parameter vector | χ = (ξ₁, b, τ_R, τ_ℓ) | θ (`fdtr_gauge` B1) | θ is the temperature rise |
| Spelling | American (PRB) | brief uses British ("linearised", "regularised") | journal style |

## 11. Warnings that remain in the build

- `Package nameref Warning: The definition of \label has changed!` — produced by REVTeX 4.2f's own `\label`
  (2022) with hyperref 7.01r/nameref 2.58 and the LaTeX 2026-06-01 kernel; it occurs for a bare
  `revtex4-2` + `hyperref` document. Cross-references are correct. Not removable without patching the class.
- `Warning--jnrlst (dependency: not reversed) set 1` from bibtex: informational control message of `apsrev4-2.bst`.
- 20 underfull boxes: two with badness 10000 inside the temporary PENDING item of the contribution list; the
  rest have badness ≤ 4266 (two-column justification and long bibliography entries). Pages inspected visually,
  including every page changed by the audit repairs (pp. 1–4, 6–10, 12, 14) and the four pages of the
  Supplemental Material.
- Supplemental Material: no warnings; script paths carry `\allowbreak` after each `/` (REVTeX conflicts with the
  `array` package, so table columns are plain `p{}` columns).

## 12. Supplemental Material (sections written so far)

| Item | Source | Status / notes |
|---|---|---|
| Sec. S2, equality networks: currents attaining K₀ = K_RTA/ς_max are j = C y with C y = ς_max diag(r) y; sufficiency (S1) when every event touches ς_max modes and ζ_{γν}y_ν is constant on each support | audit §4.6 (equality criterion); `memory_gap/review/proof-audit.md` §5.3 (sufficient condition, sign structure) | PROVED (written here); with y = 𝖶ψ the condition involves the stoichiometry only, so it holds for all positive weights and rates |
| Table S1: ς_max = 2 (5 modes, energies 1, 2, 4, 8, 16; events 2→1+1, …, 5→4+4; ψ = (1, −2, 4, −8, 16)); ς_max = 3 (8 modes, energies (6, 6, 6, 3, 3, 3, 3, 3); 7 events 6 → 3 + 3); ς_max = 4 (9 modes, energies (6, 6, 6, 2, 2, 2, 2, 2, 2); 8 events 6 → 2 + 2 + 2, 1 ↔ 3 processes); ratios 1/2, 1/3, 1/4 | audit §4.6 table; the ς = 3 and ς = 4 networks were found by `checks/audit/audit_lemma2.py` (event lists printed in `audit_lemma2.out.txt`); the chain is that of `memory_gap/13-final-report.md` §2 L3 | PROVED: energy conservation, rank n − 1, C W ψ = ς_max diag(C) W ψ for symbolic weights and rates, and K₀/K_RTA = 1/ς_max exactly for rational weights (check block "Lemma 2 sharpness": K₀ = 1504, 3707/56, 8926/63; K_RTA = 3008, 11121/56, 35704/63) |
| 2 ↔ 2 networks alone conserve the phonon number (not admissible unless all energies are equal) | audit §4.6 ("pure 2 ↔ 2 networks are never admissible"; `audit_lemma2.py`, last check) | the qualification "unless all modes have the same energy" is added here: 𝖶1 ∈ ker C, and 𝖶1 ∥ e only when ε ∝ 1 |
| Networks with events of different sizes: min_j K₀/K_RTA = 1/λ_max(G), G = diag(r)^{−1/2} C diag(r)^{−1/2}; six-mode family (energies 1, …, 6, unit weights; 2+5→3+4 with rate 1; 3→1+2, 4→1+3, 5→1+4, 5→2+3, 6→1+5, 6→2+4 with rate δ); y = (0, 1, −1, −1, 1, 0): yᵀCy = 16 + 10δ, yᵀdiag(r)y = 4 + 12δ; bound (S2); minima 1.2198, 1.0234, 1.00236, 1.000024 for δ = 10⁻¹, 10⁻², 10⁻³, 10⁻⁵ | audit §4.6 (family and minima; `audit_lemma2.py`); the identity min K₀/K_RTA = 1/λ_max(G) and the exact quotient are written here | PROVED (strictness by Lemma 2, bound (S2) exact); minima numerical (check "Lemma 2, mixed network", 50-digit eigensolver: 1.2198291, 1.0234091, 1.0023557, 1.0000236) |
| Sec. S2, self-energy linewidths | `memory_gap/13-final-report.md` §2 L3 | PENDING (not part of the audit repairs) |
| Sec. S4, setting (𝖡, \|det 𝖡\|V_cell = (2π)^d, 𝒱 = N₁^dV_cell, k₀ small) | audit §1 preamble and §1.13 | — |
| Lemma S1 (lattice sums), Eqs. (S3)–(S6), parts (a)–(c) | audit §1.13 (statement), §1.6–1.9 (dyadic shells, uniform shell bound, Stolz–Cesàro, dominated convergence, sharp constants) | PROVED (written here from the audit). Choices made in writing it out: shell bound with C₁ = (6/ρ₀)^d from a packing argument; Riemann sums on 2^{−l}ℤ^d ∩ D by Lebesgue's criterion; (a) proved by a split at radius δ instead of shell by shell; positivity of S_Φ from continuity of Φ at a point where Φ > 0 and density of the lattice directions |
| "Regularity is necessary" (two constructions) | audit §1.7 (Φ = 0 on lattice directions) and §1.10 (J set to 0 at mesh points with k ≥ k₀) | counterexamples, verified by inspection |
| Mesh laws and the four τ^RTA laws with constants | audit §1.13 | PROVED under (H1)–(H3), (H5) |
| Event-operator lower bound (S7) | audit §1.14 (only the upper side of the uniform (H1) enters) | derivation written here; part of Theorem 3(ii), DERIVED UNDER ASSUMPTIONS |
| Numerical check, Table S2 (d = 2, 𝖡 = 2πI, N₁ ≤ 3201, J² = (k̂_x² + 0.3k)e^{−k²/8}, r = (1 + 2.4k̂_x²k̂_y²)k^α(1 + 0.5k²); values 0.6291354, 2.6·10⁻³, N₁^{−0.79}, 0.3885388, 0.99657, 2.5·10⁻³, N₁^{−0.500}, 0.99785, 1.00087, 1.00158, 1.057, ratio 1.026 at α = 1.5) | audit §1.15; `checks/audit/audit_ir_lattice.py`, `.out.txt`, `.json` | numerical observation; the ratio 1.026 is `ratio/pred = 1.025713` of the output |
