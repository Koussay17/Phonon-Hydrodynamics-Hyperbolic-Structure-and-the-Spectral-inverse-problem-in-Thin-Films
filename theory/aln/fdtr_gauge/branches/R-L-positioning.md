# R-L — Prior-art verification and positioning (referee objections B-1 and A-3)

Campaign `20261111-fdtr-gauge-memory`, revision branch R-L, 12 November 2026. Answers objection B-1 (prior art and
novelty of workstream B) and the positioning part of objection A-3 (workstream A) in `12-peer-review.md`. It builds on
`08-literature.md` and `01-known-results.md` and does not repeat their searches. New verification records and searches
are in `searches/R-L/` (`verify/` for API records, `runs/` for the multi-index searches). BibTeX entries with verified
DOIs: `refs_revision.bib` (campaign root). The manuscript bibliography `paper/manuscript/refs.bib` was not modified.

**Labels.** SOURCE SAYS (text read; short quotation or paraphrase) · SOURCE (ABSTRACT) SAYS · SECONDARY (another
paper's statement about the source, not checked in the source) · OUR INTERPRETATION (derivation or reading done here).
Access level is given for every item: *full text*, *section read*, *abstract*, *metadata*.

---

## 0. Verdicts at a glance

| # | Reference (verified DOI) | Access | What it states (as verified) | What it anticipates |
|---|---|---|---|---|
| 1 | Tzou, *J. Heat Transfer* **117**, 8–16 (1995), doi:10.1115/1.2822329 | abstract | A universal lagging (DPL) constitutive law covering diffusion, waves, phonon–electron interaction and phonon scattering; identifies the parameters governing transitions between mechanisms | The DPL law itself. The equal-lag reduction is **not in the abstract**; its location in this paper is SECONDARY (unverified) |
| 2 | Zhang, Tzou, Chen, Nova chapter (2009), arXiv:1511.03566 | full text | The first-order DPL energy equation "reduces to the parabolic conduction equation … if τ_q = τ_T" | Δ ⊂ Σ (Fourier equivalence on Δ), stated in a Tzou-coauthored text |
| 3 | Tang & Araki, *Int. J. Heat Mass Transf.* **42**, 855–860 (1999), doi:10.1016/S0017-9310(98)00244-0 | metadata; citing-paper descriptions | Title and citing papers: wavy, wavelike and diffusive DPL responses of finite slabs under short-pulse heating | The "diffusion mode" label for τ_T = τ_q is attributed by a review to the **2000** paper (row 4), not to this one; primary texts unread |
| 4 | Tang & Araki, *Mater. Sci. Eng. A* **292**, 173–178 (2000), doi:10.1016/S0921-5093(00)01000-5 | metadata; SECONDARY | Four DPL modes: wave (τ_T = 0), wavelike (0 < τ_T < τ_q), diffusion (τ_T = τ_q), over-diffusion (τ_T > τ_q) (Zhmakin review, p. 16) | Δ as "diffusion mode" (naming) |
| 5 | Kitamura & Nakagiri, *SIAM J. Control Optim.* **15**, 785–802 (1977), doi:10.1137/0315050; comment Courdesses, **21**, 410–412 (1983), doi:10.1137/0321023 | abstracts | Deterministic identifiability of spatially varying and constant coefficients of a linear 1D parabolic PDE, distributed or pointwise data, boundary-condition parameters known; the comment finds the solution expression and some results incorrect and gives a necessary and sufficient condition | Problem class and method family (known domain). The technique (Laplace-domain analyticity) is not visible in the abstracts. No specific B statement |
| 6 | Pierce, *SIAM J. Control Optim.* **17**, 494–499 (1979), doi:10.1137/0317035 | abstract | Point observation of a 1D parabolic problem: eigenvalues determined, coefficients not; normal-form coefficient unique in some cases (Gel'fand–Levitan) | Spectral uniqueness technique; antecedent of the graded Fourier gauge (B-K3). Nothing on non-Fourier laws or unknown thickness |
| 7 | Lee & Cahill, *J. Appl. Phys.* **81**, 2590–2595 (1997), doi:10.1063/1.363923; companion *Microscale Thermophys. Eng.* **1**, 47–52 (1997), doi:10.1080/108939597200421 | abstracts | Thickness dependence of the apparent conductivity of < 50 nm films attributed to an interface resistance ≈ 2×10⁻⁸ m²K/W, equivalent to ≈ 20 nm of SiO₂; heater/film and film/substrate interface conductances act in series | The film/interface resistance confound in the Fourier thin-film limit and its experimental remedy (thickness series); see §2.4 for the exact relation to the Δ ∩ {R = 0} stratum |
| 8 | Zhmakin, "The Zoo of Non-Fourier Heat Conduction Models", arXiv:2212.12922 (2022) | section read (DPL, p. 16) | Tang–Araki four modes (row 4); exact-delay DPL depends only on τ_T − τ_q | Secondary pointer only. Caution: the exact-delay DPL is a different model from the first-order DPL that equals 1D GK |
| 9 | Kovács, arXiv:2302.13931 = *Phys. Rep.* **1048**, 1–75 (2024), doi:10.1016/j.physrep.2023.11.001 (already `Kovacs2024` in refs.bib) | sections read (§4.3, §6) | DPL "reproduces Fourier's solutions for any τ_q = τ_T", likened to the Fourier resonance; the GK structure is recovered in the T-representation of DPL; GK Fourier resonance named | The most direct verified statement that the DPL equal-lag case and the GK Fourier resonance are one phenomenon: anticipates Δ ⊂ Σ and B-K6 at equation level |
| 10 | Ordóñez-Miranda & Alvarado-Gil, *J. Heat Transfer* **132**, 061302 (2010), doi:10.1115/1.4000748 | abstract, full reference list, companions | Finite DPL layer on a semi-infinite layer, modulated surface flux; at large effusivity contrast, explicit formulas for both time parameters and further layer properties | Partial prior art for Theorem 1 (constructive, regime-limited). Thickness assumption and any uniqueness statement **not established** (§3) |

Workstream A (§5): Beardo et al. 2020 and 2021, the PRISM preprint (arXiv:2611.02475) and Camacho de la Rosa et al. 2025
were read in full; their claims about relaxation times are tabulated with numbers in §5.

---

## 1. Access, routes and records

| Route | Outcome |
|---|---|
| Crossref API | All DOIs in this file resolved and matched (authors, title, venue, volume, pages, date). Records `searches/R-L/verify/*.json`. The full 36-item reference list of O-M & A-G 2010 was obtained from Crossref |
| OpenAlex API | Abstracts (Tzou JHT 1995, Tzou JTHT 1995, Kitamura–Nakagiri, Courdesses, Pierce, Lee–Cahill ×2, O-M & A-G 2010, Vales-Pinzón 2012, Gentili 1991); forward citations of O-M & A-G 2010 (25 records; 2 open access, neither reproduces the method) |
| Semantic Scholar API | Records (abstracts elided by ASME/Elsevier); citation contexts for O-M & A-G 2010 (27 citing), O-M & A-G 2009 (78), Tang–Araki 1999 (143), Tang–Araki 2000 (94), Basirat Tabrizi 2009 (12) |
| Unpaywall | No open copy for any paywalled item checked (13 DOIs) |
| arXiv API and PDFs | 2212.12922, 2302.13931, 2611.02475, 2007.05487, 1511.03566, 2406.13852, 1908.07612, 2609.13536 |
| Europe PMC | Full text of O-M & A-G, *Nanoscale Res. Lett.* **6**, 327 (2011) (open access) |
| HAL API | Ordóñez-Miranda's deposits start in 2013; nothing from 2009–2012 |
| CINVESTAV institutional repository | Ordóñez-Miranda's 2011 doctoral thesis not deposited; theses of later students of Alvarado-Gil found are on other topics |
| Springer abstract page | O-M & A-G, *Granular Matter* **12**, 569–577 (2010) |
| Local full texts | Beardo et al. 2020 (`these/.build/research/beardo_2020.*`), Camacho de la Rosa et al. 2025 (`these/notes/155103_1_5.0257299.pdf`, CC BY) |
| Multi-index search (Crossref + OpenAlex; Semantic Scholar rate-limited) | 5 new queries RL01–RL05 (§8) |
| Web search | 12 targeted queries (titles, quotations, theses, PRISM) |
| Not accessible | ASME (Tzou 1995; O-M & A-G 2010; HTTP 403 to automated retrieval), SIAM (Kitamura–Nakagiri, Pierce, Courdesses; HTTP 403), Wiley (Tzou monograph; HTTP 403), ResearchGate (HTTP 403), Elsevier full texts (Tang–Araki 1999/2000, O-M & A-G 2009, Basirat Tabrizi 2009, Shen–Zhang 2008), AIP full texts (Lee–Cahill 1997, Vales-Pinzón 2012), Google Books API (daily quota exhausted). No bot-protection was circumvented. |

---

## 2. Task 1 — the classical references cited by the referee

### 2.1 Tzou (1995 papers, monograph)

- **Existence and metadata (verified).** D. Y. Tzou, "A unified field approach for heat conduction from macro- to
  micro-scales", *J. Heat Transfer* **117**(1), 8–16 (1995), doi:10.1115/1.2822329. Companions: "The generalized lagging
  response in small-scale and high-rate heating", *Int. J. Heat Mass Transf.* **38**(17), 3231–3240 (1995),
  doi:10.1016/0017-9310(95)00052-B; "Experimental support for the lagging behavior in heat propagation", *J. Thermophys.
  Heat Transf.* **9**(4), 686–693 (1995), doi:10.2514/3.725. Monograph: *Macro- to Microscale Heat Transfer: The Lagging
  Behavior* (Wiley, 2014; listed as the 2nd edition in the Zhmakin review; first edition Taylor & Francis), doi:10.1002/9781118818275.
- **SOURCE (ABSTRACT) SAYS (JHT 1995).** A universal constitutive equation between heat flux and temperature gradient is
  proposed to cover diffusion, waves, phonon–electron interaction and pure phonon scattering, generalised from the
  dual-phase-lag concept; the energy equation identifies the parameters governing transitions between mechanisms.
  JTHT 1995 (abstract): the lagging model is re-examined against femtosecond gold-film data (Qiu et al.) and superfluid
  helium pulses (Bertman–Sandiford). IJHMT 1995: no abstract available.
- **Not accessible.** Full texts of the three 1995 papers and the monograph.
- **Verified statements of the equal-lag reduction elsewhere.**
  (a) Zhang, Tzou & Chen (2009 book chapter, arXiv:1511.03566, read): the first-order DPL energy equation "reduces to the
  parabolic conduction equation (3) if τ_q = τ_T". (b) Ordóñez-Miranda & Alvarado-Gil, *Nanoscale Res. Lett.* **6**, 327
  (2011) (full text): for the exact-delay DPL law with τ_q = τ_T the law reduces to Fourier's "except for a trivial shift
  in the time scale". (c) Kovács 2024 review (row 9). (d) Zhmakin preprint (row 8).
- **VERDICT.** The Fourier equivalence on Δ = {τ_q = τ_T} is classical in the DPL literature: LITERATURE SUPPORTED by a
  Tzou-coauthored text and two reviews. Its presence in Tzou 1995 itself is SECONDARY (not verified; the abstract does not
  mention it). What it anticipates in B: that Δ is a Fourier-equivalent, hence degenerate, set (the inclusion Δ ⊂ Σ and the
  invisibility of the common lag on Δ). It does not anticipate Σ = Δ (the absence of any other degeneracy, in particular at
  matched effusivity), the ranks, injectivity off Δ, or the fibres.

### 2.2 Tang & Araki (1999, 2000)

- **Metadata (verified).** D. W. Tang, N. Araki, "Wavy, wavelike, diffusive thermal responses of finite rigid slabs to
  high-speed heating of laser-pulses", *Int. J. Heat Mass Transf.* **42**(5), 855–860 (1999),
  doi:10.1016/S0017-9310(98)00244-0 (the paper the referee cites). D. W. Tang, N. Araki, "Non-Fourier heat conduction
  behavior in finite mediums under pulse surface heating", *Mater. Sci. Eng. A* **292**(2), 173–178 (2000),
  doi:10.1016/S0921-5093(00)01000-5 (the Crossref title contains the typo "condution").
- **Not accessible.** Both full texts (Elsevier, no open copy); abstracts elided in all indices.
- **SECONDARY.** Zhmakin (arXiv:2212.12922, p. 16, read as a rendered page): "Tang & Araki [MSEA 2000] (see also … Shen &
  Zhang) distinguish four heat propagation modes of the DPL model": wave (τ_T = 0), wavelike (0 < τ_T < τ_q), diffusion
  (τ_T = τ_q), over-diffusion (τ_T > τ_q). The review dates Shen & Zhang to 2013; Crossref gives *Int. J. Heat Mass Transf.*
  **51**, 1713–1727 (2008), doi:10.1016/j.ijheatmasstransfer.2007.07.039. Citation contexts of the 1999 paper describe it as
  illustrating wavy, wavelike and diffusive behaviour of DPL slabs irradiated by short pulses, and quote lag values for
  chromium taken from it (τ_q = 0.136 ps, τ_T = 7.86 ps).
- **VERDICT.** The "diffusion mode" (τ_T = τ_q) is SECONDARY-attributed to the 2000 paper; the 1999 paper belongs to the
  same programme but the classification could not be located in it. Same bearing on B as §2.1 (naming of Δ), nothing more.
  For the manuscript: cite both Tang–Araki papers only after the 2000 text (or the 1999 text) has been checked; meanwhile the
  Zhmakin and Kovács reviews are verified pointers.

### 2.3 Kitamura & Nakagiri (1977), Courdesses (1983), Pierce (1979)

- **SOURCE (ABSTRACT) SAYS.** Kitamura & Nakagiri: identifiability (deterministic) of spatially varying and constant
  parameters of a linear 1D parabolic equation; only parameters in the equation, not in the boundary conditions, are
  unknown; distributed and pointwise measurements; the conditions depend on the state profile (distributed case) or on the
  detector position and the initial or input functions (pointwise case). Courdesses (1983): the solution expression of
  Kitamura–Nakagiri "seems to be incorrect" and consequently some results; a corrected expression and a necessary and
  sufficient condition for identifiability of constant parameters with pointwise measurement are given. Pierce: the
  solution observed at one point over a finite interval determines the eigenvalues (under conditions) but not the
  coefficients; in normal form the single coefficient is in some cases unique (spectral function or two spectra;
  Gel'fand–Levitan, Levinson).
- **Related, verified.** Nakagiri, "Review of Japanese work of the last ten years on identifiability in distributed parameter
  systems", *Inverse Problems* **9**, 143–191 (1993), doi:10.1088/0266-5611/9/2/001 (abstract only). Gentili,
  *Inverse Problems* **7**, 77–84 (1991), doi:10.1088/0266-5611/7/1/008 — SOURCE (ABSTRACT) SAYS: a uniqueness theorem for
  heat capacity and conductivity memory functions from boundary flux and temperature, proved "by using the Laplace
  transforms and the thermodynamic restrictions" (already in 08 item 30).
- **Not accessible.** SIAM full texts.
- **VERDICT.** Correct as method citations: identifiability of constant coefficients of 1D parabolic problems from point
  or boundary data on a **known domain** is classical (Kitamura–Nakagiri with the Courdesses correction; Pierce), and
  Laplace-transform uniqueness for heat conduction with memory is classical (Gentili). The specific claim that these
  papers use "Laplace-domain analyticity" is not confirmed from the abstracts (Pierce's route is spectral; the two are
  equivalent for exponential sums, OUR INTERPRETATION). None treats an unknown thickness, a hyperbolic or Jeffreys-type law,
  a semi-infinite substrate, or a scale gauge. If Kitamura–Nakagiri is cited, Courdesses 1983 must be cited with it.

### 2.4 Lee & Cahill (1997)

- **SOURCE (ABSTRACT) SAYS.** *J. Appl. Phys.* **81**, 2590: 3ω measurements on 20–300 nm PECVD SiO₂ and SiNₓ films on
  Si, 78–400 K; below ≈ 50 nm the apparent conductivity decreases with thickness, interpreted as an interface resistance
  R_I ≈ 2×10⁻⁸ K m² W⁻¹ at room temperature, equivalent to the resistance of ≈ 20 nm of SiO₂. *Microscale Thermophys. Eng.*
  **1**, 47: 8–200 nm SiO₂ and MgO films; the apparent conductivity is reduced by the heater/dielectric and
  dielectric/substrate interface conductances added in series; for SiO₂ at 250 K the series conductance is equivalent to a
  20-nm SiO₂ film.
- **OUR INTERPRETATION (one-line check).** For a Fourier film (on Δ) with a resistance R at the film/substrate interface,
  expanding Z = B⁻¹(1 + rE)/(1 − rE) to first order in A = ξ₁√p gives Z ≈ R + 1/(b_s√p) + ξ₁(1/b − b/b_s²) + O(ξ₁², Rξ₁).
  At first order in thickness, R and ξ₁ therefore trade along (dξ₁, dR) ∝ (b b_s²/(b² − b_s²), 1), which is exactly the
  kernel direction of Theorem 2′ on Δ ∩ {R = 0, b ≠ b_s}. For b ≪ b_s (SiO₂ on Si, b/b_s ≈ 0.1) the coefficient reduces to
  the pure series-resistance picture e/λ = ξ₁/b of Lee–Cahill; for comparable effusivities the film heat capacity enters.
- **VERDICT.** Lee–Cahill anticipate the physical mechanism of the R = 0 stratum in the thin-film limit and its classical
  remedy (thickness series). Theorem 2′ is a refinement (exact at every thickness and frequency, but only infinitesimal at
  R = 0; Theorem 1′ keeps global injectivity). Cite as the referee asks.

### 2.5 The two reviews used by the referee

- **arXiv:2212.12922** — A. I. Zhmakin, "The Zoo of Non-Fourier Heat Conduction Models" (v1, 25 Dec 2022; no journal
  reference; a Springer monograph by the same author, *Non-Fourier Heat Conduction* (2023), doi:10.1007/978-3-031-25973-9,
  was not read). SOURCE SAYS (p. 16): the four Tang–Araki modes (§2.2); the DPL law reduces to the single-phase-lag law for
  τ_T = 0 and to Fourier's law for τ_T = τ_q = 0; written as q(r, t) = −λ∇T(r, t + τ_T − τ_q), the exact-delay DPL depends only
  on the difference of the lags.
  *Caution (OUR INTERPRETATION).* The last statement is about the exact-delay law. The first-order (Jeffreys-type) DPL law,
  which is the 1D GK law of this paper with τ_q = τ_R, τ_T = τ_ℓ, depends on both lags separately (rank 4 off Δ). Both models
  are Fourier-equivalent on Δ. A manuscript sentence citing this review must not imply that only τ_R − τ_ℓ is observable.
- **arXiv:2302.13931** — R. Kovács, "Heat equations beyond Fourier: from heat waves to thermal metamaterials", published as
  *Phys. Rep.* **1048**, 1–75 (2024), doi:10.1016/j.physrep.2023.11.001 (Crossref; the arXiv record carries no journal
  reference; already `Kovacs2024` in refs.bib). SOURCE SAYS: §6, on Tzou's DPL law (citing Tzou 1995), the law reduces to
  Fourier's for τ_q = τ_T = 0 "but also reproduces Fourier's solutions for any τ_q = τ_T", likened to the Fourier resonance
  (citing Fehér et al., *J. Non-Equilib. Thermodyn.* **46**, 403 (2021)); the GK structure is recovered only in the
  T-representation of the DPL law; §4.3, the Fourier resonance of the GK-type T-representation is named (citing Ván et al.,
  *EPL* **118**, 50005 (2017)).
- **VERDICT.** Both are reliable secondary sources for "Δ is the DPL equal-lag (diffusion-mode) set and coincides with the
  GK Fourier resonance". Kovács 2024 is already in the bibliography and is the best single citation for that sentence.

### 2.6 What the classical Δ statement does and does not cover (OUR INTERPRETATION)

The DPL equal-lag reduction concerns the bulk law. It transfers to the film-on-substrate boundary-value problem when the
interface conditions are continuity of temperature and normal flux (or a Kapitza jump), the convention used in the DPL
multilayer literature and in Assumption M. Under the entropic or kinetic interface conditions (objection B-2), the interface
condition involves τ_ℓ and Δ is no longer Fourier-equivalent. The classical attribution of Δ (B-1) and the convention
dependence (B-2) should therefore be stated together: Δ is classical **under the local convention**, and it is exactly under
that convention that Theorem 2 adds Σ ⊂ Δ.

---

## 3. Task 2 — Ordóñez-Miranda & Alvarado-Gil, *J. Heat Transfer* 132, 061302 (2010)

### 3.1 Routes tried (all open and legal)

Publisher page (ASME, HTTP 403 to automated retrieval) · Unpaywall (no open copy) · Semantic Scholar (closed, abstract
elided) · OpenAlex (abstract, no open location) · ResearchGate (HTTP 403) · HAL (the author's deposits begin in 2013) ·
CINVESTAV repository (the 2011 doctoral thesis of J. Ordóñez-Miranda, reported in a web summary of his CV as "Heat
conduction in nanomaterials: beyond the Fourier's law" — SECONDARY; the CV file now redirects — is not deposited; theses of
later students of J. J. Alvarado-Gil are on other topics) · Europe PMC (the open 2011 paper of the same authors, read; it does
not cite the 2010 paper) · Google Books API (daily quota exhausted) · Semantic Scholar citation contexts of the 2010 and 2009
papers (read) · OpenAlex forward citations (25 records; only two open access, neither reproduces the formulas) · Springer
abstract page of the 2010 *Granular Matter* companion (read).
**The full text could not be obtained.** It must be requested through a library (inter-library loan) or from the
first author before any novelty sentence about B is finalised.

### 3.2 What is established (SOURCE (ABSTRACT) SAYS, and Crossref)

- Bibliographic data: J. Ordóñez-Miranda, J. J. Alvarado-Gil, "Determination of time-delay parameters in the dual-phase
  lagging heat conduction model", *J. Heat Transfer* **132**(6), 061302 (2010), published 25 March 2010, 36 references.
- Configuration: a semi-infinite layer in thermal contact with a finite layer, excited by a modulated heat flux; DPL model.
- Observable: the normalized amplitude and phase of the spatial component of the oscillatory surface temperature.
- Regime: when the effusivities of the layers are "quite different", a frequency range exists where the normalized
  amplitude and phase show strong oscillations.
- Result: this behaviour yields explicit formulas determining simultaneously the two time parameters "as well as additional
  thermal properties of the finite layer"; single-phase-lag limits are discussed.
- Reference list (Crossref): Tzou 1989, 1995 ×3 and monograph; Quintanilla 2002 (exponential stability); DPL multilayer papers
  (Ho et al. 2003, Ramadan 2009, Al-Huniti 2004, Lee 2007, Liu 2007, Xu 2005); photothermal two-layer methods (Salazar 1998,
  Pichardo 2001, Marín 1996, Vargas 1988, McDonald 1978); measurement-method precursors O-M & A-G 2009 (*Int. J. Therm. Sci.*
  **48**, 2053) and Basirat Tabrizi & Andarwa 2009 (*Int. Commun. Heat Mass Transf.* **36**, 186); Roetzel 2003 and Antaki
  2005 (non-Fourier experiments).
- **Not established:** whether the thickness is known; which "additional properties" are returned; whether the substrate is
  DPL or Fourier and whether its properties are known; the interface condition; which DPL version (first-order or exact
  delay); whether the formulas are exact or asymptotic in the effusivity contrast; whether a reference measurement is needed
  for the normalization; any uniqueness or multiplicity statement. The abstract contains no identifiability or uniqueness
  claim, and no citation context found describes the assumptions.

### 3.3 Indirect evidence from the same programme

| Source | Access | What it shows |
|---|---|---|
| O-M & A-G, *Int. J. Therm. Sci.* **48**, 2053–2062 (2009), doi:10.1016/j.ijthermalsci.2009.03.008 (Cattaneo precursor) | SECONDARY (citing paper, 2013) | A finite layer in contact with a semi-infinite layer under Cattaneo's law, Dirichlet and Neumann boundary conditions |
| Vales-Pinzón, Ordóñez-Miranda, Alvarado-Gil, *J. Appl. Phys.* **112**, 064909 (2012), doi:10.1063/1.4754552 (Fourier analogue) | abstract | Signal of the two-layer system **normalized by the signal of the illuminated layer alone**; four characteristic frequencies (normalized amplitude through unity and maximum, phase through zero and minimum) are independent of the effusivities and give the **thermal diffusivity** of the illuminated layer; the effusivity then follows from the full complex expression. A characteristic frequency can return a diffusivity only if the layer thickness is an input (OUR INTERPRETATION) |
| O-M & A-G, *Granular Matter* **12**, 569–577 (2010), doi:10.1007/s10035-010-0195-6 (DPL, thermal-wave resonant cavity) | abstract page | Diffusivity and both lags of the cavity medium from a cavity-length scan (controlled lengths) |
| O-M & A-G, *Mech. Res. Commun.* **37**, 276–281 (2010), doi:10.1016/j.mechrescom.2010.03.002 | SECONDARY (self-citation context) | In the Cattaneo limit, diffusivity and relaxation time of a semi-infinite layer from amplitude and phase at a position x ≠ 0 (known depth) |
| O-M & A-G, *Nanoscale Res. Lett.* **6**, 327 (2011), doi:10.1186/1556-276X-6-327 | full text | Exact-delay DPL: equal lags give Fourier up to a time shift; the exact-delay model predicts unbounded temperature growth at long times and is not a general description |

The programme consistently works with known geometry (thickness, depth, cavity length). This makes a known film thickness in
the 2010 paper likely, but it is not established.

### 3.4 What any such formula can and cannot determine (OUR INTERPRETATION)

For a homogeneous first-order DPL layer on a substrate, the front-face response depends on the layer only through
θ = (ξ₁, b, τ_q, τ_T) (Theorem G; the lags are gauge-invariant). Oscillations of the normalized signal come from the
interference factor e^{−2σl}, σ² = (iω/a)(1 + iωτ_q)/(1 + iωτ_T): for τ_T < τ_q there is a window 1/τ_q ≲ ω ≲ 1/τ_T where
Im σ ≈ ω(τ_q/a)^{1/2}, and at large effusivity contrast (|r| ≈ 1) the extremum frequencies are set by the round-trip phase
2lω(τ_q/a)^{1/2}, i.e. by ξ₁ and τ_q, while τ_T controls damping and window edges. Consequences:
(i) formulas returning τ_q, τ_T, l/√a and the effusivity are consistent with Theorem 1 and recover a point of θ in a region of
parameter space; (ii) any formula returning a, λ or l separately must use a known length (or another known property), else it
would contradict Theorem G; (iii) no oscillation window exists for τ_T ≥ τ_q, so the method, if as described, does not cover
the over-diffusive side τ_ℓ > τ_R (the side of the conserving closure is τ_ℓ < 4τ_R/5, which it may cover); (iv) if the
normalization requires a reference measurement on the layer alone (as in the 2012 Fourier analogue), the data set is richer
than the single front-face impedance assumed in Theorem 1.

### 3.5 Verdict and checklist

**Partial prior art for Theorem 1**: a constructive determination of both lags (and further layer properties) for the same
configuration, in a regime of large effusivity contrast, published 2010. It is not a uniqueness theorem, does not give the
exceptional set, the ranks, the fibres or the gauge, and its assumptions are unverified. Checklist for the reading:
(1) thickness known? (2) which properties returned? (3) substrate law and known properties; (4) interface condition;
(5) first-order or exact-delay DPL; (6) exact or asymptotic formulas, validity window (τ_T < τ_q?); (7) reference
measurement for normalization? (8) any uniqueness or multiplicity remark, any mention of τ_q = τ_T as degenerate.
The same checklist applies to the 2009 Cattaneo precursor and to Basirat Tabrizi & Andarwa 2009 (described in a citing paper
as a frequency-response method for the lag constants; full text unread).

---

## 4. What remains potentially novel for B (precise statement and scope)

**Statement (potentially novel within the literature searched).** For a homogeneous film obeying the first-order DPL law
(equivalently, the 1D GK law with τ_q = τ_R, τ_T = τ_ℓ), on a semi-infinite Fourier substrate of known effusivity b_s,
with front-face flux excitation and front-face temperature detection, zero initial state, no interior source and the
**local interface convention** (continuity of temperature and flux, or a Kapitza jump R), from exact complex data on any
frequency set with an accumulation point:
1. **Exact exceptional set.** Σ = Δ: the Jacobian of θ = (ξ₁, b, τ_R, τ_ℓ) ↦ Z has rank 4 everywhere off Δ, including the
   matched-effusivity stratum b = b_s and the boundary strata τ_R = 0, τ_ℓ = 0; rank 3 on Δ \ D₂ and 2 on D₂. The new part is
   Σ ⊂ Δ, in particular that matched effusivity, degenerate for Fourier films, is regular off Δ.
2. **Global identifiability with explicit fibres.** Injective off Δ; fibres the half-lines {(ξ₁, b, t, t)} on Δ \ D₂ and a
   two-parameter set on D₂.
3. **Interface resistance.** Rank 5 and injective off Δ; on Δ rank 4 for R > 0 and 3 for R = 0, with the explicit kernels.
4. **Finitely many frequencies.** Two generic frequencies give local rank 4 at every θ ∉ Δ but not global uniqueness
   (certified twin pair).

**Anticipated (cite).** Δ ⊂ Σ: DPL equal-lag reduction (Zhang–Tzou–Chen 2009; Kovács 2024; Tang–Araki "diffusion mode",
secondary) and GK Fourier resonance (Both et al. 2016; Kovács 2018). ξ₁ lost at b = b_s and the scale gauge for Fourier
films (Krapez–Rigollet 2017). Constructive regional determination of both lags for the same configuration (O-M & A-G 2010).
The R = 0 trade-off in the Fourier thin-film limit (Lee–Cahill 1997). The uniqueness technique on known domains
(Kitamura–Nakagiri 1977 with Courdesses 1983; Pierce 1979; Gentili 1991). The gauge itself is dimensional analysis.

**Conditional on.** Assumption M and the local convention (B-2: under the entropic or kinetic conventions, items 1 and 3
change); exact, noise-free data (B-3: the practical rank in the FDTR band is 3).

**Scope of the search behind "potentially novel".** 08-literature (66 + 10 queries, 11 forward and 6 backward citation
traces, neighbouring theorem families) plus this branch: verification of the ten references above; citation-context mining
for O-M & A-G 2009/2010 and Tang–Araki 1999/2000; 5 new queries (§8); a recent Fourier neighbour found here,
D. Wang and H. Ban, arXiv:2609.13536 (2026), which proves a global structural non-identifiability (Schur-equivalent
conductivity tensors give identical surface responses for all frequencies, delays, beam geometries and offsets) and states
that SVD-based screening cannot rule out globally equivalent solutions — Fourier only, consistent with the local-versus-global
distinction of item 4, not prior art for the non-Fourier statements. **Residual risk:** the unread full texts of
O-M & A-G 2010 and 2009, Tang–Araki 1999/2000 and Tzou's monograph.

---

## 5. Task 3 — Positioning of workstream A

### 5.1 What each paper claims about relaxation times from FDTR

| Paper | Access | Configuration | Relaxation-time claim, with numbers |
|---|---|---|---|
| Beardo et al., *Phys. Rev. B* **101**, 075303 (2020), doi:10.1103/PhysRevB.101.075303 | full text | Si with a metal transducer (5 nm Cr adhesion layer), broadband FDTR data of Regner et al., *Nat. Commun.* **4**, 1640 (2013), discussed up to and above 100 MHz, 81–416 K; KCM (GK-type) with ab initio parameters | τ is an **input**, not fitted: τ = 1002, 148, 42, 27 ps and ℓ = 3127, 1047, 185, 142 nm at 81, 154, 311, 416 K (Table I). Only a temperature-independent interface factor R/R_min = 2.71 is fitted to low-frequency data. At 311 K KCM and Fourier diverge from 5 MHz, "dominated by nonlocal effects"; ωτ ≪ 1 so resistive hydrodynamic effects are negligible at all frequencies (ωτ = 0.053 at 200 MHz for 42 ps, OUR computation). At 81 K resistive effects become comparable with nonlocal ones at the highest frequencies. No relaxation time is extracted |
| Beardo et al., *Sci. Adv.* **7**, eabg4677 (2021), doi:10.1126/sciadv.abg4677 | full text (arXiv:2007.05487) | Bare natural Ge, no transducer, 30 kHz–200 MHz | **Fitted** τ = 500 ps with α = 3×10⁻⁵ m²/s (fitted simultaneously), v ≈ 250 m/s, room temperature (3D hyperbolic model); κ = 51 W m⁻¹K⁻¹ from 3ω. With a 60-nm Au transducer wave-like effects are not observed (transducer and interface dominate). Disputed by Huberman, Zhang & Abou Haibeh, *J. Phys.: Condens. Matter* **37**, 435601 (2025), doi:10.1088/1361-648X/ae0dd6 (phase lag alone not a suitable observable; abstract) |
| Qiu et al., "PRISM", arXiv:2611.02475 (2 Nov 2026) | full text | Au/Si, Au/Ti/Si, Au/SiO₂/Si, Au/MoS₂ FDTR; optical transfer matrix + two-temperature + GK substrate; noise-whitened SVD screening | GK branch with τ = 42 ps and ℓ = 185 nm **fixed** at the 311 K values of Beardo 2020; τ is **not among the candidate parameters** (13 candidates: k_Si, C_Si, h_Au, w_x, w_y, five two-temperature parameters, ℓ, α_GK, b_NF). With phase-noise floor 0.1° and retention σ_m/σ_1 > 0.1, ℓ and α_GK (with four two-temperature parameters) are screened as unidentifiable and held at priors; the fit targets are k_Si and b_NF. Result (Au/Ti/Si stacks): k_Si = 141 W m⁻¹K⁻¹, against 116–137 for local Fourier fits across spot sizes and frequency windows; 1σ relative uncertainties 6–18 % (conductivities) and 11–65 % (interface parameters) with 5 % priors on fixed inputs. The referee's "τ held fixed as unidentifiable" should read "τ fixed at its first-principles value; ℓ and α_GK screened as unidentifiable" |
| Camacho de la Rosa et al., *J. Appl. Phys.* **137**, 155103 (2025), doi:10.1063/5.0257299 | full text (CC BY) | Simulation only: SiC half-space; Au(100 nm)/SiC(1 µm)/Si, 10 µm FWHM spot, 1 mW; Cattaneo | τ_SiC = 10⁻⁸ s ("of the order but not equal to" phonon relaxation times). Noise-free sensitivity dφ/d ln τ is significant near ωτ ≈ 1, where fitting τ "will be possible". In their band relaxation times below 10⁻⁹ s are "practically identical" (Fourier limit); at 3×10⁻⁹ s the Fourier/non-Fourier differences lie outside common FDTR ranges. Conclusion cites a phase resolution of "a few millidegrees" (attributed to Kirsch et al. 2024; not found in the authors' companion presentation, paper text unread, 08 item 36). No noise model, Fisher matrix, interface conductance or transducer nuisance |
| Saurav & Mazumder, ASME HT2020-8920 (2020), doi:10.1115/HT2020-8920; *J. Heat Transfer* **144**, 013501 (2021), doi:10.1115/1.4052655 | 08 item 51 | Au/Si FDTR 200 kHz–200 MHz, numerical | Effective τ_Si ≈ 36.5 ps; Fourier and hyperbolic phase lags "identical" |

### 5.2 Positioning statement for A (proposed)

- Beardo 2020 and Saurav–Mazumder state that a tens-of-ps relaxation time is negligible in the FDTR band because ωτ ≪ 1;
  Camacho 2025 places the Fourier limit of its own configuration at τ ≈ 1 ns; PRISM never fits τ. Result A is consistent with
  all four and supplies what they lack: a quantitative Cramér–Rao bound, with instrument, interface and transducer nuisances
  free, for a 10–100 ps Cattaneo (memory) time of an AlN film, and the leading-order law δφ ≈ w ωτ_R W(ω) that explains why
  (peak |∂φ/∂ ln τ_R| = 0.021° at 65 MHz for 500 nm, review check R8).
- Result A contradicts the abstract-level reading of Camacho 2025 ("a practical approach for characterizing the relaxation
  times") for phonon-scale times: their numbers concern τ ≥ 1 ns and noise-free sensitivity.
- Beardo 2021 is the only FDTR extraction of a relaxation time found: τ = 500 ps (≥ 5× the range of A) on a bare
  semiconductor without transducer or interface, and disputed. Its own remark that a metal transducer masks wave-like effects
  agrees with A, where the ≈ 2 ns transducer–interface time constant (Angeles 2023, 08 item 40) and the free interface
  resistances dominate the information budget.
- Novelty wording supported by this branch (scope: 08 §4.5, query RL05, full texts of PRISM and Camacho 2025): "first
  Cramér–Rao estimate of the FDTR resolution of a 10–100 ps memory time with instrument, interface and transducer nuisances,
  within the literature searched". It is a quantitative bound, not a new physical effect; the qualitative statement
  (ωτ ≪ 1 makes τ invisible) is in Beardo 2020.

---

## 6. Citation changes recommended for the manuscript (not applied)

1. Fourier resonance subsection: keep Both 2016 and Kovács 2018 for the GK naming; add the DPL equal-lag case with Tzou 1995
   for the DPL law and Kovács 2024 (already in refs.bib) and Zhang–Tzou–Chen 2009 for the reduction; add Tang–Araki ("diffusion
   mode") only after checking the 2000 text.
2. Identifiability subsection: cite O-M & A-G 2010 for the constructive regional determination of both lags (paraphrase of the
   abstract only until the text is read); cite Kitamura–Nakagiri with Courdesses 1983, Pierce 1979 and Gentili 1991 for the
   method; Lee–Cahill 1997 (JAP, optionally MTE) for the R = 0 stratum.
3. FDTR section: position A against Beardo 2020, Beardo 2021 with Huberman 2025, PRISM (arXiv:2611.02475) and Camacho 2025 as
   in §5.2; correct the description of PRISM.
4. Do not attribute the equal-lag reduction to Tzou 1995 specifically, and do not cite the exact-delay "difference only"
   statement for the first-order model.

---

## 7. Unresolved prior-art uncertainty

1. O-M & A-G 2010 (and the 2009 precursor): assumptions (thickness, normalization, substrate, DPL version) and any uniqueness
   statement — full text needed (library or author). Highest priority.
2. Tzou 1995 and the monograph: exact location of the equal-lag statement.
3. Tang–Araki 1999/2000: location of the four-mode classification.
4. Kitamura–Nakagiri 1977, Pierce 1979: techniques (eigenfunction/Dirichlet-series versus Laplace-domain arguments).
5. Lee–Cahill 1997: full-text statement of the thickness-series analysis.
6. Basirat Tabrizi & Andarwa 2009: frequency-response method for lag constants, content unknown.

---

## 8. Query and record log

- API verification records: `searches/R-L/verify/` (Crossref `*.json`, OpenAlex `oa_*.json`, Semantic Scholar `s2_*.json`
  including citation contexts `s2_cit_*.json`, arXiv `arxiv_ids.xml`, HAL `hal_*.json`).
- Multi-index searches (Crossref + OpenAlex; Semantic Scholar HTTP 429), `searches/R-L/queries.txt`, outputs in
  `searches/R-L/runs/`: RL01 dual-phase-lag equal phase lags reduces to Fourier diffusion mode · RL02 determination of phase
  lags dual-phase-lag modulated photothermal two-layer system · RL03 thermal relaxation time and thickness simultaneous
  determination thermal wave layer · RL04 identifiability of lag times dual-phase-lag model sensitivity correlation · RL05
  frequency-domain thermoreflectance relaxation time Cattaneo uncertainty estimation. Outcome: no new anticipation of B or A;
  relevant hits were already known (O-M & A-G 2009, Kovács 2024, Maillet 2019, Camacho 2025); RL01 surfaced Zhmakin's
  *Phase-Lag Models* chapter (2023) and an electric-circuit analogue of first-order DPL (*J. Heat Transfer*, 2024,
  doi:10.1115/1.4067256, not read).
- Web searches (12): title of O-M & A-G 2010 (two forms), abstract of O-M & A-G 2009, thesis of Ordóñez-Miranda,
  conference versions, InTech chapters, Tzou equal-lag wording (two forms), Tang–Araki classification, Lee–Cahill full text,
  Kitamura–Nakagiri method, PRISM.
- Full texts downloaded for reading were kept outside the campaign (temporary workspace); only records and conclusions are
  stored here.
