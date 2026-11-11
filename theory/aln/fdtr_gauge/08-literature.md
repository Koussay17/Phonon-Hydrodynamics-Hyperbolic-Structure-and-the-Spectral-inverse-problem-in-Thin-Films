# 08 — Literature trail: FDTR detectability of the film memory time (A) and rank/gauge of the 1D non-Fourier front-face response (B)

Campaign `20261111-fdtr-gauge-memory`, 11 November 2026. Branch L (literature, two passes). This trail extends, and does not
repeat, the audit of 2–3 October 2026 (`D:\ResearchLab\orchestration\campaigns\20261002-paper1-prior-art\08-literature.md`,
summarised in `C:\Users\Koussay\these\paper\PRIOR_ART.md`, claim C5) and the reading notes of the thesis repo
(`notes/05`, `06`, `07`, `15`, `16`, the Camacho reading note `notes/camacho_2025.pdf`). Sources already attributed there
(Krapez–Rigollet 2017, Krapez 2023, Krapez 2019, Both et al. 2016, Kovács 2018, Hennessy–Myers 2021, Beardo et al. 2020,
Camacho de la Rosa et al. 2025, Cahill 2004) are only re-checked where this campaign needs a new fact from them.

**Evidence labels.** SOURCE SAYS / SOURCE DERIVES / SOURCE PROVES / SOURCE OBSERVES (numerical or experimental) /
SECONDARY (a later paper's statement about an earlier one, not checked in the earlier text) / OUR INTERPRETATION (derivation or
reading done here, to be checked by the analytic, symbolic and asymptotic branches). Access level is stated for every item:
*read* (full text or a named section), *text search* (full text searched for specific statements), *abstract*, *metadata*.

**Status labels used in 01-known-results.md.** ESTABLISHED / PARTIAL / CONJECTURED / NUMERICAL / DISPUTED / OPEN.

---

## 1. Databases, access and search outcomes

| Source | Use | Outcome |
|---|---|---|
| Local multi-index search (Crossref + Semantic Scholar + OpenAlex), 15 results per provider | First pass 56 queries (B01–B30, A01–A18, S01–S08, `searches/queries.txt`); second pass 10 equivalent-term queries (P01–P10, `searches/queries2.txt`). Raw JSON/MD under `searches/<ID>/`; compiled lists `searches/compiled-hits.txt`, `searches/compiled-hits-2.txt`; logs `searches/batch-log*.txt` | Crossref 66/66, OpenAlex 66/66, Semantic Scholar 2/66 (HTTP 429). Broad queries were noisy (photothermal therapy, dielectric relaxation); no absence claim rests on a single query. |
| Crossref API (direct) | Metadata verification of every cited DOI; publication list of J.-C. Krapez (2011–2026) | Worked. Verification records `searches/verify/verify_1..5.{txt,json}`. |
| OpenAlex API (direct) | Abstracts and forward citations | Worked. Forward-citation lists `searches/forward_citations_1.txt`, `_2.txt` (11 anchor papers, 477 citing records returned; cited-by counts sum to 479). |
| Semantic Scholar API (direct) | Abstracts of paywalled papers | Records found; abstracts elided by publishers (ASME, Elsevier). |
| Unpaywall | Open-access locations | Used to locate open copies (arXiv, OSTI, publisher OA). |
| arXiv (full text) | 2007.05487, 2206.02769, 2412.02227, 2611.02475, 2610.09192, 2609.27618, 1501.03174, 1508.03864, 1709.06825, 2307.00665, 2402.04241 | Read or searched as stated per item. The usual arXiv and Semantic Scholar search interfaces were unavailable this session; arXiv was reached through its public PDF and query endpoints. |
| Publisher/repository OA PDFs | Wilson–Cahill 2014 (Nature Commun.), Angeles et al. 2023 (OSTI), Regner et al. 2013 (RSI and Nature Commun., repository copies), Saurav–Mazumder HT2020 manuscript (NSF repository), NIST-hosted presentation of Kirsch et al. | Read as stated. |
| Thesis repository local copies | `notes/155103_1_5.0257299.pdf` (Camacho 2025), `.build/research/beardo_2020.{pdf,txt}`, `.build/research/hennessy_myers_2021.*`, `notes/maillet_2000_index.md` | Read (Beardo 2020 Table I and App. F re-read for this campaign). |
| Project reference library | Checked first | 4 items (Camacho 2025, Hoque 2024, Krapez 2018/2019, Chen–Gamba–Li–Wang 2025); nothing new for this campaign. |
| General web search | 12 targeted queries (titles, critiques, noise figures, the 9/5 relation, phase calibration) | Located the Huberman critique, the PRISM preprint, the NIST presentation, the arXiv paper quoting the DPL–GK 9/5 relation. |

**Not accessible (affects verdicts, see §7).** ScienceDirect pages (HTTP 403), HAL documents (bot challenge: Maillet 2019 review,
Ezzahri–Joulain 2012 preprint), ASME full texts (Ordóñez-Miranda & Alvarado-Gil 2010 *J. Heat Transfer*; Tzou 1995), AIP
paywalled full texts (Schmidt 2009; Yang–Ziade–Schmidt 2016; Ziade 2020; Liu et al. 2013), the AIP PDF of Kirsch et al. 2024
(challenge page; the authors' NIST-hosted presentation was read instead), Tzou's monograph.

---

## 2. Query log

First pass (direct terms, A and B): B01 photothermal radiometry parametric identifiability coating thickness · B02 identifiability
thermal properties coating substrate photothermal · B03 structural identifiability heat conduction parameter estimation · B04 linear
dependence sensitivity coefficients thickness diffusivity conductivity · B05 thermal quadrupole identification thin film · B06 front
face periodic method thermal effusivity coating · B07 front face flash method thin layer diffusivity estimation · B08 dual phase lag
parameter estimation · B09 dual phase lag thermal wave photothermal · B10 estimation relaxation time hyperbolic heat conduction inverse
problem · B11 inverse problem Cattaneo relaxation time identification · B12 Guyer-Krumhansl parameter estimation heat pulse · B13
Fourier resonance Guyer-Krumhansl · B14 Ordonez-Miranda Alvarado-Gil dual phase lag · B15 non-Fourier photothermal thin film thermal
waves · B16 hyperbolic heat conduction thin film frequency domain response · B17 thermal impedance non-Fourier layer · B18 FDTR
non-Fourier · B19 TDTR hyperbolic heat conduction · B20 Tzou dual phase lag lagging behavior · B21 identifiability dual phase lag model
· B22 sensitivity analysis thermoreflectance parameter correlation · B23 uniqueness inverse problem thin film thickness diffusivity ·
B24 simultaneous estimation thickness diffusivity coating photothermal · B25 photothermal depth profiling non-uniqueness · B26
Liouville transformation effusivity identification · B27 dual phase lag equivalence Guyer-Krumhansl · B28 inverse coefficient problem
telegraph equation uniqueness · B29 Laplace transform uniqueness inverse problem layered heat conduction · B30 thermal wave
interferometry dual phase lag layered · A01 FDTR phonon mean free path spectroscopy · A02 broadband FDTR · A03 frequency dependence
thermal conductivity alloys TDTR · A04 anisotropic failure of Fourier theory TDTR · A05 hydrodynamic heat transport silicon
thermoreflectance · A06 GK thermoreflectance · A07 FDTR uncertainty sensitivity · A08 pulse accumulation TDTR/FDTR · A09 quasi-ballistic
modulation frequency thermoreflectance · A10 thermal conductivity spectroscopy · A11 relaxation time measurement thermoreflectance
non-Fourier · A12 AlN thin film thermoreflectance · A13 GaN quasi-ballistic frequency-dependent conductivity · A14 Cattaneo relaxation
time thin film picosecond measurement · A15 thermal relaxation time measurement dielectric crystal · A16 truncated Lévy thermoreflectance
· A17 FDTR phase noise lock-in · A18 nonlocal heat conduction thermoreflectance kinetic collective model.

Equivalent terms (S, P): S01 thermal impedance spectroscopy hyperbolic heat equation · S02 transfer function hyperbolic heat
conduction multilayer · S03 sensitivity matrix rank photothermal · S04 non-Fourier thermoreflectance · S05 structural identifiability
thermal model · S06 global identifiability PDE · S07 phase lag heat conduction parameter identification experiment · S08
identification of thermal relaxation time · P01 Jeffreys-type heat equation parameter identification · P02 heat conduction with memory
kernel identification boundary measurement · P03 Gurtin–Pipkin memory inverse problem · P04 transmission line parameters from input
impedance and length · P05 lossy transmission line inverse problem travel time uniqueness · P06 thermal wave interferometry coating
thickness diffusivity correlation · P07 photoacoustic coating thickness diffusivity not separable · P08 identifiability two-layer
thermal system frequency response · P09 hyperbolic heat conduction modulated thermoreflectance relaxation time sensitivity · P10 second
sound FDTR phase lag.

Backward references traced from: Krapez–Rigollet 2017 (note 05), Camacho 2025 (refs 36–53), Beardo 2020 (App. F, refs), PRISM 2026
(refs 71, 77, 79), Chen–Wang–Jiang 2026 (refs 11, 29, 55), Baratifarimani–Shomali 2023 (ref 16).
Forward citations traced (OpenAlex, titles filtered by identifiability / relaxation / non-Fourier / thermoreflectance keywords):
Krapez–Rigollet 2017 (1 citing record), Ordóñez-Miranda & Alvarado-Gil 2010 *JHT* (25), Camacho 2025 (2), Hennessy–Myers 2021 (5),
Krapez 2023 (3), Krapez 2019 (5), Beardo et al. 2021 (87), Both et al. 2016 (101), Kovács 2018 (51), Saurav–Mazumder 2021 (17),
Yang–Ziade–Schmidt 2016 (182). Crossref publication list of J.-C. Krapez, 2011–2026: no non-Fourier extension of the gauge.

---

## 3. Workstream B — identifiability of the front-face response of a film on a substrate

### 3.1 Fourier, homogeneous film: the scale gauge and what is identifiable

1. **A. Rosencwaig, A. Gersho**, "Theory of the photoacoustic effect with solids", *J. Appl. Phys.* **47**, 64–69 (1976),
   doi:10.1063/1.322296 — abstract. Earliest standard forward model of periodic front-face heating of layered solids. Historical anchor.
2. **C. A. Bennett, R. R. Patty**, "Thermal wave interferometry: a potential application of the photoacoustic effect",
   *Appl. Opt.* **21**, 49 (1982), doi:10.1364/AO.21.000049 — abstract. SOURCE SAYS: the photoacoustic signal is recast in terms of
   thermal-wave interference in the layer. OUR INTERPRETATION: the interference form 1 + Γe^{−2ξ₁√p} makes the transit time and the
   reflection coefficient the only layer quantities; earliest form of the (ξ₁, Γ) structure found in this search.
3. **S. I. Yun, H. J. Seo**, "Photoacoustic measurements of thermal diffusivity or thickness of multi-layer solids", *J. Phys.
   Colloques* **44** (C6), C6-459–462 (1983), doi:10.1051/jphyscol:1983673 — abstract. Diffusivity *or* thickness of each layer, i.e.
   one given the other: the trade-off was understood in practice.
4. **S. J. Sheard, M. G. Somekh**, *Appl. Phys. Lett.* **53**, 2715–2716 (1988), doi:10.1063/1.100161 — abstract. Coating thickness
   from amplitude and phase versus frequency (photothermal radiometry); the diffusivity assumption is not stated in the abstract.
5. **A. Bendada**, *Meas. Sci. Technol.* **13**, 1946–1951 (2002), doi:10.1088/0957-0233/13/12/319 — abstract. Sensitivity of
   thermal-wave interferometry to coating properties; practical (noise-dependent), not structural.
6. **X. Guo, A. Mandelis, J. Tolev, K. Tang**, *J. Appl. Phys.* **121**, 095101 (2017), doi:10.1063/1.4977246 — abstract.
   SOURCE SAYS: thickness, diffusivity and conductivity are "unique only when their sensitivity coefficients are linearly independent
   over the range of measurements"; identifiability maps. The claim refuted by item 7.
7. **J.-C. Krapez, F. Rigollet**, Comment, *J. Appl. Phys.* **122**, 066101 (2017), doi:10.1063/1.4997733, arXiv:1708.07362 —
   read in note 05; abstract re-verified. SOURCE PROVES (elementary, quadrupole form and the Euler identity
   2a∂_aG + λ∂_λG + e∂_eG = 0 at every frequency): the three coating parameters are fully correlated; only Q₂ = e/√a and the effusivity
   ratio are identifiable; Q₂ is lost when the effusivities match (Γ = 0). SOURCE SAYS: lifting the indeterminacy needs an independent
   measurement of a, ρc or λ, a two-dimensional (lateral) experiment, or internal sources. **Fourier only.** Response: X. Guo,
   A. Mandelis, J. Tolev, K. Tang, *J. Appl. Phys.* **122**, 066102 (2017), doi:10.1063/1.4990378 — metadata only. Forward citations
   (OpenAlex): only the Response.
8. **J.-C. Krapez**, Comment, *J. Appl. Phys.* **134**, 056101 (2023), doi:10.1063/5.0151465; Response A. Mandelis, S. Kooshki,
   A. Melnikov, *J. Appl. Phys.* **134**, 056102 (2023), doi:10.1063/5.0155055 — read in note 06. Graded Fourier media: the response
   depends only on the effusivity profile in the Liouville (travel-time) coordinate; functional gauge.
9. **A. Pierce**, "Unique identification of eigenvalues and coefficients in a parabolic problem", *SIAM J. Control Optim.* **17**,
   494–499 (1979), doi:10.1137/0317035 — abstract. SOURCE SAYS: for a 1D parabolic equation observed at one point the eigenvalues can be
   determined, "the coefficients are not uniquely determined"; in normal (Liouville) form the single coefficient is, in certain cases,
   unique (Gel'fand–Levitan). OUR INTERPRETATION: the mathematical antecedent (by 38 years) of the graded gauge of item 8; it contains
   the homogeneous gauge as a special case.
10. **S. Kitamura, S. Nakagiri**, *SIAM J. Control Optim.* **15**, 785–802 (1977), doi:10.1137/0315050 — abstract (identifiability and
    non-identifiability of parameters of 1D parabolic systems from distributed or pointwise data). **S. Gutman, J. Ha**, *Math. Comput.
    Simul.* **79**, 2192–2210 (2009), doi:10.1016/j.matcom.2008.12.002 — metadata. Both assume a known domain.
11. **D. Maillet, S. André, J.-C. Batsale, A. Degiovanni, C. Moyne**, *Thermal Quadrupoles: Solving the Heat Equation through
    Integral Transforms* (Wiley, 2000) — local copy, index `notes/maillet_2000_index.md` (§1.3.4 p. 10, §3.3.1 p. 73): the wall quadrupole
    depends on (ke, λkS) only. **W. J. Parker, R. J. Jenkins, C. P. Butler, G. L. Abbott**, *J. Appl. Phys.* **32**, 1679–1684 (1961),
    doi:10.1063/1.1728417 — abstract: rear-face flash gives the diffusivity of a specimen of known thickness; with unknown thickness only
    e²/a is measured (OUR INTERPRETATION, standard).

### 3.2 Fourier, finite beams and multilayers: which combinations enter, and does a lateral wavenumber break the gauge?

12. **A. J. Schmidt, R. Cheaito, M. Chiesa**, *Rev. Sci. Instrum.* **80**, 094901 (2009), doi:10.1063/1.3212673 — abstract: FDTR
    "can simultaneously determine the thermal conductivity and heat capacity of a sample, provided the thermal diffusivity is greater
    [than] ≈3×10⁻⁶ m²/s", plus in-plane/cross-plane conductivities and boundary conductance. **J. Liu, J. Zhu, M. Tian, X. Gu, A. Schmidt,
    R. Yang**, *Rev. Sci. Instrum.* **84**, 034902 (2013), doi:10.1063/1.4797479 — abstract: "how different combinations of κ and C
    determine the frequency-dependent TDTR signals"; κ and C for bulk Si and SiO₂ films. **C. Wei, X. Zheng, D. G. Cahill, J.-C. Zhao**,
    *Rev. Sci. Instrum.* **84**, 071301 (2013), doi:10.1063/1.4815867 — abstract (heat capacity by dual-frequency TDTR). SOURCE OBSERVES:
    radial spreading with a known spot radius separates the bulk effusivity into κ and C — the experimental form of the lateral-wavenumber
    mechanism, for isotropic media.
13. **J. P. Feser, D. G. Cahill**, *Rev. Sci. Instrum.* **83**, 104901 (2012), doi:10.1063/1.4757863; **P. Jiang, X. Qian, R. Yang**,
    *Rev. Sci. Instrum.* **88**, 074901 (2017), doi:10.1063/1.4991715; **D. Rodin, S. K. Yee**, *Rev. Sci. Instrum.* **88**, 014902
    (2017), doi:10.1063/1.4973297 — abstracts. Beam offset / variable spot size separate k_r and k_z; thickness is an input in all three.
14. **E. Ziade**, *Rev. Sci. Instrum.* **91**, 124901 (2020), doi:10.1063/5.0021917 — abstract: 1 Hz–75 MHz FDTR; "thickness and
    thermal conductivity of gold/titanium thin films" measured simultaneously. The assumption that fixes the gauge there (known C,
    isotropy, or both) is not verified (full text not obtained).
15. **T. Chen, P. Jiang**, "Unraveling intrinsic relationship of thermal properties in thermoreflectance experiments", *Acta Phys.
    Sin.* **73**, 230202 (2024), doi:10.7498/aps.73.20241369 — abstract. SVD of the FDTR sensitivity matrix; the null-space vectors give
    "the most fundamental combination parameters that determine the thermoreflectance signal … the maximum number of parameters that can
    be experimentally extracted, and the parameters that must be known beforehand". Numerical local-rank method, Fourier, Al/sapphire.
16. **T. Chen, D. Wang, P. Jiang**, "Physics-derived natural coordinates for fast uncertainty propagation in multilayer
    thermoreflectance", *Int. J. Heat Mass Transf.* **273**, 129653 (2026; issue Feb. 2027), doi:10.1016/j.ijheatmasstransfer.2026.129653,
    arXiv:2610.09192 — read §2.4 and §3.4. SOURCE DERIVES (Hankel-transformed multilayer diffusion equation): a semi-infinite substrate
    enters through one composite parameter, with limits √(k_z k_r) (low frequency) and √(k_z C) (high frequency); a finite film enters
    through γh, "controlled by the in-plane thermal diffusivity k_r/C and the cross-plane diffusion time Ch²/k_z", amplitude set by the
    areal heat capacity Ch. Used as approximately Gaussian coordinates for uncertainty propagation; no exact invariance is stated and
    Krapez–Rigollet is not cited.
    **OUR INTERPRETATION (derivation, for B-S to check).** The film quadrupole depends on γh and k_zγ with γ² = (k_rκ² + Cp)/k_z
    (κ the Hankel wavenumber). Both are invariant under (h, k_z, k_r, C) → (μh, μk_z, k_r/μ, C/μ) at fixed κ (known beam radii); the
    invariants are exactly k_r/C, Ch and Ch²/k_z. Hence a finite beam breaks the scale gauge **only if the in-plane conductivity is known
    or tied to k_z** (isotropy forces μ² = 1). The same holds with a Cattaneo factor (1+τp) common to both directions. This makes the
    qualifier in claim B3 ("a lateral wavenumber … with known radii and in-plane properties") necessary.
17. **P. Jiang, S. Chen, X. Li, R. Yang**, "Thermoreflectance-based techniques for micro- and nanoscale thermophysical property
    measurements: principles, methods, and recent advances", arXiv:2609.27618 (2026), review — text search. SOURCE SAYS: signals are
    "governed mainly by certain parameter combinations, such as the cross-plane thermal effusivity, the in-plane thermal diffusivity
    k_r/C, and the areal heat capacity of a film Ch"; separating k and C "generally requires the film thickness and thermal penetration
    depth to fall within an appropriate relative range so that … Ch becomes an independently observable quantity". Fourier.

### 3.3 Non-Fourier forward models for a film on a substrate under front-face heating

18. **J. Ordóñez-Miranda, J. J. Alvarado-Gil**, "Determination of time-delay parameters in the dual-phase lagging heat conduction
    model", *J. Heat Transfer* **132**, 061302 (2010), doi:10.1115/1.4000748 — abstract only (paywalled, no open copy). SOURCE SAYS: a
    finite layer on a semi-infinite layer, modulated surface flux, DPL law; "when the thermal effusivities of the layers are quite
    different … a frequency range can be found in which the normalized amplitude and phase … show strong oscillations. This behavior is
    used to obtain explicit formulas for determining simultaneously the time parameters as well as additional thermal properties of the
    finite layer". **Closest prior art to B found**: same configuration and, by the 1D DPL ≡ GK identity (§3.4), the same response family.
    It is a constructive determination under a large-contrast assumption, not a rank or injectivity theorem; whether the thickness is
    assumed known could not be checked. Forward citations (25): thermoelastic damping, bioheat, IC thermal modelling; none on
    identifiability.
    Same programme (metadata only unless stated): *Int. J. Therm. Sci.* **48**, 2053–2062 (2009), doi:10.1016/j.ijthermalsci.2009.03.008
    (Cattaneo, relaxation-time determination; a web summary states a semi-infinite layer with relaxation time and effusivity —
    SECONDARY); *Int. J. Therm. Sci.* **49**, 209–217 (2010), doi:10.1016/j.ijthermalsci.2009.07.005 (effective properties of layered
    hyperbolic systems); *Int. J. Thermophys.* **31**, 900–925 (2010), doi:10.1007/s10765-010-0777-x (with interface resistance);
    *Nanoscale Res. Lett.* **6**, 327 (2011), doi:10.1186/1556-276X-6-327 (abstract: single- and DPL models "mathematically equivalent"
    under causality; the exact delayed DPL "cannot provide a general description"); *Int. J. Eng. Sci.* **50**, 101–112 (2012),
    doi:10.1016/j.ijengsci.2011.08.012 (hyperbolic transport, modulated source).
19. **J. Ordóñez-Miranda, R. Yang, S. Volz, J. J. Alvarado-Gil**, *J. Appl. Phys.* **118**, 075103 (2015), doi:10.1063/1.4928770,
    arXiv:1501.03174 — abstract and text search. BTE (RTA) film on substrate under modulated heating, "the BTE extension of
    Rosencwaig–Gersho"; response depends on L/MFP, ωτ and thermal mismatch; typical average relaxation times stated as τ < 10⁻¹⁰ s.
    No identifiability analysis.
20. **M. G. Hennessy, T. G. Myers**, "Guyer–Krumhansl heat conduction in thermoreflectance experiments", in *Multidisciplinary
    Mathematical Modelling* (SEMA SIMAI Springer Series, 2021), pp. 21–34, doi:10.1007/978-3-030-64272-3_2 — read in note 15. GK
    half-space, Gaussian heating; 1D limit φ = −π/4 + ½arctan(ωτ_R) − ½arctan(3ωℓ²/a) in the project convention; resonance condition
    (factor-3 correction recorded in note 15). Forward citations (5): none on identifiability.
21. **A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, *J. Appl. Phys.* **137**, 155103 (2025), doi:10.1063/5.0257299 —
    read (project reading note). Cattaneo in an FDTR forward model; sensitivity S = dφ/dlnτ only; no noise model, Fisher matrix,
    interface conductance or inversion. Forward citations (2): none on identifiability.
22. **Y. Qiu, Z. Zhong, Z. Qian, Y. Lu, Y. J. Hong, C. Liu, Q. Zheng**, "PRISM: unifying multiphysics optothermal transport and
    identifiability-guided inference for nanoscale thermal metrology", arXiv:2611.02475 (2 Nov 2026) — text search (§2.2, §4.2.4,
    §4.2.6). SOURCE SAYS: GK branch for the Si substrate with first-principles parameters ℓ = 185 nm, τ = 42 ps taken from Beardo et al.
    2020; noise-whitened SVD screening with a 0.1° baseline phase-noise floor over a 13-parameter Au/Si FDTR model; "holding the six
    unidentifiable TTM and GK parameters at their literature priors isolates k_Si and b_NF … as the two fit targets"; τ is not among the
    candidate fit parameters. Numerical, preprint; the scale gauge is not discussed.

### 3.4 Dual-phase lag, Jeffreys and GK: equivalence, Fourier resonance, parameter estimation, contrary results

23. **D. Y. Tzou**, *J. Heat Transfer* **117**, 8–16 (1995), doi:10.1115/1.2822329 — abstract: a universal DPL law covering diffusion,
    waves, phonon–electron interaction "and pure phonon scattering". *Int. J. Heat Mass Transf.* **38**, 3231–3240 (1995),
    doi:10.1016/0017-9310(95)00052-B, and *J. Thermophys. Heat Transf.* **9**, 686–693 (1995), doi:10.2514/3.725 — metadata.
    *Macro- to Microscale Heat Transfer: The Lagging Behavior* (Taylor & Francis, 1997; 2nd ed. Wiley, 2014, doi:10.1002/9781118818275)
    — not inspected. SECONDARY: **R. Baratifarimani, Z. Shomali**, arXiv:2307.00665 (2023) — read: the GK relaxation times "are related
    to two phase lags of DPL via τ_q = τ_R and τ_T = 9τ_N/5", citing Tzou 1997.
    OUR INTERPRETATION (derivation): in 1D, with the energy balance and no interior source, q_zz = −ρc T_tz, so the GK law
    τ_R q_t + q = −λT_z + L²q_zz becomes τ_R q_t + q = −λ(T_z + τ_ℓ T_zt), τ_ℓ = L²/a: the DPL law with τ_q = τ_R, τ_T = τ_ℓ, *including
    the flux boundary and interface conditions*. With an interior source s the GK flux acquires +L²s_z and the identity fails:
    **interior sources break the 1D GK–DPL equivalence**.
24. **R. Kovács, P. Ván**, "Thermodynamical consistency of the dual-phase-lag heat conduction equation", *Contin. Mech. Thermodyn.*
    **30**, 1223–1230 (2018), doi:10.1007/s00161-017-0610-x, arXiv:1709.06825 — read §3, §5. SOURCE DERIVES: GK (l₂ = 0) and the
    Jeffreys/DPL equation (l₃ = 0) are distinct special cases of one internal-variable equation, with flux laws containing ∂_xx q and
    ∂_x∂_t T respectively. SOURCE OBSERVES (re-evaluation of a published processed-meat heat-pulse experiment): when cooling is present
    "the Fourier fitting … is not appropriate and there is no such deviation in the experimental data"; the DPL fit produced an overshoot
    instead of a cooling tail. A contrary result to a DPL validation.
25. **D. D. Joseph, L. Preziosi**, "Heat waves", *Rev. Mod. Phys.* **61**, 41–73 (1989), doi:10.1103/RevModPhys.61.41; Addendum
    **62**, 375–391 (1990), doi:10.1103/RevModPhys.62.375 — metadata (Jeffreys-type flux, historical terminology). **M. E. Gurtin,
    A. C. Pipkin**, *Arch. Ration. Mech. Anal.* **31**, 113–126 (1968), doi:10.1007/BF00281373 — metadata (memory kernels).
26. **R. Quintanilla, R. Racke**, "Qualitative aspects in dual-phase-lag heat conduction", *Proc. R. Soc. A* **463**, 659–674 (2007),
    doi:10.1098/rspa.2006.1784 — abstract. SOURCE PROVES: well-posedness for positive lags; exponential stability when τ_T > τ_q.
    OUR INTERPRETATION: the Fourier resonance τ_q = τ_T is the boundary between the two regimes; the conserving closure of note 14 gives
    τ_ℓ < τ_R (note 7), the side not covered by the stability statement.
27. **Fourier resonance.** GK: **S. Both, B. Czél, T. Fülöp, Gy. Gróf, Á. Gyenis, R. Kovács, P. Ván, J. Verhás**, *J. Non-Equilib.
    Thermodyn.* **41**, 41–48 (2016, online 2015), doi:10.1515/jnet-2015-0035 (abstract: room-temperature heat-pulse data modelled by GK
    "in its over-diffusion regime"); **P. Ván et al.**, *EPL* **118**, 50005 (2017), doi:10.1209/0295-5075/118/50005; **R. Kovács**,
    *Int. J. Heat Mass Transf.* **127**, 631–636 (2018), doi:10.1016/j.ijheatmasstransfer.2018.06.082 (already attributed). DPL form
    τ_q = τ_T: OUR INTERPRETATION (the factor (1+τp) cancels in λ_eff); not located verbatim in a primary text inspected here.
28. **Parameter estimation with lags (thickness known throughout).** **H. R. B. Orlande, M. N. Özişik**, "Simultaneous estimation of
    thermal diffusivity and relaxation time with a hyperbolic heat conduction model", *Proc. IHTC-10* (1994), 403–408,
    doi:10.1615/IHTC10.840 — metadata (earliest estimation result found). **B. Mochnacki, M. Paruch**, *J. Appl. Math. Comput. Mech.*
    **12**, 97–102 (2013), doi:10.17512/jamcm.2013.4.12 — abstract (relaxation time of a thin metal film, evolutionary algorithm).
    **M. V. França, H. R. B. Orlande**, *Int. Commun. Heat Mass Transf.* **92**, 107–111 (2018),
    doi:10.1016/j.icheatmasstransfer.2018.02.001 — metadata (DPL parameters in MOSFETs). **M. Zubert, Z. Kulesza, M. Jankowski,
    A. Napieralski**, *Energies* **14**, 4425 (2021), doi:10.3390/en14154425 — abstract (SiO₂ temperature time-lag from scattering
    parameters with a 2ω method). **L. F. Escamilla-Herrera et al.**, *J. Appl. Phys.* **136**, 175105 (2024), doi:10.1063/5.0226434 —
    abstract (τ_T as a free fit parameter in pulsed photoacoustics). **A. Fehér, R. Kovács**, *Int. J. Eng. Sci.* **169**, 103577 (2021),
    doi:10.1016/j.ijengsci.2021.103577 — abstract (fast GK parameter estimation in flash experiments). **A. Fehér et al.**, *J. Non-Equilib.
    Thermodyn.* **46**, 403–411 (2021), doi:10.1515/jnet-2021-0033 — abstract (several thicknesses per rock type; size-dependent
    behaviour; empirical relation between Fourier and non-Fourier parameters). None addresses structural identifiability.
29. **Contrary results and failed validations.** **K. Mitra, S. Kumar, A. Vedavarz, M. K. Moallemi**, *J. Heat Transfer* **117**,
    568–573 (1995), doi:10.1115/1.2822615 — metadata (claimed τ in processed meat). **H. Herwig, K. Beckert**, *J. Heat Transfer* **122**,
    363–365 (2000), doi:10.1115/1.521471 — abstract: "experimental evidence for a perfectly Fourier-like behavior" in the same class of
    materials; *Heat Mass Transf.* **36**, 387–392 (2000), doi:10.1007/s002310000081 — metadata. **D. Maillet**, "A review of the models
    using the Cattaneo and Vernotte hyperbolic heat equation and their experimental validation", *Int. J. Therm. Sci.* **139**, 424–432
    (2019), doi:10.1016/j.ijthermalsci.2019.02.021 — abstract only: the cited validations "were biased, because of a deficient
    methodological approach based on too simple assumptions or on poor data reduction techniques".

### 3.5 Mathematically and physically equivalent formulations

30. **Memory-kernel inverse problems.** **G. Gentili**, *Inverse Problems* **7**, 77–84 (1991), doi:10.1088/0266-5611/7/1/008 —
    abstract. SOURCE PROVES: uniqueness of the heat-capacity and conductivity memory functions of a rigid homogeneous conductor from
    boundary flux and temperature ("overspecified data"), by Laplace transform and thermodynamic restrictions; explicit for a half-space.
    **J. Janno, L. von Wolfersdorf**, *Inverse Problems* **15**, 1455–1467 (1999), doi:10.1088/0266-5611/15/6/305; *Z. Anal. Anwend.*
    **19**, 489–510 (2000), doi:10.4171/ZAA/963 — abstracts (existence and uniqueness of memory kernels by Laplace transform).
    **E. Pais**, *Math. Model. Anal.* **15**, 473–490 (2010), doi:10.3846/1392-6292.2010.15.473-490 — abstract (energy and flux kernels
    from flux at the ends of a bar). **S. A. Avdonin, S. A. Ivanov, J.-M. Wang**, *Inverse Probl. Imaging* **13**, 31–38 (2019),
    doi:10.3934/ipi.2019002 — abstract (1D Gurtin–Pipkin from the Dirichlet-to-Neumann map; explicit formula; local Borg–Marchenko-type
    uniqueness). All assume a known domain: uniqueness for kernels, not for the (thickness, kernel) pair of a film on a substrate.
31. **Energy-dependent Sturm–Liouville problems.** **R. Hryniv, N. Pronska**, *Inverse Problems* **28**, 085008 (2012),
    doi:10.1088/0266-5611/28/8/085008 — abstract. OUR INTERPRETATION: a non-Fourier layer is a Sturm–Liouville problem with a
    p-dependent coefficient λ_eff(p); natural class for a graded non-Fourier extension (not needed for B1).
32. **Lossy transmission lines (telegraph ≡ Cattaneo).** **M. Jaulent**, "The inverse scattering problem for LCRG transmission lines",
    *J. Math. Phys.* **23**, 2286–2290 (1982), doi:10.1063/1.525307 — abstract. SOURCE SAYS: from the data one can construct only two
    functions *of the travel time*, built from (1/4)d/dx ln(L/C) ± (1/2)(R/L − G/C) with an exponential weight. OUR INTERPRETATION: a
    uniform Cattaneo layer is an RLC line with R′ = 1/λ, L′ = τ/λ, C′ = ρc per unit length; the lumped totals R′e, L′e, C′e map exactly
    onto (ξ₁, b, τ_R): ξ₁ = √(R′e·C′e), b = √(C′e/R′e), τ_R = L′e/R′e. The scale gauge is the statement that terminal data of a uniform
    line determine lumped totals, not the length.
33. **Lossy seismology.** **B. J. Chaderjian**, *SIAM J. Appl. Math.* **54**, 1224–1249 (1994), doi:10.1137/S0036139992241157 —
    abstract: uniqueness of density and attenuation from a surface point-source response "if the wave speed is fixed". Analogue of the
    known-velocity mechanism of B3.
34. **Structural identifiability of PDE and impedance models.** **H. M. Byrne, H. A. Harrington, A. Ovchinnikov, G. Pogudin,
    H. Rahkooy, P. Soto**, *Nonlinearity* **38**, 025022 (2025), doi:10.1088/1361-6544/ada510 — abstract (differential-algebra theory
    for polynomial PDE models). **A. P. Browning, M. Taşcă, C. Falcó, R. E. Baker**, *Proc. R. Soc. A* **480**, 20230911 (2024),
    doi:10.1098/rspa.2023.0911 — abstract (linear reaction–advection–diffusion PDEs). Both treat polynomial input–output structure; the
    front-face impedance is transcendental in p, outside their setting. **F. Berthier, J.-P. Diard, L. Pronzato, E. Walter**,
    "Identifiability and distinguishability concepts in electrochemistry", *Automatica* **32**, 973–984 (1996),
    doi:10.1016/0005-1098(96)00031-3 — metadata; **F. Berthier, J.-P. Diard, R. Michel**, *J. Electroanal. Chem.* **510**, 1–11 (2001),
    doi:10.1016/S0022-0728(01)00554-X — metadata. **J. Bisquert**, *J. Phys. Chem. B* **106**, 325–333 (2002), doi:10.1021/jp011941g —
    abstract (finite-diffusion impedances, transmission-line representation; OUR INTERPRETATION: the finite-length diffusion impedance
    depends on L and D only through L²/D and a resistance — the same gauge). **Y. S. Ju**, *PNAS Nexus* **5**, pgag281 (2026),
    doi:10.1093/pnasnexus/pgag281 — abstract; state-dependent constitutive functions from sparse data; not applicable.

### 3.6 Gauge-breaking mechanisms: status in the literature

| Mechanism | Literature status | Sources |
|---|---|---|
| Thickness known | ESTABLISHED (trivial; standard practice) | Parker 1961; Krapez–Rigollet 2017; all estimation papers of §3.4 |
| λ or ρc known | ESTABLISHED (algebraic: b and ξ₁ then give the rest) | Krapez–Rigollet 2017 ("independent measurement of a, ρC or κ") |
| Lateral diffusion, known beam radii, isotropic medium | ESTABLISHED (Fourier; used experimentally to separate κ and C) | Schmidt 2009; Liu et al. 2013; Wei et al. 2013; suggested by Krapez–Rigollet 2017 |
| Lateral diffusion, anisotropic film, k_r free | Gauge persists | OUR INTERPRETATION from the Chen–Wang–Jiang 2026 coordinates (k_r/C, Ch, Ch²/k_z) |
| Lateral diffusion under isotropic GK | Broken (shear-mode length ℓ against 1/κ) | OUR INTERPRETATION from the Hennessy–Myers 2021 decay rates (note 15) |
| Interior sources | PARTIAL: suggested, not proved; the source profile must be known in absolute length (an unknown absorption coefficient β scales as 1/μ and preserves the gauge — OUR INTERPRETATION); also breaks GK ≡ DPL in 1D | Krapez–Rigollet 2017; Krapez 2023; Pailhès et al., *Int. J. Therm. Sci.* **53**, 49–55 (2012), doi:10.1016/j.ijthermalsci.2011.10.005; Krapez–Dohou, *Int. J. Therm. Sci.* **81**, 38–51 (2014), doi:10.1016/j.ijthermalsci.2014.02.007; Fuente, Apiñaniz, Mendioroz, Salazar, *J. Appl. Phys.* **110**, 033515 (2011), doi:10.1063/1.3614524 (diffusivity and absorption coefficient of slabs of known thickness) |
| Two-sided (rear-face) detection | Does not break the homogeneous gauge (quadrupole entries depend on ξ₁, b only); fixes the end of the orbit only with known thickness | Maillet et al. 2000; Parker 1961; note 06 (graded case) |
| Independent velocity with a kinetic/Debye relation | PARTIAL: project argument (note 7); analogues in seismology and transmission lines | Chaderjian 1994; Jaulent 1982 |
| Interface resistance film/substrate (extra parameter) | OPEN for the non-Fourier rank; Fourier: film and interface resistances classically separated by thickness series | S.-M. Lee, D. G. Cahill, *J. Appl. Phys.* **81**, 2590–2595 (1997), doi:10.1063/1.363923 (metadata) |

### 3.7 Verdict for B (scope stated)

- **Generic rank four / injectivity of (ξ₁, b, τ_R or τ_q, τ_ℓ or τ_T) for a film on a known Fourier substrate: no proof found.**
  Scope: 66 database queries in three indices (direct, mathematical-equivalent, physical-equivalent and historical terms), forward
  citations of 11 anchor papers, backward references of 6, the Crossref publication list of J.-C. Krapez, open-access full texts of the
  2026 Fourier-coordinate papers, and neighbouring theorem families (memory kernels, LCRG lines, lossy seismology, algebraic PDE
  identifiability). The closest work, Ordóñez-Miranda & Alvarado-Gil 2010 (DPL, finite layer on semi-infinite substrate, explicit
  determination formulas under strong effusivity contrast), was read at abstract level only; it is the main residual prior-art risk for
  B2 and must be read before any novelty statement.
- **Scale gauge under Cattaneo/GK:** not found outside the project (note 7). Implied for Fourier anisotropic multilayers by
  Chen–Wang–Jiang 2026 and, in the travel-time sense, for telegraph lines by Jaulent 1982.
- **Degenerate set Σ:** Fourier resonance (τ_R = τ_ℓ, i.e. τ_q = τ_T) as a reduction to Fourier is ESTABLISHED (Both 2016 for GK;
  elementary for DPL). Matched effusivity is ESTABLISHED as degenerate under Fourier (Krapez–Rigollet 2017). OUR INTERPRETATION: under
  the non-Fourier laws Γ(p) = (b√((1+τ_ℓp)/(1+τ_Rp)) − b_s)/(b√((1+τ_ℓp)/(1+τ_Rp)) + b_s) vanishes identically only if b = b_s **and**
  τ_R = τ_ℓ, so b = b_s alone is a candidate for poor conditioning rather than structural degeneracy; for B-S/B-X to test.

---

## 4. Workstream A — FDTR sensitivity to 10–100 ps relaxation times at MHz frequencies

### 4.1 Bandwidth, noise and instrument nuisances

35. **K. T. Regner, S. Majumdar, J. A. Malen**, *Rev. Sci. Instrum.* **84**, 064901 (2013), doi:10.1063/1.4808055 — read §II and
    Table II. SOURCE SAYS: heterodyne BB-FDTR 200 kHz–200 MHz; conventional FDTR loses signal-to-noise above ≈20 MHz (coherent and
    ambient noise); thermal phase isolated by subtracting a pump-phase measurement (φ_thermal = φ_probe,m − φ_pump,m); uncertainty budget
    of the example: transducer thickness ±5 %, transducer conductivity ±5 %, spot radius 3.2 ± 0.1 µm, G = 210 ± 10 MW m⁻²K⁻¹, and a DC
    phase shift ±0.2 applied to all frequencies (unit not legible in the extracted table; degrees presumed).
36. **D. J. Kirsch, J. Martin, R. Warzoha, M. McLean, D. Windover, I. Takeuchi**, "An instrumentation guide to measuring thermal
    conductivity using frequency domain thermoreflectance (FDTR)", *Rev. Sci. Instrum.* **95**, 103006 (2024), doi:10.1063/5.0213738 —
    abstract; companion conference presentation by the same group (NIST document page "ITCC FDTR instrumentation: common challenges and
    practical guidance") read by text extraction. SOURCE SAYS (presentation): the frequency-averaged phase standard deviation is used as the
    phase standard uncertainty "(± 0.1)" in Monte Carlo; the pump reference path is phase-matched to the pump leak "to within ± 0.1 at
    30 MHz" (degree sign lost in extraction); a standard signal cable gave "phase deviations upwards of 1 to 2" degrees when moved or
    touched; above 10 MHz noise is dominated by coherent RF pickup; sweeps 100 Hz–20 MHz with a 50 MHz lock-in.
    Camacho 2025 cites this paper for "a resolution of a few millidegrees"; that figure was not found in the presentation and the paper
    text could not be read.
    OUR INTERPRETATION: a ±0.1° reference-phase match at 30 MHz is a delay uncertainty of 0.1/(360 × 3×10⁷) s ≈ 9 ps, inside the
    10–100 ps range of the memory times of this campaign.
37. **J. Yang, E. Ziade, A. J. Schmidt**, "Uncertainty analysis of thermoreflectance measurements", *Rev. Sci. Instrum.* **87**,
    014901 (2016), doi:10.1063/1.4939671 — abstract (authors verified): precision formula for multi-parameter least squares including
    noise and the uncertainty of fixed parameters; FDTR vs TDTR. Phase-noise figures not obtained.
38. **M. Shahzadeh, M. Rahman, O. Hellwig, S. Pisana**, *Rev. Sci. Instrum.* **89**, 084905 (2018), doi:10.1063/1.5037117 —
    abstract: differential BB-FDTR; varying the spot size or offsetting the beams "to recover the thermal phase of interest and to separate
    it from the unwanted instrumental contributions". A design remedy for nuisance set N1.
39. **X. Wang, M. Jeong, A. J. H. McGaughey, J. A. Malen**, *Rev. Sci. Instrum.* **93**, 023001 (2022), doi:10.1063/5.0080119 —
    abstract (spot-radius uncertainty reduced by a reference-sample phase quotient). **J. Drew, S. Godse, Y. Liang, A. Pathak,
    J. A. Malen, R. C. Kurchin**, *Rev. Sci. Instrum.* **97**, 034903 (2026), doi:10.1063/5.0316853 — abstract (Bayesian estimation for
    FDTR; priors can shift inferred values). **Z. Wang, B. Linghu, Z. Li, J. Zhu, D. Tang**, *J. Appl. Phys.* **140**, 175302 (2026),
    doi:10.1063/5.0354480 — abstract (TDTR Al/GaN: interface conductance and substrate conductivity "not fully separable").
    **T. Liu, F. Chen, J. Li, C. Tong**, *Adv. Sci. Instrum.* **1**, 100037 (2026), doi:10.1016/j.asi.2026.100037 — abstract (FDTR review
    centred on identifiability). **J. Tu, M. A. Haque, D. Baran, W.-L. Ong**, *Fundam. Res.* **5**, 288–295 (2025),
    doi:10.1016/j.fmre.2023.01.010 — abstract (log-sensitivity-ratio criterion). All Fourier.
40. **F. Angeles et al.**, "Picosecond magneto-optic thermometry measurements of nanoscale thermal transport in AlN thin films",
    *APL Mater.* **11**, 061127 (2023), doi:10.1063/5.0149651 — read (abstract, §I–II). SOURCE SAYS: transducer response time τ ≈ hC/G;
    for h = 80 nm, C = 2.5 MJ m⁻³K⁻¹, G = 100 MW m⁻²K⁻¹ "the thermal response time is typically of the order of ≈2 ns. As a result, for a
    material with high thermal diffusivity like AlN … spatial resolution of TDTR and FDTR is limited to length scales greater than
    L ≈ 500 nm"; observed 1600 ps (80 nm Al), 200 ps (TiN), 38 ps (PtCo/TiN); AlN films of 100, 400, 1000 nm: Λ = 200 ± 80, 165 ± 35,
    300 ± 70 W m⁻¹K⁻¹. The Au/AlN/sapphire scenario of A has R = 10⁻⁸ m²K/W, hence the same ≈2 ns transducer–interface time constant,
    20–200 times the memory times at stake.

### 4.2 Frequency-dependent apparent conductivity at MHz: spatial nonlocality, mean free time, or interface

41. **Y. K. Koh, D. G. Cahill**, *Phys. Rev. B* **76**, 075207 (2007), doi:10.1103/PhysRevB.76.075207 — abstract. TDTR 0.1–10 MHz on
    InGaP, InGaAs, SiGe; reduction "consistent with a model … that phonons with mean free paths larger than the thermal penetration depth
    do not contribute". Spatial interpretation, not ωτ.
42. **K. T. Regner, D. P. Sellan, Z. Su, C. H. Amon, A. J. H. McGaughey, J. A. Malen**, *Nat. Commun.* **4**, 1640 (2013),
    doi:10.1038/ncomms2630 — read (introduction, methods). Si k decreasing over 200 kHz–200 MHz, read as an MFP accumulation function with
    L_c = L_p. SOURCE SAYS: the transducer electron relaxation time (3×10⁻¹⁴ s) is "minimally five orders of magnitude" below the
    modulation period, so G is taken frequency independent.
43. **J. P. Freedman, J. H. Leach, E. A. Preble, Z. Sitar, R. F. Davis, J. A. Malen**, *Sci. Rep.* **3**, 2963 (2013),
    doi:10.1038/srep02963 — abstract. BB-FDTR MFP spectra of GaAs, GaN, **AlN**, 4H-SiC (80–400 K): phonons with MFP > 230 ± 120,
    1000 ± 200, **2500 ± 800**, 4200 ± 850 nm carry 50 % of bulk k near room temperature. For a 500 nm AlN film most of the heat sits in
    modes whose bulk MFP exceeds the thickness (consistent with the film regularisation of assumption A1).
44. **R. B. Wilson, D. G. Cahill**, "Anisotropic failure of Fourier theory in time-domain thermoreflectance experiments", *Nat. Commun.*
    **5**, 5075 (2014), doi:10.1038/ncomms6075 — read (abstract, results by text search). Venue note: a Nature Communications paper; the
    only Wilson–Cahill PRL in Crossref is the 2012 Wiedemann–Franz paper (doi:10.1103/PhysRevLett.108.255901). SOURCE SAYS: "at room
    temperature, Fourier theory accurately predicts the thermal response of Si and Si:B up to f ≈ 17.6 MHz" although phonons with MFP
    beyond the penetration depth carry ≈50 W m⁻¹K⁻¹; failure is anisotropic; where Fourier fails, "a simple radiative boundary condition …
    cannot adequately describe interfacial thermal transport".
45. **Y. K. Koh, D. G. Cahill, B. Sun**, *Phys. Rev. B* **90**, 205412 (2014), doi:10.1103/PhysRevB.90.205412 — abstract: nonlocal BTE
    under periodic surface heating; for pure crystals the effective k from the out-of-phase signal is frequency independent; the
    "ballistic phonons do not contribute" picture applies to alloys only.
46. **F. Yang, C. Dames**, *Phys. Rev. B* **91**, 165311 (2015), doi:10.1103/PhysRevB.91.165311 — abstract. Two-flux BTE for periodic
    planar heating: "the natural variable is mean free time rather than mean free path", the cutoff criterion "involves the heater
    frequency rather than thermal diffusion length". OUR INTERPRETATION: the apparent conductivity at angular frequency ω is governed by
    modes with τ_i ≳ 1/ω (ns at 10–200 MHz), the long-time tail of the kernel; K(ω) ≈ K(0)(1 − iωτ_mem) holds only if ωτ_i ≪ 1 for all
    modes carrying weight. This links A to the infrared behaviour of τ_mem noted in the project (bulk τ_mem not mesh-converged).
47. **C. Hua, A. J. Minnich**, *Phys. Rev. B* **89**, 094302 (2014), doi:10.1103/PhysRevB.89.094302 — abstract (frequency-dependent BTE,
    transient grating). **J. Maassen, M. Lundstrom**, *J. Appl. Phys.* **117**, 135102 (2015), doi:10.1063/1.4916245, and **119**, 095102
    (2016), doi:10.1063/1.4942836, arXiv:1508.03864 — abstracts and text search: the two-flux BTE "can be recast exactly as the Cattaneo
    and hyperbolic heat equations"; "the key to capturing ballistic heat effects is to use the correct physical boundary conditions".
    **A. T. Ramu, J. E. Bowers**, *Appl. Phys. Lett.* **106**, 263102 (2015), doi:10.1063/1.4923310 — abstract: the BB-FDTR Si data "can
    adequately be explained" by a two-channel minimum-parameter model, "highlighting an important limitation of FDTR experiments".
    **T. Li, B. Jiang, Z. Chen**, *Newton* **1**, 100092 (2025), doi:10.1016/j.newton.2025.100092 — abstract (complex thermal conductivity
    for BTE-RTA under periodic heating; nonlocality in a jump boundary condition).
48. **Y. R. Koh et al.**, *Appl. Phys. Lett.* **109**, 243107 (2016), doi:10.1063/1.4972186 — abstract (Al₀.₁Ga₀.₉N films: apparent k
    drops 30 % from 0.8 to 10 MHz; "limitations in concurrent extraction of other model parameters"). **A. Zenji, G. Pernot, D. Lacroix,
    J.-M. Rampnoux, O. Bourgeois, S. Grauby, S. Dilhaire**, *Commun. Mater.* **5**, 123 (2024), doi:10.1038/s43246-024-00572-7 — abstract (TDTR + FDTR,
    few kHz to 100 GHz, InGaAs, truncated Lévy).
49. **Y. Liu, Q. Li, F. Liu, X. Wang, B. Sun**, "Boundary conditions dictate frequency dependence of thermal conductivity in silicon",
    *Appl. Phys. Rev.* **12**, 031421 (2025), doi:10.1063/5.0254248, arXiv:2412.02227 — abstract. SOURCE OBSERVES: frequency dependence
    in Al/Si with an atomically sharp interface, absent in Al/SiO₂/Si at 300 K, re-emerging at 80 K. The MHz non-Fourier signature is
    interface-controlled: nuisance set N2 is also a model-form issue.

### 4.3 Relaxation times and FDTR: what has been claimed, computed or measured

50. **A. Beardo, M. G. Hennessy, L. Sendra, J. Camacho, T. G. Myers, J. Bafaluy, F. X. Alvarez**, *Phys. Rev. B* **101**, 075303 (2020),
    doi:10.1103/PhysRevB.101.075303 — read (Table I, App. F; local copy). SOURCE SAYS: ab initio KCM parameters for Si, τ = 1002, 148,
    **42**, 27 ps and ℓ = 3127, 1047, **185**, 142 nm at 81, 154, **311**, 416 K (not fitted); at 311 K "ωτ ≪ 1 and thus hydrodynamic
    resistive effects are negligible across all frequencies"; KCM and Fourier diverge from 5 MHz, "dominated by nonlocal effects".
    OUR INTERPRETATION: in the 1D variables of note 7 these values give τ_ℓ = (1+α)ℓ²/a ≈ 0.5 ns (α = 1/3, a = κ/c = 8.9×10⁻⁵ m²/s),
    about twelve times τ_R; the FDTR-visible time in Si is the nonlocal one, not the relaxation time.
51. **S. Saurav, S. Mazumder**, "Assessment of models for extracting thermal conductivity from frequency domain thermoreflectance
    experiments", ASME 2020 Heat Transfer Summer Conference, V001T15A001 (2020), doi:10.1115/HT2020-8920 — read (results and conclusions of
    the deposited manuscript); *J. Heat Transfer* **144**, 013501 (2021), doi:10.1115/1.4052655 — abstract. SOURCE OBSERVES (numerical):
    Au/Si FDTR 200 kHz–200 MHz, effective relaxation time of Si estimated as ≈0.0365 ns; "the phase lag calculated using the two models
    [Fourier and hyperbolic] was found to be identical"; largest temperature difference 0.05 K of ≈2.1 K at 200 MHz. No noise model.
52. **Camacho de la Rosa et al. 2025** (item 21). SOURCE SAYS: "relaxation times smaller than 10⁻⁹ s are practically identical and
    represent the Fourier limit of the system" in their band; at 3×10⁻⁹ s the differences lie "outside of the frequency range of common FDTR
    setups". No Fisher analysis.
53. **A. Beardo, M. López-Suárez, L. A. Pérez, L. Sendra, M. I. Alonso, C. Melis, J. Bafaluy, J. Camacho, L. Colombo, R. Rurali,
    F. X. Alvarez, J. S. Reparaz**, "Observation of second sound in a rapidly varying temperature field in Ge", *Sci. Adv.* **7**,
    eabg4677 (2021), doi:10.1126/sciadv.abg4677, arXiv:2007.05487 — read (main text). SOURCE OBSERVES: bare natural Ge (optical
    penetration 15–17 nm), pump modulated 30 kHz–200 MHz; 3D hyperbolic fit gives τ = 500 ps, α = 3×10⁻⁵ m²/s, v ≈ 250 m/s at room
    temperature; "wave-like effects are not observed in the presence of a metallic transducer, since the thermal interface … as well as the
    transducer itself, dominate the system response in the frequency range where wave-like effects are expected". **DISPUTED** by
    **S. Huberman, C. Zhang, J. Abou Haibeh**, *J. Phys.: Condens. Matter* **37**, 435601 (2025), doi:10.1088/1361-648X/ae0dd6,
    arXiv:2206.02769 — abstract: LBTE solutions show no second sound; Hardy's driftless condition fails in Ge; RTA simulations of the 1D
    FDTR geometry show that "phase lag alone is not a suitable experimental observable for inferring second sound".
54. **PRISM 2026** (item 22): GK and two-temperature parameters unidentifiable in Au/Si FDTR at a 0.1° noise floor; τ = 42 ps fixed.
55. **Y. Ezzahri, K. Joulain**, "Dynamical thermal conductivity of bulk semiconductor crystals", *J. Appl. Phys.* **112**, 083515
    (2012), doi:10.1063/1.4759366 — abstract: κ(Ω) from the BTE (single relaxation time, modified Debye–Callaway), causal, of Gurtin–Pipkin
    convolution type; cut-off frequency sensitive to scattering parameters. Values not obtained (open copy inaccessible).
    **D. E. Crawford, Y. Zeng, J. Vidal, J. Dong**, *Phys. Rev. B* **112**, 115201 (2025), doi:10.1103/p8wg-p1j3 — abstract: time-domain
    transport function from flux correlations (Zwanzig), aimed at transient-grating and EUV/x-ray experiments; no FDTR detectability
    estimate.
56. Ultrafast DPL lags in metals: **D. Y. Tzou**, *J. Thermophys. Heat Transf.* **9**, 686–693 (1995), doi:10.2514/3.725 — metadata.
    The ps-scale lags of that literature come from femtosecond pump–probe data on metal films (electron–phonon regime), not from MHz FDTR
    on dielectrics. Cited only to prevent conflation.
57. Room-temperature second sound by transient thermal grating: **Z. Xie et al.**, *Nat. Commun.* **17**, 4806 (2026),
    doi:10.1038/s41467-026-70807-3 — abstract (isotopically purified graphite). Not FDTR; the observable for an unambiguous claim is a
    damped oscillation in a grating experiment.

### 4.4 AlN-specific measurements (context for the 500 nm scenario)

**Z. Cheng et al.**, *Phys. Rev. Mater.* **4**, 044602 (2020), doi:10.1103/PhysRevMaterials.4.044602 — abstract (intrinsic κ of thick
MOCVD AlN, 130–480 K). **Y. R. Koh et al.**, *ACS Appl. Mater. Interfaces* (2020), doi:10.1021/acsami.0c03978 — abstract (3–22 µm
films, κ ≈ 320 ± 42 W m⁻¹K⁻¹ at room temperature, TDTR and SSTR). **Angeles et al. 2023** (item 40). **M. S. B. Hoque et al.**,
*Appl. Phys. Lett.* **125**, 262201 (2024), doi:10.1063/5.0239769 — abstract (ballistic-to-diffusive transition in 1.6–2440 nm AlN
films on sapphire). **Freedman et al. 2013** (item 43). GaN context: **E. Ziade et al.**, *Appl. Phys. Lett.* **110**, 031903 (2017),
doi:10.1063/1.4974321 — metadata.

### 4.5 Verdict for A (scope stated)

- **No published estimate was found of the FDTR resolution of a 10–100 ps Cattaneo (memory) time with nuisance parameters freed**
  (instrument phase, delay, interface resistances, transducer, beam radii). Scope: queries A01–A18, S04, S07, S08, P09, P10; forward
  citations of Beardo 2021, Saurav–Mazumder 2021, Camacho 2025, Yang–Ziade–Schmidt 2016; full texts of Beardo 2020/2021,
  Saurav–Mazumder 2020, Camacho 2025, PRISM 2026, Angeles 2023, Regner 2013. The project's half-space single-parameter estimate (note 7:
  2.78 % at 0.01°, τ = 10 ps, all else known) has no published counterpart.
- What the literature does establish: at room temperature and ≤200 MHz, ωτ ≪ 1 for relaxation times of tens of ps, and published models
  report a negligible Cattaneo/resistive contribution to the phase (Beardo 2020 statement; Saurav–Mazumder numerical; Camacho threshold);
  observed MHz non-Fourier signatures are attributed to spatial nonlocality, mean-free-time tails or interface non-equilibrium, with
  disputed interpretation; metal transducers impose a ≈2 ns response time (Angeles 2023) and mask wave-like effects (Beardo 2021);
  published phase uncertainties are ≈0.1–0.2° random, with systematic reference and delay errors of the same order or larger.
- **OUR INTERPRETATION, for AB-A and A-C to check (not literature).** For the 1D half-space, φ = −π/4 + ½arctan(ωτ_R) − ½arctan(ωτ_ℓ)
  gives δφ ≈ (ω/2)(τ_R − τ_ℓ) at first order and δln|Z| ≈ (ω²/4)(τ_R² − τ_ℓ²) at second order; an instrument delay contributes −ωt_d.
  At first order in ωτ the three are collinear, so with N1 ∪ {τ_ℓ} free the information on τ_R moves to second order in ωτ (≈10⁻³ for
  50 ps at 200 MHz). Film/substrate structure and the ≈2 ns transducer time constant modify this; quantifying it is the task of the twin.

---

## 5. Failed historical approaches and contrary results relevant to both workstreams

1. Identifiability maps from determinants of finite-difference sensitivity matrices at two or three frequencies (Guo et al. 2017) —
   refuted structurally (Krapez–Rigollet 2017): rank deficiency independent of the frequency set; non-zero determinants were round-off.
2. Macroscopic "experimental validations" of Cattaneo/DPL times (Mitra et al. 1995 and successors) — contradicted by Herwig–Beckert 2000,
   re-evaluated by Kovács–Ván 2018, judged methodologically biased in Maillet's 2019 review.
3. Fourier fits with frequency-dependent conductivity cannot reproduce phase and amplitude together in Si FDTR (Beardo 2020); frequency
   dependence in Si appears or disappears with the interface structure (Liu et al. 2025); Si BB-FDTR data admit a two-channel explanation
   (Ramu–Bowers 2015).
4. Phase lag alone as evidence of second sound in Ge (Beardo 2021) — disputed (Huberman et al. 2025).
5. Sensitivity magnitude used as evidence of fittability (Camacho 2025) — not an identifiability argument (Krapez–Rigollet 2017; PRISM
   2026: "sensitivity magnitudes alone do not establish joint identifiability").

---

## 6. Prior-art position of the campaign claims (scope as in §3.7 and §4.5)

| Claim | Position |
|---|---|
| B3 gauge invariance under Cattaneo/GK in 1D | Potentially novel within the literature searched (project note 7); Fourier case and Fourier multilayer coordinates known |
| B1 generic rank four with explicit Σ | No proof found; closest constructive work Ordóñez-Miranda & Alvarado-Gil 2010 (abstract only) |
| B2 injectivity off Σ | No proof found; same residual risk; neighbouring uniqueness theorems require a known domain |
| Lateral-wavenumber breaking requires known or tied k_r | Derivable from Chen–Wang–Jiang 2026; not stated there |
| Interior sources break the 1D GK–DPL equivalence | Not found stated; elementary derivation (§3.4 item 23) |
| A: FDTR resolution of ps memory times with nuisances | No published Fisher/Cramér–Rao estimate found |
| A: delay–τ collinearity at first order | Not found stated; elementary (§4.5) |

---

## 7. Unresolved prior-art uncertainty

1. **Ordóñez-Miranda & Alvarado-Gil 2010, *J. Heat Transfer* 132, 061302** (and the 2009/2010/2012 companions) were not read beyond the
   abstract. They may contain a constructive inversion of (time parameters, layer properties) that partially anticipates B2. Obtain the
   full texts before writing any novelty sentence about B.
2. Maillet 2019 review — abstract only; it may discuss τ–diffusivity correlations in estimation.
3. Ziade 2020 — the assumption that fixes the gauge in the thickness measurement (known C, isotropy) is unverified.
4. Kirsch et al. 2024 — the "few millidegrees" figure quoted by Camacho 2025 could not be checked against the paper; the presentation gives
   ±0.1° as the standard phase uncertainty.
5. Tzou's monograph and 1995 papers — the DPL–GK 9/5 identification is quoted secondarily (arXiv:2307.00665).
6. Semantic Scholar keyword searches were rate-limited (2/66); citation-graph work relied on OpenAlex. Russian-language literature on
   telegraph-type inverse problems (Durdiev, Kozhanov, Sabitov) was seen at title level only.
7. Engineering literature on extracting per-unit-length line parameters with an unknown line length was not searched in depth (P04 returned
   power-system results).

**Alternate names to request from the mathematician and the theoretical physicist** (to repeat phases 1–3): "Weyl–Titchmarsh function of a
string with frequency-dependent density", "impedance of a Krein string", "Jeffreys-type thermal transmission line", "inverse problem for
the damped wave equation from one-port reflection", "acoustic impedance of a standard-linear-solid (Zener) layer on a half-space" — the
Zener layer has the same rational-in-p structure as GK and an analogous thickness–modulus scaling in one-port impedance measurements.
