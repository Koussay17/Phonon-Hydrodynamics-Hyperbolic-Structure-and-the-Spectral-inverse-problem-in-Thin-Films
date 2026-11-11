# Novelty protocol — film-regularised memory time: identifiability, multi-mesh bounding, mesh convergence (claims C1–C3)

Campaign `20261005-fc3-memory`, review stage, 14 October 2026.
Search records: `review/novelty-searches-film/` (multi-index query plan `queries.txt`, runner `run_queries.sh`, log
`batch-log.txt`, raw JSON/MD per claim in `C1/`–`C3/`, all returned titles `titles_all.txt`; arXiv query plan
`arxiv_queries.txt`, runner `arxiv_batch.py`, results `arxiv_results.md`, raw feeds `arxiv_raw/`; abstracts and metadata
retrieved by DOI in `abstracts/` (runner `openalex_abstracts.py`); keyword-filtered forward citations
`cites_batch1_out.md`, `cites_batch2_out.md` (runner `openalex_cites.py`); DOI verification `cited-dois.txt`,
`crossref-verification.txt`; arXiv-only records `arxiv-ids-cited.txt`). Hand-run searches: Appendix A.

## 0. Scope and conventions

**Objects.** tau_mem = N/K, K = b^T C^+ b, N = |C^+ b|^2, with C the linearised conserving three-phonon operator in
entropy coordinates and b the heat-current source vector (as in `novelty-protocol.md` §0). Film regularisation:
C_L = C + diag(|v_z|/L) (Matthiessen boundary rate of a c-oriented film of thickness L), cross-plane component
tau_z(L); at RTA level tau_mem(L) = sum c v_z^2 tau_L^2 / sum c v_z^2 tau_L with 1/tau_L = 1/tau + |v_z|/L (E26, E30).

**Claims audited** (numbers from `07-experiments.md`; AlN, 300 K):
* **C1 (identifiability).** In the force-constant-generated rate class (symmetric, sum-rule-satisfying fc3 of the 3x3x2
  supercell; per-cluster relative bound eps = 0.1, or per-cluster bounds calibrated on the difference between two
  independent first-principles workflows), static data — all mode lifetimes within ±10 %, kappa_xx and kappa_zz within
  ±1 % — leave the bulk tau_mem undetermined (admissible ratios of order 0.7 to 1.7 in these classes on the discretised
  operators, larger in wider classes; E15–E22), whereas the film value tau_z(L), L = 30–300 nm, is confined to about
  ±10 % with data on one q-mesh (100 nm, per-cluster: [0.909, 1.098]; calibrated: [0.802, 1.247]), ±5 % with data on two
  (100 nm: [0.950, 1.046]; calibrated: [0.897, 1.079]) and at most +3.8 % with data on three (100 nm maximum; the
  minimum run was still in progress) (E31, E33, E36). Reading: data imposed on several q-meshes progressively pin the
  film quantity.
* **C2 (method).** Bounding a derived transport quantity by constrained optimisation over the fc3 tensor: augmented
  Lagrangian of Powell–Hestenes–Rockafellar type with L-BFGS-B inner solves, data constraints imposed on every mesh of a
  set, objective a smooth minimum (log-sum-exp, beta = 40 by default) of signed log ratios across meshes, only iterates
  feasible for the nominal tolerances reported (inner bounds), and validation on finer meshes not used in the
  optimisation, on which single- and two-mesh optima violate the data (E19–E21, E29, E32, E36).
* **C3 (convergence).** The bulk heat-flux memory time is not mesh-converged (RTA, tetrahedron integration:
  tau_mem,x/z = 135.4/99.8, 152.8/124.1, 158.0/132.3 ps on 31x31x17, 48x48x26, 56x56x30 while kappa changes by 0.1–1.4 %;
  full operator, z: 27.4 to 70.1 ps from 5x5x3 to 12x12x8), whereas the film-regularised value is (RTA:
  tau_z(30/100 nm) = 5.295/10.728 and 5.299/10.737 ps on the two finest meshes, change < 0.1 %; full operator: within
  ±2–3 % of the mean over four meshes) (E27, E30).

**Relation to the earlier protocol.** `review/novelty-protocol.md` (11 October; claims L1–L4) established that the RTA
formula is the published heat-flux relaxation time (Sendra et al. 2021; Beardo et al. 2021, Eq. 2), that its infrared
power counting follows from Esfarjani et al. 2011 (Eq. 9) and Vermeersch et al. 2015 (Eq. 19) (L1, L2), that
size-dependent memory parameters are known phenomenologically (L3), and that identifiability of the bulk tau_mem in the
fc3 class is partially known in spirit (L4). Those verdicts and sources are used here by label (e.g. "L1-3" = item 3 of
section L1) and not repeated. New here: the film-regularised quantity as the object of identifiability (C1), the
multi-mesh bounding procedure as a method (C2), and the paired bulk/film convergence statement on finer meshes (C3).

**Evidence labels.** SOURCE SAYS / SOURCE DERIVES / SOURCE OBSERVES (numerical) / SECONDARY (one paper's statement about
another, not checked in the original) / OUR DERIVATION / OUR INTERPRETATION. Access is stated for every source: read
(text inspected), text search (targeted passages of the full text), abstract only, metadata only. All 88 DOIs cited
resolve: 86 were checked against Crossref records (title, authors, venue, volume, pages) and 2 through the DOI resolver
with OpenAlex metadata after Crossref rate limiting (`crossref-verification.txt`). arXiv identifiers were confirmed
through the arXiv API (`arxiv-ids-cited.txt`).

## 1. Verdicts

| Claim | Verdict | Closest prior work | What is not in the searched literature |
|---|---|---|---|
| **C1** bulk tau_mem not identifiable from static data in the fc3 class; film tau_z(L) constrained to ~±10 % (one mesh), ~±5 % (two), +3.8 % (three) | **PARTIALLY KNOWN** in concept; the specific result is **POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE** (numerical, inner bounds, synthetic data) | Prediction-level identifiability with unidentifiable parameters (Kreutz, Raue & Timmer 2012; Gutenkunst et al. 2007; Kurniawan et al. 2026, information matching for quantities of interest); partial identification (Tamer 2010); scale-dependent identifiability in phonon-transport inverse problems (Chen, Gamba, Li & Nair 2026: ill-posed as Kn → 0); first-principles uncertainty concentrated in long-mean-free-path acoustic phonons (Parks et al. 2020; Jiang, Lindsay & Koh 2016) and sensitivity of hydrodynamic lengths/times to the scattering model (Tur-Prats et al. 2025; Raya-Moreno, Carrete & Cartoixà 2022); static data and fc3 (Lee & Gan 2017; Xie, Gu & Bao 2017) | An identifiability analysis of any hydrodynamic closure coefficient (memory/Cattaneo time, nonlocal length, viscosity) from static data within a parametrised anharmonic class; the bulk-versus-film contrast; narrowing of an admissible range with the number of Brillouin-zone meshes carrying the data |
| **C2** multi-mesh constrained optimisation over fc3 (augmented Lagrangian, smooth minimum across meshes) to bound a derived transport quantity, with held-out finer-mesh validation | **PARTIALLY KNOWN**: framework and components KNOWN; the application is **POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE** (methodologically incremental) | Bound-to-bound data collaboration (Frenklach et al. 2004; Feeley et al. 2004; Russi et al. 2010; Edwards et al. 2014; Frenklach et al. 2016; Hegde et al. 2018); optimal UQ (Owhadi et al. 2013); prediction profile likelihood and constrained-optimisation confidence intervals (Kreutz et al. 2012; Borisov & Metelkin 2020); set-membership estimation (Milanese & Vicino 1991); multiplier methods (Hestenes 1969; Powell 1969; Rockafellar 1973); L-BFGS-B (Byrd et al. 1995); KS aggregation (Kreisselmeier & Steinhauser 1979); min–max robust formulations ensuring mesh convergence (Wang, Lazarov & Sigmund 2011); mesh-dependent optima (Sigmund & Petersson 1998); inverse crime (Wirgin 2004; Kaipio & Somersalo 2007); discretisation invariance (Lassas, Saksman & Siltanen 2009) | Optimisation over a symmetry-adapted fc3 tensor to bound another transport quantity under transport data; constraints imposed on several Brillouin-zone discretisations with a worst-case objective across them; held-out-mesh validation; any phonon-transport application of feasible-set prediction bounds (none among the forward citations of the framework papers) |
| **C3** bulk memory time not mesh-converged, film-regularised value converged | **PARTIALLY KNOWN** (bulk part known in substance; regularisation by finite size known in general; at RTA level the film convergence is guaranteed by construction); the explicit paired ab initio demonstration is **POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE** (incremental) | L1/L2 sources (Esfarjani et al. 2011; Vermeersch et al. 2015; Lepri, Livi & Politi 2003; Akhiezer-regime literature; Tur-Prats et al. 2025); 2D analogue: q-grid divergence of kappa tied to the fc3 range and long-wavelength lifetimes (Jiang, Lu, Ouyang & Chen 2022; Qin & Hu 2018); finite size as regulator of infrared-divergent coefficients (Lepri et al. 2003; Majee & Akšamija 2016; Scalfi et al. 2023; Fiorentino, Pegolo & Baroni 2023); size-dependent memory parameters (Alvarez & Jou 2007; Krewer et al. 2020); single-mesh Ge value contested (Huberman, Zhang & Abou Haibeh 2025) | A mesh sequence for the heat-flux relaxation time of a 3D crystal; its contrast with a finite-size-regularised value converged to < 0.1 % |

"Potentially novel within the search scope" always means: not found in the indexes and phases of section 2, with the
limits of section 7, together with the scope of the 11 October protocol (its section 2).

## 2. Search scope (common to C1–C3)

| Index / source | Use | Outcome |
|---|---|---|
| Local multi-index search (Crossref + Semantic Scholar + OpenAlex APIs) | 41 queries labelled by claim and phase (C1: 20, C2: 11, C3: 10; P1 direct, P2 mathematical, P3 physical, P4 historical, P5 backward, P7 neighbouring), 15 results per provider | Crossref 41/41, OpenAlex 41/41, Semantic Scholar 3/41 (HTTP 429); 1,241 records, titles in `titles_all.txt`. Noisy; no absence statement rests on it alone |
| arXiv API, scripted | 36 field-restricted queries (C1: 17, C2: 7, C3: 12), up to 50 records each | `arxiv_results.md`; exact-phrase queries for "heat flux relaxation time", film relaxation time with first principles, and q-point convergence of hydrodynamic parameters returned 0 records |
| arXiv, hand-run; full or targeted text | ~10 title/author/field queries; targeted text of 9 preprints (Appendix A) | Main source of full texts |
| Semantic Scholar, hand-run | 3 queries | 1 answered (no relevant record), 2 rate-limited; full-text snippet search unavailable |
| OpenAlex API | abstracts and metadata for 80 DOIs; keyword-filtered forward citations for 28 (DOI, keyword) pairs | worked |
| Crossref API and DOI resolver | verification of all 88 cited DOIs | all resolve |
| Publisher and repository pages | Springer (2 abstracts), PubMed Central (Jiang et al. 2022, targeted), arXiv HTML (Shang et al. 2025, targeted), AMS (Hourdin et al. 2017: statements could not be extracted; not used) | see Appendix A |
| General web search | 17 queries (3 research-filtered, 14 unfiltered) | pointers; every source used was then checked in an index or at the publisher |
| Project reference library | 2 semantic queries | same four items as on 2 and 11 October (Chen et al. 2025; Camacho de la Rosa et al. 2025; Hoque et al. 2024; Krapez 2018); none bears on C1–C3 beyond the uses below |

Terminology searched (phases 1–4): identifiability / sensitivity / uncertainty of phonon-hydrodynamic parameters,
second sound, Guyer–Krumhansl parameters; heat-flux relaxation time in films; prediction bounds consistent with data;
prediction profile likelihood; bound-to-bound data collaboration; optimal uncertainty quantification; estimable
function; set-membership estimation; partial identification; inverse determination of anharmonic force constants from
linewidths, neutron or X-ray scattering, Grüneisen parameters, third-order elastic constants; information content of
kappa about lifetimes; uniqueness and stability for inverse problems of the phonon transport equation; augmented
Lagrangian; Kreisselmeier–Steinhauser aggregation; scenario / worst-case robust optimisation; mesh dependence of optimal
designs; approximation-error approach; inverse crime; discretisation-invariant inversion; out-of-sample validation across
discretisations; differentiable / adjoint phonon BTE; q-mesh convergence of hydrodynamic parameters, thermal viscosity,
nonlocal length; boundary scattering as regulator of divergent conductivity; finite-size cutoff of long-time tails;
Burnett coefficients in finite systems; memory kernel and Mori–Zwanzig formulations of heat transport; Akhiezer damping
in films; infrared-safe observables; thickness-dependent relaxation times of electrons in films.
Neighbouring fields (phase 7): combustion-kinetics uncertainty quantification, parametrised quantum chemistry, systems
biology, econometrics, robust control, structural and topology optimisation, large-eddy simulation, statistical inverse
problems, machine-learned interatomic potentials, operator learning, molecular-dynamics finite-size corrections, thermal
metrology, conduction electrons in metal films.

---

## 3. C1 — identifiability of the memory time from static data: bulk versus film

**Claim audited.** As in section 0. The bulk part was audited as L4 (verdict: partially known in spirit; specific
statement potentially novel). The new elements are the film-regularised quantity, its much narrower admissible range,
and the narrowing with the number of meshes on which the data are imposed.

**Phases.** P1 "identifiability phonon hydrodynamic parameters first principles", "uncertainty heat flux relaxation time
first principles thin film", "uncertainty quantification second sound phonon hydrodynamics force constants", "sensitivity
Guyer-Krumhansl parameters anharmonic force constants", arXiv `abs:"phonon hydrodynamics" AND (uncertainty OR
sensitivity OR identifiability)`, `abs:"second sound" AND (uncertainty OR "exchange-correlation" OR sensitivity)`,
`abs:"force constants" AND (ensemble | Bayesian | "uncertainty quantification")`; P2 prediction bounds over data-consistent
sets, prediction profile likelihood, bound-to-bound data collaboration, optimal UQ, estimable functions, partial
identification; P3 inverse determination of fc3 (linewidths, neutron scattering, Grüneisen data), information content of
kappa about lifetimes, uniqueness/stability for the phonon transport equation, Guyer–Krumhansl parameter estimation from
heat-pulse data; P4 historical force-constant determination (least-squares Born–von Kármán fits, inelastic neutron
scattering, third-order elastic constants), set-membership estimation; P5 references of Kreutz et al. 2012, Hegde et al.
2018 and Jiang et al. 2022; P6 keyword-filtered forward citations of Parks et al. 2020 (7 citing records, all scanned),
Abou Haibeh & Huberman 2026 (1), Phonon Olympics 2025 (10), Tur-Prats et al. 2025 (3), Kurniawan et al. 2022 ("phonon" 3,
"thermal conductivity" 2), Kreutz et al. 2012 ("phonon" 0, "thermal" 2), Lee & Gan 2017 ("uncertainty" 4, "inverse" 10);
P7 combustion kinetics, systems biology, econometrics, statistical inverse problems, interatomic-potential UQ, thermal
metrology.

### Closest prior work

1. **Identifiability of predictions rather than of parameters** (neighbouring fields).
   * **C. Kreutz, A. Raue, J. Timmer**, "Likelihood based observability analysis and confidence intervals for predictions
     of dynamic models", *BMC Syst. Biol.* **6**, 120 (2012), doi:10.1186/1752-0509-6-120; arXiv:1107.0013. **Abstract
     only.** SOURCE SAYS: "reliable confidence intervals are calculated based on the prediction profile likelihood. Such
     prediction confidence intervals of the dynamic states can be utilized for a data-based observability analysis. The
     method is also applicable if there are non-identifiable parameters yielding to some insufficiently specified model
     predictions that can be interpreted as non-observability."
   * **R. N. Gutenkunst et al.**, *PLoS Comput. Biol.* **3**, e189 (2007) — L4-4 (sloppy parameters, predictions).
   * **Y. Kurniawan, T. B. Neilsen, B. L. Francis, A. M. Stanković, M. Wen, I. Nikiforov, E. B. Tadmor, V. V. Bulatov,
     V. Lordi, M. K. Transtrum**, "An information-matching approach to optimal experimental design and active learning",
     *Appl. Phys. Lett.* **128**, 064104 (2026), doi:10.1063/5.0296026. **Abstract only.** SOURCE SAYS: "Because models
     often contain many unidentifiable (sloppy) parameters, QoIs often depend on a relatively small number of parameter
     combinations"; the criterion "ensures that the selected data contain sufficient information to learn only those
     parameters that are needed to constrain downstream QoIs".
   * **E. Tamer**, "Partial Identification in Econometrics", *Annu. Rev. Econ.* **2**, 167 (2010),
     doi:10.1146/annurev.economics.050708.143401. **Abstract only.** SOURCE SAYS: identification "should not result in a
     binary answer that consists of whether the parameter is point identified"; "more assumptions beget more information".
   * **O. J. Maclaren, R. Nicholson**, "What can be estimated? Identifiability, estimability, causal inference and
     ill-posed inverse problems", arXiv:1904.02826 (2019). **Abstract only.** SOURCE SAYS: "identifiability corresponds to
     the question of uniqueness; in contrast, we take estimability to mean satisfaction of all three conditions, i.e.
     well-posedness."
   * OUR INTERPRETATION: C1 is a statement of exactly this type for a physical quantity of interest — the bulk tau_mem is
     non-observable from the static data in the fc3 class, the film tau_z(L) is (practically) observable — and its
     class dependence (eps = 0.1 versus calibrated bounds) is the dependence on assumptions that partial identification
     makes explicit. The concept is known; the phonon-transport instance is not in the searched literature.

2. **Scale dependence of identifiability in phonon-transport inverse problems.**
   * **P. Chen, I. M. Gamba, Q. Li, A. Nair**, "Stability of the reconstruction of the heat reflection coefficient in the
     phonon transport equation", *Kinet. Relat. Models* **24**, 109–136 (2026), doi:10.3934/krm.2026021; arXiv:2512.24394.
     **Abstract only.** SOURCE SAYS: "the problem becomes ill-posed as the system transitions from the ballistic to the
     diffusive regime, characterized by the Knudsen number converging to zero"; the rate of deterioration is quantified.
   * **P. Chen, I. M. Gamba, Q. Li, L. Wang**, arXiv:2502.19533 (2025) (relaxation-index reconstruction; read in the
     2 October campaign) and **M. Forghani, N. G. Hadjiconstantinou** (2016–2021; L4-5): relaxation-time reconstruction
     from thermal spectroscopy; a given F(Lambda) "does not correspond to a unique relaxation-times function".
   * OUR INTERPRETATION: there the data are dynamic temperature responses and the unknown a kinetic coefficient; in C1 the
     data are static bulk properties and the unknown a moment of the collision operator. Common thread: what the data fix
     depends on the ratio of system size to the mean free paths of the modes that carry the unknown. No source found makes
     the C1 statement.

3. **Long-wavelength phonons as the locus of first-principles disagreement and of hydrodynamic sensitivity.**
   * **H. L. Parks et al.** 2020 (L4-1: kappa ensemble differences correlated with transverse-acoustic phonons of mean
     free path 100–300 nm); **P. Jiang, L. Lindsay, Y. K. Koh** 2016 (L1-8: discrepancies "as large as 30%" for low-energy
     phonons; Akhiezer damping removes mean free paths > 20 μm at 300 K); **J. Tur-Prats et al.** 2025 (L1-3:
     low-frequency phonons dominate converged l and tau; l and tau "significantly reduced" with full anharmonicity).
   * **M. Raya-Moreno, J. Carrete, X. Cartoixà**, "Hydrodynamic signatures in thermal transport in devices based on
     two-dimensional materials: An ab initio study", *Phys. Rev. B* **106**, 014308 (2022),
     doi:10.1103/PhysRevB.106.014308. **Abstract only.** SOURCE SAYS: "boundary scattering and the relation of sample
     dimensions to the nonlocal length l are the determinant factors"; "approximations made on the scattering operator may
     have, through the value of l, qualitative consequences on the signatures of hydrodynamic behavior."
   * **J. Abou Haibeh, S. Huberman**, *Mater. Res. Express* **13**, 085901 (2026), doi:10.1088/2053-1591/ae53a0;
     arXiv:2311.05107v3. **Text search.** SOURCE SAYS: functional choice changes kappa and the hydrodynamic windows
     (L4-6). SOURCE OBSERVES: kappa converged on a 40^3 q-mesh; Appendix C shows mode-averaged normal and umklapp rates at
     20 K versus q-mesh for three functionals. No heat-flux relaxation time was found (a text search of the arXiv version
     for "relaxation time" returned no passage).
   * OUR INTERPRETATION: "first-principles uncertainty concentrates in long-mean-free-path acoustic phonons" plus "a film
     of thickness L removes contributions of mean free paths ≳ L" implies that the film quantity should be better
     constrained. C1 is consistent with these sources; none states it. Note that the Parks et al. window (100–300 nm) is
     only partly removed by 100–300 nm films.

4. **Static data and the cubic force constants.**
   * **C. H. Lee, C. K. Gan**, "Anharmonic interatomic force constants and thermal conductivity from Grüneisen
     parameters: An application to graphene", *Phys. Rev. B* **96**, 035105 (2017), doi:10.1103/PhysRevB.96.035105.
     **Abstract only.** SOURCE SAYS: "although the Grüneisen parameters cannot possibly contain all the information on the
     cubic IFCs, being derivable from spatially uniform deformations, they can still unambiguously and accurately
     determine the largest and most physically relevant ones"; the dominant cubic constants are "reverse-engineer[ed]"
     from anisotropic Grüneisen data and give an approximate kappa of graphene.
   * **H. Xie, X. Gu, H. Bao**, *Comput. Mater. Sci.* **138**, 368 (2017), doi:10.1016/j.commatsci.2017.06.047;
     arXiv:1706.08232. **Abstract only.** SOURCE SAYS: with fc3 errors at typical first-principles force accuracy "the
     calculated thermal conductivity would be quite different from the benchmark result"; "Grüneisen parameters cannot be
     used as a necessary and sufficient criterion for the accuracy of third-order IFCs in the aspect of predicting thermal
     conductivity." A counterpoint to item 4a for kappa.
   * Historical and harmonic analogues: **V. Sundara Raja, P. Jayarama Reddy**, *Pramana* **8**, 402 (1977),
     doi:10.1007/BF02862892 (**metadata only**: anharmonic force constants of GaSb from measured third-order elastic
     constants); **H. Kohl**, *phys. status solidi (b)* **130**, 151 (1985), doi:10.1002/pssb.2221300112 (**abstract
     only**: force constants "can be obtained experimentally by measuring both the frequency and a corresponding
     eigenvector of each phonon mode", whereas the double differential cross section "does not generally suffice to
     determine the eigenvector unambiguously"); **K. Suchan, S. Mohanty, H. Zhai, W. Cai**, arXiv:2603.28683 (2026)
     (**abstract only**: differentiable, symmetry-reduced harmonic force-constant fit to thermal diffuse scattering;
     recovers ground truth in benchmarks); **A. Sable et al.**, *Digital Discovery* **5**, 1545 (2026),
     doi:10.1039/D6DD00008H (**abstract only**: machine-learning inversion of force constants from inelastic neutron
     scattering); **R. J. Stoner, H. J. Maris**, in *Phonon Scattering in Condensed Matter VII*, Springer Ser. Solid-State
     Sci., 17–18 (1993), doi:10.1007/978-3-642-84888-9_5 (**metadata only**; title: "... Determination of Phonon Lifetimes
     from Thermal Conductivity"; not read).
   * OUR INTERPRETATION: these address what static data fix about the force constants or about kappa. None addresses a
     derived dynamic moment, and none maps the attainable range of a transport quantity over a constrained fc3 class.

5. **Identifiability of relaxation times from dynamic thermal measurements** (the experimental counterpart).
   * **S. Huberman, C. Zhang, J. Abou Haibeh**, "On the question of second sound in germanium: a theoretical viewpoint",
     *J. Phys.: Condens. Matter* **37**, 435601 (2025), doi:10.1088/1361-648X/ae0dd6; arXiv:2206.02769v3. **Text search.**
     SOURCE SAYS: "phase lag alone is not a suitable experimental observable for inferring second sound"; "the RTA-BTE
     contains ample flexibility in modeling deviations from diffusion in germanium".
   * **A. Camacho de la Rosa et al.** 2025 (L3-5: characterising the relaxation time by thermoreflectance "remains an open
     problem"); **T. Liu, F. Chen, J. Li, C. Tong**, *Adv. Sci. Instrum.* **1**, 100037 (2026),
     doi:10.1016/j.asi.2026.100037 (**abstract only**: thermoreflectance parameters come from "an inverse problem, where
     strong parameter coupling can render the extraction ill-posed"); **A. É. Fehér, R. Kovács**, *Int. J. Eng. Sci.*
     **169**, 103577 (2021), doi:10.1016/j.ijengsci.2021.103577 (**abstract only**: evaluation of Guyer–Krumhansl
     parameters from heat-pulse temperature histories); in the opposite direction **A. Beardo et al.**, *Phys. Rev. B*
     **105**, 165303 (2022), doi:10.1103/PhysRevB.105.165303 (**abstract only**: hydrodynamic-like equations "successfully
     predict the local system response at temperatures ranging from 100 to 300 K" in silicon).
   * OUR INTERPRETATION: complementary to C1. C1 asks whether first-principles-type static information fixes the memory
     time; these ask whether measurements do. Both lines indicate a weakly determined bulk relaxation time.

### What differs
* No source found poses identifiability of a hydrodynamic closure coefficient from static data within a parametrised
  anharmonic (fc3) class.
* No source found contrasts the identifiability of an infrared-controlled bulk quantity with that of its
  finite-size-regularised counterpart.
* No source found imposes data on several Brillouin-zone discretisations and reports how an admissible range narrows
  with their number.

### Caveats raised by the literature and by the construction
* **Monotonicity is built in.** OUR DERIVATION: for a fixed objective functional and exact global optimisation, adding
  the data of another mesh can only shrink the feasible set (F_{M1,M2} ⊆ F_{M1}) and hence the admissible range. The
  narrowing is therefore expected; its information is the rate and the limit. The reported sequence at 100 nm (1.098,
  1.046, 1.038) compares maxima evaluated on different objective meshes and uses local inner bounds (single or warm
  starts), so it is a NUMERICAL OBSERVATION, not a monotone statement; whether it tends to zero or to a finite width is
  open (E36). In data-collaboration terms this is the joint-analysis effect: Frenklach et al. 2004 (section 4, item 1)
  SAY that "combining the entire collection of data into a joint analysis extracts substantially more of the information
  content of the data".
* **Synthetic data on the inverting operator ("inverse crime").** The data are the reference tensor's own predictions on
  the same discretised operator. **A. Wirgin**, "The inverse crime", arXiv:math-ph/0401050 (2004), **abstract only**:
  "The inverse crime occurs when the same (or very nearly the same) theoretical ingredients are employed to synthesize as
  well as to invert data in an inverse problem." Real data constrain the physical crystal; discretisation and model error
  would have to be represented explicitly, e.g. as model discrepancy in feasible-set methods (Li et al. 2021, section 4)
  or misspecification-aware UQ (**D. Perez, A. P. A. Subramanyam, I. Maliyov, T. D. Swinburne**, *npj Comput. Mater.*
  **11**, 263 (2025), doi:10.1038/s41524-025-01758-4, **abstract only**: misspecification, "the inability of any one
  choice of model parameters to exactly match all training data", is "a key contributor to errors that is often
  disregarded").
* **Class dependence.** The ranges depend on the assumed class (Tamer 2010). The calibrated class rests on one workflow
  pair at one functional, while spreads between packages and functionals are larger (**Y. Srivastava, A. G. Gokhale,
  A. Jain**, *Phys. Rev. B* **112**, 195201 (2025), doi:10.1103/2n3p-hz62, abstract: ~1 % between solvers, ~10 % between
  DFT packages, > 20 % between functionals for kappa; Phonon Olympics ±15 %, L4-2).
* **Regulator.** The film quantity uses the Matthiessen regulator (section 5 caveats).

**Verdict C1: PARTIALLY KNOWN in concept** (prediction-level identifiability with unidentifiable parameters;
assumption-dependent partial identification; scale-dependent identifiability in phonon-transport inverse problems;
long-mean-free-path acoustic phonons as the locus of first-principles uncertainty). **The specific result — bulk tau_mem
not identifiable from static data in the fc3 class, film-regularised tau_z(L) constrained to about ±10 % and narrowing to
about ±5 % and +3.8 % as the data are imposed on two and three Brillouin-zone meshes — is POTENTIALLY NOVEL WITHIN THE
SEARCH SCOPE** (numerical inner bounds on discretised operators, synthetic data, one material, one temperature).

---

## 4. C2 — multi-mesh constrained optimisation over fc3 to bound a derived transport quantity

**Claim audited.** As in section 0 (implementation: `scripts_c/robust_mesh.py`; exact linearisation of the amplitudes
in scaled per-orbit coordinates; constraints on all lifetimes and on kappa_xx, kappa_zz on every constraint mesh; per-orbit
class bounds; sense max or min).

**Phases.** P1 "constrained optimization force constants thermal conductivity", "optimizing interatomic force constants
phonon properties fit", "automatic differentiation lattice thermal conductivity force constants gradient",
"differentiable Boltzmann transport equation phonon inverse design", arXiv `abs:"force constants" AND ("automatic
differentiation" OR differentiable)`, `abs:phonon AND abs:"augmented Lagrangian"` (0 records); P2 Kreisselmeier–Steinhauser
aggregation, multiplier methods, scenario / worst-case robust optimisation, worst-case quantity-of-interest bounds; P3 mesh
dependence of optimal designs, multi-fidelity robust optimisation, approximation-error approach, out-of-sample validation
across discretisations; P4 multiplier methods (1969–1973), KS function (1979), set-membership estimation (1990–1991), data
collaboration (2004); P5 references of Hegde et al. 2018 (data-collaboration lineage, Edwards et al. 2014); P6
keyword-filtered forward citations of Frenklach et al. 2016 ("phonon" 0, "thermal conductivity" 0), Hegde et al. 2018
("phonon" 0), Russi et al. 2010 ("phonon" 1, irrelevant; "thermal" 12, none in phonon transport), Frenklach et al. 2004
("phonon" 0, "interatomic" 0), Edwards et al. 2014 (all 15 citing records scanned; none in phonon transport), Owhadi et al.
2013 ("phonon" 0; "thermal conductivity" 4, none in phonon transport), Kreutz et al. 2012 ("phonon" 0), Kaipio &
Somersalo 2007 ("phonon" 1, irrelevant; "Boltzmann" 8, none in phonon transport); P7 combustion kinetics, quantum
chemistry, systems biology, topology optimisation, large-eddy simulation, machine learning, statistical inverse problems.

### Closest prior work

1. **Bounding a prediction over the set of parameters consistent with bounded data.**
   * **M. Frenklach, A. Packard, P. Seiler, R. Feeley**, "Collaborative data processing in developing predictive models
     of complex reaction systems", *Int. J. Chem. Kinet.* **36**, 57–66 (2004), doi:10.1002/kin.10172. **Abstract only.**
     SOURCE SAYS: a methodology "founded on a combination of response surface technique and robust control theory";
     "combining the entire collection of data into a joint analysis extracts substantially more of the information content
     of the data".
   * **R. Feeley, P. Seiler, A. Packard, M. Frenklach**, "Consistency of a Reaction Dataset", *J. Phys. Chem. A* **108**,
     9573–9583 (2004), doi:10.1021/jp047524w. **Abstract only.** SOURCE SAYS: "The numerical foundation of data
     collaboration lies in constrained optimization"; "Lagrange multipliers are used to identify factors that influence
     consistency".
   * **T. Russi, A. Packard, M. Frenklach**, *Chem. Phys. Lett.* **499**, 1–8 (2010), doi:10.1016/j.cplett.2010.09.009
     (**metadata only**).
   * **M. Frenklach, A. Packard, G. García-Donato, R. Paulo, J. Sacks**, "Comparison of Statistical and Deterministic
     Frameworks of Uncertainty Quantification", *SIAM/ASA J. Uncertain. Quantif.* **4**, 875–901 (2016),
     doi:10.1137/15M1019131. **Abstract only.** SOURCE SAYS: bound-to-bound data collaboration "deploys semidefinite
     programming algorithms where the initial bounds on unknowns are combined with initial bounds of experimental data to
     produce new uncertainty bounds for the unknowns that are consistent with the data and, finally, deterministic
     uncertainty bounds for prediction in new settings" (102 parameters, 76 quantities of interest; consistent with
     Bayesian calibration and prediction).
   * **A. Hegde, W. Li, J. Oreluk, A. Packard, M. Frenklach**, "Consistency Analysis for Massively Inconsistent Datasets in
     Bound-to-Bound Data Collaboration", *SIAM/ASA J. Uncertain. Quantif.* **6**, 429–456 (2018), doi:10.1137/16M1110005;
     arXiv:1701.04695v3. **Read** (introduction and the section on datasets and consistency; Appendix A by text
     search). SOURCE SAYS: feasible set
     F = {x ∈ H : L_e ≤ M_e(x) ≤ U_e, e = 1, ..., N} (Eqs. 2.1–2.2), consistency measured by the maximal uniform tightening of
     all constraints (Eq. 2.3); SDP relaxations for quadratic surrogate models (Appendix A); applications in combustion,
     atmospheric chemistry, quantum chemistry and systems biology.
   * **W. Li, A. Hegde, J. Oreluk, A. Packard, M. Frenklach**, *SIAM/ASA J. Uncertain. Quantif.* **9**, 231–259 (2021),
     doi:10.1137/19M1270185 (**abstract only**: model discrepancy added to the feasible-set framework).
   * Physical-sciences application: **D. E. Edwards, D. Yu. Zubarev, A. Packard, W. A. Lester, M. Frenklach**, "Interval
     Prediction of Molecular Properties in Parametrized Quantum Chemistry", *Phys. Rev. Lett.* **112**, 253003 (2014),
     doi:10.1103/PhysRevLett.112.253003. **Abstract only.** SOURCE SAYS: "We propose a strategy of directly predicting the
     uncertainty interval for a property of interest, based on training-data uncertainties, which sidesteps the need for
     an optimum set of parameters."
   * **H. Owhadi, C. Scovel, T. J. Sullivan, M. McKerns, M. Ortiz**, "Optimal Uncertainty Quantification", *SIAM Rev.*
     **55**, 271–345 (2013), doi:10.1137/10080782X; arXiv:1009.0679. **Abstract only.** SOURCE SAYS: optimal bounds "are
     obtained as values of well-defined optimization problems ... subject to the constraints imposed by the scenarios
     compatible with the assumptions and information".
   * **I. Borisov, E. Metelkin**, *PLoS Comput. Biol.* **16**, e1008495 (2020), doi:10.1371/journal.pcbi.1008495
     (**abstract only**: "Confidence Intervals by Constraint Optimization (CICO)" based on the profile likelihood);
     **M. J. Simpson, O. J. Maclaren**, *PLoS Comput. Biol.* **19**, e1011515 (2023), doi:10.1371/journal.pcbi.1011515
     (**abstract only**: profile-wise prediction intervals); Kreutz et al. 2012 (section 3, item 1);
     **M. Milanese, A. Vicino**, *Automatica* **27**, 997–1009 (1991), doi:10.1016/0005-1098(91)90134-N (**metadata only**:
     set-membership estimation).
   * OUR INTERPRETATION: C2's core operation — extremise a quantity of interest over the set of parameters satisfying
     interval data constraints and prior bounds — is the data-collaboration prediction problem, here with the fc3 tensor
     as parameters, an exact linearised model instead of a surrogate, and no outer bound.

2. **Numerical components** (standard). **M. R. Hestenes**, *J. Optim. Theory Appl.* **4**, 303–320 (1969),
   doi:10.1007/BF00927673 (**metadata only**); **M. J. D. Powell**, "A method for nonlinear constraints in minimization
   problems", in *Optimization*, ed. R. Fletcher, Academic Press (1969), pp. 283–298 (**not retrieved**; attribution as in
   the title of the next item); **R. T. Rockafellar**, "The multiplier method of Hestenes and Powell applied to convex
   programming", *J. Optim. Theory Appl.* **12**, 555–562 (1973), doi:10.1007/BF00934777 (**metadata only**);
   **R. H. Byrd, P. Lu, J. Nocedal, C. Zhu**, *SIAM J. Sci. Comput.* **16**, 1190–1208 (1995), doi:10.1137/0916069
   (**abstract only**: gradient projection with a limited-memory BFGS matrix for bound-constrained problems), and
   **C. Zhu, R. H. Byrd, P. Lu, J. Nocedal**, *ACM Trans. Math. Softw.* **23**, 550–560 (1997), doi:10.1145/279232.279236;
   **G. Kreisselmeier, R. Steinhauser**, "Systematic Control Design by Optimizing a Vector Performance Index", *IFAC Proc.
   Vol.* **12**, 113–117 (1979), doi:10.1016/S1474-6670(17)65584-8 (**metadata only**; the standard attribution of the
   log-sum-exp aggregation, as used in the constraint-aggregation papers returned by the searches).

3. **Worst case over several model variants; information at several resolutions** (neighbouring fields).
   * **F. Wang, B. S. Lazarov, O. Sigmund**, "On projection methods, convergence and robust formulations in topology
     optimization", *Struct. Multidisc. Optim.* **43**, 767–784 (2011), doi:10.1007/s00158-010-0602-y. **Abstract only.**
     SOURCE SAYS: "simple projection methods do not ensure local mesh-convergence"; a "modified robust topology optimization
     formulation based on erosion, intermediate and dilation projections ... ensures both global and local
     mesh-convergence" (a min–max over three realisations of the design).
   * **M. Germano, U. Piomelli, P. Moin, W. H. Cabot**, "A dynamic subgrid-scale eddy viscosity model", *Phys. Fluids A* **3**,
     1760–1765 (1991), doi:10.1063/1.857955. **Abstract only.** SOURCE SAYS: the coefficient "is computed dynamically ...
     rather than input a priori", based on "an algebraic identity between the subgrid-scale stresses at two different
     filtered levels". OUR INTERPRETATION: the closest conceptual relative of "data on several meshes pin the closure
     coefficient".

4. **Parameters fitted on one discretisation failing on another.**
   * Wirgin 2004 (section 3 caveats); **J. Kaipio, E. Somersalo**, "Statistical inverse problems: Discretization, model
     reduction and inverse crimes", *J. Comput. Appl. Math.* **198**, 493–504 (2007), doi:10.1016/j.cam.2005.09.027
     (**metadata only**).
   * **O. Sigmund, J. Petersson**, *Struct. Optim.* **16**, 68–75 (1998), doi:10.1007/BF01214002. **Abstract only.** SOURCE
     SAYS: "The mesh-dependence problem refers to obtaining qualitatively different solutions for different mesh-sizes or
     discretizations. Local minima refers to the problem of obtaining different solutions to the same discretized problem
     when choosing different algorithmic parameters."
   * **M. Lassas, E. Saksman, S. Siltanen**, *Inverse Probl. Imaging* **3**, 87–122 (2009), doi:10.3934/ipi.2009.3.87.
     **Abstract only.** SOURCE SAYS: priors are discretisation-invariant when they "represent the same a priori information
     for all n" and the reconstruction "converges to a limit estimate as k,n → ∞".
   * **R. L. S. Silva et al.**, arXiv:2604.01101 (2026). **Abstract only.** SOURCE SAYS: calibration of a temperature-dependent
     conductivity with an adaptive algorithm "that sequentially refines both the numerical discretization ... and the model
     complexity", so that "the numerical and modeling biases remain of the same order as the uncertainty imposed by the
     measurement noise".
   * **Z. Li et al.**, "Fourier Neural Operator for Parametric Partial Differential Equations", arXiv:2010.08895 (2020).
     **Abstract only.** SOURCE SAYS: "zero-shot super-resolution" (evaluation at resolutions not used in training).
   * OUR INTERPRETATION: C2's finding that coarse-mesh optima violate finer-mesh data is the phonon-transport instance of a
     feasible set that is not discretisation-invariant: constraints on a coarse mesh do not control the modes it does not
     sample, mostly Gamma-adjacent ones, and the optimiser exploits them (E19, E32: an engineered slow eigenpair among the
     objective mesh's Gamma-adjacent modes).

5. **Optimisation through phonon models** (phonon side).
   * **A. Rohskopf, H. R. Seyf, K. Gordiz, T. Tadano, A. Henry**, *npj Comput. Mater.* **3**, 27 (2017),
     doi:10.1038/s41524-017-0026-y (**abstract only**: a genetic algorithm fits empirical potential parameters to ab initio
     phonon properties); Lee & Gan 2017 (section 3, item 4); Suchan et al. 2026 (section 3, item 4).
   * **W. Shang et al.**, *npj Comput. Mater.* **11**, 129 (2025), doi:10.1038/s41524-025-01635-0; arXiv:2503.23657v2.
     **Text search.** SOURCE SAYS: a differentiable phonon-BTE solver; its inverse example determines "the unknown film
     thickness L that results in a given temperature profile"; parameters may include "material properties, geometric
     variables, or heat source intensities". No force constants.
   * **G. Romano, S. G. Johnson**, *Struct. Multidisc. Optim.* **65**, 297 (2022), doi:10.1007/s00158-022-03392-w
     (**abstract only**: density-based topology optimisation with adjoint sensitivities of the phonon BTE).
   * None constrains an anharmonic force-constant tensor by transport data to bound another transport quantity, and none
     uses several Brillouin-zone discretisations.

### What differs
* Parameter space: a symmetry-adapted fc3 tensor (thousands of parameters) with an exact linearisation of the scattering
  amplitudes, instead of rate parameters, empirical potential parameters or harmonic force constants.
* Scenarios: Brillouin-zone discretisations of one operator, instead of experiments, load cases or design realisations.
* Validation: optima tested on held-out finer meshes, which exposes data violation (a discretisation-transfer test).
* No phonon-transport application of feasible-set (data-collaboration-type) bounding was found in the forward citations of
  the five framework papers checked.

### Caveats
* All extrema are inner bounds (feasible local iterates; warm starts can select different local optima, E33). Feasible-set
  methods also compute outer bounds by relaxation (Hegde et al., Appendix A, for quadratic surrogates); C2 has no outer
  bound for the per-cluster or calibrated classes (final report, S3).
* OUR DERIVATION: for m objective meshes the smooth minimum satisfies
  min_M s_M − (log m)/beta ≤ −(1/beta) log sum_M exp(−beta s_M) ≤ min_M s_M, i.e. it underestimates the worst case by at most
  (log m)/beta = 0.017 (m = 2) and 0.027 (m = 3) in log ratio for beta = 40. The reported film ratios are evaluated mesh by
  mesh, so the smoothing affects which point is found, not the values quoted.

**Verdict C2: PARTIALLY KNOWN.** The framework — bounding a quantity of interest over the data-consistent parameter set —
is KNOWN (data collaboration since 2004, optimal UQ, prediction profile likelihood, set-membership estimation); the
numerical components are KNOWN; the failure of coarse-discretisation optima on finer discretisations is KNOWN as a
phenomenon in neighbouring fields (inverse crime, mesh-dependent optimal designs, discretisation invariance). **Its
application — constrained optimisation over the fc3 tensor of a crystal with the data imposed on several Brillouin-zone
discretisations, a smooth worst-case objective across them and held-out finer-mesh validation — is POTENTIALLY NOVEL
WITHIN THE SEARCH SCOPE as an application of known methods (methodologically incremental).**

---

## 5. C3 — bulk memory time not mesh-converged; film-regularised memory time converged

**Claim audited.** As in section 0. The bulk part extends L1 (31x31x17) by two finer RTA meshes; the film part extends L3
(convergence to 0.5–1 % on 31x31x17/48x48x26) to < 0.1 % on 56x56x30, and adds the full-operator values on four meshes.

**Phases.** P1 "q-mesh convergence phonon hydrodynamic parameters second sound", "boundary scattering regularizes
divergent thermal conductivity sample size", "size dependent thermal relaxation time nanofilm phonon Boltzmann first
principles", arXiv `abs:"heat flux relaxation time"` (0 records), `abs:"q-point" AND convergence AND (hydrodynamic |
"second sound" | viscosity)` (0), `abs:"nonlocal length" AND phonon`, `abs:"thermal viscosity" AND phonon`; P2 "finite size
Green-Kubo long-time tail cutoff", "Burnett coefficients finite system size regularization", `abs:Burnett AND
divergence`; P3 graphene length dependence and flexural phonons, memory kernel of the heat flux, Mori–Zwanzig and
mode-coupling formulations, Akhiezer damping with films or boundary scattering; P4 infrared-safe observables (one query;
no source used); P5 references of Jiang et al. 2022 (Qin & Hu 2018) and of Huberman et al. 2025; P6 forward citations of
Sendra et al. 2021 (66 records in OpenAlex, all titles scanned), Jiang, Lindsay & Koh 2016 (40, scanned), Tur-Prats et al.
2025 (3), Jiang et al. 2022 (19), Beardo et al. 2021 (OpenAlex filtered by "film" and "convergence": 0; the complete
Semantic Scholar list was scanned on 11 October); P7 molecular-dynamics finite-size corrections, glasses, conduction
electrons in metal films, kinetic theory of gases.

### Closest prior work

1. **Bulk non-convergence**: L1 and L2 — Beardo et al. 2021, Eq. 2 (single 24^3 mesh); Esfarjani et al. 2011, Eq. 9 (kappa
   extrapolation for omega^2 rates); Vermeersch et al. 2015, Eq. 19 (memory moment divergent iff n ≥ 3/2, OUR DERIVATION
   there); Lepri, Livi & Politi 2003 (t^{-3/2} tail in 3D); Akhiezer-regime damping; Tur-Prats et al. 2025 (low-frequency
   dominance of l and tau). Unchanged.

2. **New single-mesh evaluation and a contrary result for Ge.** Huberman, Zhang & Abou Haibeh 2025 (section 3, item 5).
   **Text search.** SOURCE SAYS: "The macroscopic relaxation time, tau, was calculated using the obtained ab initio inputs
   following the procedure described in Ref. [1]" (Beardo et al. 2021), with interpolation on a 48^3 q-mesh chosen "to
   yield convergence of the thermal conductivity"; by eigendecomposition of the full scattering matrix, the requirement
   "that the smallest nonzero eigenvalue must be associated with an odd eigenvector is not satisfied"; direct solutions of
   the linearised BTE "do not reveal the presence of second sound". No convergence test of tau was found in the text. A
   contrary result for the interpretation of the Ge experiment, not a convergence study.

3. **Two-dimensional analogue of the mesh-convergence problem.**
   * **J. Jiang, S. Lu, Y. Ouyang, J. Chen**, "How Hydrodynamic Phonon Transport Determines the Convergence of Thermal
     Conductivity in Two-Dimensional Materials", *Nanomaterials* **12**, 2854 (2022), doi:10.3390/nano12162854. **Read**
     (targeted full text). SOURCE OBSERVES: "When the cutoff radius is small (NNs = 10th, 11th, and 12th), the value of κ
     diverges with the increasing Q-grid"; beyond the 12th neighbour "a convergent κ with respect to the Q-grid size is
     observed in graphene. This is because the long-range interaction can affect the lifetime of low-frequency phonons";
     "this significant convergence issue with respect to rc observed in the 2D material vanishes in the bulk silicon".
   * **G. Qin, M. Hu**, *npj Comput. Mater.* **4**, 3 (2018), doi:10.1038/s41524-017-0058-3 (**abstract only**: convergence
     of kappa with the fc3 cutoff radius; "the long-debate divergence problem of thermal conductivity in low-dimensional
     systems"; the divergent-kappa statement is SECONDARY via Jiang et al.).
   * OUR INTERPRETATION: in 3D, kappa is the "safe" moment and tau_mem, one moment higher, inherits the 2D-type infrared
     sensitivity. The 2D finding that the q-grid behaviour depends on the fc3 range suggests that the infrared law of the
     AlN memory moment may also depend on the fc3 truncation and on sum-rule enforcement; no source tests this for a 3D
     memory time (cf. campaign E4 on real-space truncation).

4. **Finite size as regulator of infrared-divergent transport coefficients.**
   * Lepri, Livi & Politi 2003 (L2-3: divergent kappa(L) in low-dimensional lattices); **A. K. Majee, Z. Akšamija**,
     *Phys. Rev. B* **93**, 235423 (2016), doi:10.1103/PhysRevB.93.235423 (**abstract only**: for ribbon lengths below
     100 μm "a logarithmic divergence of thermal conductivity", arrested by renormalisation of the flexural dispersion).
   * **L. Scalfi, D. Vitali, H. Kiefer, R. R. Netz**, *J. Chem. Phys.* **158**, 191101 (2023), doi:10.1063/5.0151406.
     **Abstract only.** SOURCE SAYS: periodic boundary conditions significantly affect the long-time behaviour of the
     velocity autocorrelation and memory kernel; an analytical correction allows them "to unmask the power-law long-time
     tails of the memory kernel".
   * **A. Fiorentino, P. Pegolo, S. Baroni**, *npj Comput. Mater.* **9**, 157 (2023), doi:10.1038/s41524-023-01116-2.
     **Abstract only.** SOURCE SAYS: "Leveraging hydrodynamic arguments, we show how this issue can be overcome through a
     simple formula to extrapolate a reliable estimate of the bulk thermal conductivity of glasses from finite models".
   * Electrons: Krewer et al. 2020 (L3-3: macroscopic response time set by the mean and variance of the microscopic
     relaxation times, thickness-dependent); phenomenological memory/nonlocal parameters depending on the Knudsen number:
     Alvarez & Jou 2007, 2008 (L3-1).

5. **First-principles hydrodynamic modelling with size** (no film-dependent memory time with a mesh study).
   * **K. Ghosh, A. Kusiak, J.-L. Battaglia**, *Phys. Rev. Mater.* **5**, 073605 (2021), doi:10.1103/PhysRevMaterials.5.073605;
     arXiv:2107.04731. **Text search.** SOURCE SAYS: kinetic-collective model with boundary scattering "included via the
     Matthiessen's rule" in the kinetic term and "a form factor F ... to incorporate boundary scattering in the collective
     term"; a nonlocal length and a "heat wave propagation length" are extracted versus size. No size-dependent heat-flux
     relaxation time or q-mesh study was found in the passages inspected.
   * **Y. Guo et al.**, *Phys. Rev. B* **104**, 075450 (2021), doi:10.1103/PhysRevB.104.075450 (**abstract only**: size
     scaling of kappa in graphite ribbons); Raya-Moreno et al. 2022 (section 3, item 3); **U. Yadav et al.**, *Phys. Rev.
     Appl.* **24**, 034024 (2025), doi:10.1103/xvbk-75bf (**abstract only**: thickness-dependent effective kappa and heat flux
     from a generalised transport equation compared with Monte Carlo BTE); Beardo et al. 2019, 2020 (L3-2).
   * **J. Dragašević, B. Rajkov, M. Simoncelli**, *Phys. Rev. Lett.* **136**, 186302 (2026), doi:10.1103/nbbn-56hr;
     arXiv:2303.12777v7. **Text search** ("converge"): space-dependent simulations of graphite at 70 K use parameters from a 17x17x3
     grid whose bulk kappa is "compatible within 15 %" with the converged 49x49x3 value; no infrared discussion of the
     thermal viscosity was found.

6. **Boundary treatment for transient film responses.**
   * **N. K. Ravichandran, A. J. Minnich**, "Role of thermalizing and nonthermalizing walls in phonon heat conduction along
     thin films", *Phys. Rev. B* **93**, 035314 (2016), doi:10.1103/PhysRevB.93.035314. **Abstract only.** SOURCE SAYS: "in
     steady state, the thermal transport is governed by the Fuchs-Sondheimer theory and is insensitive to whether the
     boundaries are thermalizing or not. In contrast, under transient conditions, the thermal decay rates are significantly
     different for thermalizing and nonthermalizing walls"; "the thermalizing boundary condition is unphysical due to
     violation of heat flux conservation at the boundaries". (In-plane transport along films.)
   * **A. J. Minnich**, *Appl. Phys. Lett.* **107**, 183106 (2015), doi:10.1063/1.4935160 (**abstract only**:
     Fuchs–Sondheimer generalised to arbitrary dispersion; the isotropic equation "vastly overestimates the boundary
     scattering" in graphite). Cross-plane: Vermeersch, Carrete & Mingo 2016; Hua & Minnich 2015; Maassen & Lundstrom 2015,
     2016 (L3-4).

### Interpretation
* OUR DERIVATION (Rayleigh-quotient form). With C symmetric positive semidefinite, b orthogonal to its null space and
  y = (C^+)^{1/2} b: K = y^T y and N = y^T C^+ y, so tau_mem = y^T C^+ y / y^T y lies between 1/lambda_max and 1/lambda_min of
  C restricted to the complement of the conserved subspace. Bulk divergence requires lambda_min → 0 under mesh refinement
  together with enough spectral weight of b near it (the n ≥ 3/2 condition of L2). This links C3 to the smallest-eigenvalue
  analyses of the collision operator: the slowest relaxons of Cepellotti & Marzari 2016 (L1-6: in graphene they "converge —
  from above"), Hardy's criterion as applied by Huberman et al. 2025, and the temporal Knudsen number of Malviya &
  Ravichandran (L1-7).
* OUR DERIVATION (RTA film bound). With 1/tau_L = 1/tau + |v_z|/L, every mode satisfies v_z^2 tau_L^2 ≤ L^2 and
  v_z^2 tau_L ≤ |v_z| L. Each mode's contribution to N_z and K_z is therefore bounded uniformly in q, the region
  |q| < q_c contributes O(q_c^3 L^2) to N_z, and tau_z(L) converges at the rate of ordinary Brillouin-zone quadrature
  whatever the infrared law of tau. At RTA level the convergence half of C3 holds by construction; the numbers are the
  content. The full operator C_L = C + diag(|v_z|/L) has no such elementary bound (no regularisation of modes with
  v_z = 0, which C couples to the others), so its film convergence (E30) remains a NUMERICAL OBSERVATION.
* OUR INTERPRETATION: for L ≤ 1 μm the regulator also removes the Akhiezer-regime modes (L3), so the regularised quantity
  avoids the beyond-golden-rule physics that makes the bulk quantity ill-defined for the physical crystal (L2).

### What differs
* No source found reports the heat-flux relaxation time of a 3D crystal on a mesh sequence, nor contrasts it with a
  finite-size-regularised counterpart converged to < 0.1 %.
* Published first-principles values of this time are single-mesh (Beardo et al. 2021, Ge, 24^3; Huberman et al. 2025, Ge,
  48^3; kinetic-collective-model values for Si, L1-7).

### Caveats raised by the literature
* The |v_z|/L Matthiessen regulator is heuristic. Transient film responses depend on the wall condition (Ravichandran &
  Minnich 2016, in-plane geometry; OUR INTERPRETATION: the cross-plane memory time, a transient quantity, is exposed to the
  same dependence), and cross-plane kappa(L) is reproduced accurately only by calibrated suppression functions or full BTE
  solutions (L3-4). A film memory time should ultimately be the low-frequency (first-moment) coefficient of the film's
  response function under stated boundary conditions (L3).
* The bulk divergence is model-level: for the three-phonon golden-rule operator it remains UNRESOLVED; for the physical
  crystal divergence is expected from long-time-tail and Akhiezer physics (L2).

**Verdict C3: PARTIALLY KNOWN.** Bulk part: KNOWN in substance (L1, L2: follows from published power counting and from
the long-time-tail and Akhiezer literature; an analogous q-grid divergence of kappa is documented in 2D, Jiang et al.
2022). Film part: finite-size regularisation of infrared-divergent transport coefficients is KNOWN in general (low-
dimensional kappa(L), finite-box long-time tails, hydrodynamic finite-size scaling), and at RTA level the convergence is
guaranteed by construction (OUR DERIVATION). **The explicit ab initio paired demonstration for the heat-flux memory time of
a 3D crystal — bulk value still growing at 56x56x30 while the film value is converged to < 0.1 % (RTA) and to ±2–3 %
(full operator, four meshes) — is POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE (incremental convergence result).** No contrary
result was found.

---

## 6. Consequences for paper 1 (positioning)

1. **C3: convergence warning plus well-posed replacement.** Cite the L1/L2 sources for the bulk behaviour and Jiang et al.
   2022 for the documented two-dimensional analogue (q-grid divergence of kappa tied to the fc3 range). State the RTA bound
   v_z^2 tau_L^2 ≤ L^2 so that the convergence of the film quantity is seen to follow by construction, and present the
   full-operator convergence as numerical. Cite Ravichandran & Minnich 2016 together with the L3 film-solution sources for
   the boundary-condition caveat.
2. **C1: identifiability of a quantity of interest.** Frame with prediction-level identifiability (Kreutz et al. 2012;
   Gutenkunst et al. 2007; Kurniawan et al. 2026) and partial identification (Tamer 2010); cite Chen, Gamba, Li & Nair 2026
   for scale-dependent identifiability in phonon-transport inverse problems. Qualify every number as an inner bound on a
   discretised operator with synthetic data generated by that operator; report the multi-mesh sequence as a numerical
   observation whose limit is open; state that monotone narrowing is expected for a fixed objective.
3. **C2: application of established bounding methods.** Cite data collaboration (Frenklach et al. 2004; Feeley et al.
   2004; Frenklach et al. 2016; Hegde et al. 2018; Edwards et al. 2014), optimal UQ, profile-likelihood prediction
   intervals, the multiplier methods, L-BFGS-B and KS aggregation, robust formulations (Wang, Lazarov & Sigmund 2011),
   mesh-dependent optima (Sigmund & Petersson 1998) and the inverse crime (Wirgin 2004; Kaipio & Somersalo 2007). A
   defensible statement is "not previously applied to first-principles phonon transport, within the literature searched";
   "a new method" is not.

## 7. Unresolved prior-art uncertainty and limits of the search

1. Semantic Scholar was rate-limited (3/41 batch queries, 1 of 3 hand-run queries; no full-text snippet search). Full-text
   coverage relied on arXiv versions and open repositories.
2. Read at abstract or metadata level only: (a) every neighbouring-field source of sections 3–4 except Hegde et al. 2018
   (read in part); (b) the phonon and thermal-transport sources Lee & Gan 2017; Rohskopf et al. 2017; Sable et al. 2026;
   Kohl 1985; Sundara Raja & Jayarama Reddy 1977; Stoner & Maris 1993; Raya-Moreno et al. 2022; Beardo et al. 2022; Liu et
   al. 2026; Fehér & Kovács 2021; Chen, Gamba, Li & Nair 2026; Qin & Hu 2018; Majee & Akšamija 2016; Scalfi et al. 2023;
   Fiorentino et al. 2023; Guo et al. 2021; Yadav et al. 2025; Minnich 2015; Ravichandran & Minnich 2016. Powell 1969 was
   not retrieved. Full or targeted text was inspected for Huberman et al. 2025, Jiang et al. 2022, Ghosh et al. 2021, Abou
   Haibeh & Huberman 2026, Dragašević et al. 2026 and Shang et al. 2025 (arXiv or open-access versions).
3. Not searched or not extracted: econometric partial identification beyond Tamer 2010; robust-control model-invalidation
   literature beyond the data-collaboration lineage; climate-model tuning across resolutions (Hourdin et al., Bull. Am.
   Meteorol. Soc. 98, 589 (2017), doi:10.1175/BAMS-D-15-00135.1, retrieved but no statement could be extracted; not used);
   lattice-QCD continuum-limit fitting; Russian-language phonon kinetics and the monographs listed in the earlier protocol.
4. Alternate terminology was generated within this review, not obtained from separate mathematical or theoretical-physics
   consultation. Recommended for a follow-up search: "identified set" of a functional; "observability of predictions";
   Chebyshev–Markov–Krein and Stieltjes moment bounds, and the Kantorovich inequality, for the ratio N/K of inverse moments
   of the spectral measure of C relative to b (Rayleigh-quotient form, section 5); "infrared-safe observable"; "scale-aware
   parametrisation"; "discretisation-robust calibration".
5. Coverage ends on 14 October 2026; arXiv was searched through abstract fields, so works using unusual vocabulary may be
   missed.
6. Highest-consequence possible misses: (a) a materials-science application of feasible-set prediction bounds to phonon or
   thermal transport (none among forward citations of five framework papers); (b) a mesh study of the
   kinetic-collective-model relaxation time in the theses of the Barcelona group (Sendra thesis still not retrieved, as on
   11 October); (c) supplementary materials of Huberman et al. 2025 and Abou Haibeh & Huberman 2026 (arXiv versions only).

## 8. Bibliography (DOIs verified; access level as stated in the text; L-labels refer to `novelty-protocol.md`)

Abou Haibeh, J., Huberman, S. Mater. Res. Express 13, 085901 (2026), doi:10.1088/2053-1591/ae53a0; arXiv:2311.05107.
Alvarez, F. X., Jou, D. Appl. Phys. Lett. 90, 083109 (2007), doi:10.1063/1.2645110; J. Appl. Phys. 103, 094321 (2008), doi:10.1063/1.2913057.
Beardo, A. et al. Phys. Rev. Appl. 11, 034003 (2019), doi:10.1103/PhysRevApplied.11.034003.
Beardo, A. et al. Phys. Rev. B 101, 075303 (2020), doi:10.1103/PhysRevB.101.075303.
Beardo, A. et al. Sci. Adv. 7, eabg4677 (2021), doi:10.1126/sciadv.abg4677.
Beardo, A., Alajlouni, S., Sendra, L., Bafaluy, J., Ziabari, A., Xuan, Y., Camacho, J., Shakouri, A., Alvarez, F. X. Phys. Rev. B 105, 165303 (2022), doi:10.1103/PhysRevB.105.165303.
Borisov, I., Metelkin, E. PLoS Comput. Biol. 16, e1008495 (2020), doi:10.1371/journal.pcbi.1008495.
Byrd, R. H., Lu, P., Nocedal, J., Zhu, C. SIAM J. Sci. Comput. 16, 1190–1208 (1995), doi:10.1137/0916069.
Camacho de la Rosa, A., Esquivel-Sirvent, R., Becerril, D. J. Appl. Phys. 137, 155103 (2025), doi:10.1063/5.0257299.
Cepellotti, A., Marzari, N. Phys. Rev. X 6, 041013 (2016), doi:10.1103/PhysRevX.6.041013.
Chen, P., Gamba, I. M., Li, Q., Nair, A. Kinet. Relat. Models 24, 109–136 (2026), doi:10.3934/krm.2026021; arXiv:2512.24394.
Chen, P., Gamba, I. M., Li, Q., Wang, L. arXiv:2502.19533 (2025).
Dragašević, J., Rajkov, B., Simoncelli, M. Phys. Rev. Lett. 136, 186302 (2026), doi:10.1103/nbbn-56hr; arXiv:2303.12777.
Edwards, D. E., Zubarev, D. Yu., Packard, A., Lester, W. A., Frenklach, M. Phys. Rev. Lett. 112, 253003 (2014), doi:10.1103/PhysRevLett.112.253003.
Esfarjani, K., Chen, G., Stokes, H. T. Phys. Rev. B 84, 085204 (2011), doi:10.1103/PhysRevB.84.085204.
Feeley, R., Seiler, P., Packard, A., Frenklach, M. J. Phys. Chem. A 108, 9573–9583 (2004), doi:10.1021/jp047524w.
Fehér, A. É., Kovács, R. Int. J. Eng. Sci. 169, 103577 (2021), doi:10.1016/j.ijengsci.2021.103577.
Fiorentino, A., Pegolo, P., Baroni, S. npj Comput. Mater. 9, 157 (2023), doi:10.1038/s41524-023-01116-2.
Forghani, M., Hadjiconstantinou, N. G., Péraud, J.-P. M. Phys. Rev. B 94, 155439 (2016), doi:10.1103/PhysRevB.94.155439.
Forghani, M., Hadjiconstantinou, N. G. Phys. Rev. B 97, 195440 (2018), doi:10.1103/PhysRevB.97.195440; Int. J. Heat Mass Transf. 175, 121119 (2021), doi:10.1016/j.ijheatmasstransfer.2021.121119.
Frenklach, M., Packard, A., Seiler, P., Feeley, R. Int. J. Chem. Kinet. 36, 57–66 (2004), doi:10.1002/kin.10172.
Frenklach, M., Packard, A., García-Donato, G., Paulo, R., Sacks, J. SIAM/ASA J. Uncertain. Quantif. 4, 875–901 (2016), doi:10.1137/15M1019131.
Germano, M., Piomelli, U., Moin, P., Cabot, W. H. Phys. Fluids A 3, 1760–1765 (1991), doi:10.1063/1.857955.
Ghosh, K., Kusiak, A., Battaglia, J.-L. Phys. Rev. Mater. 5, 073605 (2021), doi:10.1103/PhysRevMaterials.5.073605; arXiv:2107.04731.
Guo, Y., Zhang, Z., Bescond, M., Xiong, S., Wang, M., Nomura, M., Volz, S. Phys. Rev. B 104, 075450 (2021), doi:10.1103/PhysRevB.104.075450.
Gutenkunst, R. N. et al. PLoS Comput. Biol. 3, e189 (2007), doi:10.1371/journal.pcbi.0030189.
Hegde, A., Li, W., Oreluk, J., Packard, A., Frenklach, M. SIAM/ASA J. Uncertain. Quantif. 6, 429–456 (2018), doi:10.1137/16M1110005; arXiv:1701.04695.
Hestenes, M. R. J. Optim. Theory Appl. 4, 303–320 (1969), doi:10.1007/BF00927673.
Hoque, M. S. B. et al. Appl. Phys. Lett. 125, 262201 (2024), doi:10.1063/5.0239769 (reference library; context only).
Hua, C., Minnich, A. J. J. Appl. Phys. 117, 175306 (2015), doi:10.1063/1.4919432.
Huberman, S., Zhang, C., Abou Haibeh, J. J. Phys.: Condens. Matter 37, 435601 (2025), doi:10.1088/1361-648X/ae0dd6; arXiv:2206.02769.
Jiang, J., Lu, S., Ouyang, Y., Chen, J. Nanomaterials 12, 2854 (2022), doi:10.3390/nano12162854.
Jiang, P., Lindsay, L., Koh, Y. K. J. Appl. Phys. 119, 245705 (2016), doi:10.1063/1.4954674.
Jin, J. S. J. Appl. Phys. 135, 045103 (2024), doi:10.1063/5.0187032 (Akhiezer damping; via L2).
Kaipio, J., Somersalo, E. J. Comput. Appl. Math. 198, 493–504 (2007), doi:10.1016/j.cam.2005.09.027.
Kohl, H. phys. status solidi (b) 130, 151–159 (1985), doi:10.1002/pssb.2221300112.
Kreisselmeier, G., Steinhauser, R. IFAC Proc. Vol. 12, 113–117 (1979), doi:10.1016/S1474-6670(17)65584-8.
Kreutz, C., Raue, A., Timmer, J. BMC Syst. Biol. 6, 120 (2012), doi:10.1186/1752-0509-6-120; arXiv:1107.0013.
Krewer, K. L. et al. Appl. Phys. Lett. 116, 102406 (2020), doi:10.1063/1.5142479.
Kurniawan, Y. et al. J. Chem. Phys. 156, 214103 (2022), doi:10.1063/5.0084988.
Kurniawan, Y., Neilsen, T. B., Francis, B. L., Stanković, A. M., Wen, M., Nikiforov, I., Tadmor, E. B., Bulatov, V. V., Lordi, V., Transtrum, M. K. Appl. Phys. Lett. 128, 064104 (2026), doi:10.1063/5.0296026.
Lassas, M., Saksman, E., Siltanen, S. Inverse Probl. Imaging 3, 87–122 (2009), doi:10.3934/ipi.2009.3.87.
Lee, C. H., Gan, C. K. Phys. Rev. B 96, 035105 (2017), doi:10.1103/PhysRevB.96.035105.
Lepri, S., Livi, R., Politi, A. Phys. Rep. 377, 1–80 (2003), doi:10.1016/S0370-1573(02)00558-6.
Li, W., Hegde, A., Oreluk, J., Packard, A., Frenklach, M. SIAM/ASA J. Uncertain. Quantif. 9, 231–259 (2021), doi:10.1137/19M1270185.
Li, Z., Kovachki, N., Azizzadenesheli, K., Liu, B., Bhattacharya, K., Stuart, A., Anandkumar, A. arXiv:2010.08895 (2020).
Lindsay, L., Broido, D. A., Reinecke, T. L. Phys. Rev. B 87, 165201 (2013), doi:10.1103/PhysRevB.87.165201 (ab initio kappa of AlN; context).
Liu, T., Chen, F., Li, J., Tong, C. Adv. Sci. Instrum. 1, 100037 (2026), doi:10.1016/j.asi.2026.100037.
Maassen, J., Lundstrom, M. J. Appl. Phys. 117, 135102 (2015), doi:10.1063/1.4916245; 119, 095102 (2016), doi:10.1063/1.4942836.
Maclaren, O. J., Nicholson, R. arXiv:1904.02826 (2019).
Majee, A. K., Akšamija, Z. Phys. Rev. B 93, 235423 (2016), doi:10.1103/PhysRevB.93.235423.
McGaughey, A. J. H. et al. (Phonon Olympics) J. Appl. Phys. 138, 135108 (2025), doi:10.1063/5.0289819.
Milanese, M., Vicino, A. Automatica 27, 997–1009 (1991), doi:10.1016/0005-1098(91)90134-N.
Minnich, A. J. Phys. Rev. Lett. 109, 205901 (2012), doi:10.1103/PhysRevLett.109.205901.
Minnich, A. J. Appl. Phys. Lett. 107, 183106 (2015), doi:10.1063/1.4935160.
Owhadi, H., Scovel, C., Sullivan, T. J., McKerns, M., Ortiz, M. SIAM Rev. 55, 271–345 (2013), doi:10.1137/10080782X; arXiv:1009.0679.
Parks, H. L., Kim, H.-Y., Viswanathan, V., McGaughey, A. J. H. Phys. Rev. Mater. 4, 083805 (2020), doi:10.1103/PhysRevMaterials.4.083805.
Perez, D., Subramanyam, A. P. A., Maliyov, I., Swinburne, T. D. npj Comput. Mater. 11, 263 (2025), doi:10.1038/s41524-025-01758-4.
Powell, M. J. D. "A method for nonlinear constraints in minimization problems", in Optimization, ed. R. Fletcher, Academic Press, London (1969), pp. 283–298 (not retrieved).
Qin, G., Hu, M. npj Comput. Mater. 4, 3 (2018), doi:10.1038/s41524-017-0058-3.
Raue, A. et al. Bioinformatics 25, 1923 (2009), doi:10.1093/bioinformatics/btp358.
Ravichandran, N. K., Minnich, A. J. Phys. Rev. B 93, 035314 (2016), doi:10.1103/PhysRevB.93.035314.
Raya-Moreno, M., Carrete, J., Cartoixà, X. Phys. Rev. B 106, 014308 (2022), doi:10.1103/PhysRevB.106.014308.
Rockafellar, R. T. J. Optim. Theory Appl. 12, 555–562 (1973), doi:10.1007/BF00934777.
Rohskopf, A., Seyf, H. R., Gordiz, K., Tadano, T., Henry, A. npj Comput. Mater. 3, 27 (2017), doi:10.1038/s41524-017-0026-y.
Romano, G., Johnson, S. G. Struct. Multidisc. Optim. 65, 297 (2022), doi:10.1007/s00158-022-03392-w.
Russi, T., Packard, A., Frenklach, M. Chem. Phys. Lett. 499, 1–8 (2010), doi:10.1016/j.cplett.2010.09.009.
Sable, A., Linjawi, B., Bradbury, K., Malof, J., Delaire, O. Digital Discovery 5, 1545–1557 (2026), doi:10.1039/D6DD00008H.
Scalfi, L., Vitali, D., Kiefer, H., Netz, R. R. J. Chem. Phys. 158, 191101 (2023), doi:10.1063/5.0151406.
Sendra, L., Beardo, A., Torres, P., Bafaluy, J., Alvarez, F. X., Camacho, J. Phys. Rev. B 103, L140301 (2021), doi:10.1103/PhysRevB.103.L140301.
Sendra, L. et al. Phys. Rev. B 106, 155301 (2022), doi:10.1103/PhysRevB.106.155301 (2D divergences of collective-limit coefficients; via L2).
Shang, W., Zhou, J., Panda, J. P., Xu, Z., Liu, Y., Du, P., Wang, J.-X., Luo, T. "JAX-BTE: a GPU-accelerated differentiable solver for phonon Boltzmann transport equations", npj Comput. Mater. 11, 129 (2025), doi:10.1038/s41524-025-01635-0; arXiv:2503.23657.
Sigmund, O., Petersson, J. Struct. Optim. 16, 68–75 (1998), doi:10.1007/BF01214002.
Silva, R. L. S. et al. arXiv:2604.01101 (2026).
Simoncelli, M., Marzari, N., Cepellotti, A. Phys. Rev. X 10, 011019 (2020), doi:10.1103/PhysRevX.10.011019 (thermal viscosity; via L1/L2).
Simpson, M. J., Maclaren, O. J. PLoS Comput. Biol. 19, e1011515 (2023), doi:10.1371/journal.pcbi.1011515.
Srivastava, Y., Gokhale, A. G., Jain, A. Phys. Rev. B 112, 195201 (2025), doi:10.1103/2n3p-hz62.
Stoner, R. J., Maris, H. J. Springer Ser. Solid-State Sci. (Phonon Scattering in Condensed Matter VII), 17–18 (1993), doi:10.1007/978-3-642-84888-9_5.
Suchan, K., Mohanty, S., Zhai, H., Cai, W. arXiv:2603.28683 (2026).
Sundara Raja, V., Jayarama Reddy, P. Pramana 8, 402–407 (1977), doi:10.1007/BF02862892.
Tamer, E. Annu. Rev. Econ. 2, 167 (2010), doi:10.1146/annurev.economics.050708.143401.
Torres, P. et al. Phys. Rev. B 95, 165407 (2017), doi:10.1103/PhysRevB.95.165407 (low-frequency treatment in first-principles KCM; via L1).
Tur-Prats, J., Han, Z., Beardo, A., Ruan, X., Alvarez, F. X. Nano Lett. 25, 11203–11209 (2025), doi:10.1021/acs.nanolett.5c00855.
Vermeersch, B., Carrete, J., Mingo, N., Shakouri, A. Phys. Rev. B 91, 085202 (2015), doi:10.1103/PhysRevB.91.085202.
Vermeersch, B., Carrete, J., Mingo, N. Appl. Phys. Lett. 108, 193104 (2016), doi:10.1063/1.4948968.
Wang, F., Lazarov, B. S., Sigmund, O. Struct. Multidisc. Optim. 43, 767–784 (2011), doi:10.1007/s00158-010-0602-y.
Ward, A., Broido, D. A. Phys. Rev. B 81, 085205 (2010), doi:10.1103/PhysRevB.81.085205 (first-principles relaxation-time forms; context).
Wirgin, A. "The inverse crime", arXiv:math-ph/0401050 (2004).
Xie, H., Gu, X., Bao, H. Comput. Mater. Sci. 138, 368–376 (2017), doi:10.1016/j.commatsci.2017.06.047; arXiv:1706.08232.
Yadav, U., Goutham, C. N. S., Meshram, K., Pathak, A., Bansal, D., Agrawal, A. Phys. Rev. Appl. 24, 034024 (2025), doi:10.1103/xvbk-75bf.
Zhu, C., Byrd, R. H., Lu, P., Nocedal, J. ACM Trans. Math. Softw. 23, 550–560 (1997), doi:10.1145/279232.279236.

---

## Appendix A. Hand-run search log (14 October 2026)

"Hit" = record judged relevant and followed up. Scripted batches are documented by their query files and logs.

**A.1 arXiv (interactive and API).**

| query or target | outcome |
|---|---|
| `abs:"phonon hydrodynamics" AND (uncertainty OR sensitivity OR identifiability)`; `abs:"force constants" AND (UQ OR identifiability OR Bayesian) AND abs:"thermal conductivity"`; automatic differentiation of kappa | timeouts / HTTP 429; re-run in the scripted batch (C1 P1, C2 P1) |
| `abs:"force constants" AND abs:"thermal conductivity" AND (inverse OR identifiability OR "constrained optimization")` | 11 records, none relevant |
| Guyer–Krumhansl / hydrodynamics / second sound AND film AND first principles | timeout, then empty API response; covered by the batch query of C1 P1 (0 records) and by L3 |
| `ti:"second sound in germanium"` | hit: 2206.02769 (text search: "mesh", "relaxation time", "macroscopic relaxation time") |
| Beardo et al. 2022 (silicon, 100–300 K) by title and authors | not on arXiv; abstract via OpenAlex |
| 2107.04731 (Ghosh et al.) | text search: "relaxation time", "grid", "propagation length" |
| 2311.05107 (Abou Haibeh & Huberman) | text search: "relaxation time" (no occurrence), "mesh", "q-point mesh" |
| 2303.12777 (Dragašević et al.) | text search: "converge" |
| 1701.04695 (Hegde et al.) | read sections 1–2; text search "prediction", "feasible set" |
| 2503.23657v2 (Shang et al.) | targeted reading of the HTML (inverse problems, parameters) |
| abstracts: 2604.01101, 1904.02826, 1706.08232, math-ph/0401050, 2010.08895, 2603.28683 | used as stated |

**A.2 Semantic Scholar.** "thickness-dependent heat flux relaxation time thin film first principles phonon Boltzmann
boundary scattering" (2015–): 4 records, none relevant. Two further queries (uncertainty of phonon hydrodynamics;
identifiability of relaxation times / fc3 from kappa data): HTTP 429.

**A.3 OpenAlex.** Abstracts and metadata for 80 DOIs (`abstracts/batch1`–`batch7`). Keyword-filtered forward citations
for 28 (DOI, keyword) pairs (`cites_batch1_out.md`, `cites_batch2_out.md`); hits from Sendra et al. 2021 citing records:
Huberman et al. 2025, Yadav et al. 2025, Beardo et al. 2022, Ghosh et al. 2021, Guo et al. 2021, Raya-Moreno et al. 2022;
from Kurniawan et al. 2022: Kurniawan et al. 2026, Perez et al. 2025.

**A.4 General web search** (research-filtered and unfiltered). UQ of phonon-hydrodynamic parameters (pointers: Parks et al.
2020, Tur-Prats et al. 2025); data collaboration (Frenklach, Packard, Seiler); optimal UQ; prediction profile likelihood;
Kaipio & Somersalo; Lassas, Saksman & Siltanen; Sigmund & Petersson; inverse determination of fc3 from kappa or linewidths
(forward methods only; no inverse method found); automatic differentiation of kappa with respect to force constants (none
found); differentiable phonon BTE (hits: Shang et al. 2025, Suchan et al. 2026); bounds on a derived transport property over
force constants (no relevant hit); "identifiability" with anharmonic force constants (hit: Lee & Gan 2017); heat-flux /
Cattaneo relaxation time with q-mesh convergence (pointer: Jiang et al. 2022); Wang, Lazarov & Sigmund; Germano et al.;
identification of the Cattaneo–Vernotte relaxation time from temperature data (Mochnacki & Paruch, J. Appl. Math. Comput.
Mech. 2013; peripheral, DOI not verified, not used).

**A.5 Publisher and repository pages.** Springer abstract pages of Sigmund & Petersson 1998 and Wang, Lazarov & Sigmund
2011; PubMed Central full text of Jiang et al. 2022 (targeted passages on the cutoff radius and q-grid); AMS page of
Hourdin et al. 2017 (statements on resolution dependence could not be extracted; not used).

**A.6 Crossref and DOI resolver.** 88 DOIs (`cited-dois.txt`, `crossref-verification.txt`): 86 Crossref records, 2 via the
DOI resolver with OpenAlex metadata.

**A.7 Project reference library.** Two semantic queries ("identifiability of thermal transport parameters ... constrained
optimization bounds"; "phonon hydrodynamics relaxation time thin film boundary scattering second sound first principles"):
the same four items as on 2 and 11 October.
