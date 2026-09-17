# L — Independent literature frontier and AlN data availability

Status: first-pass findings persisted 2026-09-16; bounded verification of recent preprints remains in progress. No other branch output has been read. This report makes no novelty claim.

## Research question and conclusion supported so far

The campaign asks what a conserving, anisotropic hydrodynamic closure for wurtzite AlN can determine from frequencies, velocities and total diagonal rates. Existing theory already derives collective relaxation modes, nonlocal conductivity, and viscosity from the **full collision operator**. The available Rao diagonal data do not provide that operator or an N/U split. Separate N/U lifetimes would improve a relaxation model, but published counterexamples show that a Callaway model is not automatically a quantitatively accurate substitute for the full operator. A unique AlN hydrodynamic window is not established by these sources or by the supplied 300 K arrays.

The strongest actionable result of this search is that one can obtain a **consistent replacement calculation**, without starting a new force-constant DFT calculation: the Phonon Olympics author repository contains independent AlN harmonic/cubic force constants and polar corrections. These inputs must not be mixed with Rao rates.

## Read scope and evidence categories

- SOURCE DERIVES means equations inspected in the actual paper/preprint; it does not mean this branch has audited the proof.
- SOURCE OBSERVES means a result from the source's numerical examples.
- OUR INTERPRETATION identifies consequences for this campaign.
- An abstract-only source is explicitly labelled and is not treated as an inspected derivation.
- Local context read: `00-question.md`, `these/theory/aln/README.md`, `these/theory/aln/sources.json`, and `these/notes/17_AlN_spectral_temperature.tex`.

## A. Closest theoretical prior art

### L1. Collision eigenmodes and the limits of diagonal lifetimes

Andrea Cepellotti and Nicola Marzari, **Thermal Transport in Crystals as a Kinetic Theory of Relaxons**, *Physical Review X* **6**, 041013 (2016), [DOI](https://doi.org/10.1103/PhysRevX.6.041013), [arXiv:1603.02608v3](https://arxiv.org/abs/1603.02608v3).

Read: actual arXiv HTML, introduction and Sections II–III, including Eqs. (1)–(7), with additional text through the beginning of the relaxon construction. SOURCE DERIVES: Bose-weight symmetrization of the linearized collision operator and its eigenmode representation. Off-diagonal couplings make the equilibration of a single phonon population generally non-exponential. A relaxon lifetime is the inverse collision eigenvalue, rather than the reciprocal of a diagonal phonon out-scattering rate. Its homogeneous steady conductivity representation has explicit assumptions; surfaces and localized sources invalidate the spatially homogeneous reduction. The paper credits Hardy for the earlier formal eigenmode construction.

**Correction alert:** an [erratum, PRX 10, 049901 (2020)](https://doi.org/10.1103/PhysRevX.10.049901) is linked by the publisher. Its full text has not yet been obtained. The 2016 arXiv text contains a statement about the equilibrium-associated zero eigenvector that should not be reused uncritically. The correct energy-null-vector construction is explicit in L2 Appendix A.1. No claim about the exact content of the erratum is made here.

### L2. Tensor viscosity and viscous heat equations

Michele Simoncelli, Nicola Marzari and Andrea Cepellotti, **Generalization of Fourier's Law into Viscous Heat Equations**, *Physical Review X* **10**, 011019 (2020), [DOI](https://doi.org/10.1103/PhysRevX.10.011019), [full arXiv text](https://arxiv.org/html/1906.09743v2).

Read: Sections I–III and targeted Appendix A/B passages, Eqs. (4)–(9), (27)–(53), and boundary discussion. SOURCE DERIVES: energy/drift variables, a fourth-rank thermal viscosity, collision-eigenvector formulas, and coupled temperature/drift equations for general dispersion. Odd/even collision sectors enter conductivity/viscosity respectively. SOURCE ASSUMES: near equilibrium, a phonon quasiparticle regime, gradient expansion, and an approximation discarding collision action on the local drift in Eq. (37), justified there by limiting high- and low-temperature arguments. SOURCE USES: no-slip drift conditions for examples, while acknowledging more general surface reflection is outside scope. OUR INTERPRETATION: this is decisive prior art for an anisotropic closure, but not a universal error theorem for arbitrary wavevector/frequency. It does not supply AlN coefficients.

### L3. Historical eigenmode and hydrodynamic formulations

- R. A. Guyer and J. A. Krumhansl, **Solution of the Linearized Phonon Boltzmann Equation**, *Physical Review* **148**, 766–778 (1966), [publisher](https://journals.aps.org/pr/abstract/10.1103/PhysRev.148.766). Actual publisher abstract inspected: normal-collision eigenvectors, wavevector/frequency-dependent conductivity, damped second sound and Poiseuille equations. Its simple current/zero-mode relation is stated for an isotropic dispersionless setting; full derivation not inspected here.
- Robert J. Hardy, **Phonon Boltzmann Equation and Second Sound in Solids**, *Physical Review B* **2**, 1193 (1970), [publisher](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.2.1193). Actual publisher abstract inspected: exact linearized collision eigenmode formulation and distinct drifting/driftless second-sound conditions. Thus neither spectral phonon dynamics nor a single rate-ordering criterion should be presented as new.

### L4. Space/time response and reduced-order dynamics

Chengyun Hua and Lucas Lindsay, **Space-time dependent thermal conductivity in nonlocal thermal transport**, *Physical Review B* **102**, 104310 (2020), [arXiv:2008.07596](https://arxiv.org/abs/2008.07596), [author manuscript](https://www.osti.gov/servlets/purl/1669763). Abstract inspected and full text retrieved for targeted reading. The source describes eigendecomposition of the full linearized collision matrix to obtain nonlocal conductivity and transient-grating response. This is direct prior art for finite-wavevector/frequency validation; detailed equations still to inspect in this pass.

Nikhil Malviya and Navaneetha K. Ravichandran, **Efficient calculation of phonon dynamics through a low-rank solution of the Boltzmann equation**, [arXiv:2502.00337v1](https://arxiv.org/abs/2502.00337v1) (2025 preprint; journal publication not verified here). Read: introduction and targeted Appendix B. SOURCE OBSERVES: about 3% of collision eigenmodes account for 99% of bulk conductivity in its natural-diamond 100 K example; proposes a low-rank dynamic method. OUR INTERPRETATION: a possible computational route after obtaining the operator, not evidence that AlN has the same low-rank structure.

### L5. Projection memory is established general machinery

- Robert Zwanzig, **Memory Effects in Irreversible Thermodynamics**, *Physical Review* **124**, 983–992 (1961), [publisher](https://journals.aps.org/pr/abstract/10.1103/PhysRev.124.983). Publisher abstract inspected: exact non-Markovian reduction, causal memory convolution, and low-frequency transport limit.
- Hazime Mori, **Transport, Collective Motion, and Brownian Motion**, *Progress of Theoretical Physics* **33**, 423–455 (1965), [publisher](https://academic.oup.com/ptp/article/33/3/423/1925580), DOI 10.1143/PTP.33.423. Publisher abstract and equation excerpt inspected; attempted full PDF access failed. The reduced equation contains a memory integral and unresolved-variable forcing.

OUR INTERPRETATION: finite-dimensional block elimination/Schur complements provide an elementary exact realization of this established idea. A memory-free local closure still requires a justified approximation. This search has not located an exact prior theorem with this campaign's particular fixed-diagonal, fixed-invariant identifiability hypotheses; a constructive counterexample would need its own audited mathematics, and absence from this bounded search is not novelty evidence.

## B. Contrary results and AlN-specific cautions

### L6. Direct wurtzite AlN benchmark of Callaway approximations

Jinlong Ma, Wu Li and Xiaobing Luo, **Examining the Callaway model for lattice thermal conductivity**, *Physical Review B* **90**, 035203 (2014), [publisher](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.035203), [author-uploaded full text](https://www.researchgate.net/publication/275027119_Examining_the_Callaway_model_for_lattice_thermal_conductivity).

Read: actual author manuscript, introductory/model discussion, Eq. (2), first-BZ convention, and conclusions; publisher abstract cross-checked. SOURCE OBSERVES: the calculated room-temperature AlN anisotropy is 5%, compared with 7% for RTA, 19% for original Callaway and 29% for Allen's modification. These percentages are properties of that calculation, not estimates for Rao's inputs. SOURCE EXPLAINS: N/U labels depend on the reciprocal-cell convention; it uses the first Wigner–Seitz Brillouin zone. Absorption/emission are distinct from N/U. Its low-frequency normal/Umklapp asymptotics are branch-specific results, not a licence to extrapolate missing Rao data or attach phenomenological constants. OUR INTERPRETATION: this is concrete AlN evidence against assuming that adding N/U scalar times validates anisotropic closure.

### L7. More recent Callaway failure tests

Nikhil Malviya and Navaneetha K. Ravichandran, **Callaway approximation to the Boltzmann equation fails to predict phonon hydrodynamics**, *Journal of Physics: Condensed Matter* **37**, 265701 (2025), [DOI](https://doi.org/10.1088/1361-648X/ade219), [OpenAlex record](https://openalex.org/W4411081276).

Access limitation: publisher fetch blocked; only the indexed abstract has been inspected. It reports first-principles comparisons for 20 semiconductors below 200 K and increasing disagreement in the hydrodynamic regime as U processes weaken, including transient-grating dynamics. No AlN-specific result or universal mathematical impossibility theorem is attributed to it. This is a high-priority full-text follow-up.

## C. Current preprints — provisional evidence only

### L8. Rate-only hydrodynamic indicators

Nikhil Malviya and Navaneetha K. Ravichandran, **Examining the microscopic origin of a computationally inexpensive thermal-conductivity-based indicator for phonon hydrodynamics**, [arXiv:2605.17947v2](https://arxiv.org/abs/2605.17947v2) (2026 preprint).

Read: actual arXiv HTML v2, abstract, introduction, Eqs. (1)–(2), concluding limitations and references. Version number confirmed by arXiv tool; exact version date still to verify. Its LPBE is driven by the full entropy-symmetrized collision matrix, whereas its RTA sets off-diagonals to zero. The authors report numerical correlations between full-BTE/RTA conductivity ratios and drift eigenmode signatures, and sensitivity to Brillouin-zone sampling. They explicitly limit the ratio to screening: it does not establish a particular hydrodynamic flow profile or second-sound response. Exact inspected statement from the conclusion: “approaches that utilize only the diagonal part — the phonon scattering rates, cannot reliably predict the strength of hydrodynamic transport.”

OUR INTERPRETATION: closely relevant caution, not a substitute for a general identifiability proof. No wurtzite-AlN prediction is claimed.

### L9. Generalized spectral Knudsen conditions

Nikhil Malviya and Navaneetha K. Ravichandran, **Bridging the continuum and the kinetic-Boltzmann theories of heat flow through generalized Knudsen numbers**, [arXiv:2606.17829v1](https://arxiv.org/abs/2606.17829v1) (submitted 16 June 2026, per arXiv discovery record).

Read so far: actual v1 abstract/introduction, Eqs. (1)–(4), and targeted later discussion. The authors propose collision-eigenspectrum-based conditions for reducing the LPBE to Fourier, weakly quasiballistic and hydrodynamic equations. Their computational examples comprise **22 cubic semiconductors and graphene**, as stated in the introduction. A whole-text search for `AlN` returned no match. This cannot be cited as a wurtzite-AlN validation. The detailed reduction equations remain to inspect before using its precise validity claims. Publication/review status not established; treat as a preprint.

## D. Reproducible AlN data and software route

### L10. Phonon Olympics benchmark and exact input locations

Alan J. H. McGaughey et al., **Phonon Olympics: Phonon property and lattice thermal conductivity benchmarking from open-source packages**, *Journal of Applied Physics* **138**, 135108 (2025), [DOI](https://doi.org/10.1063/5.0289819), [full author/publication PDF](https://mdr.nims.go.jp/filesets/3113ca3b-841a-4fe0-bc12-d2b28c9eea1c/download).

Read: actual PDF abstract, Sections II–IV and repository/data-availability discussion, pp. 1–5. Its model uses zero-temperature phonons with three-phonon/isotope scattering. Independent teams chose force-constant settings; their AlN inputs are not one shared calculation. The common-force-constant solver comparison concerns Ge. The paper stresses cubic-force-constant sensitivity and incomplete experimental agreement. Therefore cross-package agreement is a useful check, not validation of every higher moment or low-temperature eigenvalue.

The following remote contents were verified through GitHub MCP at commit **0640f07735059be9717a7565c2a0f22dc0da7a17**:

- [ShengBTE input directory](https://github.com/McGaughey-Lab/Phonon-Olympics/tree/0640f07735059be9717a7565c2a0f22dc0da7a17/Aluminum%20Nitride/ShengBTE/AlN_final%20result/5.thermal%20conductivity/input): `CONTROL` (2,152 bytes), `espresso.ifc2` (1,030,119 bytes), `FORCE_CONSTANTS_3RD` (1,359,454 bytes).
- `CONTROL` content actually read: four atoms; `ngrid=45 45 45`; `scell=6 6 6`; T=300 K; `scalebroad=1`; `espresso=.true.`; `isotopes=.false.`; **`convergence=.false.`**. It includes the dielectric tensor and Born effective charges. The supplied configuration is thus not itself an iterative-conductivity/convergence certificate; flags must be deliberately set and documented for a rerun.
- [phono3py input directory](https://github.com/McGaughey-Lab/Phonon-Olympics/tree/0640f07735059be9717a7565c2a0f22dc0da7a17/Aluminum%20Nitride/phono3py/AlN_kappa_input_files): `POSCAR` (744 bytes), `BORN` (411 bytes), `fc2.hdf5` (40,530 bytes), `fc3.hdf5` (25,581,870 bytes), `README` (100 bytes). README actually read: harmonic 5×5×3 and cubic 3×3×2 supercells.
- The phono3py parent directory also contains `Convergence_tests`, mode-scattering spreadsheets, and conductivity/accumulation outputs. Contents of those spreadsheets and convergence results were **not** audited here. Existence is not verification of convergence.

No precomputed AlN full collision operator or N/U arrays were identified in the specific inspected directories. This is a bounded inventory, not proof that none exists anywhere in the repository or literature.

### L11. Software export capabilities and normalization traps

Official current [phono3py command options](https://phonopy.github.io/phono3py/command-options.html), [input/output documentation](https://phonopy.github.io/phono3py/input-output-files.html), [direct LBTE solution](https://phonopy.github.io/phono3py/direct-solution.html), and [API reference](https://phonopy.github.io/phono3py/api-reference.html) were inspected.

- `--nu` writes separate `gamma_N`, `gamma_U`; the API describes this split for RTA output. `gamma` is a half-linewidth in ordinary THz; rates require both the factor two and conversion to angular frequency: **r = 4π gamma × 10^12 s^-1**. This differs from directly importing ShengBTE rates in ps^-1.
- `--lbte --write-collision` / `--read-collision` support collision-matrix files and distributed grid-point computation. The documented default symmetry-reduced object has dimension `N_ir × N_bands × 3`, not simply the scalar full-zone population dimension.
- OUR INTERPRETATION: before using an exported matrix for temperature, drift and even-parity stress dynamics, document its basis, weighting and symmetry sector and verify that it represents the required full scalar population operator. The conductivity-oriented reduced matrix cannot silently be relabelled as that operator. The appropriate full/reducible export/API path remains to verify.

### Proposed minimum next dataset (our specification)

Use either complete ShengBTE inputs or complete phono3py inputs from one calculation. At each temperature and converged mesh, export frequencies, polarizations, full-zone group velocities, full-zone weights, the reciprocal lattice and first-zone mapping, the separate normal/Umklapp/isotope/defect collision components, and either a fully documented matrix action or event-level data sufficient to reproduce it. Preserve both parity sectors. Record energy/momentum residuals, detailed-balance symmetry, positive-semidefinite checks, rate conventions, and all physical-model omissions. Converge conductivity, relevant inverse-operator moments, drift decay, even-sector viscosity, and the smallest relevant spectral rates independently. A scalar conductivity match is an insufficient validation target.

## E. Search and access trail

Local multi-index tool used: `C:/Users/Koussay/ResearchLab/software/research-tools/research-search.cmd`. Completed queries and raw JSON/Markdown records are under `C:/Users/Koussay/ResearchLab/literature/searches/20260916-aln-spectral-closure/`:

1. `phonon relaxons viscous heat equations` — Crossref and OpenAlex returned results; Semantic Scholar HTTP 429.
2. `aluminum nitride normal umklapp phonon hydrodynamic` — Crossref, Semantic Scholar and OpenAlex returned results.
3. `phonon collision operator inverse relaxation times memory projection` — Crossref/OpenAlex results; Semantic Scholar HTTP 429. Broad query was poorly selective and does not support an absence claim.
4. `Generalization of Fourier Law into Viscous Heat Equations` — Crossref and OpenAlex independently matched DOI 10.1103/PhysRevX.10.011019; Semantic Scholar HTTP 429.

Live web searches additionally used exact titles and synonyms: relaxons/collective phonon excitations; projection memory/Mori–Zwanzig; collision inverse/diagonal lifetime identification; normal/Umklapp; AlN hydrodynamics; Callaway accuracy/failure; first-zone convention; full collision export. Source pages and full-text read scope are recorded above.

Citation tracing: Semantic Scholar retrieved 30 citing and 30 referenced works for the viscous-heat-equations paper. This located L8–L9 and the low-rank/dynamic papers. This was a bounded page, not exhaustive forward/backward citation coverage. The arXiv papers' reference sections were used for further tracing.

Access limits: APS full text/harvest fetches for the relaxon erratum and Ma paper failed; the latter was recovered from an author-uploaded manuscript. The Mori PDF fetch failed. IOP's 2025 Callaway article was blocked. A fifth local scholarly query was rejected by automatic approval review because the service reported a usage limit; it did not run. Earlier shell access required escalation because sandbox initialization failed with `SetNamedSecurityInfoW ... failed: 5`. These tool failures must not be misreported as negative literature results.

## F. Unresolved prior-art and data questions

1. Read the relaxon erratum before repeating any disputed zero-mode statement from 2016.
2. Complete bounded equation-level inspection of L8/L9, check version dates, and distinguish numerical correlations from necessary/sufficient theorems.
3. Obtain the full 2025 Callaway paper and determine its exact material list and closure assumptions.
4. Verify a full-zone scalar collision export path, including even parity, for a fixed code version.
5. Inspect and reproduce the AlN force-constant convergence records; assess low-temperature acoustic resolution and higher-order scattering separately.
6. The exact inverse-identification proposition in this campaign remains a mathematical task, not a claim established merely by this literature search.
