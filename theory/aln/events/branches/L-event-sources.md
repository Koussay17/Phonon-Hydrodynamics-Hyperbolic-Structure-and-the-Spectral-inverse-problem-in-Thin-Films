# L — Sources for physical three-phonon events and the AlN import barrier

Campaign: 20260917-232509-aln-physical-events. Report saved 2026-09-18 after an interrupted bounded search. Scope read: `00-question.md`, the existing repository `theory/source_access.md`, primary papers and the specified public input directories. No other new branch report was read. No repository edits, BTE calculation or large material-data download was performed. No novelty claim is made.

## 1. Main primary source and exact inspected locations

Giorgia Fugallo, Michele Lazzeri, Lorenzo Paulatto and Francesco Mauri, **Ab initio variational approach for evaluating lattice thermal conductivity**, *Physical Review B* **88**, 045430 (2013), [DOI](https://doi.org/10.1103/PhysRevB.88.045430), [arXiv:1212.0470v2 PDF](https://arxiv.org/pdf/1212.0470v2), [HTML](https://arxiv.org/html/1212.0470v2).

SOURCE DERIVES: Eq. (3) uses population perturbations proportional to `nbar*(1+nbar)*f`; Eq. (4) has signed combinations `f_i+f_j-f_k`, with a factor 1/2 on the ordered decay sum. Eqs. (5)-(6) give equilibrium absorption/decay weights using cubic matrix elements, momentum selection, Bose factors and an energy delta. The paragraph after Eq. (13) invokes detailed balance **under exact resonance**. Appendix A, Eqs. (A1)-(A2) in the PDF (renumbered (31)-(32) in HTML), exhibits the positive contribution with block `(1,1,-1)*(1,1,-1)^T`. Appendix B, Eq. (B2) in the PDF / (36) in HTML, relates a modal lifetime to a sum over transitions; it is not an event list.

SOURCE ASSUMES / MODIFIES: Section IV, immediately after Eq. (30), states that replacing the energy delta by a Gaussian makes detailed balance approximate. Its symmetric matrix definition retains positivity at finite smearing by replacing one direction's equilibrium weight with the reverse direction's weight. This is an explicit discretization choice; it is not an exact nonlinear Bose Jacobian for arbitrary nonresonant triples. Appendix D studies mesh/smearing dependence for diamond, not AlN.

Access: actual HTML and the 10-page arXiv v2 PDF were obtained and inspected. The local Crossref search independently matched the exact title and DOI. A later shell request to independently retrieve additional Crossref metadata was rejected by automatic approval review due to a reported usage limit; it did not run. No unseen publisher-PDF equation is claimed as inspected.

### Additional primary coordinate-convention source

Andrea Cepellotti and Nicola Marzari, **Thermal Transport in Crystals as a Kinetic Theory of Relaxons**, *Physical Review X* **6**, 041013 (2016), [DOI](https://doi.org/10.1103/PhysRevX.6.041013), [arXiv:1603.02608v3](https://arxiv.org/html/1603.02608v3). Actual Eq. (5)-(6) and Appendix A Eq. (15)-(17) were inspected. SOURCE DEFINES: the entropy-weighted population as `delta n/sqrt(nbar*(1+nbar))`; the population-space collision matrix is `A*S^(-2)`, and the symmetric matrix is `S^(-1)*A*S^(-1)`. These formulas directly support the coordinate distinction below. Appendix A uses Fugallo's event construction, so this is an additional explicit convention source, not an independent derivation of all event weights. Its finite-Gaussian implementation is also explicitly discussed. The disputed zero-mode sentence elsewhere in the 2016 text is not used as evidence here.

## 2. Translation to this campaign's coordinates

The following is OUR TRANSLATION of the signed source structure into a finite event model, not a literal equation copied from the paper. Assume positive mode energies, finite positive temperature, an equal-weight full grid and event `i+j <-> k` with exact resonance `epsilon_i+epsilon_j=epsilon_k`.

Define `a_l=nbar_l*(1+nbar_l)`, `S=diag(sqrt(a_l))`, entropy coordinates `h=S^(-1)*delta n`, and the stoichiometric vector `s_e=e_i+e_j-e_k`. Let a **single reversible event** have the explicitly declared nonlinear convention

\[
 J_e(n)=\alpha_e\{n_i n_j(1+n_k)-(1+n_i)(1+n_j)n_k\},
 \qquad \dot n|_e=-s_e J_e(n),\qquad \alpha_e\geq0.
\]

All quantum-matrix-element, cell/grid normalization, quadrature, multiplicity and unit factors must already be contained in the declared alpha. This report assigns no absolute material value to it.

At exact Bose equilibrium, `F_e=nbar_i*nbar_j*(1+nbar_k)=(1+nbar_i)*(1+nbar_j)*nbar_k`. Direct differentiation gives

\[
 \delta J_e=\alpha_e F_e\,s_e^T S^{-2}\delta n,
 \quad b_e=S^{-1}s_e,
 \quad \dot h|_e=-\alpha_e F_e\,b_e b_e^T h.
\]

Thus the source's positive signed block becomes an entropy-coordinate outer product after congruence by `S^(-1)`. Fugallo's `f` variable and this campaign's `h` must not be identified: apart from the driving-gradient scalar, `delta n=S^2*f`, whereas `h=S*f`.

This also identifies a normalization test: a population-space Jacobian is `-S*Omega*S^(-1)`, not generally `-Omega`. A diagonal modal rate alone does not provide any of the off-diagonal signed products.

## 3. Counting and repeated indices

### 3.1 Ordered-source convention versus one canonical event

OUR COUNTING TRANSLATION of Fugallo Eq. (4): if `P_ij^k=P_ji^k` denotes the **ordered** equilibrium coefficient appearing there, an equivalent full ordered-pair sum is

\[
 A_3=\frac12\sum_{i,j,k}P_{ij}^k\,
 (e_i+e_j-e_k)(e_i+e_j-e_k)^T.
\]

For one canonical channel with the two inputs unordered, this becomes an event weight

\[
 w_{\{i,j\}\to k}=\frac{P_{ij}^k}{1+\delta_{ij}}.
\]

For `i != j`, the two ordered terms cancel the global half. For `i=j`, there is only one ordered term, so the half remains and the vector is `2*e_i-e_k`. This can be checked directly against the source row equations: the repeated input row receives `P_ii^k*(2*f_i-f_k)` and the output row receives `-(P_ii^k/2)*(2*f_i-f_k)`. Equivalently, `w=P/2` multiplying `(2,-1)*(2,-1)^T` gives these two rows.

The paper's Appendix A displays a block with three labelled rows and does not separately state this repeated-index convention. The above is a finite-index translation to be independently checked in the campaign, not a source theorem about a specific exporter. **Do not divide an exported weight by two again if its definition already counts identical inputs.** A second-quantized occupation-number transition probability, an ordered Boltzmann coefficient, and a canonical event weight are different objects until their conventions are mapped.

### 3.2 Other duplications that must remain distinct

- Forward and reverse scattering are already the two terms of one `J_e`; storing both as independent reversible events with the same weight doubles the operator.
- Exchange of the two input labels is a permutation, not another canonical event.
- Wavevector reversal produces a generally different time-reversed event. Pair it once with equal physical weight to preserve parity, except when it is already the same canonical event.
- If indices coincide, **add stoichiometric entries**. Never assign three array entries in a way that overwrites a repeated index.
- With all mode energies strictly positive and exact resonance, the output cannot equal either input; such equality would require the remaining input energy to vanish. Zero-frequency modes are outside this campaign's stated finite-positive-frequency model.

## 4. Exact resonance and broadened material grids

OUR ALGEBRAIC CHECK: for a triple with mismatch `Delta_e=epsilon_i+epsilon_j-epsilon_k`, the forward/reverse Bose product ratio is `exp(-Delta_e/(kB*T))`. Multiplying both channels by the same positive Gaussian weight therefore does not restore exact detailed balance when `Delta_e != 0`.

For an independently constructed positive outer-product model, let `u_E=S*epsilon`. Then `b_e^T*u_E=Delta_e`, and

\[
 u_E^T\Omega u_E=\sum_e w_e\,\Delta_e^2.
\]

Consequently a positive matrix built from nonresonant triples does not have the exact physical energy null vector merely because it is symmetric positive semidefinite. Pairing q with -q also leaves the energy mismatch unchanged. These observations distinguish three separate claims: positivity, parity, and exact energy conservation.

The source's finite-smearing symmetrization is a useful published example of this distinction. It must not be described as the unmodified nonlinear Bose collision Jacobian. Tetrahedron integration similarly represents continuum energy-surface integration; a nonzero quadrature weight attached to a grid triple is not, by itself, proof that the three stored grid frequencies are exactly resonant. No detailed tetrahedron conserving-event reconstruction has been verified here.

## 5. Concrete public AlN inputs rechecked without large downloads

GitHub MCP rechecked the author repository at fixed commit **0640f07735059be9717a7565c2a0f22dc0da7a17**. Related benchmark: Alan J. H. McGaughey et al., **Phonon Olympics: Phonon property and lattice thermal conductivity benchmarking from open-source packages**, *Journal of Applied Physics* **138**, 135108 (2025), [DOI](https://doi.org/10.1063/5.0289819), [public paper](https://mdr.nims.go.jp/filesets/3113ca3b-841a-4fe0-bc12-d2b28c9eea1c/download). This follow-up relies on the existing source-access record for the paper's detailed benchmark reading and independently rechecks the input inventory below.

### 5.1 phono3py input family

[Exact directory](https://github.com/McGaughey-Lab/Phonon-Olympics/tree/0640f07735059be9717a7565c2a0f22dc0da7a17/Aluminum%20Nitride/phono3py/AlN_kappa_input_files):

| File | Verified size | What it supplies |
|---|---:|---|
| POSCAR | 744 bytes | Crystal input |
| BORN | 411 bytes | Polar correction input |
| fc2.hdf5 | 40,530 bytes | Harmonic force constants |
| fc3.hdf5 | 25,581,870 bytes | Cubic force constants |
| README | 100 bytes | Harmonic 5x5x3 and cubic 3x3x2 supercells |

README text was actually read. File sizes and Git blob hashes were returned by the directory listing; the HDF5 payloads were not downloaded. Force constants are inputs from which event couplings can be recomputed; they are not a precomputed signed event table. The Git blob hashes are identifiers, not substitute SHA-256 file checksums.

### 5.2 ShengBTE input family

[Exact directory](https://github.com/McGaughey-Lab/Phonon-Olympics/tree/0640f07735059be9717a7565c2a0f22dc0da7a17/Aluminum%20Nitride/ShengBTE/AlN_final%20result/5.thermal%20conductivity/input): `CONTROL` (2,152 bytes), `espresso.ifc2` (1,030,119 bytes), and `FORCE_CONSTANTS_3RD` (1,359,454 bytes).

`CONTROL` was actually read. It specifies four atoms; 45x45x45 grid; 6x6x6 supercell; T=300 K; `scalebroad=1`; `espresso=.true.`; `isotopes=.false.`; **`convergence=.false.`**, with dielectric/Born data. It is a concrete alternative complete force-constant family, not a completed iterative BTE or event-export certificate. Do not run that production-size grid as the first importer pilot.

The two code-specific input families are independent calculations. Do not mix either family's couplings with Rao's diagonal rates. The existing source-access record identifies the Rao T300K dataset as diagonal spectral input with no N/U split, not an event table. No precomputed full signed event dataset was found in these two inspected input directories; that is a bounded inventory, not a claim about every file in the repository.

### 5.3 Concrete lower-level export candidate: interactions, not final event rates

Official [phono3py output documentation](https://phonopy.github.io/phono3py/input-output-files.html#pp-hdf5) identifies `pp-*.hdf5` as squared three-phonon interaction strengths in eV^2, with triplet and three band indices; its `--full-pp` option uses a simpler full array. This is a better starting object for an event importer than the conductivity matrix. The same documentation identifies `gamma_detail` as triplet contributions to **half-linewidths**, and shows that multiplying by stored triplet weights and summing recovers the modal self-energy. Those are distinct exported quantities.

A fixed source check at phono3py **v4.5.0**, commit **21fa8f3817fbcc603254656f525bb5aec113afb6**, confirms [`file_IO.py`, lines 1220-1262](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/file_IO.py#L1220): `write_pp_to_hdf5` supports `pp`, `triplet`, `weight`, `triplet_map`, and `triplet_all` fields in its full-array branch. Its alternative sparse branch stores `nonzero_pp`, shape and integration-zero masks, making recovery more involved. [`conductivity/utils.py`, lines 476-503](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/utils.py#L476) is the writer call site. Optional fields must be inspected in an actual pilot file; supported writer arguments are not a guarantee that every caller supplies every map.

The documented phonon file may include multiple translationally equivalent Brillouin-zone boundary points. Therefore its leading dimension is not necessarily the number of independent mesh modes. A future importer must use explicit grid maps and avoid counting boundary copies as distinct physical modes. No example file or backend was executed in this task, and no `pp` output for AlN was found in the two input directories inspected here.

## 6. Export barriers and minimum next step

An import requires an explicit contract for: full-grid mode identifiers and time-reversal map; positive frequencies and the handling of acoustic zero modes; all three event indices; absorption/decay orientation; reciprocal momentum vector; channel-resolved cubic coupling; quadrature/mesh and star multiplicities; energy mismatch; units and Bose factors; input-exchange/identical-mode multiplicity; whether reverse scattering is already combined; and whether the scalar is a bare coupling, equilibrium event flux or modal linewidth contribution.

A phono3py `gamma` or `gamma_detail` object must not simply be renamed `w_e`: linewidth contributions already combine occupation factors and channel/counting conventions. The inspected standard phono3py conductivity-matrix path is also not an established physical even-sector event export, as recorded in the pre-existing source-access record. No new general export capability is assumed here.

**Minimal practical next step:** choose the complete phono3py AlN force-constant family above, pin a code version, and create a tiny interaction-export pilot based on `pp` plus the full grid/phonon maps. Emit separately labelled energy channels before conductivity-specific rearrangement. First validate its counting and units against an explicitly enumerated synthetic resonant event and its modal self-energy contribution. For actual grid triples report mismatches and the chosen integration prescription; do not label broadened triples as exact conservative events. Only after that contract is validated should a small AlN slice be imported. This is a software/data validation task, not a converged material prediction.

Absolute normalization remains unresolved until the cubic-amplitude convention, hbar placement, ordinary versus angular frequency, delta-function units, grid normalization and combinatorial factors are reconciled with the pinned source. The common phono3py relation `1/tau_ps=4*pi*gamma_THz` converts modal half-linewidths; it does not by itself normalize a new event coefficient.

## 7. Search/access trail and remaining scope

Local multi-index command ran one exact-title-oriented query, `Ab initio variational approach evaluating lattice thermal conductivity Fugallo`, limit 3, with JSON/Markdown under `C:/Users/Koussay/ResearchLab/literature/searches/20260917-aln-physical-events/`. Crossref and OpenAlex responded; Semantic Scholar returned HTTP 429. The Crossref result matched Fugallo's title/DOI; the broad OpenAlex results were not independent exact matches and are not represented as such.

Live access inspected the actual Fugallo v2 HTML/PDF, the relaxon coordinate equations and Appendix A, official phono3py output documentation and two pinned writer files, official ShengBTE documentation landing page, and the versioned author-repository inventory/README/CONTROL. The direct DOI fetch for Fugallo failed in the web tool. A later shell metadata query was rejected by auto-review's usage-limit response before execution. None of those failures is negative scientific evidence.

Remaining: independent audit of the campaign's event model and repeated-index translation; a fully pinned source-to-event normalization/export mapping; actual AlN event generation; conservation-compatible treatment of continuum energy integration; and all material convergence tests. No rate magnitude, hydrodynamic window or experimental claim follows from this report.
