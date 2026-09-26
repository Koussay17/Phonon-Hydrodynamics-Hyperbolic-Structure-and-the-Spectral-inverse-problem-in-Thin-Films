# L: degenerate phonon collision action - primary literature first pass

2026-09-25. Independent bounded review of `00-question.md` and `02-assumptions.md`; no other campaign reports read. **Finding:** phonon matrix kinetics, exact-degeneracy secular kinetics, and exact aggregation have established but distinct assumptions. None of the inspected sources establishes that block sums of squared cubic amplitudes alone determine a closed block-population collision operator.

## 1. Phonon basis freedom and the collision approximation

**SOURCE SAYS.** M. Simoncelli, N. Marzari and F. Mauri, *Wigner Formulation of Thermal Transport in Solids*, Phys. Rev. X **12**, 041011 (2022), [DOI](https://doi.org/10.1103/PhysRevX.12.041011); inspected [arXiv:2112.06897v2](https://arxiv.org/html/2112.06897v2) and the [published PDF](https://iris.uniroma1.it/retrieve/f025069e-421a-49cd-b717-9c5aeda52fd7/Simoncelli_Wigner-formulation_2022.pdf). Sections IV B and VI, Eq. (35), Eq. (50), and Appendix A2 discuss matrix distributions, phase conventions and eigenvector freedom. Inside an exactly degenerate eigenspace, one Cartesian velocity component can be diagonalized by rotating eigenvectors; different components generally cannot be diagonalized simultaneously. The smooth/steplike Fourier-convention issue is distinct from arbitrary rotations within a degenerate eigenspace.

**SOURCE ASSUMES.** Their collision Eq. (38), pp. 10-11, uses cubic anharmonicity, perturbative mass disorder, linearization around equilibrium and neglect of energy renormalization. It has population repumping/depopulation and off-diagonal decay `-(Gamma_s+Gamma_s') n_ss'/2`, without population-coherence scattering transfer. The stated quasiparticle condition compares linewidth with mode frequency. This particular formula is not evidence for a universally covariant collision closure under arbitrary degenerate rotations. Here `Gamma=1/tau`; do not identify it directly with an exported half-linewidth.

The predecessor, the same authors' [*Unified theory of thermal transport in crystals and disordered solids*, arXiv:1901.01964v1 (2019)](https://arxiv.org/html/1901.01964v1), explicitly describes Eq. (11) as a closed-form approximation of many-body collisions using one-body quantities. Its Wigner Eq. (9) evolves a matrix, including coherences. A scalar population BTE is therefore an additional representation/closure choice, not the complete matrix dynamics merely written differently.

## 2. Exact degeneracy does not remove intra-block coherences

**SOURCE SAYS AND DERIVES.** A. Trushechkin, *Unified Gorini-Kossakowski-Lindblad-Sudarshan quantum master equation beyond the secular approximation*, Phys. Rev. A **103**, 062226 (2021), [DOI](https://doi.org/10.1103/PhysRevA.103.062226); inspected [arXiv:2103.12042v3](https://arxiv.org/html/2103.12042v3), Sections II-III. Equation (4) constructs transition operators with full energy projectors. Secular averaging removes terms of unequal Bohr frequencies, not distinct basis states having the same energy. Section III.2 explicitly retains transfers between populations and zero-Bohr-frequency coherences inside degenerate eigenspaces. Nearly equal frequencies require an additional timescale argument: Eq. (7) uses `H_S=H_S^(0)+lambda^2 delta H_S`, and Eqs. (9)-(13) retain within-cluster terms.

**SCOPE LIMIT.** This is a weak-coupling system-bath construction with specified interaction operators and bath correlations. It does not supply a microscopic three-phonon collision generator for the present finite toy domain. Its relevance is precise: neither exact degeneracy nor the word "secular" justifies deleting all intra-block coherence. The phonon approximation above and this general secular structure must not be conflated.

## 3. Exact closure has a separate algebraic criterion

Known terminology includes **exact linear lumping**, **invariant observable subspaces**, and, for Markov chains, **ordinary lumpability**. G. Li and H. Rabitz, *A general analysis of exact lumping in chemical kinetics*, Chem. Eng. Sci. **44**, 1413-1430 (1989), [DOI/publisher abstract](https://doi.org/10.1016/0009-2509(89)85014-6), identifies invariant Jacobian subspaces. Its full text was not inspected. A directly inspected proof is in A. Ovchinnikov, I. Perez Verona, G. Pogudin and M. Tribastone, [*CLUE: Exact maximal reduction of kinetic models by constrained lumping of differential equations*, arXiv:2004.11961v2](https://arxiv.org/html/2004.11961v2), Section 2 and Supplementary Proposition II.1: the row space of the lumping map must be invariant under every relevant Jacobian.

**APPLICATION, not a new theorem:** for finite constant `dot(n)=-L n`, block totals `y=M n` obey `dot(y)=-Lbar y` for every initial `n` precisely when

\[
M L=\bar L M.
\]

With any lifting `R` satisfying `M R=I`, this is equivalent to `M L (I-R M)=0`; only then is `Lbar=M L R` exact for all initial states. Merely choosing equal occupations within each block supplies an ansatz, not this closure condition. A block total and an average differ by block dimension, which must be included in `M` and `R`.

P. Buchholz, *Exact and ordinary lumpability in finite Markov chains*, J. Appl. Probab. **31**, 59-75 (1994), [DOI](https://doi.org/10.2307/3215235), is a terminology/reference pointer only here: metadata and abstract inspected, not its theorem text. A phonon population collision matrix need not itself be a Markov probability generator, so the general linear criterion is the safer formulation.

## 4. What remains unresolved and the next importer requirement

**OUR INTERPRETATION:** `N'=U^dagger N U` preserves a degenerate block trace but generally mixes diagonal populations and off-diagonal coherences. For a complete Cartesian tensor block, `sum |V_ijk|^2` is its unitary-invariant squared norm. That scalar invariance neither identifies all rotated component rates nor tests dynamical closure. Generic rotated squared amplitudes require interference cross-products missing from componentwise `|V|^2`.

The importer should preserve basis/projector metadata and distinguish fixed-basis populations, block totals/averages, and any matrix kinetic model. Before calling a block operator exact, specify the collision generator and test the intertwining condition above. A matrix model additionally needs complex amplitudes or sufficient cross-product data and a stated secular/bath approximation. Which such data suffice for the actual AlN exporter remains unresolved; no material conclusion follows from this review.

## Search scope and access

`research-search.cmd` queries: "phonon Wigner transport degenerate gauge covariance density matrix"; "secular approximation degenerate energy levels Davies master equation"; "exact ordinary lumpability finite Markov chains Buchholz"; and the exact Wigner title. Crossref and OpenAlex returned results; Semantic Scholar returned HTTP 429 throughout. A later direct Crossref DOI batch also hit 429. Local query outputs persist under `C:/Users/Koussay/ResearchLab/literature/searches/20260925-degenerate-collision-action/`. Live search plus primary-paper inspection supplied the equation checks above; title/DOI records were checked against available Crossref or publisher/author records. Several publisher/PDF fetches failed, with accessible arXiv versions used as stated. No exhaustive citation-forward, earliest-priority or novelty search was performed; no novelty claim is made.


## Supplement: alternative formulations and their scope

Second targeted pass, 2026-09-25. The following maps the supplied alternative terms to established constructions; it adds no novelty or material claim.

**Conditional expectation onto the center.** D. Ranard, M. Walter and F. Witteveen, *A Converse to Lieb-Robinson Bounds in One Dimension Using Index Theory*, Ann. Henri Poincare (2022), [DOI](https://doi.org/10.1007/s00023-022-01193-x), [inspected paper, Eq. (2.2) and text preceding Eq. (4.5)](https://d-nb.info/1269764101/34), gives Haar averaging over an algebra's unitaries as a conditional expectation onto its commutant. Its standard specialization here is

\[
\mathcal A=\bigoplus_a B(H_a),\quad
E_Z(X)=\int_{\prod_a U(H_a)}UXU^\dagger\,dU
      =\sum_a\frac{\operatorname{Tr}(P_aX)}{d_a}P_a,
\quad d_a=\dim H_a.
\]

The range is `Z(A)=span{P_a}`. This retains block traces and replaces each block by its average. In contrast, the pinching map `X -> sum_a P_a X P_a` retains full matrices inside blocks, including their coherences. Identifying these two projections would discard extra information. The source's quantum-lattice application is not being imported into phonon kinetics.

**Passive covariance versus physical isotropy.** In representation-independent notation, if `T_U(X)=U^dagger X U`, a passive change of coordinates sends a generator to `L'=T_U L T_U^(-1)`. Requiring the *same* generator to obey `L T_U=T_U L` is an additional physical symmetry. If that symmetry holds for the entire block-unitary group, Haar averaging gives `L E_Z=E_Z L`, a sufficient closure condition. Degeneracy of the free Hamiltonian alone does not impose this symmetry on interactions. This is an elementary specialization, not a source claim about AlN. [Trushechkin, Sec. III.5, Eq. (15)](https://arxiv.org/html/2103.12042v3) verifies covariance under the one-parameter group generated by `H_S^(0)`; that is weaker than arbitrary block-unitary isotropy. His Sec. III.2 retains zero-frequency coherences, as recorded above.

**Lumpability/intertwining.** The [CLUE criterion, Supplementary Proposition II.1](https://arxiv.org/html/2004.11961v2) already supplies the required general linear formulation. In finite dimensions, `ML=Lbar M` is equivalently the semigroup intertwining `M exp(-tL)=exp(-tLbar) M` for every `t`. This specifies closure for every initial condition, the property intended by the campaign's "strong lumpability" wording. Markov-chain nomenclature depends on convention; Buchholz's distinction between ordinary and exact lumpability should not be replaced by an unchecked synonym. No Markov-generator property is assumed for the phonon operator.

**Mori-Zwanzig and Schur-complement memory.** Q. Du, B. Engquist and X. Tian's author review, [*Multiscale Modeling, Homogenization and Nonlocal Effects: Mathematical and Computational Issues*, arXiv:1909.00708v1, Secs. 3.1 and 3.4, Eq. (3.12)](https://arxiv.org/html/1909.00708v1), explicitly connects Schur elimination to the exact projected memory equation. In the present sign convention, split `dot(x)=-Lx`, `L=[[A,B],[C,D]]`, `x=(u,v)`. The familiar finite-dimensional identity reads

\[
\dot u=-Au+\int_0^t B e^{-D(t-s)}C u(s)\,ds-B e^{-Dt}v(0).
\]

Its resolvent Schur operator is `zI+A-B(zI+D)^(-1)C`, when the inverse exists. Thus a projected instantaneous matrix alone omits both memory and unresolved initial data; a vanishing memory kernel alone need not remove the latter. This is the cited identity in our notation, without a Markov or timescale approximation.

Historical primary pointer: R. Zwanzig, [*Memory Effects in Irreversible Thermodynamics*, Phys. Rev. **124**, 983 (1961)](https://doi.org/10.1103/PhysRev.124.983). Its publisher abstract supports the convolution-memory lineage. The [1972 correction by Zwanzig, K. S. J. Nordholm and W. C. Mitchell, Phys. Rev. A **5**, 2680](https://doi.org/10.1103/PhysRevA.5.2680) concerns neglected fluctuations in the nonlinear transport-equation derivation. Both publisher abstracts were inspected; their full derivations were not. The explicit finite linear identity above was checked in Du et al., rather than inferred from the historical abstracts.

**Added search limits.** Two paired scholarly/live searches covered "conditional expectation twirling commutant quantum symmetry" and "Mori Zwanzig Schur complement memory kernel projection"; one live exact-title search checked Zwanzig's 1961 source and revealed its correction. OpenAlex worked; Semantic Scholar returned 429, and Crossref returned 429 for the memory query. New local search records are in the existing trail's `supplement-conditional-expectation/` and `supplement-memory/` folders. No broader priority search was performed. This pass preserves the unresolved need to specify and validate the actual phonon collision generator.
