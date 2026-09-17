# L2 — Bounded literature and phono3py source follow-up

Date: 2026-09-17. Read scope: frozen `branches/L-literature.md`, primary papers, official source code and bibliographic indexes. No red-team or discovery reports were read. This supplement makes no novelty claim and does not modify L.

## Findings established so far

1. The relaxon erratum has been obtained and read. It corrects an inequality in Appendix C, not a zero-mode statement.
2. The two exact arXiv version records and their submission timestamps have been verified.
3. phono3py v4.5.0 contains a full-mesh scalar, reducible collision-matrix route, but the inspected kernel implements Chaput's operator that is equivalent to the physical collision operator on odd populations. It does not supply the correct even-sector action merely by being reducible. The source also distinguishes raw rows, assembly and subsequent in-place eigensolver mutation. The exact formula/code comparison is recorded in Section 3.3.
4. No AlN collision calculation, convergence study, matrix-free implementation, or large allocation was run.

## 1. Relaxon erratum: actual content

Andrea Cepellotti and Nicola Marzari, **Erratum: Thermal Transport in Crystals as a Kinetic Theory of Relaxons [Phys. Rev. X 6, 041013 (2016)]**, *Physical Review X* **10**, 049901(E) (2020), published 26 October 2020, [DOI](https://doi.org/10.1103/PhysRevX.10.049901), [publisher full text](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevX.10.049901/fulltext).

SOURCE SAYS: Appendix C Eq. (C7) had its inequality reversed. The corrected expression is `1/k >= 1/k1 + 1/k2`. The authors state that the paper's discussion of relevance to Matthiessen's rule already follows this corrected inequality. The one-page erratum contains no correction concerning the equilibrium-associated zero vector.

Access: web-tool APS/Harvest fetches failed. A direct HTTP request to the publisher's Harvest endpoint succeeded (HTTP 200, application/pdf, 75,363 bytes); the complete one-page PDF was extracted and inspected. Crossref independently verified DOI, title, authors, publication date and publisher links. Consequently the caution in L about the 2016 zero-mode discussion remains a separate issue and must not be attributed to this erratum.

## 2. Exact recent-preprint versions and limited claim scope

### 2.1 arXiv:2605.17947v2

Nikhil Malviya and Navaneetha K. Ravichandran, **Examining the microscopic origin of a computationally inexpensive thermal-conductivity-based indicator for phonon hydrodynamics**. [Version record](https://arxiv.org/abs/2605.17947v2), [actual HTML](https://arxiv.org/html/2605.17947v2).

The arXiv submission history explicitly records v2 as **19 August 2026, 14:08:48 UTC**, following v1 on 18 May 2026. The HTML title block also displays 24 August 2026; that displayed document date must not replace the authoritative version-history timestamp.

SOURCE OBSERVES: conductivity-ratio correlations with collision-eigenmode hydrodynamic signatures and sensitivity to Brillouin-zone sampling. SOURCE USES: the full collision matrix for the LPBE numerator; this is not a diagonal-rate-only diagnostic. Main-text Eq. (1), Eq. (2), the effective RTA/Callaway matrices, and the numerical BZ-refinement discussion were inspected. The main text says its collision construction includes three-phonon, four-phonon and isotope processes, with details delegated to references and Supplementary Note 1. Those implementation details have not been independently audited here. Its correlations are numerical evidence for the studied systems, not a theorem that every fixed-diagonal collision family is identifiable or unidentifiable. No wurtzite-AlN window follows from this paper.

### 2.2 arXiv:2606.17829v1

The same authors, **Bridging the continuum and the kinetic-Boltzmann theories of heat flow through generalized Knudsen numbers**. [Version record](https://arxiv.org/abs/2606.17829v1), [actual HTML](https://arxiv.org/html/2606.17829v1).

The arXiv submission history records v1 as **16 June 2026, 11:56:23 UTC**. Its HTML title block also displays 24 August 2026; again use the submission history for version dating.

SOURCE DERIVES: a collision-eigenbasis equation, a split into odd/even sectors, elimination of the nonzero even modes, and a frequency-dependent odd-sector matrix Gamma [Eqs. (5)-(9)] preceding a temperature-response equation [Eq. (10)]. SOURCE ASSUMES for this construction: a symmetric positive-semidefinite collision operator, an energy null vector, parity symmetry, a source proportional to each mode's heat capacity, and then one-dimensional heating/transport. The numerical examples are 22 cubic semiconductors and graphene, as stated in the introduction; this is not wurtzite-AlN validation.

Exact algebraic elimination must be separated from the later asymptotic reductions. Eqs. (11)-(13) introduce an entrywise smallness criterion used for a diagonal approximation. This follow-up has not audited that criterion as a dimension-uniform matrix-inverse error bound, nor established necessity/sufficiency of every advertised regime condition. Treat the work as a preprint; peer-review/publication status was not established. Important equation-level caveats should be checked before importing its precise error or regime claims.

## 3. phono3py: fixed source snapshot and verified plumbing

Repository: [phonopy/phono3py](https://github.com/phonopy/phono3py). GitHub MCP verified tag **v4.5.0**, commit **21fa8f3817fbcc603254656f525bb5aec113afb6**. All numbered code pointers below are at that commit. The separately inspected default branch was c64b97ab0e7e42b29e3d363588ebc9363ec4e98b (15 September 2026); the `phonon3/collision_matrix.py` blob matched the tag snapshot. No installed local version was assumed.

### 3.1 Confirmed scalar route

- [`cui/phono3py_argparse.py`, lines 753-759](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/cui/phono3py_argparse.py#L753) registers `--reducible-colmat`; lines 933-939 register `--write-collision`.
- [`api_phono3py.py`, lines 2493-2496](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/api_phono3py.py#L2493) documents `is_LBTE=True, is_reducible_collision_matrix=True` as experimental construction/solution of the full collision matrix.
- [`phonon3/collision_matrix.py`](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/collision_matrix.py): `rot_grid_points=None` chooses the reducible path; each row has shape `(num_band0, num_mesh_points, num_band)`, without Cartesian-vector indices. The class explicitly stores the imaginary-self-energy/main diagonal separately from its other collision contribution.
- [`conductivity/collision_matrix_kernel.py`, lines 1141-1175](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/collision_matrix_kernel.py#L1141) allocates `(n_sigma, n_T, N_q, N_b0, N_q, N_b)` for the reducible kernel. With all bands, each temperature/broadening slice is an `(N_q N_b)` square scalar array. Lines 1363-1364 use unit full-mesh collision weights.

This verifies an actual scalar source path beyond the usual Cartesian, irreducible conductivity representation. **The formula check below finds that it remains an odd-equivalent conductivity operator, so changing to full-mesh scalar indices does not restore the physical even-sector action.**

### 3.2 Raw writes, assembly and mutation

- [`conductivity/lbte_collision_solver.py`, lines 123-178](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/lbte_collision_solver.py#L123): per grid point, `compute` returns `gamma` and copied `collision_row` separately, for each temperature and broadening.
- [`conductivity/lbte_init.py`, lines 425-446](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/lbte_init.py#L425) attaches a per-grid-point writer callback before finalization. [`conductivity/calculators.py`, lines 212-273](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/calculators.py#L212) orders the grid-point loop/callbacks before post-loop processing and finalization.
- [`conductivity/output.py`, lines 297-375](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/output.py#L297) writes separate `gamma`, optional `gamma_isotope`, and `collision_matrix` datasets. Reducible row addressing uses `bzg2grg`. [`file_IO.py`, lines 793-835](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/file_IO.py#L793) shows the HDF5 fields actually written. Grid maps, all basis conventions, and a declaration of raw/assembled state are not comprehensively encoded there and should be saved separately by a future exporter.
- [`collision_matrix_kernel.py`, lines 1228-1234](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/collision_matrix_kernel.py#L1228) expands rows by symmetry when requested, adds the main diagonal, then symmetrizes. Lines 1263-1303 expand by simultaneous permutations of both scalar grid indices. Lines 1336-1359 add the diagonal once. Lines 426-444 show that optional isotope and boundary rates enter this diagonal sum. Their existence does not verify a conserving full off-diagonal isotope operator.
- Lines 494-521 symmetrize by `(Omega + Omega.T)/2`; lines 539-563 then pass the same collision array into `diagonalize_collision_matrix` before solving for conductivity. Downstream solve routines use this storage as eigenvectors or an inverse (lines 648-708 and 1446-1466). A snapshot intended as a physical generator must be captured **after the required assembly, but before eigensolver mutation**. Re-running assembly on an already assembled snapshot risks adding the diagonal twice.
- [`lbte_init.py`, lines 448-456 and 545-567](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/lbte_init.py#L448) places the additional full-array write after `lbte.run` in the fresh-computation path. Therefore a filename or the `collision_matrix` dataset name alone does not establish raw-generator content; the exact lifecycle must be recorded. No ready-to-run exporter relying on this path is certified here.

### 3.3 Physical collision operator versus the odd-equivalent operator

Primary source inspected: Laurent Chaput, **A direct solution to the phonon Boltzmann equation**, [arXiv:1303.4062v1](https://arxiv.org/pdf/1303.4062), submitted 17 March 2013. The published work is **Direct Solution to the Linearized Phonon Boltzmann Equation**, *Physical Review Letters* **110**, 265506 (2013), [DOI](https://doi.org/10.1103/PhysRevLett.110.265506). The actual arXiv PDF, not an abstract, was inspected: printed p.3 defines physical Omega-prime immediately after Eq. (2); printed p.4, after Eq. (3), invokes oddness of the driven population and replaces it by Omega. The web PDF extraction positions are lines 104-173 and 215-277 respectively. The journal publisher abstract verifies publication/title/DOI. A subsequent attempt to fetch the published full PDF was rejected by automatic approval review due to a reported usage limit, so the formula comparison is explicitly against the inspected arXiv v1, not an unseen published-page equation.

**SOURCE SAYS:** because the heat-current response is odd under wavevector reversal, another collision matrix can give the same action on that response. The paper introduces Omega by changing the dummy wavevector in the collision sum. This replacement occurs **before** the later reduction to the irreducible Brillouin zone.

**OUR FORMULA TRANSCRIPTION AND ALGEBRA:** let lambda=(q,p), bar(lambda)=(-q,p) with a consistent time-reversed branch mapping, and let P act as `(Pf)_lambda=f_bar(lambda)`. Set D to the lifetime diagonal. Define the nonnegative channel weight

\[
 W_{\lambda\lambda'b}
 =\frac{\pi}{\hbar^2}
   \frac{|F_{\mathbf q,\mathbf q',\mathbf q_b}^{p p' p_b}|^2
   \Delta(\mathbf q+\mathbf q'+\mathbf q_b)}
   {\sinh(\hbar\omega_b/(2k_BT))},
\]

where Delta is one when its argument is a reciprocal vector and zero otherwise. With exact energy delta functions, define

\[
 A_{\lambda\lambda'}=\sum_b W_{\lambda\lambda'b}
          \delta(\omega_{\lambda'}+\omega_\lambda-\omega_b),
\]

\[
 B_{\lambda\lambda'}=\sum_b W_{\lambda\lambda'b}
  \left[\delta(\omega_{\lambda'}-\omega_\lambda+\omega_b)
       +\delta(\omega_{\lambda'}-\omega_\lambda-\omega_b)\right].
\]

In this notation the two displayed source formulas are

\[
 \boxed{\Omega'=D+A-BP},\qquad
 \boxed{\Omega=D+A+B}.
\]

The factor BP expresses the physical first term's wavevector `q-q'+q_b` and interaction `F(q,-q',q_b)`. Frequency equality under time reversal leaves its energy delta functions unchanged. Therefore

\[
 (\Omega-\Omega')f=B(I+P)f
 =\begin{cases}0,&Pf=-f,\\2Bf,&Pf=f.\end{cases}
\]

This is a direct algebraic comparison of the two source expressions. It establishes equality on the odd sector and shows why equality on a general even perturbation does not follow. It does not assume that every individual even vector has a nonzero discrepancy; special null directions or absent channels can be degenerate exceptions.

**SOURCE CODE MATCH:** at the fixed v4.5.0 commit, [`phonon3/collision_matrix.py`, lines 275-306](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/collision_matrix.py#L275), `_run_py_reducible_collision_matrix` explicitly identifies its result with the second term of Chaput's Omega. It multiplies squared interaction strengths by inverse sinh and `g[2]` with a positive sign. [`phonon3/triplets.py`, lines 240-252](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/triplets.py#L240) forms `g[2]=g0+g1+g2`, where the three Gaussian arguments are `omega0-omega1-omega2`, `omega0+omega1-omega2`, and `omega0-omega1+omega2`. The tetrahedron Python path does the same at lines 509-528. These are the three plus-sign channels of A+B. The separately added lifetime diagonal gives D+A+B. The C reducible kernel also accumulates the supplied collision integration weights with positive inverse-sinh/interaction factors ([`c/collision_matrix.c`, lines 249-284](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/c/collision_matrix.c#L249)). Rust dispatch is present, but the separate `phonors` backend source was not audited here; no claim of independent Rust verification is made.

**BOUNDED CONCLUSION:** the inspected phono3py reducible kernel is not a verified export of the physical population generator for energy/even-stress dynamics. Merely enabling `--reducible-colmat`, disabling the later symmetry reduction, or removing Cartesian indices does not undo the earlier odd-response replacement. Retain the separate physical channel signs and reversal maps in any future event-level construction. A conserving even-sector action must be independently implemented/obtained and checked; it is not established by a conductivity match.

This conclusion concerns the specific standard three-phonon collision path inspected at v4.5.0. It is **not** a repository-wide proof that no other full operator, custom implementation, or external matrix-free route exists. No general impossibility of reconstructing a physical operator from richer event-level inputs is asserted.

### 3.4 Normalization and lifecycle conclusions

Chaput's population variable is `f_lambda=sinh(hbar*omega_lambda/(2*kB*T))*delta n_lambda`. The entropy-weighted variable is `h_lambda=delta n_lambda/sqrt(nbar_lambda*(1+nbar_lambda))=2*f_lambda`; this uniform factor of two does not change the operator. It must be distinguished from the **matrix normalization** below.

- [`phonon3/imag_self_energy.py`, lines 200-207](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L200) defines the conversion to Gamma in ordinary THz; the reducible kernel uses that same conversion factor.
- [`collision_matrix_kernel.py`, lines 721-737](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/collision_matrix_kernel.py#L721) explicitly states that the collision matrix is half of Chaput's and divides the solved Y by two. [`conductivity/utils.py`, lines 319-328](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/utils.py#L319) separately contains the ordinary-to-angular-frequency factor of 2*pi in the conductivity conversion.
- The [official gamma documentation](https://phonopy.github.io/phono3py/input-output-files.html#gamma) gives `tau_ps=1/(4*pi*gamma)`. Together these source checks support multiplying the **assembled code matrix on its justified odd subspace** by `4*pi*10^12` to obtain s^-1. This scaling does not repair its even-sector action. No numerical normalization test was run in this follow-up.
- [`conductivity/utils.py`, lines 372-392](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/utils.py#L372) explicitly documents that `diagonalize_collision_matrix` overwrites the supplied matrix with eigenvectors. This confirms the lifecycle warning in Section 3.2 directly. Raw per-grid-point files, assembled matrices, eigenvectors, and pseudoinverses must be distinguished in a future exporter.
- Optional isotope and boundary contributions enter the inspected assembly as diagonal rates, not as verified energy-conserving full scattering operators. They must be tracked separately.

No public matrix-free operator API, runtime-tested exporter, or converged AlN dataset has been established here. The lower-level row/event machinery may support future development, but streaming the present odd-equivalent rows alone does not solve the even-sector issue. A small-mesh pilot would test implementation, normalization and invariant residuals; it would not establish a converged hydrodynamic window.

Resource limit: the PI reports 27.86 GiB RAM and 12 logical CPUs. A float64 scalar matrix for 12 bands on a 24^3 mesh requires about 205 GiB for one array, before working copies/eigensolver storage. It is not a runnable dense target on this machine. Authorized additional storage on D: does not remove the RAM limitation.

## 4. Prior-art terminology and search trail

Search vocabulary added: Schur complement; Feshbach or Feshbach-Schur projection/map; transfer-function or resolvent realization; matrix Stieltjes/spectral moment problem; positive-semidefinite matrix completion; correlation-matrix elliptope; prescribed nullspace; diagonal-plus-low-rank decomposition. These are terminology bridges, not novelty judgments.

A particularly close candidate was located: James Saunderson, Venkat Chandrasekaran, Pablo A. Parrilo and Alan S. Willsky, **Diagonal and Low-Rank Matrix Decompositions, Correlation Matrices, and Ellipsoid Fitting**, *SIAM Journal on Matrix Analysis and Applications* (2012), DOI [10.1137/120872516](https://doi.org/10.1137/120872516). Crossref matched the exact title and DOI via the local search. Its full theorem scope has not yet been inspected in this follow-up; do not cite it as a proof of this campaign's identifiability proposition.

Gene H. Golub and Gerard Meurant's **Matrices, moments and quadrature** is available as an [author-hosted primary manuscript](https://web.stanford.edu/class/cme335/sccm93-07.pdf). The indexed abstract describes bounds/approximations to bilinear matrix functions via quadrature and Lanczos. Full derivations have not yet been read here. This establishes a relevant search direction for resolvent moments, not a claim that the campaign's exact inverse problem has been settled.

Local tool: `C:/Users/Koussay/ResearchLab/software/research-tools/research-search.cmd`. Two bounded searches were saved under `C:/Users/Koussay/ResearchLab/literature/searches/20260917-aln-spectral-closure-followup/`:

1. `positive semidefinite matrix completion prescribed diagonal moments`, limit 5: Crossref and OpenAlex returned results; Semantic Scholar HTTP 429. Broad results were not treated as an absence test.
2. `Diagonal and low-rank matrix decompositions correlation matrices ellipsoid fitting`, limit 3: Crossref returned the close exact-title match. Other service outcomes are recorded in the raw search trail.

Additional live searches used exact erratum DOI/title, full collision/reducible code terminology, prescribed correlation-matrix nullspaces, Schur/Feshbach resolvents, and the Golub-Meurant title. No exhaustive forward-citation search or novelty audit was performed. Windows sandbox initialization failed with `SetNamedSecurityInfoW ... failed: 5`; successful shell reads/searches used escalation. Those failures are access limitations, not scientific evidence.

## 5. Final scope

Completed: erratum content; exact arXiv version timestamps; fixed phono3py tag/commit and lifecycle source trace; physical-versus-odd-equivalent formula comparison; source-level rate normalization. Unresolved: independent validation of that formula/code interpretation, alternative full physical operator paths, a conserving event-level implementation, full backend equivalence, and every AlN convergence/data question. The PI requested an independent physical/source check before adopting the code-path conclusion. No novelty claim, BTE run, large allocation, or unsupported command was introduced.
