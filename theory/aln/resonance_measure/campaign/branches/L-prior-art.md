# Prior art: conserving resonance integration

Independent first pass; sources inspected 2026-09-27, persisted 2026-09-28. Read campaign question/assumptions only; no new peer reports. Bounded methodological search, no novelty or material claim.

## 1. Phonon quadrature and finite-width detailed balance

**SOURCE SAYS.** A. Togo, L. Chaput and I. Tanaka, *Distributions of phonon lifetimes in Brillouin zones*, Phys. Rev. B **91**, 094306 (2015), [DOI](https://doi.org/10.1103/PhysRevB.91.094306), inspected [1501.00691v3, Appendix C](https://arxiv.org/html/1501.00691v3). Momentum elimination reduces the double Brillouin-zone sum. Linear tetrahedron integration replaces delta functions by weights calculated from mesh frequencies and assigned to vertices. It removes the freely selected smearing width; it is not presented as a theorem of exact nonlinear finite-mesh Bose equilibrium or energy conservation.

**SOURCE SAYS.** G. Fugallo, M. Lazzeri, L. Paulatto and F. Mauri, *Ab initio variational approach for evaluating lattice thermal conductivity*, Phys. Rev. B **88**, 045430 (2013), [DOI](https://doi.org/10.1103/PhysRevB.88.045430), inspected [1212.0470v2, section IV after (30)](https://arxiv.org/html/1212.0470v2). The authors explicitly state that Gaussian replacement makes detailed balance approximate. Their symmetrized rate construction retains a symmetric nonnegative linearized matrix at finite width, while differing from the Omini-Sparavigna construction. This is direct prior evidence that matrix positivity and exact detailed balance are separate requirements.

**PINNED CODE.** Installed phono3py 4.5.0, upstream commit `21fa8f3817fbcc603254656f525bb5aec113afb6`: [triplets.py L218-L271](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/triplets.py#L218) selects Gaussian or tetrahedron integration. L509-L526 passes daughter-vertex frequency sums/differences to the tetrahedron engine. [imag_self_energy.py L775-L793](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L775) multiplies those weights by Bose factors and interaction strengths at the labelled mesh triplets. Thus a nonzero tetrahedron weight is not a certificate that that mesh tuple itself is resonant.

Normalization matters: [phonon/func.py L44-L46](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon/func.py#L44) uses $e^{-x^2/(2\sigma^2)}/(\sqrt{2\pi}\sigma)$, whereas Fugallo (30) uses $e^{-(\Delta E/\sigma)^2}/(\sqrt\pi\sigma)$. The widths differ even after converting ordinary THz to energy. Installed phonopy 4.5.0 `phonon/tetrahedron_method.py:256-278` confirms `run(..., value='I')` denotes delta-function integration. Native backend equivalence and full-operator conservation were not audited here.

## 2. Coarea, root quadrature, and singular resonance measures

**SOURCE ANALYZES.** Y.-K. Shi and G. L. Eyink, *Resonance Van Hove singularities in wave kinetics*, Physica D **332**, 55-72 (2016), [DOI](https://doi.org/10.1016/j.physd.2016.05.014), inspected [1507.08320v1](https://arxiv.org/html/1507.08320v1). Equation (3) writes the regular resonance measure as surface area divided by the difference of daughter group velocities. Section 3 uses the Morse normal form to analyze nondegenerate saddle points: local radial behavior is proportional to $\int_0^\eta r^{D-3}dr$, finite for $D>2$ and logarithmically divergent for $D=2$. Degenerate critical sets and cancellations by interaction/occupation factors require additional analysis. These are geometric finiteness results and physical examples, not a universal quadrature convergence theorem.

Equations (38)-(41) distinguish finite-time sinc-squared and Lorentzian regularizations from their limiting measures; the detailed measure construction is deferred to another reference. Do not apply the regular coarea denominator at a critical root or infer that every critical point has the same behavior. In particular, isolated extrema should be checked with the declared limiting prescription.

**SOURCE CONSTRUCTS/TESTS.** R. I. Saye, *High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles*, SIAM J. Sci. Comput. **37**, A993-A1019 (2015), [DOI](https://doi.org/10.1137/140966290), inspected [author's full PDF](https://math.lbl.gov/~saye/96629.pdf), sections 1 and 3. It converts implicit surfaces into local height graphs, recursively combining one-dimensional root finding and Gaussian quadrature. Weights are positive and surface nodes lie on the implicit surface; numerical tests demonstrate high orders on smooth examples. This establishes relevant quadrature prior art. It does not establish conservation for a phonon discretization or regularize a divergent resonance density.

**OUR MAPPING.** For a regular energy mismatch $\Delta(q)$,

\[
\int F(q)\delta(\Delta(q))\,dq
=\int_{\Delta=0}\frac{F}{|\nabla\Delta|}\,dS.
\]

A surface rule must integrate $F/|\nabla\Delta|$, not merely $F$. Alternatively, a graph/root rule carries $1/|\partial_{q_j}\Delta|$. Positivity survives this division only where the denominator is finite and nonzero. Root tolerances, missing components, chart changes, periodic representatives, and critical points remain numerical obligations.

## 3. Entropic and conservative finite representations

**SOURCE PROVES A LIMITED PACKAGE.** Z. Cai, Y. Fan and L. Ying, *An Entropic Fourier Method for the Boltzmann Equation*, SIAM J. Sci. Comput. **40**, A2858-A2882 (2018), [DOI](https://doi.org/10.1137/17M1127041), inspected [1704.07369v2](https://arxiv.org/html/1704.07369v2). Theorem 2 proves semidiscrete mass conservation, nonnegativity and a discrete H-theorem. Section 5 explicitly states that momentum and energy conservation are lost. Section 1.1, equations (12)-(16), shows why reversible nonnegative discrete-velocity collisions that preserve the invariants have an H-theorem. These are classical binary-collision results, not Bose three-phonon results.

**SOURCE DERIVES THE GALERKIN PRINCIPLE.** M. R. A. Abdelmalik and E. H. van Brummelen, *Moment Closure Approximations of the Boltzmann Equation Based on phi-Divergences*, J. Stat. Phys. **164**, 77-104 (2016), [full primary paper](https://doi.org/10.1007/s10955-016-1529-5). Sections 3.2-3.3, especially (18)-(23), express entropy-based closure as a Galerkin method in a transformed distribution variable. Including collision invariants in the test space yields conservation; placing the entropy derivative in that space allows testing the equation by it to obtain dissipation. The paper states admissibility/realizability conditions. Its generalized phi-entropy claims assume the collision operator dissipates that entropy; they are not automatic for an arbitrary operator or quadrature.

**INDEPENDENT CONSERVATION CONSTRUCTION.** I. M. Gamba and S. Rjasanow, *Galerkin-Petrov approach for the Boltzmann equation*, J. Comput. Phys. **366**, 341-365 (2018), [DOI](https://doi.org/10.1016/j.jcp.2018.04.017), inspected [1710.05903v1, test-space construction and section 6](https://arxiv.org/html/1710.05903v1). Polynomial test functions contain all five classical collision invariants, producing conservation to machine accuracy. This is useful prior art for preserving invariants by the weak form, not evidence of a Bose H-theorem.

## 4. Consequence for this campaign

**INTERPRETATION, NOT A LITERATURE THEOREM FOR THE PROPOSED SCHEME.** Combine resonance quadrature with a compatible weak representation; they are separate design choices. For a positive quadratic event form, the energy test contributes $\sum_e w_e\Delta_e^2$. Positivity alone cannot make this vanish when weighted events are off shell. At Bose equilibrium, the absorption/decay occupation-factor ratio is $e^{\beta\Delta_e}$, so an off-shell tuple is not individually balanced.

A candidate should declare one energy representation, evaluate occupations and entropy variables consistently at its resonant quadrature points, and preserve its invariant through the same weak form. Conservation of interpolated energy is not automatically conservation of the original dispersion. Classical $f\log f$ and $|v|^2$ must be replaced by Bose entropy and the actual crystal energy; classical particle-number and momentum invariants cannot simply be imported into three-phonon Umklapp kinetics. The inspected literature supports these ingredients but does not prove convergence of their proposed combination. No universal relation between width and mesh size was established.

## Search and access limits

Local `research-search.cmd` queries: `phonon tetrahedron integration energy conservation detailed balance Gaussian collision`; `entropic Fourier method Boltzmann conservation Cai Fan Ying`. For the first, Crossref and Semantic Scholar returned HTTP 429 and OpenAlex returned five mostly irrelevant candidates. For the second, Crossref/OpenAlex returned five each; Semantic Scholar returned 429. Subsequent DOI-specific Crossref requests succeeded for all seven cited papers.

Live searches used tetrahedron/Gaussian detailed balance, entropic Fourier conservation, resonance Van Hove singularities, entropy/Galerkin collision invariants, Galerkin-Petrov, and implicit-surface root quadrature. Full primary HTML/PDF was inspected for every technical attribution above; Saye's implementation was not run. No exhaustive history/novelty search, no peer reports, and no material computation.

Reproducibility trail: `C:\Users\Koussay\ResearchLab\literature\searches\20260927-resonance-measure\`, including `crossref-metadata.json` and `installed-source-lf-sha256.json`. A first local tetrahedron-helper read used the obsolete `structure/` path and failed; the actual `phonon/` helper was then found and inspected. No unresolved access failure affects the cited passages.
