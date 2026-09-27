# Closure terminology and primary evidence

Bounded second literature pass, 2026-09-27. First-pass L report preserved. Read first-pass B/D for the precise question; no material calculation. Scope: factorial moments, quantum molecular chaos, quasifree closure, and population-versus-spectral relaxation. No novelty claim.

## Finding

Primary literature distinguishes an occupation-probability master equation from a closed equation for mean populations. Quasifree/Wick or molecular-chaos factorization is an additional approximation. Inspected sources do **not** directly establish the campaign's repeated-daughter factor-of-two comparison, nor a universal correction to phono3py. In particular, a continuum phonon equation does not by itself prescribe finite coincident-mode kinetics.

## 1. Quantum Boltzmann master equation versus quantum Boltzmann equation

**SOURCE SAYS.** C. W. Gardiner and P. Zoller, *Quantum kinetic theory: A quantum kinetic master equation for condensation of a weakly interacting Bose gas without a trapping potential*, Phys. Rev. A **55**, 2902-2921 (1997), [DOI](https://doi.org/10.1103/PhysRevA.55.2902), inspected [quant-ph/9611043v1](https://arxiv.org/html/quant-ph/9611043v1). Section V.1, equations (100)-(103), evolves probabilities of occupation configurations. Equation (110) gives mean evolution containing joint occupation moments; equation (111) follows only after factorizing them. Section V.2 calls the resulting kinetic equation the Uehling-Uhlenbeck equation. Section V explicitly relates the approximation to molecular chaos.

**COVERAGE LIMIT.** This is bosonic two-body scattering, not a cubic phonon channel. It discusses discrete energy levels (V.1.1), but the displayed collision factors use separately indexed occupations. I did not locate an explicit coincident-index replacement by falling factorials. Those formulas are therefore not verified counting rules for our repeated discrete daughter.

**INDEPENDENT FOLLOW-UP SOURCE.** D. Jaksch, C. W. Gardiner and P. Zoller, *Quantum kinetic theory. II. Simulation of the quantum Boltzmann master equation*, Phys. Rev. A **56**, 575-586 (1997), [DOI](https://doi.org/10.1103/PhysRevA.56.575); inspected [quant-ph/9701008v2](https://arxiv.org/html/quant-ph/9701008v2), whose preprint title says **III**, not the published **II**. Section II.1.2 requires an effectively continuous spectrum for its Markov approximation despite simulating finite systems. Section II.4, equation (18), explicitly factorizes occupation averages and notes that the mean equation omits occupation fluctuations. Section II.3 separately introduces an ergodic approximation across degenerate levels. Thus equal energy, statistical independence, and discarding coherences are different assumptions. This source also does not explicitly resolve our coincident-index factors.

## 2. Locally quasifree phonon kinetic limit

**SOURCE ASSUMES/DERIVES.** H. Spohn, *The Phonon Boltzmann Equation, Properties and Link to Weakly Anharmonic Lattice Dynamics*, J. Stat. Phys. **124**, 1041-1104 (2006), [DOI](https://doi.org/10.1007/s10955-005-8088-5), inspected [math-ph/0505025v2](https://arxiv.org/html/math-ph/0505025v2). Equation (9.6) defines gauge-invariant quasifree moments through a permanent of two-point covariances. Around (9.10)-(9.12), a finite-box diagonal quasifree density matrix is exponential in number operators. The paragraph after (9.16) assumes local quasifreeness remains a good approximation at kinetic times. Equation (10.1) performs Wick factorization; (10.5) is the continuum-momentum phonon collision integral. Equation (13.7) gives its linearized event quadratic form.

**COVERAGE LIMIT.** The kinetic theory uses an infinite spatial lattice and continuous Brillouin-zone momentum, with scaling limits; the finite box in section 9 is used to compute entropy. No exact invariant product-geometric manifold for a finite interacting channel is established. A [2006 erratum](https://doi.org/10.1007/s10955-006-9144-5), J. Stat. Phys. **123**, 707, was found, but its correction text was inaccessible behind the publisher's access page. I therefore use the inspected ansatz/setup, not its stationary-state classification or collision prefactors.

**OUR INTERPRETATION.** A finite diagonal gauge-invariant quasifree state has geometric mode occupations; its same-mode Wick contractions give $\langle (a^\dagger)^2a^2\rangle=2\langle a^\dagger a\rangle^2$. This algebraic statement does not establish propagation of that state family. Nor may the coincident set be dropped merely because it is Lebesgue-null before enforcing resonance: one must check its mass under the actual delta-constrained collision measure. Critical resonances, singular occupations, and flat dispersions need separate treatment. No such measure analysis for the AlN channel was performed here.

## 3. Explicit repeated discrete-mode source

**SOURCE SAYS.** T. E. Lee and H. R. Sadeghpour, *Quantum Synchronization of Quantum van der Pol Oscillators with Trapped Ions*, Phys. Rev. Lett. **111**, 234101 (2013), [DOI](https://doi.org/10.1103/PhysRevLett.111.234101); inspected [1306.6359v2](https://arxiv.org/html/1306.6359v2). Equation (2) contains the single-mode jump operator $a^2$. The following paragraph gives its two-phonon event rate as $2\kappa_2\langle a^{\dagger 2}a^2\rangle$. Footnote 1 connects this dissipator to a nonlinear system-bath coupling $a^{\dagger 2}b+a^2b^\dagger$. This genuinely addresses two quanta of the **same discrete oscillator**, not two continuum labels.

**DERIVED ILLUSTRATION, NOT A CLAIM QUOTED FROM THE PAPER.** Isolate its loss term and write $\dot\rho=\gamma\mathcal D[a^2]\rho$, with $\mathcal D[L]\rho=L\rho L^\dagger-\{L^\dagger L,\rho\}/2$ and $\gamma=2\kappa_2$. Let $F_r=\langle N(N-1)\cdots(N-r+1)\rangle$. The jump $N\to N-2$ gives

\[
\dot{\langle N\rangle}=-2\gamma F_2,\qquad
\dot F_2=-4\gamma F_3-2\gamma F_2.
\]

This is a factorial-moment hierarchy. At a geometric state, $F_2=2n^2$, $F_3=6n^3$; the actual $\dot F_2=-24\gamma n^3-4\gamma n^2$ differs from the geometric tangent $4n\dot n=-16\gamma n^3$ for $n>0$. Hence even this simpler nonlinear loss model does not preserve geometric statistics. It is an open-system example, not a derivation of rates for an isolated resonant $p\leftrightarrow a+a$ Hamiltonian.

## 4. Spectral lifetime versus population evolution

**SOURCE SAYS.** A. F. Kemper, O. Abdurazakov and J. K. Freericks, *General Principles for the Nonequilibrium Relaxation of Populations in Quantum Materials*, Phys. Rev. X **8**, 041009 (2018), [DOI](https://doi.org/10.1103/PhysRevX.8.041009); inspected [1708.05725v2](https://arxiv.org/html/1708.05725v2). Sections IV.2-IV.3 distinguish relative-time propagator damping from average-time population evolution. Equations (15)-(18) contain both Green-function and self-energy distributions. Equation (21) retains both $\delta f^G$ and $\delta f^\Sigma$; dropping the latter is an additional limit discussed after (23).

**COVERAGE LIMIT.** These are electronic pump-probe models. They establish the general danger of equating a linewidth with population relaxation; they neither analyze repeated discrete phonons nor prove our factor two. In our application, differentiating every occurrence of a coincident mode in a kinetic closure is a separate operation from evaluating an equilibrium self-energy.

## What this licenses, and what remains open

Useful terms are **factorial-moment hierarchy**, **occupation-probability/Pauli master equation**, **quantum molecular chaos**, **quasifree (Wick) closure**, and **spectral versus population relaxation**. A product-geometric closure can be a declared reduced model; the sources do not make it exact for the finite channel. The first-pass factor comparison must therefore remain conditional, not a code-defect claim. Next: retain the finite-channel probability distribution or explicit additional moments and test whether the proposed reduced state family is tangent to its generator. A continuum reduction additionally needs a scaling and resonance-measure argument.

## Search/access record

Local `research-search.cmd` queries: `Spohn phonon Boltzmann equation locally quasifree Gaussian decoupling`; `Quantum Kinetic Theory III Simulation Quantum Boltzmann Master Equation Jaksch Gardiner Zoller`. Crossref returned five results for each; OpenAlex returned zero and five, respectively; Semantic Scholar returned HTTP 429. Live searches added `two-photon loss factorial moments master equation`, `quantum van der Pol Sadeghpour 2013`, and `Kemper population relaxation quasiparticle lifetime self energy collision integral general principles nonequilibrium populations`. Backward tracing inspected Gardiner-Zoller from Jaksch et al.; the Spohn erratum was checked but not obtained. All cited article metadata were independently retrieved from Crossref. Primary HTML, not abstracts alone, supports the equation pointers above.

Logs and metadata: `C:\Users\Koussay\ResearchLab\literature\searches\20260927-channel-closure\` (including `crossref-metadata.json`). This is a focused terminology/source pass, not an exhaustive historical search. No inspected paper directly addresses the complete finite repeated-daughter/linewidth comparison.
