# Research programme — theory, simulation and digital twin (October 2026 – January 2027)

Scope: everything that can be established theoretically and computationally. Laboratory measurements are out of
scope. Wherever a result would need them, the programme produces a **digital twin**: a high-fidelity synthetic
experiment with declared inputs and uncertainty. These synthetic results are never presented as measurements.

Two publications are targeted:
- **Paper 1 (theory/methods):** what static transport data determine about a conserving dynamic closure, and
  what additional operator and measurement information is needed.
- **Paper 2 (quantitative AlN/FDTR, computational):** first-principles AlN phonon transport → conserving kinetic
  closure → thin-film response → FDTR digital twin → identifiability of non-Fourier parameters, with the
  laboratory step left as an explicit, prepared interface.

## Workstreams

| ID | Workstream | Main deliverables | Depends on |
|---|---|---|---|
| **A** | Central claim and prior art (paper 1) | Precise theorem statements; targeted novelty audit; manuscript outline | Notes 14, 18–23 |
| **B** | First-principles AlN transport | Validated phono3py pipeline (RTA + full LBTE) against Phonon Olympics; mesh, smearing and tetrahedron convergence; κ(T) for 100–1000 K | Olympics fc2/fc3 |
| **C** | Collision operator and spectral structure | Full collision matrices; N/U/isotope split; relaxon spectrum; slow modes; energy and momentum invariants; Guyer window vs T | B |
| **D** | Kinetic → hydrodynamic coefficients | Spectral (not grey) τ_R, τ_N, ℓ, viscosity; validity of GK vs T; comparison of conserving and historical closures on real data | C, note 14 |
| **E** | Conserving resonance measure on the real dispersion | Globally enumerated resonance-surface quadrature (note 23 on AlN); Gaussian vs tetrahedron vs surface rule; coincident-label weight | B, note 23 |
| **F** | Thin-film kinetic transport | Cross-plane spectral BTE (deterministic and Monte Carlo) with boundary scattering; effective conductivity and memory vs thickness; GK boundary conditions tested against BTE | C, D |
| **G** | FDTR digital twin | Au/AlN/sapphire axisymmetric model: Fourier, Cattaneo, GK and BTE-informed responses; beams, interfaces, calibration, noise; synthetic data generation | F, src/fdtr.py |
| **H** | Inverse problem on the twin | Identifiability and Fisher/Bayesian analysis; model selection (Fourier vs GK vs kinetic); bias from model mismatch; optimal experiment design; detectability thresholds | G |
| **I** | Formal and symbolic verification | Lean proofs of central finite lemmas (paper 1); symbolic checks of closure algebra | A |
| **J** | Manuscripts | Paper 1 draft → review → final; paper 2 draft → review → final; unified notation; figures; supplements | all |

## Order of execution
1. **Wave 1 (now):** A (prior art and central claim), B (pipeline validation), G0 (twin architecture: Fourier/GK
   multilayer with full nuisance model, built on `src/fdtr.py`).
2. **Wave 2:** C and E on the validated pipeline; I on the claims fixed by A.
3. **Wave 3:** D and F; G with BTE-informed responses.
4. **Wave 4:** H on the full twin; J drafting both papers.

## Standards (all workstreams)
- Every campaign follows the repository protocol: question, assumptions, independent branches, experiments,
  failed approaches, four independent reviews, final report (`theory/aln/<campaign>/`).
- Validation against external references comes before any new result: Phonon Olympics κ, analytic limits,
  manufactured solutions.
- Convergence (mesh, smearing, time step, particle count) is reported with numbers for every material result.
- Status vocabulary as in the README. No novelty claim before the targeted prior-art audit.
- Bulky data on D:; the repository keeps scripts, hashes and summarised results.

## Compute budget
Laptop: 6 cores / 12 threads, 28 GB RAM, RTX 3050 Ti (4 GB). Feasible: phono3py LBTE up to about 24³ q-points on
wurtzite AlN (irreducible collision matrix of order 10⁴); 1D spectral BTE; axisymmetric FDTR models. Larger meshes
use the matrix-free and irreducible formulations, or are reported as convergence trends.
