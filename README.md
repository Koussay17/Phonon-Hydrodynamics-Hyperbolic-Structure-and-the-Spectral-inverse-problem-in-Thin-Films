# Phonon Hydrodynamics, Hyperbolic Structure, and the Spectral Inverse Problem in Thin Films

**From kinetic closure to the limits of thermal depth profiling**

> Research in progress — CRTEn internship, September 2026 – January 2027.
> No experimental validation and no general novelty claim. Every result below carries an explicit status.

---

## Overview

**Question.** When can a surface thermal measurement, such as frequency-domain thermoreflectance (FDTR),
determine the properties of a thin film separately, and what does non-Fourier heat transport change?

The project has two connected strands.

1. **Macroscopic models and the inverse problem.** Fourier, Cattaneo and Guyer–Krumhansl (GK) conduction
   in thermal-quadrupole form: their identifiability, entropy admissibility and transport regimes, plus a
   synthetic FDTR design study. Material: AlN films on sapphire.
2. **Microscopic foundation.** From first-principles AlN phonon data towards a conserving kinetic closure:
   collision operators, three-phonon event counting, degeneracies and energy-conserving resonance integration.

Notes are written in French (LaTeX); code, documentation and reports are in English.

---

## Status at a glance

| Part | Topic | Status |
|---|---|---|
| I | Kinetic closure (Boltzmann → GK) | Grey conserving closure derived (note 14). Quantitative AlN closure open. |
| II | Admissibility | Entropy production and propagation speeds established for Fourier/Cattaneo/GK (note 11). |
| III | AlN transport regimes | Regime map and spectral RTA analysis (notes 12, 17); no validated N/U hydrodynamic classification. |
| IV | Inverse problem | Scale invariance, conditional sensitivities, Fourier resonance (note 7); synthetic FDTR design (note 16). |
| Microscopic | AlN collision operator | Audited building blocks (notes 18–23); a complete conserving material operator is still missing. |
| Laboratory | Real measurements | Not yet available: sample, beams, interfaces, calibration and data are required. |

---

## Main results

Status vocabulary: **proved** (audited proof, not machine-checked), **derived under assumptions**,
**numerically demonstrated**, **prior art** (reproduced, attributed), **open**.

**Inverse problem and macroscopic models**
- **Scale invariance survives non-Fourier laws** *(derived, numerically checked)*. In a homogeneous 1D film,
  the response depends only on `b, ξ₁, τ_R, τ_ℓ`. Independent thickness or constitutive constraints can break the gauge. (note 7)
- **Relaxation times degrade identifiability** *(numerically demonstrated, one design)*. Freeing τ_R inflates
  the conductivity/capacity standard errors by 3.13 and 1.61. (note 7)
- **Fourier resonance** τ_R = τ_ℓ *(prior art: Kovács 2018)*. Its position on the regime map depends on the
  closure: x = 1.8 only when τ_c ≈ τ_N; x = 0.8 for the historical closure with τ_c kept; no resonance in the
  conserving closure. (notes 7, 12, 14)
- **ωτ = 1 maximizes single-frequency phase sensitivity**; it is not a detection threshold. (notes 7, Camacho reading note)

**Kinetic theory**
- **Conserving grey closure** *(derived under assumptions)*: nonlocal coefficient α = 1/3, giving
  L² = 4ℓ²/3 and τ_ℓ = 4τ_N/5, against the historical α = 2 (L² = 3ℓ², τ_ℓ = 9τ_N/5). (note 14)
- **GK is second-law admissible but has infinite propagation speed** for ℓ > 0. (note 11)
- **A 500 nm AlN film is boundary-dominated**, not hydrodynamic. Its extracted conductivity is an apparent,
  thickness-dependent quantity. (notes 12, 17)

**Microscopic AlN programme**
- **The operator class decides identifiability** *(proved in class)*. Diagonal rates and static response
  do not fix memory in the broad positive class; in a fixed finite event cone, bounded DC response bounds the memory. (note 18)
- **Exactly resonant Bose-event reference operator**, validated by four independent approaches *(numerically
  demonstrated)*; a hidden-detuning defect was found and fixed. (note 19)
- **First real AlN interaction export** with independent linewidth reconstruction, and a diagnosed polar-cutoff
  defect with its remedy *(numerically demonstrated, coarse mesh)*. (note 20)
- **Invariant block compression is not closed kinetics** *(proved, finite counterexamples)*. (note 21)
- **Hamiltonian, Fock, statistical and linewidth counting factors differ** *(derived under assumptions)*. (note 22)
- **Conserving resonance integration** *(proved locally under finite hypotheses)*. A finite entropy-variable weak
  form conserves quadrature energy and produces entropy. Positive reweighting of off-shell events cannot restore
  conservation, invariants do not constrain rates, and the domain can be left in finite time. (note 23)

Full claim-to-evidence mapping: [paper/CLAIM_MAP.md](paper/CLAIM_MAP.md).

---

## What remains

**Theory/methods paper** — can proceed without new measurements:
1. Choose one central contribution and run a targeted prior-art assessment for that exact statement.
2. Write the integrated manuscript with unified notation; move audits and failed approaches to a supplement.
3. Trace every theorem to its proof and every number to its script; make final figures; run an independent final review.

**Quantitative AlN/FDTR paper** — needs further physics and data:
1. Declare and justify the kinetic approximation, including the treatment of degenerate-mode coherences.
2. Build the full conserving AlN collision operator: globally enumerated resonance surfaces on the real
   dispersion, normal/umklapp/isotope/defect separation, and mesh and broadening convergence.
3. Produce the temperature and thickness map, and model film boundaries and interfaces for the actual geometry.
4. Obtain laboratory data (sample, beams, interfaces, calibration), fit competing models and quantify uncertainty.

Plan and completion criteria: [paper/RESEARCH_PLAN.md](paper/RESEARCH_PLAN.md).

---

## Reference scenarios

For d = 500 nm, C = 2.41×10⁶ J m⁻³ K⁻¹, v = 6000 m/s and an illustrative sapphire substrate:

| Quantity | λ = 60 W m⁻¹ K⁻¹ | λ = 321 W m⁻¹ K⁻¹ |
|---|---|---|
| Characteristic frequency a/(2πd²) | 15.85 MHz | 84.80 MHz |
| Fourier perfect-contact reflection Γ | 0.077 | 0.460 |
| Grey conductivity-derived length 3λ/(Cv) | 12.45 nm | 66.60 nm |
| Grey Knudsen number | 0.0249 | 0.1332 |

The value 60 is illustrative. The value 321 comes from 18 and 22.5 µm films (Cheng et al. 2020), not from a 500 nm
film. The grey length is not an all-collision mean free path. In the synthetic FDTR design (note 16), the local
uncertainty on λ grows from 0.58 % with two free thermal parameters to 7.33 % with ten.

---

## Repository layout

| Folder | Contents |
|---|---|
| [`notes/`](notes) | Authored notes (LaTeX sources and PDFs), project [conventions](notes/conventions.md) |
| [`src/`](src) | Thermal quadrupoles, non-Fourier forward model, FDTR, Stehfest inversion, Fisher analysis, fitting, Bayesian sampling, spectral and collision-event tools |
| [`tests/`](tests) | 241 tests: analytical benchmarks, independent ODE and high-precision checks, regressions |
| [`notebooks/`](notebooks) | Figure scripts (energy trajectory, ωτ threshold, design map, GK blindness, regime map) |
| [`scripts/`](scripts) | Analyses, data preparation, PDF/figure rebuild, campaign reproduction runners |
| [`theory/`](theory) | Synthetic results, AlN data and the archived research campaigns (`theory/aln/*/`) |
| [`paper/`](paper) | Manuscript plan and claim-to-evidence map |
| [`figures/`](figures) | Generated figures |

### Notes

| # | Note | Part |
|---|---|---|
| 00 | Corrections, proofs and the Schrödinger analogy | I–IV support |
| 01–04 | Krapez (2019): reading note, understanding, reminders, contribution and limits | Prior art |
| 05 | Krapez & Rigollet (2017): identifiability comment | Prior art |
| 06 | 2023 identifiability controversy (Mandelis et al., Krapez) | Prior art |
| 07 | Identifiability under non-Fourier laws | IV |
| 08 | The phonon gas | I |
| 09 | From Boltzmann to Guyer–Krumhansl | I |
| 11 | Admissibility: entropy and propagation | II |
| 12 | AlN transport regimes | III |
| 13 | Exercises, Part I | I |
| 14 | Conserving kinetic closure derivation | I |
| 15 | Comparison with prior work | I, IV |
| 16 | Synthetic FDTR experimental design | IV |
| 17 | AlN spectral data, temperature and boundaries | III |
| 18 | Audited spectral closure and operator classes | Microscopic |
| 19 | Validation of the collision-event operator | Microscopic |
| 20 | First microscopic AlN export | Microscopic |
| 21 | Degeneracies and closure | Microscopic |
| 22 | Orientation and counting of cubic channels | Microscopic |
| 23 | Resonant integration and conservation | Microscopic |
| — | Camacho de la Rosa et al. (2025): reading note | Prior art |

---

## Reproducibility

```powershell
python -m pip install -r requirements-dev.txt      # Python >= 3.10
python -B -m pytest tests/ -q -p no:cacheprovider  # 241 tests
python scripts/rebuild.py                           # figures and PDFs (XeLaTeX, pdfLaTeX, DejaVu Sans)
```

| Command | Reproduces |
|---|---|
| `python -X utf8 -B scripts/analyse_experiment.py` | Synthetic FDTR design study (note 16) |
| `python -X utf8 -B scripts/analyse_aln_spectrum.py` | AlN spectral analysis (note 17), offline |
| `python -B scripts/reproduce_resonance_measure.py --output-dir <new dir>` | Resonance-measure checks (note 23), 14 checks |
| `python -B scripts/reproduce_degenerate_action.py --output-dir <new dir>` | Degenerate-action checks (note 21) |
| `python -B scripts/reproduce_channel_counting.py --output-dir <new dir>` | Channel-counting checks (note 22) |

The material-level reproductions (notes 20–22) also need the pinned phono3py 4.5.0 environment and the public
Phonon Olympics force constants; their hashes are recorded in each campaign directory. Laboratory data are not distributed.

Each research campaign in `theory/aln/*/` keeps its question, assumptions, independent branches, experiments,
failed approaches, reviews and final report (`13-final-report.md`). The change history is in [CORRECTIONS.md](CORRECTIONS.md).

---

## Prior work this builds on

The boundary between inherited and new results is kept explicit.

- **J.-C. Krapez, F. Rigollet**, comment on photothermal identifiability, [arXiv:1708.07362](https://arxiv.org/abs/1708.07362) (2017) — structural correlation of coating thickness, diffusivity and conductivity.
- **J.-C. Krapez**, comment on depth-profile reconstructions, J. Appl. Phys. **134**, 056101 (2023) — the graded, functional counterpart.
- **J.-C. Krapez**, Int. J. Therm. Sci. **136**, 182–199 (2019) — Liouville transformation and solvable effusivity profiles.
- **A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, J. Appl. Phys. **137**, 155103 (2025) — closed-form Cattaneo FDTR response, used as external validation.
- **R. Kovács**, [arXiv:1804.05225](https://arxiv.org/abs/1804.05225) (2018) — Fourier resonance of the GK equation.
- **M. G. Hennessy, T. G. Myers**, [GK heat conduction in thermoreflectance experiments](https://doi.org/10.1007/978-3-030-64272-3_2) (2021).
- **J. Ma, W. Li, X. Luo**, Phys. Rev. B **90**, 035203 (2014) — Callaway model for AlN.
- **Z. Cheng et al.**, Phys. Rev. Materials **4**, 044602 (2020) — measured AlN conductivity.
- **M. S. B. Hoque et al.**, Appl. Phys. Lett. **125**, 262201 (2024) — ballistic–diffusive transition in AlN films.
- **A. J. H. McGaughey et al.**, *Phonon Olympics*, J. Appl. Phys. **138**, 135108 (2025) — AlN force constants and inter-code benchmark.
- **G. Fugallo et al.**, Phys. Rev. B **88**, 045430 (2013) — variational ab initio solution of the phonon BTE.
- **G. Lebon, P. C. Dauby**, [Phys. Rev. A 42, 4710](https://doi.org/10.1103/PhysRevA.42.4710) (1990); **L. Sendra et al.**, [Phys. Rev. B 106, 155301](https://doi.org/10.1103/PhysRevB.106.155301) (2022) — kinetic GK coefficients.
- **D. Maillet et al.**, *Thermal Quadrupoles*, Wiley (2000) — transfer-matrix formalism.
- **R. A. Guyer, J. A. Krumhansl**, Phys. Rev. **148**, 766 and 778 (1966); **J. Callaway**, Phys. Rev. **113**, 1046 (1959).

---

## Licence and contact

Code in `src/`, `tests/` and `notebooks/` is released under the MIT licence. Notes and figures are not covered by it.
The converted Rao et al. AlN data retain CC BY 4.0 (see [theory/aln/README.md](theory/aln/README.md)).

**Koussay Mansouri** — internship supervised at CRTEn. Issues and corrections are welcome, including on the derivations.
