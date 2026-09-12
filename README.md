# Phonon Hydrodynamics, Hyperbolic Structure, and the Spectral Inverse Problem in Thin Films

**From Kinetic Closure to the Limits of Thermal Depth Profiling**

> **Status — work in progress.** Research internship, September 2026 to January 2027, CRTEn.


---

## Scope


The physical setting is aluminium nitride thin films, where the phonon mean free path and
the phonon relaxation time are not negligible compared with the film thickness and the
modulation period. In that regime Fourier's law is no longer the obvious starting point,
and the inverse problem changes character.

The work has four strands:

- **Part I — kinetics.** Deriving the Guyer–Krumhansl equation from the Boltzmann transport
  equation through a moment hierarchy and a closure, using the Callaway separation between
  normal and umklapp scattering.
- **Part II — admissibility.** Which non-Fourier models are compatible with a convex entropy
  and with finite characteristic speeds, and what that imposes on the relaxation time and the
  non-local length.
- **Part III — regimes.** A transport regime map for AlN as a function of temperature and film
  thickness, anchored to published measurements rather than to estimates.
- **Part IV — the inverse problem.** After a Liouville transformation the spatial operator takes
  the form of a stationary Schrödinger operator. The question is what the position of the
  associated spectral parameter in the complex plane implies for the conditioning of parameter
  recovery.

A single quantity runs through all four: the **relaxation time τ**. It arises from the kinetics,
it is constrained by thermodynamics, and it governs the conditioning of the inverse problem.

---

## Repository layout

| Folder | Contents |
|---|---|
| `notes/` | Reading notes, one file per paper. Novelty report. Out-of-scope ideas, dated. |
| `theory/` | Derivations in LaTeX. One file per step of the chain, so that a broken step can be found without reading everything. |
| `src/` | Reusable code: quadrupole assembly, numerical Laplace inversion, forward model, inversion routines. |
| `notebooks/` | Exploratory and figure-producing notebooks, numbered in the order they were written. |
| `tests/` | Unit tests. Every forward-model routine must pass an analytical benchmark before it touches real data. |
| `data/` | **Never committed.** Laboratory measurements live here locally and are excluded by `.gitignore`. |
| `figures/` | Generated figures. Regenerable from `notebooks/`; committed for convenience only. |
| `paper/` | Manuscript sources. |

---

## Naming conventions

**Notebooks** — `NN_short_description.ipynb`, numbered in creation order:

```
notebooks/01_cattaneo_energy_trajectory.ipynb
notebooks/02_quadrupole_homogeneous_benchmark.ipynb
notebooks/03_aln_regime_map.ipynb
```

The number records the order in which the work was done, which matters more than alphabetical
order when the notebooks are read six months later.

**Source files** — `snake_case.py`, lower case, underscores only. No hyphens: a file named
`thermal-quadrupoles.py` cannot be imported.

**Theory files** — `NN_topic.tex`, matching the step numbering in `notes/`.

**Reading notes** — `firstauthor_year.md`, e.g. `krapez_2018.md`.

---

## Reproducing the figures

Every figure in `figures/` is produced by exactly one notebook, named in the figure caption.
Notebooks run top to bottom with no hidden state; restart the kernel and run all before
committing.

Requirements: `numpy`, `scipy`, `matplotlib`, and `emcee` for the posterior sampling in the
inversion notebooks. Google Colab is sufficient for everything except the longer sampling runs.

---

## Data policy

**Laboratory data is not published in this repository, and never will be.**

The `data/` folder is excluded by `.gitignore`, which was added before the first commit rather
than after. Notebooks that require real measurements will fail on a fresh clone; that is
intended. Synthetic data generators are provided in `src/` so that every result in the
theoretical parts can be reproduced without laboratory access.

Papers in PDF form are not committed either. Bibliographic records live in a Zotero library;
this repository holds citations, not copies.

---

## Prior work this builds on

The following establish results that this work uses rather than reproves. They are listed here
so that the boundary between what is inherited and what is new is visible from the outset.

- **J.-C. Krapez**, *Linear, trigonometric and hyperbolic profiles of thermal effusivity in the
  Liouville space and related quadrupoles*, International Journal of Thermal Sciences **136**
  (2018) 182–199. Establishes the Liouville transformation of the heat equation, identifies the
  Schrödinger potential as the relative curvature of the square root of the thermal effusivity,
  and gives the corresponding quadrupole formulations. Explicitly declines the spectral reading
  of the resulting operator.
- **J.-C. Krapez**, *Heat diffusion in inhomogeneous graded media: chains of exact solutions by
  joint Property and Field Darboux Transformations*, International Journal of Heat and Mass
  Transfer **99** (2016) 485–503. The foundational paper of the above.
- **P. Chen, I. M. Gamba, Q. Li, L. Wang**, *Reconstruction of heat relaxation index in phonon
  transport equation*, arXiv:2502.19533 (2025), accepted in SIAM Journal on Applied Mathematics.
  Numerical reconstruction of the relaxation time from surface temperature by PDE-constrained
  optimisation.
- **A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, *Relaxation times of non-Fourier
  materials using frequency-domain thermoreflectance*, Journal of Applied Physics **137** (2025)
  155103.
- **M. S. B. Hoque et al.**, *Experimental observation of ballistic to diffusive transition in
  phonon thermal transport of AlN thin films*, Applied Physics Letters **125** (2024) 262201.
- **D. Maillet, S. André, J.-C. Batsale, A. Degiovanni, C. Moyne**, *Thermal Quadrupoles: Solving
  the Heat Equation through Integral Transforms*, Wiley (2000). Reference text for the transfer
  matrix formalism used throughout `src/`.

---



## Licence

Code in `src/`, `tests/` and `notebooks/` is released under the MIT licence. Manuscript text and
figures are not covered by it. Laboratory data is not distributed..

---

## Contact

Koussay Mansouri — internship supervised at CRTEn.
Issues and corrections are welcome, including on the derivations.
