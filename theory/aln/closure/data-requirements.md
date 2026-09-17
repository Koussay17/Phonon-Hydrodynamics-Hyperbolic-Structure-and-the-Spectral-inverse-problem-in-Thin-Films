# Collision data contract and computational feasibility

## Required from one consistent calculation
- Crystal cell, atom/isotope composition, harmonic/cubic force constants with hashes and units; temperature and thermal-expansion convention.
- Full-zone q coordinates, reciprocal basis including 2pi convention, q weights, symmetry operations and branch eigenvectors if symmetry reconstruction is used.
- Frequencies, group velocities and Bose capacity weights at each temperature.
- Separate normal, Umklapp, isotope and defect rates, each with a precise event-counting and reciprocal-vector convention. Absorption/emission is an independent classification.
- Collision operator action in an explicitly documented normalization, including detailed balance and energy conservation residuals. Prefer a matrix-free action or event factorization to a dense matrix.
- Mesh sequences for each target temperature, with convergence of DC transport, response at nonzero frequencies, current-weighted inverse moments, and slow/fast projection quantities. A stable total conductivity alone is insufficient.
- Interface geometry, orientation, specular/diffuse reflection and transmission model, and physical sample parameters; no double-counted boundary rate.

## Available locally
Compact Rao data are irreducible and single-temperature. No N/U or full collision action is present. Python numerical tools exist; phono3py, phonopy, spglib and h5py were not importable in the active interpreter on campaign initialization. This is an environment observation, not a claim that no external installation exists.

## Storage arithmetic
See experiments/operator_storage.json. At 24^3 q points and 12 branches, a dense float64 collision matrix is about 205 GiB. At 45^3 it is about 8909 GiB. These are arithmetic storage estimates, not measured run costs. Do not attempt a dense full-zone operator on this PC without a resource assessment. Matrix-free projected solves and event representations may avoid this obstruction.

## Unresolved external calculation
Public force constants can enable a new calculation, but their availability is not equivalent to having validated rates or an operator. No first-principles rerun has been performed by this campaign. An implementation choice must be based on verified software capability, convergence and resource cost.

## Machine inventory on 2026-09-17
Windows reports 29,909,643,264 bytes physical RAM (about 27.86 GiB), 12 logical processors, and about 32.03 GiB free on the campaign volume. The 205.03 GiB dense 24^3 scalar collision matrix alone exceeds both RAM and free disk. This rules out that naive representation on the present machine, not all possible matrix-free or distributed calculations. Export basis and parity coverage must be verified before selecting a smaller representation. No resource-heavy calculation was launched.
