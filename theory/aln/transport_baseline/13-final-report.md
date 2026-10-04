# 13 — Final report: AlN lattice thermal-conductivity baseline (workstream B)

Campaign `D:\ResearchLab\orchestration\campaigns\20261002-aln-transport-baseline\`; details in
`07-experiments.md`, negative results in `09-failed-approaches.md`, tables in `results/`.
Inputs: Phonon Olympics AlN phono3py files (commit 0640f077; POSCAR, BORN, fc2 5x5x3, fc3 3x3x2),
SHA-256-checked at every run. Software: phono3py/phonopy 4.5.0, Python 3.14.7.

## 1. Validation against the benchmark (300 K, W/(m K), in-plane / cross-plane)
| | RTA | LBTE |
|---|---|---|
| Olympics phono3py team (2.1.0, 31x31x17) | 252.992 / 231.947 | 285.048 / 271.262 |
| this work, 2.1.0-compatible settings, 31x31x17 | 252.992 / 231.947 | 285.048 / 271.262 |
| this work, phono3py 4.5.0 defaults, 31x31x17 | 270.47 / 247.09 | 302.08 / 285.56 |
| Olympics ShengBTE (own force constants) | 271 / 251 | 298 / 291 |
| Olympics ALAMODE | 282 / 263 | – |

* **NUMERICALLY DEMONSTRATED:** with C kernels, `make_r0_average=False` and finite-difference
  group velocities, the team's 300 K values are reproduced at every mesh 9x9x5–31x31x17 to
  <= 0.02 % (31x31x17: < 5e-4 % for RTA and LBTE), and the RTA kappa(T) to 0.001 % for
  100–1000 K (+1.25 % at 20 K).
* phono3py 4.5.0 defaults give +6.0/+5.3 % (LBTE) and +6.9/+6.5 % (RTA) at 31x31x17. The cause
  is the fc3 real-to-reciprocal phase convention (`make_r0_average`, default changed after 2.x);
  velocity method and kernel language contribute <= 4e-4. This sensitivity reflects the
  truncation of fc3 in the 3x3x2 supercell and is the largest model-level uncertainty found.

## 2. Convergence (production setting = phono3py 4.5.0 defaults, tetrahedron, 300 K)
LBTE 273.96/268.04 (11x11x7), 286.38/275.51, 300.66/277.68, 298.25/281.09, 299.05/283.58,
302.08/285.56 (31x31x17). Last refinement +1.0 % (xx), +0.7 % (zz). In-plane is non-monotone;
cross-plane still rises. Extrapolation over the four finest meshes (kappa_inf + b N^-p):
xx 300.4–302.1, zz 287.5–298.6 for p = 1 ... 1/3. **Estimate: LBTE 302 +/- 3 (xx) and
286 with an upward mesh correction of +0.7 to +4.6 % (zz); RTA 270 +/- 3 (xx) and 247
(+0.7 to +4.5 %) (zz).** The
LBTE/RTA ratio is converged (1.117 xx, 1.155 zz, +/-0.1 % from 19x19x11). Gaussian smearing:
sigma = 0.05 THz within 0.6 % of the tetrahedron; sigma = 0.1 / 0.2 THz lower by up to
0.9 % / 3.2 % (RTA); Gaussian LBTE sigma = 0.1 at 19x19x11 within -0.5/+0.3 %.
Status: NUMERICALLY SUPPORTED to ~1 % in-plane; cross-plane convergence below 1 % is not
demonstrated (monotone drift).

## 3. kappa(T), production setting, 31x31x17 (xx / zz)
| T (K) | RTA | LBTE |
|---|---|---|
| 100 | 2565.9 / 2315.2 | 3422.1 / 3153.3 |
| 200 | 549.2 / 499.2 | 624.0 / 586.2 |
| 300 | 270.5 / 247.1 | 302.1 / 285.6 |
| 500 | 135.3 / 124.3 | 150.2 / 142.6 |
| 1000 | 62.5 / 57.6 | 69.2 / 65.9 |
RTA at 16 temperatures in `results/tables_production.md`. Relative to the Olympics phono3py
LBTE: +9.8/+10.9 % (100 K) decreasing to +4.7/+4.1 % (1000 K). Relative to ShengBTE iterative
(independent force constants): -1.7/-5.3 % (100 K), +1.2/-1.9 % (300 K), +1.3/-1.5 % (1000 K).

## 4. Isotope scattering (natural abundance)
LBTE 300 K: -0.14 % (xx and zz). RTA: -1.2/-1.0 % (100 K), -0.11/-0.09 % (300 K),
-0.05/-0.04 % (1000 K). Negligible above 200 K compared with the mesh uncertainty.

## 5. Saved data (production, 31x31x17, D:, SHA-256 in `results/saved_data_inventory.json`)
| item | path (under runs/) | size | SHA-256 (first 16) |
|---|---|---|---|
| assembled matrix 100 K | ooc-m313117-prod/A_T0.npy | 7.74 GB | 7f2499986582da17 |
| assembled matrix 200 K | ooc-m313117-prod/A_T1.npy | 7.74 GB | 8dbe20048988ad03 |
| assembled matrix 300 K | ooc-m313117-prod/A_T2.npy | 7.74 GB | 1c43a6982268ef80 |
| assembled matrix 500 K | ooc-m313117-prod/A_T3.npy | 7.74 GB | 7b10f25817dd37f1 |
| assembled matrix 1000 K | ooc-m313117-prod/A_T4.npy | 7.74 GB | 688140811ac07422 |
| 300 K matrix before averaging/symmetrization | ooc-m313117-prod/A_T2_prephono3pyavg_presym.npy | 7.74 GB | 9c0d80939b7d6eb8 |
| meta (ir points, weights, rotations, frequencies, velocities) | ooc-m313117-prod/meta.npz | 0.5 MB | 7cf63be64ec25f3e |
| gamma (5 T) | ooc-m313117-prod/gamma.npy | 0.4 MB | baf89c9c81a8eb26 |
| LBTE solutions Y (5 T) | ooc-m313117-prod/Y_T*.npy | 0.25 MB each | see inventory |
| RTA all T incl. gamma_N, gamma_U, cv, gv, weights | prod-rta-m313117-NU/kappa-m313117.hdf5 | 10.4 MB | 6e1694f9afb023f8 |
| RTA all T with isotopes (gamma_isotope) | prod-rta-m313117-iso/kappa-m313117.hdf5 | 8.0 MB | eb051e367f1b44bf |
Reproduction-setting matrix (300 K): `ooc-m313117-v2C/A_T0.npy`.
Conventions (irreducible kappa-star representation, n = 864 x 12 x 3, THz, weights
sqrt(star/24), degeneracy averaging, (A+A^T)/2, kappa formula) are in the inventory file and
the docstring of `scripts/lbte_ooc.py`. gamma_N + gamma_U = gamma to 1e-16; at 300 K normal
processes carry 69 % (xx) / 61 % (zz) of the kappa-weighted linewidth.

**Collision-matrix checks.** Symmetric after assembly; pre-symmetrization asymmetry 0.13–0.57 %
of max|A| (100–1000 K); positive on the symmetry-allowed subspace (PCG curvature >= 3.8e-4 THz,
preconditioned Ritz values 0.59–4.0). Energy is not a null vector of phono3py's matrix: on the
scalar grid ||Omega e||/||D e|| = 2.3–2.4 at all meshes, as required by Omega - Omega' =
(C0 + C2)(I + J). The physical operator Omega' rebuilt from the three triplet channels
annihilates e to 3.7e-3 (tetrahedron, 9x9x5) and 5.5e-4 (Gaussian 0.1 THz, 7x7x5); the two
operators coincide on odd vectors (6e-16). **The exported matrix is valid for the odd
(current-carrying) sector only.**

## 6. Independent checks
RTA re-summation from saved arrays: <= 8.5e-15 relative (all T, with and without isotopes);
LBTE by own assembly + projected PCG vs phono3py dense pseudo-inverse: 1.2e-12 (15x15x9),
assembled matrix bit-identical (9x9x5). The 31x31x17 reproduction LBTE was obtained only by
this second route.

## 7. Failures (see 09-failed-approaches.md)
Default-setting mismatch (r0 convention); C path lacks analytic NAC velocities; C/Rust phonons
differ by <= 6e-7 THz, amplified by the tetrahedron method (single linewidths up to 1.6 % at
9x9x5); unprojected PCG diverged; dense LBTE infeasible at 31x31x17 on this machine; D: drive
loss and memory-pressure termination of the first queues (work resumed, nothing lost).

## 8. Recommended production setting for downstream work
phono3py 4.5.0, defaults (Rust kernels, `make_r0_average=True`, analytic velocities), linear
tetrahedron, 31x31x17, LBTE through `scripts/lbte_ooc.py`; operator files above. Quote
kappa(300 K) = 302 / 286 W/(m K) with mesh uncertainty +/-1 % (xx) and +1 to +5 % (zz), and
the fc3 phase-convention spread (6 %) as a model uncertainty. For energy/viscosity or film
problems needing the even sector, rebuild Omega' from channel-resolved weights (Gaussian
smearing preferred: symmetric to 1e-8 and time-reversal invariant), not phono3py's matrix.
Next discriminating calculation: fc3 from a larger supercell (e.g. 4x4x3) to resolve the
r0-convention spread.
