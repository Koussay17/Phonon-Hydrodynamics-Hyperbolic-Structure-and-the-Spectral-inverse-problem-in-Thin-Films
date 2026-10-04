# 07 — Experiments (workstream B, AlN transport baseline)

All runs: phono3py/phonopy 4.5.0 (`D:\ResearchLab\envs\aln-phono3py-4.5.0`, Python 3.14.7),
pinned inputs verified by SHA-256 at every start, 300 K unless stated, linear tetrahedron
unless stated, NAC (Gonze-Lee, default cutoff), no isotope unless stated, 6 OpenMP/Rayon
threads on a Ryzen 5 5600H (6 cores), 28 GB RAM (2–4 GB free during the campaign).
Each run directory holds `run.json` (settings, versions, input and output hashes, timings,
peak memory) and phono3py outputs; logs in `runs/logs/`. Tables: `results/tables.md`.

Two settings are used throughout:
* **v2C (reproduction)** — emulates the phono3py-2.1.0 run of the Olympics team: C kernels,
  `make_r0_average=False`, finite-difference group velocities (`gv_delta_q=1e-5`).
* **prod (production)** — phono3py 4.5.0 defaults: Rust kernels, `make_r0_average=True`,
  analytic group velocities including the Gonze-Lee NAC derivative.

## E1. Version switches at 9x9x5 (RTA, W/(m K), xx / zz)
| run | kernels | r0 average | gv | kappa |
|---|---|---|---|---|
| t1-rta-m995-default | Rust | on | analytic | 231.182 / 212.418 |
| t1-rta-m995-gvfd | Rust | on | FD 1e-5 | 231.182 / 212.418 (diff 5e-9) |
| t1-rta-m995-lang-C-gvfd | C | on | FD 1e-5 | 231.152 / 212.343 |
| t1-rta-m995-r0off | Rust | off | analytic | 224.950 / 206.254 |
| t1-rta-m995-v2compat-C | C | off | FD 1e-5 | 224.918 / 206.184 |
| Olympics phono3py team | (2.1.0) | (off) | (FD) | 224.932 / 206.227 |
The replica interpreter reproduced the D: interpreter bit-for-bit (231.18181689782404).

## E2. Reproduction of the Olympics phono3py team (v2C settings), 300 K
Table from `results/tables.md` (Olympics values from their spreadsheets; ours from phono3py's
dense LBTE for 9x9x5 and 15x15x9 and from the rows+PCG route otherwise).

| mesh | Olympics LBTE xx/zz | ours LBTE | diff % | Olympics RTA | ours RTA | diff % |
|---|---|---|---|---|---|---|
| 9x9x5 | 255.558/243.287 | 255.528/243.287 | -0.012/-0.000 | 224.932/206.227 | 224.918/206.184 | -0.006/-0.021 |
| 15x15x9 | 275.638/264.698 | 275.621/264.745 | -0.006/+0.018 | 244.249/227.031 | 244.238/227.053 | -0.005/+0.010 |
| 19x19x11 | 285.019/266.578 | 285.016/266.586 | -0.001/+0.003 | 253.240/228.162 | 253.237/228.163 | -0.001/+0.000 |
| 23x23x13 | 282.998/268.195 | 282.998/268.194 | +0.000/-0.000 | 251.052/229.670 | 251.052/229.669 | +0.000/-0.000 |
| 27x27x15 | 283.146/269.548 | 283.147/269.551 | +0.000/+0.001 | 251.300/230.664 | 251.301/230.666 | +0.000/+0.001 |
| 31x31x17 | 285.048/271.262 | 285.048/271.262 | -0.000/+0.000 | 252.992/231.947 | 252.992/231.947 | +0.000/+0.000 |

31x31x17 LBTE: 285.0478/271.2622 W/(m K), 19 PCG iterations, true residual 8.0e-11.
RTA kappa(T) at 31x31x17: within 0.001 % of the team's table for 100–1000 K; +0.02 % at 50 K,
+0.27 % at 30 K, +1.25 % at 20 K (low-temperature tetrahedron sensitivity).

## E3. Production convergence, 300 K
Mesh series (tetrahedron), W/(m K):

| mesh | N | LBTE xx / zz | RTA xx / zz | LBTE/RTA xx / zz |
|---|---|---|---|---|
| 11x11x7 | 847 | 273.96 / 268.04 | 243.81 / 231.06 | 1.124 / 1.160 |
| 15x15x9 | 2025 | 286.38 / 275.51 | 255.46 / 239.06 | 1.121 / 1.153 |
| 19x19x11 | 3971 | 300.66 / 277.68 | 269.30 / 240.52 | 1.117 / 1.155 |
| 23x23x13 | 6877 | 298.25 / 281.09 | 266.77 / 243.65 | 1.118 / 1.154 |
| 27x27x15 | 10935 | 299.05 / 283.58 | 267.65 / 245.65 | 1.117 / 1.154 |
| 31x31x17 | 16337 | 302.08 / 285.56 | 270.47 / 247.09 | 1.117 / 1.156 |

* Last refinement (27x27x15 -> 31x31x17): +1.0 % (xx), +0.7 % (zz); spread over the four finest
  meshes 1.3 % (xx), 2.8 % (zz). xx is non-monotone (19x19x11 > 23x23x13), so power-law fits
  kappa = k_inf + b N^-p over the four finest meshes are poor for xx (rms 1.4 W/(m K)):
  k_inf(xx) = 300.4–302.1 for p = 1 ... 1/3. zz rises monotonically: k_inf(zz) = 287.5 (p = 1),
  290.2 (p = 2/3), 298.6 (p = 1/3). RTA behaves the same (xx 268.8–270.2; zz 248.9–258.1).
* The LBTE/RTA ratio is converged to ~0.1 % from 19x19x11 on: 1.117 (xx), 1.155 (zz).
* Gaussian vs tetrahedron (RTA): sigma = 0.05 THz within +0.4/+0.4 % (15x15x9) and 0.0/+0.6 %
  (19x19x11); sigma = 0.1: -0.5/-0.5 % and -0.9/-0.2 %; sigma = 0.2: -1.5/-1.8 % and -3.2/-1.3 %.
  Gaussian LBTE (sigma = 0.1, 19x19x11): 299.01/278.39 vs tetrahedron 300.66/277.68 (-0.5/+0.3 %).
* Production vs reproduction settings at 31x31x17: +6.0 %/+5.3 % (LBTE), +6.9 %/+6.5 % (RTA);
  9x9x5: +2.8 %/+3.0 % (RTA). The fc3 phase convention (make_r0_average) is the dominant
  model-level sensitivity of this force-constant set.

## E4. kappa(T), production setting, 31x31x17 (W/(m K), xx / zz)
| T | RTA | LBTE | Olympics phono3py LBTE | ShengBTE iterative |
|---|---|---|---|---|
| 100 | 2565.9 / 2315.2 | 3422.1 / 3153.3 | 3118.0 / 2843.2 | 3481.1 / 3329.5 |
| 200 | 549.2 / 499.2 | 624.0 / 586.2 | 581.7 / 549.2 | 621.5 / 603.1 |
| 300 | 270.5 / 247.1 | 302.1 / 285.6 | 285.0 / 271.3 | 298.4 / 291.0 |
| 500 | 135.3 / 124.3 | 150.2 / 142.6 | 142.9 / 136.6 | 148.3 / 145.0 |
| 1000 | 62.5 / 57.6 | 69.2 / 65.9 | 66.1 / 63.3 | 68.3 / 66.9 |
RTA at all 16 temperatures 100–1000 K and the Olympics/ShengBTE RTA columns: `results/tables_production.md`.
Production LBTE vs ShengBTE iterative: -1.7/-5.3 % (100 K), +0.4/-2.8 % (200 K), +1.2/-1.9 % (300 K),
+1.3/-1.7 % (500 K), +1.3/-1.5 % (1000 K). The ShengBTE team used its own force constants
(QE, own FORCE_CONSTANTS_3RD), so this is a comparison of workflows, not of solvers.

## E5. Collision-matrix structure
* Irreducible matrix (31x31x17, n = 31104): asymmetry before (A+A^T)/2 = 0.13, 0.32, 0.43, 0.52,
  0.57 % of max|A| at 100, 200, 300, 500, 1000 K; exactly symmetric after. PCG on the
  symmetry-allowed subspace: 19 iterations at every T; preconditioned Ritz values 0.59–4.0;
  minimum curvature p.Ap/p.p = 3.8e-4 THz (> 0).
* Energy vector (scalar grid, `results/physop_*.json`): ||Omega e||/||D e|| = 2.3–2.4 for phono3py's
  matrix at every mesh/integration; physical Omega' = D + C1 - (C0+C2)J gives 8.4e-3, 4.9e-3,
  3.7e-3 (tetrahedron 5x5x3, 7x7x5, 9x9x5; 0.11, 0.047, 0.037 after symmetrization) and
  8.7e-4, 5.5e-4 (Gaussian 0.1 THz, 5x5x3, 7x7x5). (Omega - Omega')e = 2(C0+C2)e to 1e-15;
  Omega = Omega' on odd vectors to 6e-16; they differ by 12–55 % on even vectors.

## E6. Isotope scattering (natural abundance; g2(Al) = 0, g2(N) = 1.84e-5)
LBTE 300 K: 302.08 -> 301.65 (xx), 285.56 -> 285.17 (zz), i.e. -0.14 %/-0.14 %.
RTA: -1.24/-1.05 % (100 K), -0.20/-0.16 % (200 K), -0.11/-0.09 % (300 K), -0.05/-0.04 % (1000 K).

## E7. Independent checks
* RTA re-summation (own heat capacities, POSCAR volume, 6mm star average): <= 8.5e-15 relative
  (no-isotope and isotope runs, 100–1000 K), mode by mode <= 8e-16; with exact SI-2019
  constants -1.0e-6 to -1.6e-6 (phonopy constant set).
* LBTE rows+PCG vs phono3py dense pseudo-inverse: 4e-13 (9x9x5), 1.2e-12 (15x15x9); assembled
  matrix bit-identical at 9x9x5.

## E8. Cost (6 threads)
RTA 31x31x17 all T: 38–54 min, 0.4 GB. Rows+PCG LBTE 31x31x17: 47 min (1 T), 2 h 38 min (5 T;
rows 70 min, assembly 8–12 min/T, PCG 5.4 min/T), RAM < 1 GB, 7.2 GiB disk per T. Rows at 300 K:
13 s (11x11x7), 58 s, 211 s, 614 s, 1524 s (27x27x15). phono3py dense LBTE 15x15x9: 233 s
(dsyev 183 s); 31x31x17 dense not feasible here (7.2 GiB, ~13 h estimated).
