# Second grid check: actual gp1-m333 export

2026-09-24. **PASS for the full q-grid orbit partition and all eight weights**, under the explicitly reconstructed unshifted, time-reversal-enabled, swappable grid. This is stronger than checking only the weight sum. No interaction values or rates were calculated.

## Independent method and provenance

Script: ../scripts/B_check_grid_orbits.py. Detailed machine-readable results, source hashes, input hashes, and reconstructed maps: ../runs/gp1-m333/B_grid_orbits.json.

Read the actual POSCAR, phonon grid addresses, reduced/full triplets and weights. Construct direct-cell symmetry from spglib at symprec=1e-5, convert each integer rotation R to R^(-T), add time reversal, and build permutations of Z_3^3. Group closure and orbit enumeration use independent integer permutation composition and breadth-first search. Neither exported weights nor native map_triplets are inputs to that enumeration.

Only afterward reconstruct the C-backend BZGrid and native triplet maps for comparison. Their addresses, triplet rows, weights, and partitions agree exactly with the independent result and actual HDF5 arrays.

Versions: phono3py 4.5.0, phonopy 4.5.0, spglib 2.7.0, NumPy 2.5.3, h5py 3.16.0. Reconstructed structure: P6_3mc, number 186. The check uses the C backend; no Rust equivalence is inferred.

## Exact result

Anchor BZ 1 is regular-grid index 1, address (1,0,0), q=(1/3,0,0). The reconstructed Q and P are identity and PS=0. Regular indexing is r=a_x+3a_y+9a_z modulo 3.

There are 12 direct rotations and 24 distinct reciprocal permutations after adding time reversal. The anchor little group has order 4 and 12 orbits. Adjoining the swap x -> -r-x gives order 8 and the following eight classes:

| Export row | BZ triplet | Complete orbit of regular q1 indices | Export weight | Enumerated |
|---:|---|---|---:|---:|
| 0 | (1,0,2) | 0,2 | 2 | 2 |
| 1 | (1,1,1) | 1 | 1 | 1 |
| 2 | (1,3,10) | 3,4,7,8 | 4 | 4 |
| 3 | (1,7,8) | 5,6 | 2 | 2 |
| 4 | (1,13,28) | 9,11,18,20 | 4 | 4 |
| 5 | (1,14,27) | 10,19 | 2 | 2 |
| 6 | (1,16,36) | 12,13,16,17,21,22,25,26 | 8 | 8 |
| 7 | (1,20,34) | 14,15,23,24 | 4 | 4 |

These disjoint sets cover all 27 indices. Every row passes orbit-stabilizer. As another exact count, the eight group elements have sorted fixed-point counts (1,3,3,3,9,9,9,27); their sum 64 divided by group order 8 gives eight orbits.

Both reconstructed native maps are idempotent. map_q reproduces the independent little-group partition and map_triplets reproduces the independent little-group-plus-swap partition. The reduced/full triplets and weights are identical between pp and gamma_detail exports.

## Grid and momentum checks

- The 39 BZ rows form 27 regular-grid fibers: 21 singleton fibers and six fibers of size three.
- Address-derived bzg2grg agrees with the reconstructed library map. bzg2grg[grg2bzg[r]]=r holds for every regular index. Dense gp_map differences equal fiber sizes.
- The full 27-row triplet array contains each regular q1 exactly once and has q2=-r-q1 modulo 3.
- Every reduced and full row satisfies exact all-plus momentum closure. In reduced rows 1 and 5, the displayed BZ addresses sum to reciprocal vector (1,0,0); the other reduced rows sum to zero. This is not an energy-resonance test.
- Time reversal is an involution with one regular fixed point, Gamma. The anchor maps to regular/BZ 2, not itself. Negating every independently enumerated class produces exactly the independently enumerated partition at anchor 2.
- BZ-alias frequency arrays differ by at most 6.818098086114333e-8 THz. This is a measured discrepancy, not exact identity or a band-gauge validation.

## What remains unsupported

The HDF5 files do not contain dedicated P, Q, PS, bzg2grg, grg2bzg, rotation, map_q, or map_triplets datasets. Their present values were reconstructed from the cell and matching grid configuration. Persisting those maps and flags is still necessary for a self-contained export.

No interaction tensor at anchor 2 was inspected: the negative-anchor result concerns its grid partition only. No symmetry sewing matrix, band permutation through degeneracies, eigenvector phase relation, band-resolved pp reciprocity, absolute normalization, or energy integration was established. Orbit weights count q-grid partners, not complete physical collision events.

## Reproduction

Run the B-prefixed script with the pinned environment:

    D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe scripts\B_check_grid_orbits.py --run-dir runs\gp1-m333 --output runs\gp1-m333\B_grid_orbits.json

Working directory: the current pilot campaign. The script reads existing HDF5 inputs and writes only its B-prefixed JSON. No PI script, repository file, or HDF5 input was modified.
