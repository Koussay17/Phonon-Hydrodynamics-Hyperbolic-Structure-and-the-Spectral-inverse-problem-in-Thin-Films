# Bounded independent review of the clustering fix

**Status: COMPLETE, 2026-09-26. The reported overlapping-block defect is resolved in the portable routine for the checked cases.** This follow-up inspected `frequency_blocks.py`, `check_frequency_blocks.py` and the modified `reconstruct_complex_triplet.py`. It executed only the small helper tests and compared saved arrays; no material calculation was rerun.

The helper sorts frequencies into disjoint connected groups and rejects a group whose total span reaches the tolerance. The original adversarial chain (1, 1+0.75e-10, 1+1.5e-10) THz, embedded in the same 12-mode example, now raises `ValueError: ambiguous chained near-degeneracy; specify physical blocks`. All eight supplied regressions pass. Source inspection confirms that the portable script also compares A/B partitions, checks overlap singular values, and checks the assembled block transformation's unitarity.

For the actual saved triplet, the new groups equal the prior groups on all three legs, and the A/B partitions match. Recomputed block rotations are bitwise identical to the legacy saved rotations. The minimum overlap singular value is 0.9999999999999982; the largest block-unitarity Frobenius residual is 9.355067438917448e-16.

I independently compared the saved hardened replay at `D:\ResearchLab\scratch\complex-triplet-hardened-20260926` with the original artifact. All NPZ arrays are bitwise identical, and all common reported metrics are identical, including full transport 2.147819025539619e-15 and block-only transport 2.457103956627176e-14. The replay's helper/script hashes match the inspected files; the legacy script hash remains unchanged.

Evidence is persisted in `experiments/R-clustering-fix.json`. This closes the specific numerical overlap objection for the portable workflow. The partition remains a numerical convention, not a certificate of physical exact degeneracy. The original audit's source-independence, normalization, channel-precision and material-scope limitations are unchanged.
