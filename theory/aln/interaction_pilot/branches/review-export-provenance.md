# Export provenance review — completed 2026-09-25

**Verdict: VALID for the saved 1.0/1.5/2.0 cutoff comparison.** No stale-export, changed-input, or weakened-inspection defect was found. This bounded review addresses provenance and executed checks only; it does not assess the physical interpretation of degenerate-mode interaction entries.

## Audited evidence

- **Fixed Lambda:** each generated script is exactly the original pilot plus one NAC assignment before mesh/interaction initialization. All specify `Lambda=0.26615773936469767`; only `G_cutoff` changes, from 1.1964811966806446 to factors 1.5 and 2.0. The injection anchor occurs exactly once. Generated-script hashes match their recorded hashes. The installed API stores the supplied NAC dictionary for initialization; this is not a post-export metadata change.
- **Identical inputs:** independently rehashed POSCAR, BORN, fc2.hdf5 and fc3.hdf5 in all three run directories. Every hash matches the source input and its per-run record. Other pilot settings remain unchanged, including no force-constant symmetrization.
- **Fresh executions:** `mkdir(exist_ok=False)` prevents reuse of each work directory. Only the four input files are copied. Each subprocess starts a new interpreter, loads a new Phono3py object and writes under its own script-relative directory. Nonzero exit codes abort the comparison. Saved output hashes all match the recorded hashes. No old phonon, pp or gamma export is copied into the new runs.
- **Original inspections really pass:** at factors 1.5 and 2.0, both `inspect_pilot.py` and `reconstruct_gamma.py` are byte-identical to their campaign originals. The driver executes them and checks their exit status. Both saved inspection JSONs have every original check true; both reconstruction JSONs pass. The original strict duplicate-frequency tolerance remains `1e-10 THz`.
- **Independent read-only recomputation:** all 18 duplicate pairs give defects 6.818098086114333e-8 THz at factor 1 and 2.4868995751603507e-14 THz at both larger cutoffs. Triplet indices match. From the HDF5 arrays, the 1.5-to-2.0 pp maximum scaled change is 0.1097100507683089 and the maximum branch-relative gamma change is 1.384998383719554e-13, reproducing the summary.

## Minor gaps; no effect on these results

The driver does not assert unique text replacement; this audit verified it explicitly. Per-run settings omit the effective NAC values, so provenance requires retaining the generated script and comparison JSON. The summary uses `<=1e-10`, whereas the original inspection uses `<1e-10`; the present results satisfy both. Original inspections are rerun only for factors greater than one, so the claim must remain specific to those corrected exports.

**Disposition:** accept the stated identical-input, fixed-Lambda, fresh-run comparison and corrected-export inspection passes. Stable gamma does not by itself certify entrywise pp convergence; that interpretation is outside this review. No new material runs were performed. Only this report was written.
