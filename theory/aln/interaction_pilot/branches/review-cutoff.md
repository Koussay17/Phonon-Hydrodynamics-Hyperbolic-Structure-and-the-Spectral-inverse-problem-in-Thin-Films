# Independent numerical review of the reciprocal-cutoff control

**Complete bounded follow-up, 25 September 2026.** Inspected `cutoff_control.py`, `cutoff-control.json`, all four saved `cutoff-control-*.npz` files, the original saved failed matrices, and the relevant installed parameter-handling source. No material calculation was run and no other files were edited.

## 1. Parameter handling: the recorded control varies the intended cutoff

The script first obtains the baseline dataset, then creates a fresh phono3py object for each factor. It copies the loaded NAC parameters, explicitly assigns both `G_cutoff` and `Lambda`, and does so before initializing the interaction/dynamical matrix. It retains NAC, the input paths, backend, q coordinates and other displayed settings.

Installed `phonopy/harmonic/dynamical_matrix.py` confirms:

- `Gonze_nac_dataset` returns `(Gonze_force_constants, dd_q0, G_cutoff, G_list, Lambda)`, so the script unpacks the intended fields;
- `_set_nac_params` consumes explicit `G_cutoff` and `Lambda` values;
- when Lambda is omitted, its default depends on the cutoff, so explicitly holding it fixed matters;
- the G list is selected by the reciprocal-vector norm cutoff.

Installed `phono3py/api_phono3py.py` passes the stored NAC parameters into dynamical-matrix initialization. The JSON records **actual dataset values**, not just requested labels: Lambda is exactly `0.26615773936469767` in all four records; actual cutoff ratios are 1, 1.25, 1.5 and 2; the G-vector counts rise from 311 to 579, 1013 and 2381.

No parameter-override or stale-object error was identified in this bounded source inspection. The fresh-object pattern and returned values make this a substantially better controlled experiment than changing a cutoff while allowing the default splitting parameter to move.

## 2. Independent checks of saved arrays

The factor-1 matrices are **bitwise identical** to the original `matrix-control-nac1.npz` matrices that exhibited the failed periodicity check. Thus the sequence starts from the actual failed case.

All saved matrices are 12x12 and exactly Hermitian in their stored entries. The phase vectors are identical across factors. Independently forming an explicit diagonal phase matrix and aligning the two BZ-equivalent matrices gives:

| Cutoff factor | G count | Relative periodic-matrix defect | Relative change of d0 to factor 2 | Relative change of d1 to factor 2 |
|---|---:|---:|---:|---:|
| 1 | 311 | 1.3539770654e-8 | 4.9469e-9 | 9.9773e-9 |
| 1.25 | 579 | 1.5466773666e-12 | 5.7133e-13 | 1.1275e-12 |
| 1.5 | 1013 | 4.7802e-16 | 1.6876e-17 | 7.0759e-18 |
| 2 | 2381 | 4.7896e-16 | Reference | Reference |

Norms are Frobenius norms. The approximately 5e-16 floor differs slightly from the script's elementwise phase multiplication because of floating-point operation ordering; this is not a change of tolerance. Scaled eigenpair residuals for the eight saved matrices range from **5.90e-16 to 9.57e-16**.

The important additional check is the fixed-q comparison, not only the periodicity defect. Merely making two matrices equal could leave a common error; here each matrix also stabilizes between factors 1.5 and 2 to differences below the overall binary64 matrix scale. Eigenvalue comparisons at that point are limited by eigensolver roundoff, approximately 1e-15 relative.

The recorded frequency differences are 6.8180993e-8, 4.9436011e-12, 1.7763568e-14 and 1.4210855e-14 THz. The **original unchanged 1e-10 THz threshold** therefore still fails for the original factor-1 case and passes for this selected pair at the larger tested cutoffs. No statement about other q pairs follows from that local pass.

## 3. What the experiment establishes, and where it stops

**Supported numerical interpretation:** with the splitting parameter fixed, enlarging the finite reciprocal cutoff removes the selected pair's defect, while the assembled matrices stabilize. Combined with the earlier small eigenpair residuals, this strongly supports **finite reciprocal-cutoff truncation in the assembled NAC calculation** as the source of the original local periodicity discrepancy. The data do not point to an eigenvalue-solver failure.

This is more specific evidence than the previous NAC-on/off association. It is still a result for these two matrices, this parameter sequence and these inputs. In particular:

- Four cutoffs and a floating-point plateau do not give a rigorous infinite-cutoff remainder bound or a uniform Brillouin-zone error bound.
- The irregular growth in discrete reciprocal shells and the roundoff floor do not justify a power-law observed order or Richardson extrapolation.
- Fixed-Lambda stabilization does not establish independence from Lambda, convergence of other real-space/reciprocal approximations, or convergence of the force constants themselves.
- The control uses the requested Rust configuration. It does not separately certify a larger-cutoff C export, although the earlier default-cutoff spectra showed a shared defect.
- It does not resolve the Gamma acoustic signs; those modes were not part of these two non-Gamma matrices.
- Small relative changes in a dynamical matrix are not automatic relative-error bounds on individual eigenvectors, small/vanishing cubic couplings, or narrowly integrated rates. Near degeneracies and cancellations matter.

The saved NPZ files contain matrices and phases, while q coordinates and parameter provenance reside in the supplied script/JSON. The control does not include fresh per-case input hashes. The original-baseline bitwise match is a strong consistency check, but this review has not independently reconstructed every run's complete provenance.

## 4. Original export remains failed; pp/gamma impact is unmeasured

The original phonon/pp/gamma export was produced using the default cutoff and still fails its recorded duplicate-frequency threshold. These harmonic matrix controls do **not** replace or retroactively validate those files. No larger-cutoff pp or gamma calculation was supplied or run in this review.

Consequently the previously observed scalar reciprocity and gamma-reconstruction checks remain evidence about the **original stipulated coarse calculation**. They do not establish that its interaction strengths or rates are unchanged when the cutoff defect is removed. That impact is explicitly **unresolved**.

### Final numerical classification

- **Observed convergence with respect to the tested cutoff, at the selected pair:** matrix changes become **floating-point limited** at factors 1.5-2.
- **Original default-cutoff consistency check:** still **failed**; its local mechanism is now strongly supported as reciprocal-cutoff truncation.
- **Eigenvalue solves on the saved matrices:** residuals are at binary64 roundoff scale; no evidence that they caused the original defect.
- **Whole-grid harmonic, cubic-interaction and linewidth accuracy:** **inconclusive** at enlarged cutoff because those validations were not performed.
- **Material convergence, mesh/broadening limits and physical uncertainty:** **unestablished**.

This review does not authorize substituting the larger-cutoff control matrices into the old export, merging incompatible settings, loosening the original threshold, or claiming an AlN transport result. It concludes the requested bounded audit without further experiments.
