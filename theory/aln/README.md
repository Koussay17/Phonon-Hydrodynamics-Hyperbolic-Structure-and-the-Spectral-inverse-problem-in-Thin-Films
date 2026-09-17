# AlN spectral inputs and results

## Sources and attribution

`rao_300K_modes.npz` is a converted subset of the dataset by **Xixin Rao, Yipeng Wu, Songcheng Li, Haitao Zhang and Chengdi Xiao (2025)**, *Multispeed lattice Boltzmann method with ab initio scattering rates for phonon non-equilibrium thermal analysis in GaN heterostructures*, Mendeley Data V1, [doi:10.17632/w9hg2mnnwy.1](https://doi.org/10.17632/w9hg2mnnwy.1).

The source data are licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: selected the raw AlN ShengBTE arrays, converted to SI, reshaped by q index and branch, and calculated the primitive-cell volume from the reciprocal lattice. No rates were fitted, extrapolated or relabelled as normal/resistive. Attribution does not imply endorsement. This converted data file retains CC BY 4.0; it is not covered by the project's code licence.

`olympics_temperature.json` transcribes numerical values from the author repository for **McGaughey et al. (2025), Phonon Olympics**, [doi:10.1063/5.0289819](https://doi.org/10.1063/5.0289819). The repository commit and exact source URL/checksum are recorded in `sources.json`. This table is an independent bulk calculation; do not combine it with Rao lifetimes as if both described the same calculation.

## Arrays

| Name | Shape | Unit / meaning |
|---|---|---|
| omega_rad_s | (793,12) | Angular frequency, rad/s |
| velocity_m_s | (793,12,3) | Cartesian group velocity at irreducible representatives, m/s |
| degeneracy | (793,) | q-point multiplicities; sum 13824 |
| q_fractional | (793,3) | Original fractional reciprocal coordinates; not a verified Wigner–Seitz-folded momentum basis |
| cell_volume_m3 | scalar | Primitive-cell volume |
| total_rate_s | (793,12) | Total diagonal RTA rate at 300 K, s^-1 |
| anharmonic_rate_s | (793,12) | Three-phonon diagonal rate, normal + umklapp combined |
| isotope_rate_s | (793,12) | Isotope rate |

Only symmetry-invariant scalar quantities can be summed directly on the irreducible representatives. This study uses v_z^2 and (v_x^2+v_y^2)/2. A general tensor or momentum-conserving collision model requires full-zone symmetry reconstruction and a consistent first Brillouin-zone convention.

Branch indices are retained. Acoustic character near Gamma does not establish a global LA/TA assignment at crossings. The three zero-frequency Gamma translations contribute their kB capacity limits and zero discrete transport.

## Reproduce

From the repository root:

```powershell
python -X utf8 -B scripts/prepare_aln_sources.py
python -X utf8 -B scripts/analyse_aln_spectrum.py
python -B -m pytest tests/test_spectral.py -q -p no:cacheprovider
```

The importer downloads only when its ignored `.build/research` cache is absent. The analysis and tests use the compact repository inputs and require no network or DFT installation. Source archives stay in `.build/`. `spectral_results.json` and figures 07–08 are generated from the compact inputs.

## Limits

- The Rao archive has only T300K, no separate N/U rates, no force constants or full collision matrix. Absorption/emission (`plus/minus`) is not N/U.
- The RTA reconstruction checks normalization and units, not convergence of the DFT or q mesh. Higher lifetime moments are especially sensitive.
- The lowest positive frequency is about 0.499 THz. About 12.97% of the basal transport weight has effectively zero normal velocity on this grid; the thin-film limit requires refinement.
- The film suppression is a steady **in-plane** RTA calculation between identical reflecting surfaces. It is not a cross-plane FDTR or hydrodynamic boundary model.
- The temperature table includes a missing 20 K row, preserved as null. Its 30–75 K values are archived but not used in the transport plot or certified here. The displayed range is 100–1000 K; higher-order scattering and thermal expansion are not added.
- Harmonic C(T) changes Bose weights at fixed frequencies. No temperature law is inferred for the 300 K lifetimes.
- Static and memory times are observable-specific RTA moments, not identified GK tau_N/tau_R.

## Next quantitative calculation

Use one consistent set of harmonic and cubic force constants (available in the Phonon Olympics author repository), converge q grids at each temperature, and output normal and umklapp rates using a documented first-zone convention. Keep isotopes/defects separate. Export frequencies, full-zone velocities, weights, rates and, for a dynamic closure, the collision operator or enough information to validate its projection. Do not mix the two existing datasets to fill missing fields. No such first-principles rerun has been performed in this project.

## Audited closure investigation

See closure/13-final-report.md, closure/errata.md and note 18. Four independent mathematical methods and four hostile reviews distinguish broad operator ambiguity from the tighter fixed finite event cone, for which an audited memory bound holds. Its geometric constant is not a converged AlN estimate.

scripts/analyse_aln_response.py writes bulk_response.json with input/software fingerprints and cancellation-free errors. It is a uniform bulk RTA diagnostic, not an FDTR simulation. Higher inverse moments are strongly concentrated in the lowest resolved modes; no new physical mesh sequence was generated.

The inspected phono3py v4.5.0 reducible Python/C formula is conductivity-equivalent on odd populations. Its scalar full-grid shape does not certify the physical even action. A signed event operator or separately validated even action is needed before energy/viscosity closure. This does not rule out other implementations. Exact source scope and export mutation/normalization caveats are in the reports.
