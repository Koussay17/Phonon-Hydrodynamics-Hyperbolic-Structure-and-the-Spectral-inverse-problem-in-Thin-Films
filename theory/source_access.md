# Primary-source access record — 16 September 2026

This records what was actually consulted. Downloaded papers remain in the ignored local .build/research folder; no new third-party PDFs are added to the repository. The later spectral extension adds an attributed CC BY numerical subset, not a paper PDF.

| Source | Access and use |
|---|---|
| [Hennessy & Myers, 2021](https://doi.org/10.1007/978-3-030-64272-3_2) | Accepted manuscript obtained from [Oxford](https://ora.ox.ac.uk/objects/uuid%3A4a96e55f-51ae-4c33-a41f-bd6b5dafc9d8/files/sk06988134?download=1). Equations (2), (11)–(18), Gaussian averaging and boundary conditions checked. Page 9 rendered to check the factor-three discrepancy in the prose. Claims about that discrepancy concern this manuscript, not a checked final publisher version. |
| [Beardo et al., PRB 101, 075303 (2020)](https://doi.org/10.1103/PhysRevB.101.075303) | [Full paper](https://ddd.uab.cat/pub/artpub/2020/42018670bd68/physrevb_a2020v101p075303.pdf) and boundary-condition appendices consulted. A silicon study; not evidence validating the illustrative AlN sample. |
| [Sendra et al., PRB 106, 155301 (2022)](https://doi.org/10.1103/PhysRevB.106.155301) | [Full manuscript](https://ddd.uab.cat/pub/artpub/2022/3f5b5b33661a/Sendra2022_Preprint.pdf) consulted. Equation (21) gives the Debye 3D coefficient 1/3 and discusses the historical coefficient 2. |
| [Struchtrup, 2025](https://doi.org/10.1007/978-3-031-93918-1_27) | [Institutional record and abstract](https://dspace.library.uvic.ca/items/f16d29ee-f219-4c08-ad1a-4bcdfcbb7895) consulted. The file link leads to a copy-request form. No request was sent. Not used as a substitute for reading the full chapter. |
| [Lebon & Dauby, 1990](https://doi.org/10.1103/PhysRevA.42.4710) | Publisher abstract and existing project reference. [Institutional PDF](https://orbi.uliege.be/bitstream/2268/59807/1/1990_PhysRevA_42_4710.pdf) downloaded but scanned; text extraction is empty. No new detailed claim is based on an unread scan. |
| [Krapez, IJHMT 99, 485–503 (2016)](https://doi.org/10.1016/j.ijheatmasstransfer.2016.03.122) | Publisher abstract/introduction accessible. Full article not obtained. Detailed article-specific comparison remains incomplete. The Darboux calculation in note 15 is independently derived and tested, not represented as a transcription of this article. |
| [Krapez, JPCS 745, 032059 (2016)](https://doi.org/10.1088/1742-6596/745/3/032059) | [HAL record](https://hal.science/hal-01394669) located through HAL API; author-deposited indexed excerpts consulted. HAL/ONERA/IOP PDF requests did not return a usable PDF. This is a separate conference paper. |
| [Camacho de la Rosa et al., 2025](https://doi.org/10.1063/5.0257299) | Main paper already in repository and reading note. Supplement remains unavailable. No claim of exact reproduction of all material tables or multilayer figures. |
| [Cahill, RSI 75, 5119–5122 (2004)](https://doi.org/10.1063/1.1819431) | [Author's institutional record](https://experts.illinois.edu/en/publications/analysis-of-heat-flow-in-layered-structures-for-time-domain-therm/) confirms the established layered Gaussian-beam framework. The implemented normalization is explicitly derived in note 16 and tested against an independent analytic half-space solution. |

## What the new work establishes

- Note 14 completes the linear grey conserving moment closure under its stated asymptotic assumptions.
- Note 15 compares the accessible sources and identifies exactly which full-text comparisons remain blocked.
- Note 16 and the FDTR code provide a synthetic preparation study. Actual sample geometry, calibration, noise and measurements are still needed.
- Note 17 now adds a published 300 K spectral dataset and a separate bulk conductivity temperature table; the quantitative normal/resistive dynamic closure remains open.

## Spectral AlN extension — 16 September 2026

- [Rao et al. dataset, DOI 10.17632/w9hg2mnnwy.1](https://data.mendeley.com/datasets/w9hg2mnnwy/1): archive retrieved through its public file endpoint, published SHA-256 verified. Raw AlN ShengBTE files and MATLAB sampling script inspected. Only T300K is present; no N/U split. A converted CC BY subset is included under theory/aln with attribution; the full archive stays ignored.
- [McGaughey et al., Phonon Olympics](https://doi.org/10.1063/5.0289819): accessible [NIMS paper](https://mdr.nims.go.jp/filesets/3113ca3b-841a-4fe0-bc12-d2b28c9eea1c/download), AlN section, tables XII–XIII, convergence discussion and data-availability statement consulted. The [author repository](https://github.com/McGaughey-Lab/Phonon-Olympics) was inspected at commit 0640f07735059be9717a7565c2a0f22dc0da7a17. The numerical temperature spreadsheet and associated CONTROL were read; rate-file header inspected only. Force constants were located but no BTE calculation was rerun.
- [Morelli and Slack chapter](https://djena.engineering.cornell.edu/hws/slack_high_thermal_conductivity_solids.pdf): public chapter consulted, notably table 2.3 and discussion of AlN defects. Its global parameters were not substituted for branch-resolved N/U rates.
- [Ravichandran and Minnich (2016)](https://doi.org/10.1103/PhysRevB.93.035314): [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevB.93.035314) accessible; steady-state film formulation and boundary conditions consulted. Our modal stationary suppression is independently derived and checked against transport ODE integration; the paper's transient calculation is not reproduced.
- [ShengBTE official documentation page](https://www.shengbte.org/documentation) links to the Bitbucket README; direct retrieval of that README was unsuccessful. Output-unit descriptions were read in a [public code mirror](https://github.com/wxmwy/ShengBTE), then independently checked by reconstructing the original dataset's heat capacity and RTA conductivities. No software version-specific N/U output capability is assumed.

Note 17 provides a 300 K spectral baseline and a sourced bulk conductivity trajectory. It does not claim a normal/resistive temperature trajectory, a converged AlN GK closure, or a cross-plane film prediction.

## Closure audit — 17 September 2026

Detailed read scope, versioned code pointers and access failures are in [L](aln/closure/reports/L-literature.md), [L2](aln/closure/reports/L2-literature.md), [peer sources](aln/closure/reports/red-peer-sources.md), and [independent parity check](aln/closure/reports/operator-parity-check.md).

- [Relaxon erratum](https://doi.org/10.1103/PhysRevX.10.049901): complete one-page publisher PDF inspected. It reverses a Matthiessen inequality, not the zero-mode statement.
- [Viscous heat equations](https://doi.org/10.1103/PhysRevX.10.011019): primary derivations and assumptions inspected.
- [Chaput author manuscript](https://arxiv.org/abs/1303.4062v1): physical versus odd-equivalent collision operators compared independently with the [phono3py v4.5.0 source snapshot](https://github.com/phonopy/phono3py/tree/21fa8f3817fbcc603254656f525bb5aec113afb6). The inspected Python/C path does not certify even-sector physics. Alternate backends and an actual exported AlN operator were not validated.
- [Fugallo et al.](https://arxiv.org/html/1212.0470v2): Appendix A supports positive event contributions; the campaign's finite-cone bound has its own audited proof.
- [arXiv:2605.17947v2](https://arxiv.org/abs/2605.17947v2) and [arXiv:2606.17829v1](https://arxiv.org/abs/2606.17829v1): exact version dates and targeted equations inspected. These remain preprints and are not wurtzite-AlN validation.
- Current research did not resolve the previously recorded Krapez full-text or Camacho supplement gaps.

No novelty claim, full physical event export, converged temperature/rate calculation or experimental validation follows from this source access.

## Event-level continuation ? 24 September 2026
[The source report](aln/events/branches/L-event-sources.md) records checked
Fugallo equations, exact-resonance/broadening scope, phono3py interaction-export
fields and pinned Phonon Olympics input directories.
Five public AlN input files (POSCAR, BORN, fc2.hdf5, fc3.hdf5, README) were
downloaded to D from commit 0640f07735059be9717a7565c2a0f22dc0da7a17.
Sizes and SHA-256 hashes were verified again after D was reconnected.
Only the manifest is published. No material calculation or source-to-event
unit/counting conversion was performed.
