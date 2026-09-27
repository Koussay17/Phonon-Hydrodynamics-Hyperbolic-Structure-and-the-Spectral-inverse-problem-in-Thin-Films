# Channel orientation and counting

**WORK IN PROGRESS: hostile review incomplete. See [publication status](PUBLICATION_STATUS.md).**

Campaign 20260926-105330, continued 27 September 2026. Working copy on D; compact review archive in this repository. Raw force constants and amplitude arrays remain outside Git.

## Scope
Separate the ordered cubic Hamiltonian coefficient, collected monomial, Fock transition amplitude, spectral golden-rule density, geometric population closure and linearized entropy event weight.
Read 13-final-report.md for the final status and review resolutions; original branches remain unmodified historical evidence.

## Reproduce finite checks
From the repository root:
~~~powershell
python -B scripts/reproduce_channel_counting.py --output-dir D:/ResearchLab/scratch/channel-new-finite
~~~
Requirements: Python >=3.10, NumPy, SciPy and SymPy. This runs finite Fock matrices, exact symbolic/rational moment checks and complex sewing with a negative control. Choose a new output directory.

## Add pinned source and selected material checks
~~~powershell
python -B scripts/reproduce_channel_counting.py --output-dir D:/ResearchLab/scratch/channel-new-full --source-python D:/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe --pilot-root D:/ResearchLab/orchestration/campaigns/20260924-105051-aln-interaction-pilot
~~~
The source interpreter needs phono3py=phonopy=4.5.0 and the recorded source hashes. The pilot-root must contain the retained corrected-cutoff run at runs/export-cutoff-2.0/runs/gp1-m333 plus pilot-run.json. Its input and output hashes are checked before material contraction. No inputs are downloaded automatically.

Fresh reproduction passed all 13 top-level checks; reproduction-status.json records the scope. Source checks may be requested without the material pilot. Without optional flags, no material validation is claimed.

## Evidence
- branches/A-D: genuinely different isolated first passes.
- branches/A2, B2, C2: second-generation sewing and closure investigations.
- branches/L and L2: primary sources, inspected equations and coverage limits.
- experiments/: portable scripts, compact outputs and provenance.
- D-source-factor-check.py/json: pinned implementation comparison.
- REVIEW_CANDIDATE.md: frozen propositions sent to hostile reviewers.
- branches/red-*: independent attacks; final report resolves their objections.
- material_permutations_initial.py is retained as a failed general normalization shortcut; the reproduction runner uses material_permutations.py.

## Limitations
One selected tensor is not a full-mesh or continuum collision action. Finite geometric closure is conditional and is not invariant under the assumed exact Fock jump process. Numerical frequency clusters do not certify exact degeneracy. No experimentally validated AlN transport classification or novelty claim follows.
