# AlN interaction-export pilot
Question: can one tiny pinned phono3py calculation export squared cubic interactions, phonons and full reciprocal-grid maps with enough verified normalization and counting metadata for subsequent physical-event reconstruction?
Scope: a coarse pilot and one or a few grid points, not a converged BTE calculation.
Inputs: complete Phonon Olympics phono3py family at commit 0640f07735059be9717a7565c2a0f22dc0da7a17, locally on D.
Previous finite-event reference and hostile reviews: C:/Users/Koussay/these/theory/aln/events.
Success: reproducible executable pilot with hashes, actual array shapes/units/maps, independently checked counting identities, explicit energetic mismatch/integration limits and source-to-rate conversion status.
Falsification: incompatible force-constant shapes or atom maps, inconsistent triplet multiplicity, missing units, unphysical energy-conservation claims or unsupported backend equivalence.
Do not call raw pp or gamma_detail a complete physical collision operator.
