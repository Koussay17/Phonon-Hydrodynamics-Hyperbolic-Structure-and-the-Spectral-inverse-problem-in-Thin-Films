# Second-generation importer/kinetic contract

Cross-examination of A/C/D and `experiments/PI-amplitude-scope.md` plus `complex-triplet.json`; saved after interruption. Finite positive-frequency modes, positive temperature, declared linearized model. No implementation or new full kinetic model.

## Four distinct products

| Product | Meaning and required evidence |
|---|---|
| **Recorded-basis population events** | Incoherent Bose kinetics in the saved eigenbasis. Require mode identity, signed incidence, rate/counting convention, resonance/integration provenance, and an explicit population/secular approximation. |
| **Invariant block compression** | Projected instantaneous response or entropy Galerkin model. Require complete block sums, common channel/Bose prefactors, block dimensions, susceptibilities, and reconstruction. This does not establish autonomous dynamics. |
| **Closed block kinetics** | Autonomous block totals for the declared initial-data class. Require a supplied generator and a closure certificate; distinguish exact algebra from numerical residual evidence. |
| **Density-matrix action** | Dynamics of matrix occupations and matrix observables. Require complex amplitudes or sufficient contraction maps AND a specified kinetic approximation, including weak coupling, statistical/bath assumptions, secular treatment, units, and conservation/positivity scope. Amplitudes alone are not a generator. |

## Mandatory guards

Preserve source/input hashes, backend, units/overrides, q/branch and reciprocal maps, eigenvectors/projectors, cutoffs/omissions, integration metadata, and exact versus numerically clustered degeneracy. Polarization gauge blocks in a homogeneous crystal preserve q modulo a reciprocal vector as well as energy. Mixing different q requires transforming translation/momentum operators and reconsidering the homogeneous state space. A clustering tolerance is not an exact-degeneracy certificate.

**Passive covariance is not physical symmetry.** A basis change transforms amplitudes, states, observables, and the generator together; it need not commute with the fixed interaction. General degenerate-block unitaries generate coherences. Only permutations and phases preserve the whole diagonal population algebra. Rotating squared entries alone or discarding coherences after rotation changes the information/model. Symmetry-based closure needs the actual event set and state covariance, with representation multiplicities checked.

**Amplitude tuples are not yet events.** Record annihilation/creation orientation, reciprocal daughter labels, conjugation conventions, and the physical resonance/momentum condition. An all-incoming tuple q0+q1+q2=G cannot be assigned population incidence merely from its array order. Keep forward/reverse reactions, reciprocal partners, daughter permutations, repeated blocks, and identical modes separate. Declare ordered tensors versus normalized symmetric daughter-pair spaces. Accumulate repeated indices before the incidence outer product; derive their prefactors explicitly.

## Closure and normalization certificate

For dot x=-Lx and block totals N=Bx, the reversible event model has W_ii=n_i^0(1+n_i^0). Define

    Sigma=BWB^T,  H=WB^T Sigma^(-1),  G=BLH.

Then BH=I. G is a compression. Autonomy for every initial x additionally requires

    BL(I-HB)=0, equivalently BL=GB.

For an exact-degenerate block A, Sigma_AA=d_A n_A^0(1+n_A^0) and (HN)_i=N_A/d_A. Mean occupations require the correspondingly rescaled generator. In reversible entropy coordinates, closure is vanishing retained/discarded coupling. For general nonreversible matrix dynamics, arbitrary-state trace closure and preservation of initially isotropic matrices are different one-sided conditions. A small defect alone gives no long-time error certificate.

If closure fails, retain a labeled entropy/variational compression, polarization or bright/dark components, or a controlled fast-mixing approximation. Exact degeneracy supplies no secular phase averaging; discarded splittings must be resolved against collision rates over an admissible coarse-graining interval. Reaction generators need not have Markov-chain signs. Conserved energy/momentum and amplitude-intertwining charges provide checks, not closure automatically.

## Status of the recovered AlN triplet

The PI records complex-amplitude recovery for one nonzero-frequency triplet, squared-interaction agreement below 8e-16 on its global scale, and block-restricted transport error 2.46e-14 with 1e-10 THz clustering. Assign **amplitudes available; contraction checked for this scope**. Full-unitary transport allowing nondegenerate mixing is not an exact-degeneracy test. Shared reciprocal-force input, the r0 convention, selected-triplet coverage, and exclusion of Gamma limit the evidence.

No promotion to Hamiltonian-channel-normalized, full-population, closed-block, or density-matrix-generator status follows. Next: declare one channel and kinetic approximation; independently check orientation/counting and normalization, then conservation, passive covariance, and closure separately. No complete AlN collision action is certified.