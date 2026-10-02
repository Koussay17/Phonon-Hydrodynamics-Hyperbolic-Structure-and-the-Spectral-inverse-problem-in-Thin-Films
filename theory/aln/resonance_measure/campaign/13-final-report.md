# Conserving resonance measure — final scientific synthesis (2 October 2026)

**Question.** How can energy conservation of three-phonon events be integrated without turning broadened weights into discretely resonant events, and which finite structure survives?

**Strongest justified conclusion.** Take a declared scalar Bose kinetic model with exact resonant quadrature nodes. Its finite entropy-variable weak form has a proved local structure:
- positive capacity;
- conservation of quadrature energy;
- nonnegative entropy production;
- Bose stationarity;
- a PSD linearization with the energy null vector.

Positive reweighting of unchanged off-shell events cannot restore conservation.

These identities do **not** constrain the rates or certify the measure. They do not preserve the domain globally: an explicit example leaves it in finite time. Nor do they establish an AlN operator. The implementation uses log-space evaluation, which an independent hostile numerical review verified.

**Not established.**
- global realizability;
- convergence of evolving solutions or spectra;
- a converged, globally enumerated AlN resonance measure;
- the excluded acoustic, critical and coherent sectors;
- the kinetic regime;
- novelty.

**Bridge to the material problem** (repository audit, findings F2 and F6).
- On the 3³ mesh, real pilot triplets have detuning ≥ 1e-3 k_BT, about 1e5 times the note-19 resonance cap. The material operator therefore requires actual surface integration.
- Coincident labels form a zero-measure set under a regular measure. A finite repeated-oscillator factor must therefore not be attached at a quadrature node.

**Next.** Declare the kinetic approximation for the AlN observable. Then build a globally enumerated resonance-surface quadrature on the real phonon dispersion (Olympics force constants), refining the mesh, the surface rule and the basis independently.

Reviews: 12-peer-review.md. Proof status: 10-proof-status.md. Claims: CLAIM_LEDGER.md.
Reproduction (14 checks): `python -B scripts/reproduce_resonance_measure.py --output-dir <new directory>` from the repository root.
