# Principal-investigator amplitude experiment
Pinned phono3py source retains complex normal-mode amplitudes immediately before taking absolute squares. The present script invokes the package real-to-reciprocal prototype with its actual r0-average setting and independently contracts atom/cartesian indices with saved eigenvectors by numpy.einsum. It reproduces the Rust-exported squared interactions for triplet4 at both corrected cutoffs to below 8e-16 globally scaled.

Two different checks must not be confused:
1. Full-unitary basis transport is a contraction consistency identity, allowing nondegenerate mixing.
2. Projection of overlaps to frequency-degenerate blocks (1e-10 THz clustering), followed by block polar factors, reproduces the complex tensor to 2.46e-14 relative Frobenius error. This is the more restrictive observation supporting degeneracy-basis redistribution for this triplet.

No full collision generator, canonical channel orientation or independent absolute normalization is inferred. Near-degenerate numerical grouping is not a theorem of exact crystal degeneracy. The same reciprocal force tensor and physical model are shared; this is a cross-implementation contraction test, not independent material physics. The selected nonzero-frequency triplet avoids Gamma; the script has not been validated for arbitrary zero-mode triplets.

Script, saved tensors, source and input hashes are recorded. This result is withheld from independent first-pass branches until all reports are complete.
