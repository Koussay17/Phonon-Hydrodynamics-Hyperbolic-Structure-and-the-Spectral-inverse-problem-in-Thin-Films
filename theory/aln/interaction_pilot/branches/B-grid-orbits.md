# B — Finite-group checks for the grid export

Independent first pass, 2026-09-24. Read only 00-question.md and upstream sources; no peer/prototype results. No repository edits. **These are acceptance checks, not validation of an unseen export.** Scope below: unshifted grid, PS=0; retain and inspect shift metadata before applying them.

## 1. Quotient and BZ representatives

Let D=diag(d1,d2,d3), N=d1*d2*d3, and G=Z^3/DZ^3. For integer address a use the reference index

    I(a)=(a1 mod d1)+d1(a2 mod d2)+d1*d2(a3 mod d3).

With unimodular Q, fractional reciprocal coordinates are q=Q D^(-1)a. BZ boundary rows can represent the same element of G. The grid source distinguishes these rows from regular-grid points. [phonopy v4.5.0 grid.py](https://github.com/phonopy/phonopy/blob/v4.5.0/phonopy/phonon/grid.py)

Require integer positive d, integer unimodular Q, map bounds, and exact identities:

- I(addresses[b])=bzg2grg[b] for every BZ row b.
- bzg2grg[grg2bzg[r]]=r for every 0<=r<N.
- Every regular index occurs; each fiber contains congruent addresses modulo D.
- Do **not** require grg2bzg[bzg2grg[b]]=b: the return selects one boundary representative.
- In dense-map mode, check gp_map endpoints and successive differences against fiber sizes. Do not impose that layout on sparse mode.

Thus N, not len(addresses), counts independent grid populations. Duplicate representatives should have matching frequency multisets; eigenvectors require gauge/subspace-aware comparison.

## 2. Time reversal and its fixed points

On G, tau(a)=-a mod D. Construct it independently from addresses, then require a permutation and tau(tau(r))=r. Its fixed-point count is

    F=product_i gcd(2,d_i),
    number of point orbits=(N+F)/2.

A canonical BZ return need only reverse q modulo an integer reciprocal vector, not reproduce the literal negative BZ row.

Two independently enumerated arithmetic controls, with spatial point-group reduction disabled:

| D_diag | N | Time-reversal fixed points | Gamma triplets after partner swap |
|---|---:|---:|---|
| (2,2,2) | 8 | 8 | 8 weights, all 1 |
| (3,3,1) | 9 | 1 | 5 weights: 1,2,2,2,2 |

The first grid cannot exercise distinct q versus -q. For a shifted grid the map becomes a -> -a-PS modulo D; the above fixed-point formula and three-leg closure must be rederived.

## 3. Triplets: exact momentum and orbit weights

The pinned triplet source uses q0+q1+q2=G_recip and BZ indices; its maps use regular-grid indices. It constructs weights as map-fiber counts and checks their sum equals N. [phono3py v4.5.0 triplets.py](https://github.com/phonopy/phono3py/blob/v4.5.0/phono3py/phonon3/triplets.py)

For every row (b0,b1,b2), require the prescribed regular anchor r=I(a0), and exactly

    (a0+a1+a2) mod D=0,
    G_recip=Q[(a0+a1+a2)/D] in Z^3.

Keep the BZ representatives and this reciprocal vector. This is momentum bookkeeping, not energy resonance. A physical decay convention needs an explicit sign/time-reversal conversion; do not substitute q0=q1+q2 silently.

At fixed r, each x=q1 determines q2=-r-x: there are exactly N ordered choices. Reconstruct the effective group of grid permutations from the exported address rotations and time-reversal flag. Use address rotations, not fractional-coordinate matrices without conversion. Verify bijectivity and group closure. Its little group H_r fixes r. If swapping is enabled, add

    s_r(x)=-r-x.

Enumerate orbits under K_r=<H_r,s_r>; omit s_r when disabled. Do not multiply little-group counts blindly by two: s_r may already act as an element of H_r, notably at Gamma.

Acceptance is the **entire partition**, not only a sum:

- map_q equals the H_r partition; map_triplets equals the K_r partition.
- Both maps have valid representative labels and are idempotent.
- Each exported triplet row represents exactly one map_triplets fiber.
- Its weight equals that fiber's independently enumerated cardinality.
- All weights are positive integers and sum to N.
- Orbit-stabilizer, |orbit(x)|=|K_r|/|Stab(x)|, gives a second count.

With swapping alone, the number of rows is (N+f_r)/2, where f_r counts solutions of 2x=-r modulo D. Each coordinate contributes gcd(2,d_i) solutions if that gcd divides r_i, otherwise none. With all reductions disabled, require N rows and unit weights.

## 4. Self-reciprocal triplets and counting limits

Negation maps anchor r to -r; it is an internal operation at fixed anchor only when 2r=0. A single generic anchor therefore cannot test its companion without exporting -r.

For any declared triplet equivalence, canonicalize complete keys first, then compare key(t) with key(-t). On a negation-closed set of M keys with F_t fixed keys, there are (M+F_t)/2 reciprocal orbits. Expanding one representative from each of J orbits gives 2J-F_t keys. Never double self-partners.

At Gamma, an unordered **q-only** pair {q,-q} is self-reciprocal. This need not hold for full mode keys with different band labels. Degenerate bands require a supplied symmetry/subspace convention; grid orbit weights are neither band multiplicities nor forward/reverse-event factors.

## Handoff

Export D_diag, P, Q, PS, addresses, both grid maps, gp_map/layout flag, address and reciprocal rotations, triplets, weights, both reduction maps, anchor, symmetry/swap flags, backend, package versions and source hashes. Phono3py 4.5.0 imports its grid machinery from phonopy and permits phonopy 4.5.x; pin the actual dependency too. [Pinned dependency declaration](https://github.com/phonopy/phono3py/blob/v4.5.0/pyproject.toml)

**Earliest failure to reject:** treating BZ rows as distinct regular populations or treating a total triplet weight as proof that the orbit partition is correct. These tests establish indexing/counting only; they do not establish interaction normalization, energy integration, or physical collision reconstruction.
