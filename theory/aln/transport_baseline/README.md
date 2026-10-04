# AlN transport baseline (phono3py 4.5.0) — October 2026

Validated first-principles pipeline for wurtzite AlN with the Phonon Olympics force constants
(commit 0640f077, hash-pinned). Start with `13-final-report.md`; tables in `results/tables.md`
and `results/tables_production.md`; data inventory with SHA-256 hashes in `results/saved_data_inventory.json`.

Key results (300 K, W/(m K), in-plane / cross-plane):
- Exact reproduction of the Olympics phono3py team at 31×31×17 with 2.x-compatible settings:
  RTA 252.992 / 231.947, LBTE 285.048 / 271.262 (difference below 5e-4 %).
- Production (phono3py 4.5.0 defaults): LBTE 302 / 286, RTA 270 / 247; the fc3 phase convention
  (`make_r0_average`) alone shifts κ by about 6 %, the largest material-model uncertainty found.
- Natural isotopes: −0.14 % at 300 K. Normal processes carry 69 % / 61 % of the κ-weighted linewidth.
- phono3py's exported collision matrix is valid on odd (current) populations only; the physical operator
  Ω′ = D + C1 − (C0+C2)J conserves energy and must be used for even-sector and film problems.

Large arrays (about 52 GB: 31×31×17 collision matrices at five temperatures) are kept outside the repository;
their paths and hashes are listed in the inventory. `scripts/lbte_ooc.py` is the out-of-core LBTE route that
reproduces phono3py's dense solution to 1e-12.
