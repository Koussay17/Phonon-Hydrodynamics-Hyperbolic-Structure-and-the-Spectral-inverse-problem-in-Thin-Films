# 09 — Failed approaches, bugs and negative results (kept on purpose)

1. **Single-run ascent for tau_max is unreliable.** The maximisation landscape on F(r, K0) has
   many local maxima (different sets of suppressed events). For 1D N = 9 the multiplicative ascent
   from the reference stops at 509.8 (a boundary KKT point with 8 events switched off), SLSQP from
   another start reaches 627.5. In the scans the best maximum is typically reached by only 1 of
   about 50 starts once n >= 28. Consequence: every reported tau_max is a LOWER bound on the
   supremum (inner interval); no claim of global optimality is made for the maximum.

2. **Scipy sparse overhead.** The first implementation spent most of its time constructing sparse
   temporaries (78k evaluations at ~1 ms of overhead each for a 16-mode problem). Replaced by
   bincount assembly (fastops.py), about 3x faster with identical results (checked on 1D N = 9,
   2D N = 3: same extrema to 1e-9).

3. **Slow convergence at boundary optima.** Minimisation runs crept to the 1500-iteration cap
   (multiplicative steps approach g_alpha = 0 slowly). A stall criterion (relative change of N
   below 1e-6 over 30 accepted steps) was added; the minima changed by < 1e-5 relative.

4. **Augmented-Lagrangian L-BFGS on log-rates (alb.py) is worse than the multiplicative
   projected gradient** on the AlN tied problem: from the reference, tau_x max 105 ps and min
   22.85 ps after 6 outer iterations, versus 540 ps and 19.3 ps after only 200 multiplicative
   steps. The log parametrisation stalls when rates must go to zero. Not used for results.

5. **AlN event export bug (found and fixed).** phono3py's `CollisionMatrix._get_gp2tp_map`
   returns, for q1 outside the irreducible set of the little group of q0, the third phonon of the
   REPRESENTATIVE triplet (`tp2s[gp1] = tp2s[map_q[gp1]]`). phono3py only uses its frequency
   (rotation invariant), so the baseline Omega' rebuild was unaffected. Using it as the identity of
   the third daughter assigned some events to rotated modes: entry-wise comparison with Omega'
   showed ratios clustered at exactly 0, 1, 2, 3 and sign flips (off-diagonal relative difference
   0.75). Fix: q2 = -q0 - q1 (mod G) from the grid addresses. After the fix: off-diagonal relative
   difference 0.073 (remaining: band-gauge in degenerate subspaces, conserving stoichiometry,
   4-sigma truncation), kappa and tau_mem within 0.01-2 % of Omega'.

6. **Degeneracy was not the explanation of the pre-fix mismatch.** Restricting the comparison to
   non-degenerate pairs (thresholds 1e-4 to 0.2 THz) left the 72 % mismatch unchanged; the integer
   ratios pointed to a counting/identity error, which item 5 confirmed.

7. **The 1D exact model cannot test the infrared-exponent law.** Collinear normal events are
   exactly resonant for all same-sign pairs (a 2D region of pair space) while umklapp resonances
   lie on lines, so the umklapp fraction vanishes like 1/N; crystal momentum becomes a
   near-invariant, and K0/N and tau_ref grow linearly in N even for alpha = 0.25 (where tau_RTA
   converges). This is the 1D anomalous-conduction mechanism, not an infrared-rate effect.

8. **KL projection onto prescribed diagonals occasionally fails** (1D alpha = 0.5, N = 99;
   alpha = 1.5, N = 51; 2D alpha = 0.5, N = 17; alpha = 1.0, N = 7; 3D alpha = 1.0, N = 9): the dual
   Newton does not reach 1e-13 within 200 iterations (prescribed r close to the boundary of the
   cone W R_+^m, or poor scaling). These points are reported as missing, not imputed.

9. **Exact enumeration of M does not scale.** 1D N = 9 needs 9.66e6 subsets (331 s); the next
   sizes need C(32,14) = 4.7e8 (2D N = 3) and C(64,26) ~ 1e17 (1D N = 15). Only lower bounds by
   basis pivoting are available beyond 1D N = 9 (the pivoting search recovers the exact 1D N = 9
   value).

10. **Operational.** A broad process termination used to stop the first scan also terminated
    unrelated Python processes on the machine; later jobs were managed by PID only. An import of
    the read-only baseline helper module regenerated its bytecode cache
    (`theory/aln/transport_baseline/scripts/__pycache__/common.cpython-314.pyc`, 04:07) in the
    thesis tree; no source file was modified, and all later runs use a local copy (`aln_common.py`,
    identical SHA-256) with bytecode writing disabled.
