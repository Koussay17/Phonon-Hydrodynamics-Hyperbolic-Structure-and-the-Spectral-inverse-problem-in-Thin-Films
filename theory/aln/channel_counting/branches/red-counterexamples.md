# Independent counterexample review — channel counting

**Status: completed bounded audit, 27 September 2026.** Read the frozen candidate, note22 draft, A/B/C/D and A2/B2/C2 derivations, synthesis, and the cited sewing/material/source scripts and results. No other red report was read. Only this report and the uniquely named small sector probe were written.

## Verdict and severity

No in-scope counterexample was found to the conditional 6V/3V coefficients, 1/(1+delta_ab) event factor, repeated geometric Jacobian, or frequency-compatible sewing identities. The counterexamples below attack the proposed extensions and stronger interpretations, not the candidate's explicit scope restrictions. Failing to falsify the algebra is not proof.

## 1. B2: a reversible finite sector does not give autonomous means — HIGH extension risk

For c=1, the sector ell=2 consists of (m,n)=(1,0),(0,2). Each transition has rate 2, so their equal mixture is exactly stationary. Its means are (p,a)=(1/2,1), and its exact net flux is zero. The geometric mean equation at those same means instead predicts

    J_geo = 2[p(1+2a)-a^2] = 1.

Thus even a stationary, reversible, finite-sector distribution defeats that scalar closure. This does not refute the full sector process: it demonstrates what is lost when replacing the process by its means.

A stronger near-equilibrium obstruction survives removal of the mean-energy direction. Use B2's Bose equilibrium with daughter ratio 1/2, parent ratio 1/4, means (P,A)=(1/3,1), and define

    phi = m-1/3 -(4/9)(n-1).

Under the product equilibrium pi, phi is centered and orthogonal to ell-<ell>, but its sector-conditioned mean is nonzero:

    E(phi|ell=2k)   = (k+2)/18,
    E(phi|ell=2k+1) = (k-6)/18.

The sector weights are (3/8)(k+1)4^-k and (3/16)(k+1)4^-k, respectively. Each sector is irreducible and reversible. Therefore the equilibrium autocorrelation has the exact long-time plateau

    <phi^2> = 68/81,
    lim_(t->infinity) <phi,exp(t G)phi> = 34/729,
    normalized plateau = 1/18.

These follow by summing the displayed geometric series. The rank-one geometric compression instead predicts complete decay of its non-energy mode. An exponential tilt of pi by phi provides positive nearby initial laws; this is not an invalid negative-probability perturbation.

**Consequence:** retaining an exact sector process is consistent, but mean energy alone does not specify its stationary information. A number-response relaxation time, unique thermalization, or decaying reduced response cannot be inferred from the positive compressed matrix. Exact projection memory must carry the missing sector information; calling it memory does not remove that requirement.

## 2. A2: the channel/Sym² object is not an event generator — HIGH extension risk

The symmetric pair representation passes the elementary multiplicity attack. For a repeated-daughter monomial g p(a†)^2, the exact-degenerate change a†=(b†+c†)/sqrt(2) gives

    (g/2) p[(b†)^2+2 b†c†+(c†)^2].

The normalized vacuum-pair amplitudes are g/sqrt(2), g, g/sqrt(2); their squared sum is 2|g|^2, equal to the original. Hence the normalized Sym² construction does preserve this strength.

What fails is promoting those three entries to autonomous scalar events: the same rotation transforms an occupied a mode into a state with b/c coherence. Dropping that coherence changes the statistical model. The channel operator preserves amplitude information; it does not supply thermal factorization, geometric factorial moments, dephasing, or a Markov limit.

Likewise, a finite isolated resonant channel restricted to the one-parent/two-daughter subspace has coherent population oscillations cos²(|T|t/hbar), not irreversible exponential decay. T† as the reverse map establishes Hermiticity, not detailed-balance Markov rates.

**Material gate:** in the saved selected tensor, all 1472 retained first-leg entries fail the stated resonance threshold. Strict resonant projection removes those entries in that orientation. Their excellent conjugation residual does not make them physical resonant events; below-floor entries and other orientations remain unclassified.

## 3. Repeated-label and singular-limit interpretation — MEDIUM

The per-channel ratio (1,2) is not a demonstrated finite correction to continuum AlN transport. For a fixed parent, identical full daughter labels require 2q_a=q_p+G. Under a regular smooth continuum resonance measure these isolated momentum solutions have zero surface measure. A finite repeated-channel contribution in a thermodynamic limit requires a separately justified singular measure, occupation scaling, or finite-size regime. This is a conditional limiting objection, not an AlN result.

At zero temperature W becomes singular, so the displayed entropy coordinates and compression cannot be extended by substituting W^(-1) at the vacuum boundary. B2 explicitly excludes that limit. Approximate frequency clusters likewise do not establish the exact intertwining needed by A2.

## 4. Failed attacks and unresolved scope

The sewing contraction cancels by unitarity with the stated conjugations; replacing the conjugation is already caught by the negative control. The unordered-pair norm identity correctly includes the repeated Fock normalization. The parent/daughter derivative mismatch cannot be removed by one extra common factor without spoiling the parent comparison. These attacks supplied no defect within the declared assumptions.

The numerical material checks share reconstructed arrays and source conventions. Their relative covariance/permutation agreement does not independently establish absolute normalization, global event enumeration, exact degeneracy, or an energy-conserving material integration. No actual failed saved numerical check is alleged here.

## Executed provenance and unfinished work

experiments/red_counter_sector_probe.py ran with exit code 0. Exact Fraction arithmetic verified the two-state stationary counterexample and enumerated complete sectors ell=0,...,120. Errors versus the analytic variance and plateau were 1.99e-32 and 1.78e-34; omitted probability was 1.74e-35. The probe used no material input or repository code.

All work required for this bounded review is complete. Material collision construction, continuum limits, microscopic Markov derivation, and an independent proof audit of the additional plateau calculation remain unestablished; they were not performed or credited. No novelty is claimed.

