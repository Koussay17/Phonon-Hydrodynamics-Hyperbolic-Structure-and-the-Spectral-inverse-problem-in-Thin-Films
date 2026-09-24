# Independent red-team proof audit: physical events

**Status: COMPLETED bounded review, 2026-09-24.** The essential reproduction and exact algebra checks executed successfully on 2026-09-23; a usage-limit rejection interrupted only the final report write. This save completes that write. No further exploration was performed during finalization.

Frozen sources reviewed: review-target.md; review-inputs/collision_events.py; review-inputs/test_collision_events.py; review-inputs/19_Validation_evenements_collision.tex; branch C's first-pass parity report. No other red-team report or discussion was used. No repository source or frozen input was changed.

## Verdict

**The exact resonant rank-one linearization is VALID under its stated hypotheses. The parity counterexample is also VALID.** No missing sign or counting factor was found in the explicitly defined distinct-mode nonlinear flux.

**The numerical admission contract is FALSE on an explicit accepted input.** The constructor reports zero thermal detuning for an exact mismatch of one in energy/(k_B T), and admission depends on daughter ordering. Its returned symmetric matrix then differs substantially from the true entropy-coordinate nonlinear Jacobian. This is an implementation-validation defect, not a counterexample to the exactly resonant theorem.

| Finding | Classification | Severity |
|---|---|---|
| Thermal-detuning check hides an order-one mismatch and admits a nonsymmetric true Jacobian | FALSE numerical guarantee | High within the advertised extreme-energy input contract |
| Module states wave-vector reversal is a separate event without the fixed-point exception | FALSE at self-reciprocal triplets; the note is correctly qualified | Moderate documentation/counting inconsistency |
| Detuned thermal energy production is written strictly positive while Gamma=0 is allowed | VALID only for Gamma>0; otherwise nonnegative | Minor domain omission |
| Exactly resonant C=A(D^-1 s)(D^-1 s)^T, invariants and entropy Hessian | VALID | No mathematical defect found |
| Odd-block invisibility and finite-k parity Schur argument | VALID | No mathematical defect found |

## 1. Earliest substantive failure: resonance admission

Target: frozen collision_events.py lines 49-54 and the stored diagnostic at line 74; the note's assertion that an additional absolute detuning check controls thermal mismatch.

At the default tolerance,

    DecayEventOperator([1e16, 1., 1e16], [[0,1,2]], [1.])

is accepted and records

    relative_energy_residual = -1.1102230246251565e-16
    detailed_balance_log_residual = 0.0.

All indices are distinct and all energies and the coefficient are positive and finite. Exact arithmetic on those represented floating-point inputs gives

    Delta = 1e16 - 1 - 1e16 = -1.

Thus the forward/reverse Bose-factor ratio is exp(1), approximately 2.718281828459045, not a roundoff deviation from one.

**Earliest failed inference:** line 51 computes `(x_i-x_j)-x_k`. Binary64 rounds the first subtraction `1e16-1` back to `1e16`, so the second subtraction produces zero. The relative check also passes. A cancellation-prone computed zero is therefore treated as certification of the absolute thermal mismatch.

Swapping the daughters changes the result:

    DecayEventOperator([1e16, 1., 1e16], [[0,2,1]], [1.])

is rejected because `(1e16-1e16)-1=-1`. The admission decision is not invariant under an exchange that leaves both the stated flux and the physical unordered-daughter event unchanged.

### Consequence for the claimed linearization

Write c=1/(e-1), approximately 0.5819767068693265. The returned binary64 matrix is

    [[1+c, 0, -(1+c)],
     [0,   0,  0],
     [-(1+c), 0, 1+c]].

An independent 100-digit evaluation of the true nonlinear-flux derivative, transformed with mathematically nonzero entropy scales, gives to this precision

    -J_entropy = [[1+c, 0, -c],
                  [0,   0,  0],
                  [-(1+c), 0, c]].

Omitted entries are exponentially small. The maximum entry error is 1.0000000000000002; the relative Frobenius error is 0.5932501380835192. The true derivative is nonsymmetric because detailed balance fails.

This changes representable order-one entries; it is not merely the loss of a tiny occupation. The frozen tests themselves exercise energies of order 1e16, so the input lies within the advertised numerical stress regime. This does not demonstrate a failure for ordinary moderate inputs or disprove the exact resonant theorem.

**Unestablished obligation:** the code must actually certify its declared absolute thermal mismatch on the supplied represented values, or state when it cannot. No fix was implemented during this review.

## 2. Exact event formula: audited derivation

Assumptions: finite positive energies and temperature; distinct modes; equal full-grid weights; an explicitly defined reversible event counted once; Gamma>=0; and exact x_i=x_j+x_k.

1. Since 1+N_mu=exp(x_mu)N_mu, the direct/reverse Bose-product ratio is exp(-x_i+x_j+x_k)=1. The two equilibrium activities A in the note therefore agree.
2. Cubic terms cancel exactly: F=Gamma[n_i(1+n_j+n_k)-n_j n_k]. Its three derivatives are Gamma(1+N_j+N_k), Gamma(N_i-N_k), and Gamma(N_i-N_j). Detailed balance converts these into -A s^T D^-2. Independent symbolic substitution N_j=u/(1-u), N_k=v/(1-v), N_i=uv/(1-uv) gives zero residual for all three entries.
3. With delta_n=D y and dot(n)=s F, the linearized equation is `dot(y)=-A(D^-1 s)(D^-1 s)^T y`. The damping sign and inverse-time units are correct.
4. For A>=0 the quadratic form is A[(D^-1 s)^T y]^2. If s^T h=0 then C(Dh)=0. For a sum of events, only active affinity vectors constrain the kernel. The note correctly allows accidental invariants and does not infer finite DC from this sparse model.
5. The convex negative Bose entropy has Hessian 1/[N(1+N)]. Relative entropy has the same Hessian; its quadratic term is one half of ||y||^2. There is no missing factor of two in C.
6. F already includes forward minus reverse reaction. A separately listed daughter permutation adds the same event again. There is no missing universal factor of one half because Gamma is defined for this precise counted event; conversion from a material convention is not claimed.
7. With u_mu=1-exp(-x_mu), the row magnitude is `sqrt(A)/d_l = sqrt(Gamma) exp[-(x_i-x_l)/2] u_l/sqrt(u_i u_j u_k)`. The logarithmic implementation is algebraically this expression. Its identification as the nonlinear Jacobian requires the detailed balance invalidated by the admission counterexample.

Verdict: **VALID for exactly resonant inputs in exact arithmetic.** A tolerated detuning is an approximation whose diagnostic must be trustworthy; the discovered case defeats that numerical guarantee.

## 3. Detuned energy sign and its zero-rate exception

The formula

    epsilon^T dot(n)
      = -k_B T Delta Gamma (1+N_i)N_jN_k [exp(-Delta)-1]

has the correct sign. For every nonzero real Delta, the product Delta[exp(-Delta)-1] is strictly negative. Therefore an active event Gamma>0 produces strictly positive thermal energy change. Opposite signs of detuning do not cancel these contributions when weights are nonnegative.

The note globally permits Gamma=0, however. Such an inactive event has zero energy production and leaves the thermal state stationary. Consequently the unconditional inequality is nonnegative; strictness requires an active detuned event. This is a minor domain omission, not a sign failure for positive rates.

## 4. Parity and event counting

Let J preserve energies and commute with D. Equal coefficients for distinct reciprocal partners make their summed operator commute with J. In a self-reciprocal event, positive resonant energies make the parent uniquely highest in energy. Reversal cannot exchange that parent with a daughter; it can swap the daughters. Hence the stoichiometric and entropy-affinity vectors are even, and the event outer product vanishes exactly on the odd subspace.

With JVJ=-V, streaming has only off-diagonal parity blocks. Exact elimination gives

    S_odd=z+C_odd+k^2 V_oe(z+C_even)^(-1)V_eo.

The positive k^2 sign follows from the two factors i k in the subtracted Schur product. Re z>0 ensures invertibility despite accidental invariants. Thus equal C_odd does not in general determine the finite-k response.

For the three-mode example the even block is [[a,-sqrt(a b)],[-sqrt(a b),b]], with a=2g/w_q and b=g/w_0. Its relevant inverse entry is `(z+b)/[z(z+a+b)]`, giving the reported exact scalar response. At real z>0 and k!=0, its dependence on event strength is nonzero.

An independent solve with the frozen operator at z=0.25,k=0.3 gives current responses 3.5673008347861406 and 3.593303644699138 after doubling the self-reciprocal event. The odd-action change is 1.60e-17 numerically and zero analytically. This establishes the specific counterexample, not a claim that every even-sector error must affect every current.

### Source-contract inconsistency

The module docstring, lines 5-6, says wave-vector reversal is a separate event without the fixed-point exception correctly stated in the note, lines 49-52. For energies [1,1,2], reversal exchanging modes 0 and 1 maps event [2,0,1] to [2,1,0], the same unordered-daughter event. Listing both doubles C. The constructor sums both rows and does not identify reciprocal orbits.

This is a documentation inconsistency, not a claim that the constructor automatically generates wrong events. The note's counting statement is correct; the unqualified module statement can direct a reader to introduce precisely the hidden even-sector error discussed in the note.

## 5. Executed evidence, omissions and disposition

Saved independent artifacts:

- experiments/red_proof_event_checks.py
- experiments/red_proof_event_checks.json

The script imports only the frozen source, checks the exact represented-input mismatch, independently differentiates the resonant flux symbolically, computes the counterexample Jacobian at 100 digits, and checks the self-reciprocal three-mode example. **All assertions executed and passed.** Environment: Python 3.14.7, NumPy 2.5.3, SymPy 1.14.0, mpmath 1.3.0. The JSON records the source hash and outputs.

Reproduction:

    py -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260917-232509-aln-physical-events\experiments\red_proof_event_checks.py

**Unexecuted extras:** the entire repository suite was not rerun; no seven-mode campaign-wide numerical revalidation, literature attribution audit, material import, or broader physical closure investigation was performed. None is required to reproduce the specific admission defect or assess the stated elementary proof.

**Final disposition:** the exact event theorem and parity example survive. The advertised absolute-detuning admission guarantee does not survive its explicit counterexample and remains unvalidated until corrected and checked. The two smaller domain/documentation issues remain as stated. No repository or frozen-target edits were made.
