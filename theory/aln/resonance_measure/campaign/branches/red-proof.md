# Independent red-team proof audit — completed 2026-09-28

**Verdict:** B2's finite local theorem is VALID under its stated hypotheses. No counterexample to that theorem survived this audit. The compressed coarea statement in note 23 has a trace-regularity gap; a secondary tube-mass rate in A2 needs a density bound. Neither issue refutes the finite ODE theorem. Global dynamics, correct physical measure, operator convergence, and material kinetics remain unestablished.

**Scope:** REVIEW_CANDIDATE.md; B/B2, A/A2, C2; note23-draft.tex; C2 script, saved JSON and failure record. No other red reports were read. No material calculation or replay of the candidate experiments was performed.

## 1. Finite theorem: exact audit trail

**Domains and local existence — VALID.** The finite evaluation matrices and weights must be fixed and finite, volume weights strictly positive, Phi full column rank, and every required entropy evaluation strictly positive. Then occupations and capacities are finite and positive. For v!=0,

    v^T M v = sum_l m_l n_l(1+n_l) (Phi_l v)^2 > 0.

The coefficient domain is open and nonempty because beta e lies in it when the represented energies are strictly positive. On this domain M is smooth and invertible. The logarithmic mean is smooth even at equality: for A,B>0,

    Lambda(A,B) = integral_0^1 A^t B^(1-t) dt.

Consequently -M^(-1)K alpha is locally Lipschitz and has the asserted unique maximal solution while it remains in the domain. This argument provides no uniform conditioning, boundary control, positivity between unsampled evaluations, or global lifetime.

**Flux and entropy signs — VALID.** Bose identities give A/B=exp(-b alpha), so F=A-B=-Lambda b alpha. Parent loss therefore yields dot U=-sum omega b F=K alpha, while DU=-M gives M dot alpha=-K alpha. Since ds/dn=xi, the chain rule gives dot S_Q=alpha^T dot U=alpha^T K alpha>=0. The sign is consistent throughout A, B and note 23.

**Represented energy and support obstruction — VALID.** The SAME fixed coefficient vector e must reproduce the chosen energy at every volume and reaction evaluation. Then b e=Delta and

    e^T K e = sum_active omega Lambda Delta^2.

Strict positivity of each active omega Lambda makes K e=0 equivalent to every active Delta=0. At alpha=beta e, dot E_Q=beta sum omega Lambda Delta^2. Positive and negative detunings cannot cancel. This is an obstruction for unchanged event rows and the common-rate Bose law, not for arbitrary modified mobilities or alternative energies. Exact thermal examples with beta=log 2 and energies (1,1,1) and (3,1,1) give Delta=-1,+1 and heating 2,4/7, respectively.

**Equilibrium derivative and kernel — VALID.** Because K(alpha)e=0 identically throughout the domain, DK[h]e=0. Thus at alpha_*=beta e the derivative of K(alpha)alpha is K_*h; the allegedly missing derivative-of-K term really vanishes. The derivative of M^(-1) contributes nothing because the equilibrium RHS is zero. Therefore the coefficient Jacobian is -M_*^(-1)K_*, and conjugation by M_*^(1/2) gives the claimed symmetric PSD generator. Positive active weights imply ker K_*=ker b; invertible metric congruence preserves nullity. Its energy null vector is M_*^(1/2)e. Empty reaction sets and extra invariants are admitted; no spectral gap or unique equilibrium follows without the stated additional rank condition. Under ker b=span(e), stationary states are exactly beta e with beta>0. The finite positive-energy sum E_Q(beta) decreases strictly from infinity to zero, so a fixed positive energy selects one beta.

**Quadrature versus physical energy — VALID with the stated distinction.** Conservation concerns E_Q=e^T U. It is not conservation of the continuum integral unless that integral is used consistently or its quadrature error is controlled along the whole trajectory. The surrogate defect bound 3a and the bound 2a E_h/epsilon_min follow by the triangle inequality and positivity. Agreement of surrogate and physical energies at volume nodes can make their quadrature functionals identical while their reaction evaluations and Bose fields differ. B2 explicitly preserves this distinction.

## 2. Concrete objections and counterexamples

**INCOMPLETE AS WRITTEN — note 23, section 1, “integrable amplitude f.”** Ordinary volume L1 integrability alone does not define a trace on the prescribed resonance surface. For Delta(x,y)=y on [-1,1]^2, f=0 and f=1_{y=0} are the same volume L1 element but give surface integrals 0 and 2. A finite surface integral is also not guaranteed by volume integrability. The displayed fixed-level delta/coarea identity needs a specified trace integrable against dS/|grad Delta|, or suitable continuous/smooth amplitude assumptions. A's continuous/smooth regular setting supplies stronger hypotheses, but the compressed note omits this qualification. This is a regularity gap, not a counterexample to the finite algebra.

**VALID UNDER AN EXTRA DENSITY BOUND — A2 section 2 tube rate.** Smooth transverse geometry gives codimension d, but local integrability of the weighted surface density gives vanishing tube mass, not necessarily O(ell^d). In regular coordinates normal to the daughter diagonal, density |z|^(-a), 0<a<d, gives tube mass proportional to ell^(d-a). A locally bounded density, together with uniformly regular compact geometry, supports the displayed O(ell^d) rate. Do not export that rate to merely integrable singular weights. The zero-mass diagonal conclusion itself is valid under local integrability.

**FALSE extension — energy preservation certifies occupation interpolation.** An exact independent counterexample uses b=(-3/2,1/2), e=(1,3), t=1001/1000 and

    alpha=log(2)e+2log(t)b.

The volume Bose occupations are N1=1/(2/t^3-1)>0 and N2=1/(8t-1)>0. Wrong occupation interpolation gives parent (N1+N2)/2 and two daughters N1, hence F_bad=(N1+N2)/2+N1*N2>0. Yet b alpha=5log(t)>0 and e^T b=0. Therefore dot E_Q=0 while dot S_Q=-5log(t)F_bad<0, approximately -0.0035877701352. This supports the candidate's exclusion and uses exact signs, not a numerical threshold.

The regular coarea Jacobians, chart-coverage qualification, and the normalized Gaussian sampling formula in note 23 have the correct factors. For fixed nonzero v, the Gaussian formula's zero-shift error contains 2exp(-2pi^2 r^2)/|v|, proving necessity of r->infinity for uniform-in-shift convergence in that specific model; the convergent Fourier tail gives sufficiency. This does not establish any universal mesh-width law, critical integrability, or nonlinear evolution limit.

## 3. C2 correction and exact checks

The revised diagnostic follows by subtracting its two heating expressions:

    H-beta sum omega Lambda Delta^2
      = sum omega Lambda Delta (affinity-beta Delta).

Its absolute-value term is therefore justified algebraically. The added 32 machine-epsilon allowance is a numerical diagnostic, not a proved floating-point bound. Reading the saved JSON reproduces the first original-predicate failure at rho=-0.0001 and confirms that all recorded revised bounds pass. The final script's SHA-256 matches its saved JSON. No original pre-correction source snapshot was independently audited, so the historical assertion that only the predicate/diagnostics changed is not independently certified by this branch.

Independent `experiments/red_proof_resonance_checks.py` completed with **PASS: 21 exact zero-residual checks**, including the capacity/dual identities, forward/reverse ratio, zero-affinity limit, two-node equilibrium Jacobian, coefficient linearization, kernel and eigenvalue 85/24. It also records the exact negative-entropy example and opposite-detuning heating examples. Results: `experiments/red_proof_resonance_checks.json`. Candidate code was inspected but not imported or rerun.

**Completed scope:** bounded finite proof audit and small exact checks. Not performed: material calculations, global-domain proof/counterexample search, acoustic/critical classification, operator/time-evolution convergence proof, independent literature verification, or formal proof-assistant verification. None is silently counted as established.
