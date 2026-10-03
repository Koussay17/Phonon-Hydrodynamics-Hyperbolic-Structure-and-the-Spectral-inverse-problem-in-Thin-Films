# Paper 1 — prior-art audit and central contribution (3 October 2026)

Targeted audit of the candidate claims of the theory/methods paper. The seven-phase protocol was applied
per claim: direct terms, mathematical and physical equivalents, historical terms, backward and forward
citations, neighbouring fields. Searches covered Crossref, OpenAlex, arXiv full text and publisher and
repository PDFs; Semantic Scholar was largely rate-limited. The full literature trail is kept with the
campaign record. "Potentially novel" always means *within the literature searched*.

| Claim | Verdict | Closest prior work |
|---|---|---|
| C1 — conserving closure, coefficient 1/3 (note 14) | **KNOWN** | Guo & Wang, Int. J. Therm. Sci. 171, 107178 (2022), Eq. 7, including the τ_c form; Sendra et al., PRB 103, L140301 (2021) and PRB 106, 155301 (2022); Hardy & Albers, PRB 10, 3546 (1974) |
| C2a — lifetimes, invariants and DC response do not fix the memory moment (note 18) | **PARTIALLY KNOWN**: the mathematics is elementary (Stieltjes indeterminacy); an explicit phonon identifiability statement was not found | Saunderson et al., SIAM J. Matrix Anal. Appl. 33, 1395 (2012); Cepellotti & Marzari, PRX 6, 041013 (2016); Simoncelli et al., PRX 10, 011019 (2020) |
| C2b — fixed event cone: τ_mem ≤ M²K(0) (note 18) | **Inequality KNOWN** (Ben-Tal & Teboulle, Linear Algebra Appl. 139, 165 (1990)); **application to kinetic memory potentially novel** | Stewart (1989), O'Leary (1990), Todd (1990), Vavasis & Ye (1996); computing such constants is NP-hard in general (Tunçel 1999) |
| C3 — compression vs closure (note 21) | **KNOWN** (exact lumping) | Wei & Kuo (1969); Aoki (1968); Li & Rabitz (1989); Buchholz (1994) |
| C4 — conserving resonance weak form (notes 19, 23) | **PARTIALLY KNOWN**; the specific construction is incremental | Mielke (2011); Maas & Mielke, J. Stat. Phys. 181, 2257 (2020); Abdelmalik & van Brummelen (2016); Tran et al., J. Differ. Equ. 269, 4332 (2020) |
| C5 — non-Fourier extension of the thermoreflectance gauge (note 7) | **PARTIALLY KNOWN**; the extension to Cattaneo/GK and the full chain are potentially novel | Krapez & Rigollet, J. Appl. Phys. 122, 066101 (2017); Hennessy & Myers (2021); Beardo et al. (2020); Camacho de la Rosa et al. (2025) |

Fourier resonance predates Kovács (2018): the condition is named in Both et al., J. Non-Equilib. Thermodyn. 41, 41 (2016).

## Central contribution (recommended)

For linear, energy-conserving phonon kinetics in entropy coordinates:
1. Static transport data (DC conductivity tensor, mode lifetimes, collision invariants) do not determine the
   dynamic part of the closure. The first memory moment τ_mem = −K′(0)/K(0), which fixes the Cattaneo
   relaxation time, is unbounded in the positive-semidefinite class.
2. If the collision events are fixed by the harmonic spectrum and only their nonnegative rates are unknown, then
   τ_mem ≤ M²K(0), with a sharp, rate-independent geometric constant M. This is presented as an application of
   Ben-Tal & Teboulle.
3. A one-dimensional front-face thermoreflectance measurement identifies the closure only through the
   gauge-invariant times (τ_R, τ_ℓ), not through the nonlocal length or the closure coefficient.

C1 and C3 become attributed background; C4 moves to the quantitative AlN paper.

## Results still needed before submission
1. The identifiability gap in a real event geometry: certified M and the feasible range of τ_mem at fixed lifetimes
   and K(0), first for a Debye event set, then for the AlN event set, with mesh and temperature dependence.
2. The continuum statement: memory diverges when the infrared rate exponent satisfies 2α ≥ d, so M diverges under refinement.
3. The FDTR consequence: translate that range into phase differences for the 500 nm AlN scenario
   (10 kHz–200 MHz) against 0.01–0.1° noise, using the twin (`src/twin_fdtr.py`).
4. A written proof that the 1D Cattaneo/GK response has generic rank four in (ξ₁, b, τ_R, τ_ℓ), and of what breaks the gauge.
5. An attribution rewrite of the background sections, then proof audit, numerical audit and hostile review.

Target journals: Physical Review B, if results 1 and 3 show a quantitative gap; Journal of Statistical Physics,
for the operator-theoretic framing; Journal of Applied Physics, if the measurement side leads.
