# Independent proof/unit review — completed 2026-09-25

**Verdict: VALID for the recorded single-point Gaussian reconstruction; no unit, mesh-factor or channel-sign defect found.** This is a source-formula/accounting check, not validation of canonical events. Read C's source trace, reconstruction script/JSON, pilot metadata, A's conditional coefficient statement, and the relevant installed 4.5.0 source sections. No other review was read; no material calculation or repository edit was performed.

## Audited implications

1. **Mesh normalization — VALID under the default conversion.** `interaction.py:177–190` contains `1/prod(mesh_numbers)` in the conversion to eV². `imag_self_energy.py:767–795` subsequently multiplies by triplet multiplicity, with no extra mesh divisor. Thus using exported `pp` and weights exactly once is correct. The recorded weight sum is 27. Another division by 27 would duplicate the normalization; the JSON's 0.962963 scaled error is the expected algebraic consequence, not independent physical evidence. Custom interaction-unit overrides would change this conclusion.

2. **Channels and units — VALID.** `triplets.py:240–252` and the self-energy accumulation give exactly

       (n1+n2+1) d(f0-f1-f2)
       +(n1-n2)[d(f0+f1-f2)-d(f0-f1+f2)].

   The script preserves both absorption signs. At positive frequencies and temperature, Bose occupation decreases with frequency; for an even Gaussian decreasing with absolute argument, the absorption difference has the same sign as n1-n2. There is no hidden negative-rate channel here. The numerical conversion `18*pi/(hbar*2*pi*1e12)^2` converts eV² times reciprocal-THz integration weights to THz: its mixed-unit coefficient has units THz²/eV². No extra universal factor 1/2 belongs in this ordered linewidth sum.

3. **Lifetime conversion — VALID in the linewidth convention.** Gamma is a half-linewidth expressed in ordinary THz. Half-width to population lifetime contributes 2, ordinary to angular frequency contributes 2*pi, and THz to s^-1 contributes 1e12. Hence `tau^-1=4*pi*1e12*gamma`, consistent with the [official gamma documentation](https://phonopy.github.io/phono3py/input-output-files.html#gamma). This does not identify collective collision eigenvalues or a bare reversible-event coefficient.

## Remaining limits / severity

- **Moderate if generalized:** cutoff `1e-4 THz` is hardcoded; it matches the inspected API default and the source's strict daughter-frequency mask. The script assumes all 12 external bands, the first temperature, default conversion, and this Gaussian recipe. It is not a general importer for altered cutoffs, selected bands, integration truncation, tetrahedra or unit overrides.
- **Major if promoted to event validation:** the recorded maximum error, 6.94e-18 THz (scaled 2.96e-16), verifies the stipulated pp-to-linewidth sum using shared frequencies/constants/quadrature. Gaussian-weighted detuned tuples need not satisfy detailed balance or eventwise energy conservation. Multiplicities do not establish full-mode permutations or unordered-event counting.
- A's Eq. (2) is explicitly conditional on unfolded, consistently weighted events. This executed check neither proves nor refutes that conditional coefficient; its applicability and counting remain unvalidated. C correctly leaves canonical export conversion unresolved.

**Disposition:** accept the narrow reconstruction/convention check. No canonical-event normalization, conserving material operator, cutoff/mesh convergence, or material lifetime accuracy is established. No additional numerical run was needed for this bounded source audit.
