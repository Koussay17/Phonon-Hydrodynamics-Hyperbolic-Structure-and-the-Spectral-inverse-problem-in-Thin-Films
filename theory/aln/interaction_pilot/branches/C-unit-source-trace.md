# C: phono3py 4.5.0 interaction-to-linewidth source trace

2026-09-24; independent first pass. Scope: dimensional/source tracing, without reading other campaign reports. **Established:** standard `pp` is an energy-squared interaction including the mesh normalization; the source-to-half-linewidth formula below is explicit. **Unresolved:** a complete, independently checked conversion from symmetry-reduced exports into canonical signed physical events.

## 1. Version, interaction and mesh normalization

Inspected installed phono3py **4.5.0** and phonopy **4.5.0**, under `D:/ResearchLab/envs/aln-phono3py-4.5.0/Lib/site-packages/`. Upstream pin: phono3py v4.5.0, commit `21fa8f3817fbcc603254656f525bb5aec113afb6`. Installed `interaction.py` and `imag_self_energy.py` match the previously fetched pinned sources after CRLF-to-LF normalization. Their normalized SHA256 values are respectively `08fd28b07cd4a73af24c5dc4ce3348edea847bf1fdfdb496bb7d2e1d6b0d8b19` and `5c47e9953bf42b06aa266c0e82148a3421540b8918baa734e315e4691df9988c`.

**SOURCE SAYS:** [interaction.py:177-190](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/interaction.py#L177) explicitly converts to **eV^2**, with factors `1/36`, `1/8`, `(2*pi*THz)^(-3)`, `AMU^(-3)`, and **`1/prod(mesh_numbers)`**. This is the standard default conversion, assuming compatible force-constant, mass and length units; an explicit conversion override changes this contract. The factor is applied after the C and Rust interaction calls (lines 1000 and 1038); the Python path computes `abs(fc3_normal)**2 * unit_conversion` (1094-1096).

**PRIMARY THEORY:** Togo, Chaput and Tanaka, *Distributions of phonon lifetimes in Brillouin zones*, Phys. Rev. B **91**, 094306 (2015), [DOI](https://doi.org/10.1103/PhysRevB.91.094306), [arXiv:1501.00691v3, Eqs. (10)-(13)](https://arxiv.org/html/1501.00691v3). Their energy-valued cubic interaction contains `1/(3! sqrt(Nq))` and three oscillator factors `sqrt(hbar/(2m omega))`; squaring explains `1/36`, `1/8`, and `1/Nq`. Thus **do not divide exported `pp` by Nq again**. The array axes are `(triplet, selected_band0, band1, band2)` (`interaction.py:254-260`). "Full" retains all those interaction entries subject to frequency cutoffs; it does not mean all external grid points were computed.

## 2. Explicit linewidth formula and lifetime conversion

Let nu denote numerical ordinary frequency in THz; let p denote numerical `pp` in eV^2, and m_t the exported triplet multiplicity. Define quadrature weights d0 for `nu0-nu1-nu2`, d+ for `nu0+nu1-nu2`, and d- for `nu0-nu1+nu2`, using the code's reciprocal-THz convention. The inspected Python accumulation is

\[
\gamma_{0b_0}=C_\Gamma\sum_{t,b_1,b_2}m_t p_{t b_0b_1b_2}
\{(n_1+n_2+1)d_0+(n_1-n_2)(d_+-d_-)\},
\]

\[
C_\Gamma=\frac{18\pi}{[\hbar_{\mathrm{eV\,s}}\,2\pi\,10^{12}]^2}
=3306215.697005045,\qquad
\tau^{-1}[\mathrm{s}^{-1}]=4\pi\,10^{12}\gamma[\mathrm{THz}].
\]

Sources: [imag_self_energy.py:200-207 and 767-795](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/imag_self_energy.py#L200); primary paper Eqs. (11), (13); [official gamma documentation](https://phonopy.github.io/phono3py/input-output-files.html#gamma). The two factors in `4*pi` are **2 for half-linewidth** and **2*pi for ordinary versus angular frequency**. Equivalently, tau in ps is `1/(4*pi*gamma)`. This is the phonon lifetime relation; identifying it with an RTA transport relaxation time is the paper's approximation, Eq. (18).

A constant-only import in the installed environment gave hbar = `6.582118985531608e-16` eV s, THzToEv = `0.00413566733` eV, and `4*pi*1e12*Cgamma = 4.1547131779577225e19`. No interaction/BTE calculation was run by this branch.

[triplets.py:130-133](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/triplets.py#L130) asserts `sum(m_t)=prod(D_diag)`; lines 240-252 verify the Gaussian channel signs above. There is **no additional 1/2** in this source linewidth accumulation. An unordered-event representation requires a separate input-exchange/repeated-index counting derivation; do not insert a universal 1/2 into this formula.

## 3. Concrete tiny export path and written fields

For an already validated `Phono3py` object `ph3` containing matching FC2, FC3, structure, masses and NAC settings, the official [API, lines 2414-2447](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/api_phono3py.py#L2414) supports this one-external-point pilot:

```python
ph3.mesh_numbers = [3, 3, 3]
ph3.init_phph_interaction()
gp = int(ph3.grid.grg2bzg[1])
ph3.run_thermal_conductivity(
    temperatures=[300], grid_points=[gp],
    is_full_pp=True, write_pp=True, write_gamma_detail=True)
```

This is source-verified syntax, **not an executed pilot in this branch**. `imag_self_energy.py:239-245` takes the full path without passing an integration-zero mask. The [writer, file_IO.py:1220-1262](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/file_IO.py#L1220) then writes `version`, `pp`, and supplied triplet metadata. The [conductivity caller, utils.py:476-507](https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/conductivity/utils.py#L476) supplies **`triplet`, `weight`, `triplet_all`**, but **does not supply `triplet_map`**, although the writer accepts it. The sparse-mask branch instead packs selected values and must not be treated as the same full export.

## 4. Remaining barrier and minimal next step

Persist phonons/eigenvectors and explicit BZ/regular-grid maps separately, alongside settings, cutoffs and hashes: the pp file is not a complete event archive. Boundary-equivalent BZ entries require mapping to unique modes. Symmetry multiplicities reproduce the scalar sum; they do not by themselves specify the mode permutations needed for a full event matrix. Orient all-plus momentum triples consistently with physical decay/absorption, including wavevector reversal of the output mode. Repeated input indices require their own counting check.

**Minimal next check:** reconstruct one exported mode's gamma from pp, the actual integration weights, Bose factors and multiplicities using Section 2; compare to the program's gamma before constructing event weights. The pp-to-canonical-event prefactor, permutation expansion and repeated-index rule remain unresolved here. Gaussian/tetrahedron quadrature weights do not establish exact resonance of each stored mode triple; finite off-resonance events cannot simply inherit exact energy conservation or detailed balance. Raw pp and gamma_detail are not a complete physical collision operator, and this coarse pilot cannot establish converged AlN transport.

Access limit: GitHub MCP failed with an invalid-session response earlier; pinned raw source, installed source, the primary preprint and official documentation were accessible. No broad novelty search or backend-equivalence claim was made.
