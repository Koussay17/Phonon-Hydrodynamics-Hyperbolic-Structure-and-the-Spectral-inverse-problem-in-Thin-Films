# Novelty protocol — memory-time claims L1–L4 (AlN, conserving phonon closure)

Campaign `20261005-fc3-memory`, review stage, 11 October 2026.
Search records: `review/novelty-searches/` (query plan `queries.txt`, runner `run_queries.sh`, log `batch-log.txt`,
raw JSON/MD per claim in `L1/`–`L4/`, keyword-filtered titles `titles_filtered.txt`, hand-run searches
`manual-search-log.md`, DOI verification `crossref-verification.txt`).

## 0. Scope and conventions

**Object.** tau_mem = N/K, K = b^T C^+ b, N = |C^+ b|^2 (b the heat-current source vector, C the linearised
energy-conserving three-phonon operator in entropy coordinates). In the relaxation-time approximation (RTA)
tau_mem = sum c v^2 tau^2 / sum c v^2 tau, i.e. the conductivity-weighted mean lifetime, the first moment of the
heat-current autocorrelation divided by its integral, and -kappa'(0)/kappa(0) of the frequency-dependent conductivity.

**Relation to the earlier audit.** The prior-art campaign `20261002-paper1-prior-art` (`08-literature.md`,
`13-final-report.md`) classified the identifiability statement in the free-rate positive-semidefinite class (C2(a):
partially known, explicit phonon statement potentially novel) and the event-cone inequality (C2(b): corollary of
Ben-Tal & Teboulle 1990). Those verdicts are not repeated. This protocol covers what was not searched there: mesh and
infrared behaviour of the memory time (L1, L2), its film regularisation (L3), and identifiability in the
force-constant-generated rate class (L4).

**Evidence labels.** SOURCE SAYS / SOURCE DERIVES / SOURCE OBSERVES (numerical) / SECONDARY (one paper's statement
about another, not checked in the original) / OUR DERIVATION (a short step derived here from a source's equation) /
OUR INTERPRETATION. Access is stated for every source: read (text inspected), text search (targeted search in full
text), abstract only, metadata only. Every DOI below resolves in Crossref (`crossref-verification.txt`); titles,
authors, venues, volumes and pages were checked against the citations.

## 1. Verdicts

| Claim | Verdict | Closest prior work | What is not in the searched literature |
|---|---|---|---|
| **L1** tau_mem carried by acoustic modes below ~3 THz, not converged on meshes where kappa is nearly converged; part above 3 THz converged | **PARTIALLY KNOWN.** Mechanism known; the specific demonstration is **potentially novel within the search scope** (modest: a convergence finding) | Beardo et al., Sci. Adv. 2021, Eq. 2 (this exact RTA formula as the heat-flux relaxation time, ab initio, one 24^3 mesh, Ge: 0.460 ns at 300 K, no convergence test); Sendra et al., PRB 2021 (same expression, via Beardo thesis Eq. 73); Tur-Prats et al., Nano Lett. 2025 ("precise sampling of low-frequency phonons is crucial" for l and tau, graphene); Esfarjani et al., PRB 2011, Eq. 9 (kappa(L) = kappa(inf) - E/L for omega^2 rates) | A mesh sequence for this memory time; its non-convergence while kappa converges; the frequency-cutoff split; the full conserving operator; any AlN memory-time study |
| **L2** finite bulk limit iff n < 3/2; AlN n ~ 0.8–1.3 below 1 THz, so plausibly finite | Criterion **KNOWN** (one-line corollary of published results); AlN exponents are new data (standard method); the conclusion "plausibly finite" is **CONTRADICTED for the physical crystal** by beyond-golden-rule physics and holds at most for the three-phonon Fermi-golden-rule operator | Vermeersch et al., PRB 91, 085202 (2015), Eqs. 19–20 (kappa_Sigma(tau) ~ 1 - (tau_min/tau)^(1/2) for n = 2; the conductivity-weighted mean lifetime then diverges iff n >= 3/2); Lepri, Livi & Politi, Phys. Rep. 2003, Eqs. 134 and 141 (3D: sound damping ~ q^2, <J(t)J(0)> ~ t^(-3/2)); Brey 1983 and Ernst, Hauge & van Leeuwen 1971 (long-time tails, divergent Burnett coefficients in 3D); Akhiezer-regime damping of long-wavelength phonons (Simons 1964; Woodruff & Ehrenreich 1961), attenuation ∝ omega^2 (Shvaika et al. 2021), used as an infrared cut-off or correction to kappa (Brito-Orta & Klemens 1981; Jiang, Lindsay & Koh 2016; Jin 2024) | Converged fixed-q low-frequency exponents for AlN; the explicit statement for the conductivity-weighted lifetime / Cattaneo time of a crystal |
| **L3** film regulator (\|v_z\|/L) makes the cross-plane memory time mesh-converged up to ~1 um (5.3 ps at 30 nm to 35.6 ps at 1 um) | **PARTIALLY KNOWN** (size-dependent memory/relaxation parameters are known); ab initio film-regularised memory time with a convergence study: **potentially novel within the search scope** (incremental); regulator is heuristic | Alvarez & Jou, APL 2007 / JAP 2008 (Knudsen-dependent effective coefficients); Beardo et al. 2019/2020 (ab initio hydrodynamic parameters in films and FDTR); Krewer et al., APL 2020 (electrons: thickness-dependent macroscopic response time set by mean and variance of microscopic relaxation times); Vermeersch, Carrete & Mingo, APL 2016 (ab initio cross-plane kappa(L), ln L in single crystals) | An ab initio heat-flux memory time versus film thickness, for AlN or any material, with Brillouin-zone convergence |
| **L4** static data (lifetimes within 10 %, kappa within 1 %) leave tau_mem undetermined by large factors in the symmetric-fc3 class; freedom in long-wavelength acoustic modes; witnesses mesh-specific | **PARTIALLY KNOWN** in spirit (sloppiness; kappa uncertainty from force constants; non-unique relaxation-time reconstruction); the specific statement is **potentially novel within the search scope** | Parks et al., PRMaterials 2020 (functional-ensemble uncertainty of kappa concentrated in TA phonons, MFP 100–300 nm); McGaughey et al. (Phonon Olympics), JAP 2025 (+-15 % spread driven by cubic-IFC lifetimes); Xie, Gu & Bao 2017 (fc3 errors -> kappa); Forghani & Hadjiconstantinou, IJHMT 2021 (F(Lambda) does not fix tau(omega)); Kurniawan et al., JCP 2022 (sloppy, unidentifiable potentials); Tur-Prats et al. 2025 and Abou Haibeh & Huberman 2026 (hydrodynamic quantities sensitive to scattering model / functional) | Constraining a symmetry-respecting fc3 class by static data and mapping the attainable range of a dynamic moment; locating the freedom in long-wavelength acoustic modes; testing transfer across meshes |

"Potentially novel within the search scope" always means: not found in the indexes and phases listed in section 2,
with the limits of section 8.

## 2. Search scope (common to all claims)

| Index / tool | Use | Outcome |
|---|---|---|
| Local multi-index tool (Crossref + Semantic Scholar + OpenAlex) | 41 queries labelled by claim and phase (P1 direct, P2 mathematical, P3 physical, P4 historical, P7 neighbouring) | Crossref 41/41, OpenAlex 41/41, Semantic Scholar 2/41 (HTTP 429). Noisy; no absence claim rests on it alone |
| arXiv (search API, HTML full text, targeted text search; PDFs via pdftotext) | ~35 searches; full or targeted reading of 20 preprints | Main source of full texts (section headings in the claim sections) |
| Semantic Scholar | title matches; complete forward-citation lists of Beardo et al. 2021 (78) and Sendra et al. 2021 (39) | snippet (full-text) search rate-limited |
| OpenAlex API | abstracts (~45 DOIs); keyword-filtered forward citations of Esfarjani et al. 2011 (765 citing), Vermeersch et al. 2015 I (98), Parks et al. 2020 (7), Phonon Olympics 2025 (10) | worked |
| Crossref API | verification of all 82 cited DOIs | all resolve |
| Publisher / repository full text | Sci. Adv. HTML (Beardo 2021); Nano Lett. author PDF (Tur-Prats 2025); Beardo PhD thesis (UAB 2021, 192 pp., §2.3 and §3.5 read); local copies of Beardo 2020 and Sendra 2022 | UAB record of Sendra 2021: HTTP 504; Sendra thesis: blocked (Cloudflare) |
| Firecrawl (research search, scrape), web search | pointers and publisher pages until the free-tier limit | see manual log |
| Zotero (user library) | checked first | same four items as on 2 October; two relevant (Hoque 2024, Camacho de la Rosa 2025) |

Terminology searched (phases 1–4): conductivity-weighted / kappa-weighted mean relaxation time; heat-flux relaxation
time; Cattaneo / Maxwell–Cattaneo–Vernotte relaxation time from first principles; memory time; first moment of the
heat-current autocorrelation; -kappa'(0)/kappa(0), low-frequency expansion of kappa(omega); second moment of the
lifetime or mean-free-path distribution; nonlocal (Guyer–Krumhansl) length from kinetic moments; temporal Knudsen
number; longest-lived relaxons, transport waves; long-time tails; Burnett coefficients; Akhiezer and Landau–Rumer
damping; Herring processes; Pomeranchuk; role of low-frequency phonons; extrapolation of kappa to infinite mesh.
Neighbouring fields (phase 7): kinetic theory of fluids, anomalous heat conduction in lattices, electron transport in
thin metal films, polymer rheology (weight-average relaxation time, relaxation-spectrum inversion), parameter
identifiability and sloppiness (systems biology, interatomic potentials).

---

## 3. L1 — the memory time is an infrared quantity that does not converge with the mesh

**Claim audited.** For AlN at 300 K, tau_mem (RTA and full conserving operator) is controlled by long-wavelength
acoustic phonons below ~3 THz and is not converged on Brillouin-zone meshes on which kappa is close to converged
(RTA, tetrahedron, 9x9x5 -> 15x15x9 -> 31x31x17: tau_mem,x 60 -> 70 -> 135 ps while kappa changes by 11 %); the
contribution of modes above 3 THz converges to ~2 %.

**Phases.** P1 "thermal conductivity weighted average relaxation time first principles Cattaneo", "phonon
relaxation time second moment convergence Brillouin zone mesh", "aluminum nitride ... phonon lifetimes",
`abs:"heat flux relaxation time"`; P2 "first moment heat current autocorrelation", "frequency-dependent thermal
conductivity first principles", time-domain transport function, smallest collision eigenvalue; P3 nonlocal length
from kinetic moments, MFP accumulation tail, longest-lived relaxons; P4 Herring 1954, Klemens 1981, extrapolation to
infinite mesh; P5 references of Beardo 2021 (expansion taken from Struchtrup 2005; ShengBTE; Vermeersch II), of
Tur-Prats 2025 and of the Beardo thesis; P6 all citing records of Beardo 2021 and Sendra 2021, keyword-filtered citing
records of Esfarjani 2011 and Vermeersch 2015; P7 rheology, electrons in films, Green–Kubo finite-size scaling.

### Closest prior work

1. **A. Beardo, M. López-Suárez, L. A. Pérez, L. Sendra, M. I. Alonso, C. Melis, J. Bafaluy, J. Camacho, L. Colombo,
   R. Rurali, F. X. Alvarez, J. S. Reparaz**, "Observation of second sound in a rapidly varying temperature field in
   Ge", *Sci. Adv.* **7**, eabg4677 (2021), doi:10.1126/sciadv.abg4677; arXiv:2007.05487v3. **Read** (published main
   text; arXiv v3 including supplementary sections SM4, SM6, SM8).
   SOURCE DERIVES (main text, Eq. 2; derivation announced for section S6, which in arXiv v3 reads "This section will be
   updated in a future version"): the relaxation time of the hyperbolic heat equation is
   tau_ss = sum_l hbar omega_l nu_l^2 tau_l^2 (d f_l/dT) / sum_l hbar omega_l nu_l^2 tau_l (d f_l/dT),
   i.e. exactly the RTA tau_mem. SOURCE SAYS: omega, tau, nu from DFT force constants (VASP, LDA; ShengBTE); RTA used
   "after verifying that the full iterative BTE picture does not alter the prediction"; "The agreement of the predicted
   values with those obtained from the experiments is remarkable for T > 100 K". SOURCE OBSERVES (arXiv v3, SM4 and
   Table SM4): one 24x24x24 q-grid; tau_ss = 0.460 ns at 300 K, 10.4 ns at 100 K. A search of the full arXiv text for
   "converg", "grid", "mesh" finds only the 24^3 statement (for kappa); no convergence test of tau_ss. The published
   supplementary file was not inspected.
2. **L. Sendra, A. Beardo, P. Torres, J. Bafaluy, F. X. Alvarez, J. Camacho**, "Derivation of a hydrodynamic heat
   equation from the phonon Boltzmann equation for general semiconductors", *Phys. Rev. B* **103**, L140301 (2021),
   doi:10.1103/PhysRevB.103.L140301. **Journal text not obtained** (repository unreachable). Its RTA expressions are
   reproduced, with attribution, in **A. Beardo**, *Generalized Hydrodynamic Heat Transport in Semiconductors*, PhD
   thesis, Universitat Autònoma de Barcelona (2021), https://ddd.uab.cat/pub/tesis/2021/hdl_10803_673590/abr1de1.pdf —
   **read** §2.3 (Eqs. 68–73) and §3.5. SOURCE SAYS (thesis, SECONDARY for the journal paper): with
   <x> = integral hbar omega d_T f x dk / c, kappa = (1/3) c <v^2 tau> (Eq. 71),
   l^2 = (1/5) <tau^2 v^3 / v_p> / <v / v_p> (Eq. 72), and the heat-flux relaxation time
   tau = <v^2 tau^2> / <v^2 tau> (Eq. 73). §3.5: "the heat flux relaxation time is one order of magnitude larger than the
   usual relaxation time <tau_k> obtained using the specific heat as the averaging weight"; in the Si and Ge experiments
   modelled in Chapters 4–7 of that thesis "memory effects are negligible and, hence, the heat flux relaxation time
   values can not be simultaneously validated"; Chapter 9 compares fitted and ab initio heat-flux relaxation times for
   Ge at several temperatures (the setting of item 1). No convergence or infrared discussion was found in these sections.
3. **J. Tur-Prats, Z. Han, A. Beardo, X. Ruan, F. X. Alvarez**, "High-Order Anharmonicities Shape Phonon Hydrodynamic
   Effects in Graphene", *Nano Lett.* **25**, 11203–11209 (2025), doi:10.1021/acs.nanolett.5c00855. **Read** (author
   PDF; supporting information not read). SOURCE SAYS (p. 11207): "In particular, precise sampling of low-frequency
   phonons is crucial, since they display a dominant contribution on the converged value of l and tau whereas
   displaying a relatively smaller contribution on kappa." Abstract: the nonlocal length and heat-flux relaxation time
   "are significantly reduced if considering full anharmonicity" (three- plus four-phonon). The values are reported as
   converged; 2D material.
4. **K. Esfarjani, G. Chen, H. T. Stokes**, "Heat transport in silicon from first-principles calculations", *Phys.
   Rev. B* **84**, 085204 (2011), doi:10.1103/PhysRevB.84.085204; arXiv:1107.5288. **Read** (text search: §III.2,
   lattice-dynamics section, conclusions). SOURCE DERIVES (Eq. 9, Klemens rates 1/tau ∝ omega^2, classical capacity,
   constant velocity, DOS ∝ omega^2): kappa(L) = kappa(inf) - D omega_cut = kappa(inf) - E/L = kappa(inf) - F/sqrt(Lambda_c).
   SOURCE SAYS: "To get an accurate estimate of the thermal conductivity, one needs to extrapolate the data obtained
   from a finite number of k-mesh points, according to Eq. (9)"; "Such extrapolation is justified for relaxation rates
   which are quadratic in frequency at low frequencies"; for large supercells "relaxation times scale as
   1/omega_cut^2 ∝ L^2". This is the first-moment (kappa) version of the power counting behind L1.
5. **P. Torres, A. Torelló, J. Bafaluy, J. Camacho, X. Cartoixà, F. X. Alvarez**, "First principles kinetic-collective
   thermal conductivity of semiconductors", *Phys. Rev. B* **95**, 165407 (2017), doi:10.1103/PhysRevB.95.165407;
   arXiv:1606.01149. **Read** (text search). SOURCE SAYS: low-frequency rates replaced by analytical forms (Han for
   umklapp, Herring 1/tau_N ∝ omega^2 T^3 for normal processes); "Using analitical expressions for low frequencies,
   accurate results can be achieved even with a coarse k-point grid" (for kappa).
6. **A. Cepellotti, N. Marzari**, "Thermal Transport in Crystals as a Kinetic Theory of Relaxons", *Phys. Rev. X*
   **6**, 041013 (2016), doi:10.1103/PhysRevX.6.041013; arXiv:1603.02608. **Read** (text search, graphene section).
   SOURCE SAYS: "an improvement in the integration mesh has new phonon modes appear in the long-lifetime region;
   instead the longest relaxon relaxation times converge - from above - to the discretized values". Contrasting
   precedent: a strongly hydrodynamic 2D material whose slowest relaxons are drift-like and mesh-stable. Companion
   preprint A. Cepellotti, N. Marzari, "Transport waves as crystal excitations", arXiv:1612.04317: "computational costs
   prevent a systematic convergence study of these isolated modes".
7. Single-mesh first-principles dynamic quantities without an infrared analysis (all **read** by text search unless
   stated): C. Hua, L. Lindsay, *Phys. Rev. B* **102**, 104310 (2020), doi:10.1103/PhysRevB.102.104310 (kappa(omega)
   from the full collision matrix on reduced integration grids, "sufficiently converged for our purposes"; "At low
   frequency ... the real part remains constant"); D. E. Crawford, Y. Zeng, J. Vidal, J. Dong, *Phys. Rev. B* **112**,
   115201 (2025), doi:10.1103/p8wg-p1j3, arXiv:2501.00679 (time-domain transport function Z(t) for Si on a 21^3 grid;
   remarks only that Z(t) becomes non-integrable without umklapp); Y. Zeng, J. Dong, *Phys. Rev. B* **113**, 235204
   (2026), doi:10.1103/hx28-4ksc (gradient hierarchy "asymptotic ... rather than a uniformly convergent series");
   N. Malviya, N. K. Ravichandran, arXiv:2502.00337 (2025) (number of collision eigenmodes "required to achieve
   convergence on the steady-state kappa" chosen from the kappa-accumulation curve; meshes chosen by kappa convergence)
   and arXiv:2606.17829 (2026) (temporal Knudsen number eta/min(sigma ≠ 0) from the smallest nonzero collision
   eigenvalue; mesh data in a supplement not inspected); A. Beardo et al., *Phys. Rev. B* **101**, 075303 (2020),
   doi:10.1103/PhysRevB.101.075303 (local copy read: KCM ab initio tau = 1002 / 148 / 42 / 27 ps and l = 3127 / 1047 /
   185 / 142 nm for Si at 81 / 154 / 311 / 416 K, Table I; KCM averaging, not Eq. 73).
8. **P. Jiang, L. Lindsay, Y. K. Koh**, "Role of low-energy phonons with mean-free-paths >0.8 μm in heat conduction in
   silicon", *J. Appl. Phys.* **119**, 245705 (2016), doi:10.1063/1.4954674. **Abstract only.** SOURCE SAYS: the
   contribution of low-energy phonons is "still inconclusive, as exemplified by the discrepancies as large as 30%
   between different first-principles calculations"; macroscopic (Akhiezer) damping "eliminates the contribution of
   phonons with mean-free-paths >20 μm at 300 K", 15 W m^-1 K^-1 (10 %) of the calculated kappa.
9. Historical: **P. G. Klemens**, *Int. J. Thermophys.* **2**, 55–62 (1981), doi:10.1007/BF00503574 (**abstract
   only**): with l(omega) ∝ omega^-2 T^-1 "the spectral thermal conductivity is independent of frequency, and low
   frequencies play a larger role than they do in the heat content". OUR DERIVATION: a flat spectral conductivity
   with tau ∝ omega^-2 makes integral (spectral kappa) tau d omega ∝ integral omega^-2 d omega divergent, i.e. Klemens'
   law already implies an infinite conductivity-weighted mean lifetime. **C. Herring**, *Phys. Rev.* **95**, 954 (1954),
   doi:10.1103/PhysRev.95.954 (**metadata only**; low-energy phonon laws, see L2).
10. Context: **J. A. Johnson et al.**, *Phys. Rev. Lett.* **110**, 025901 (2013), doi:10.1103/PhysRevLett.110.025901;
   arXiv:1204.4735 (**text search**): "the commonly cited textbook values of an 'average' phonon MFP are of little
   relevance ... Perhaps a more useful parameter would be the 'median thermal conductivity MFP'". The "average" there
   is the gray kinetic MFP, not the conductivity-weighted mean; the remark is consistent with, but does not state, L1.

### What differs
* The RTA memory-time formula is published (items 1–2) and evaluated ab initio on a single mesh. No source found
  reports it on a mesh sequence, reports non-convergence while kappa converges, or splits it at a frequency cutoff.
* No source evaluates the conserving-operator memory time N/K on a mesh sequence (the event-operator sequence to
  12x12x8 of E24 is new within the search scope).
* No memory-time, Cattaneo-time or nonlocal-length study of AlN was found.
* Contrasting result: for graphene the slowest relaxons converge from above (item 6); the AlN behaviour is the
  opposite, which is expected for a 3D crystal without a drift-dominated slow sector (OUR INTERPRETATION).

### Interpretation
* OUR DERIVATION (from item 4, same assumptions): one moment higher, the infrared deficit of N for omega^2 rates is
  ∝ integral_{omega_cut} omega^-2 d omega ∝ 1/omega_cut, so tau_mem grows linearly in 1/f_min while kappa converges as
  f_min; this matches the campaign's fit tau_x = 4.6 + 110/f_min ps (review). L1 is therefore an expected, now
  quantified, consequence of known power counting rather than a new mechanism.
* OUR INTERPRETATION: the same applies to published single-mesh values — the Ge tau_ss = 0.46 ns of item 1 and any
  l^2 from Eq. 72 of item 2 (for acoustic branches v_p ≈ v, so the numerator of Eq. 72 is the same infrared moment as
  N). A Ge mesh sequence for Eq. 2 would test this directly and cheaply.

**Verdict L1: PARTIALLY KNOWN.** Formula, low-frequency dominance and power counting: LITERATURE SUPPORTED (items 1–5,
9). The explicit non-convergence of this memory time under mesh refinement at near-converged kappa, its cutoff split,
and the conserving-operator version: POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE (a convergence finding of modest
novelty but direct consequence for published ab initio Cattaneo times).

---

## 4. L2 — finite bulk limit iff n < 3/2; AlN exponents; plausibly finite

**Claim audited.** Whether tau_mem has a finite bulk limit is decided by the low-frequency linewidth exponent n
(gamma ~ omega^n): divergent iff n >= 3/2 (constant group velocity, classical heat capacity). Converged first-principles
AlN linewidths at fixed q give n ~ 0.8–1.3 below ~1 THz, so the limit is plausibly finite but approached very slowly.

**Phases.** P1 "infrared divergence phonon relaxation time moments", "low-frequency phonon scattering rate power law
exponent first principles", "Cattaneo relaxation time divergence low-frequency phonons"; P2 "long-time tail heat
current autocorrelation crystal", "nonanalytic frequency dependence thermal conductivity square root", "Burnett
coefficients divergence three dimensions"; P3 Lévy superdiffusion and scattering exponent, thermal viscosity,
nonlocal length from kinetic moments, MFP accumulation tails; P4 Herring, Pomeranchuk, Akhiezer, Klemens; P5
references of Vermeersch I (on n) and of Lepri–Livi–Politi (Pomeau–Résibois, Ernst et al.); P6 keyword-filtered
citing records of Vermeersch I and Esfarjani 2011; P7 kinetic theory of fluids, anomalous lattice heat conduction,
2D conductivity divergence, glasses (Akhiezer-type damping), rheology.

### Closest prior work

1. **B. Vermeersch, J. Carrete, N. Mingo, A. Shakouri**, "Superdiffusive heat conduction in semiconductor alloys.
   I. Theoretical foundations", *Phys. Rev. B* **91**, 085202 (2015), doi:10.1103/PhysRevB.91.085202;
   arXiv:1406.7341v2. **Read** (section "Physical origin of Lévy dynamics", Eqs. 16–21; page rendered). SOURCE ASSUMES:
   single branch, constant group velocity v0, constant mode capacity, isotropic crystal, tau_k = tau_min (ka/pi)^-n
   (Eq. 16). SOURCE DERIVES: kappa_Sigma(tau) ~ integral_{(tau_min/tau)^(1/n)}^1 dk / k^(n-2) (Eq. 19); for n = 2,
   kappa_Sigma ~ 1 - sqrt(tau_min/tau); n = 3, ln(tau/tau_min); n > 3, (tau/tau_min)^(beta_n - 1) - 1 (Eq. 20); abstract:
   alpha = 3 - beta = 1 + 3/n = 2 - gamma for n > 3. The paper treats n > 3 (superdiffusion) and does not discuss the
   first moment of the accumulation. Note: the arXiv text states that the single-branch conductivity "diverges for
   n > 2", whereas its Eq. 19 converges for n < 3 (logarithmic at n = 3); the published text was not checked.
   OUR DERIVATION from Eq. 19: d kappa_Sigma/d tau = (1/n)(tau_min/tau)^((3-n)/n)/tau, so the conductivity-weighted mean
   lifetime integral tau d kappa_Sigma / kappa ∝ integral (tau_min/tau)^((3-n)/n) d tau is finite iff (3-n)/n > 1, i.e.
   **n < 3/2**, logarithmic at n = 3/2, ∝ tau_max^(1/2) ∝ 1/f_min for n = 2 — L2's criterion under L2's assumptions.
2. **K. Esfarjani, G. Chen, H. T. Stokes**, PRB **84**, 085204 (2011) — item L1-4: the same power counting for kappa.
3. **S. Lepri, R. Livi, A. Politi**, "Thermal conduction in classical low-dimensional lattices", *Phys. Rep.* **377**,
   1–80 (2003), doi:10.1016/S0370-1573(02)00558-6; arXiv:cond-mat/0112193. **Read** (§ on mode-coupling theory, Eqs.
   133–142). SOURCE DERIVES (mode-coupling theory, momentum-conserving lattices): renormalised dispersion
   Omega(q) ≃ c q - i c' q^2 in d = 3 (Eq. 134); heat-current autocorrelation as a mode sum
   sum_k v_k^2 <E_k^2> e^(-Gamma_k t), "dominated by the low-q contribution at large times" (Eq. 139); hence
   <J(t)J(0)> ~ t^(-3/2) in d = 3 (Eq. 141) and finite kappa in 3D (Eq. 142). OUR DERIVATION: the first moment
   integral t <J(t)J(0)> dt then diverges (∝ T^(1/2) for a time cutoff T), i.e. kappa(omega) - kappa(0) ∝ omega^(1/2):
   the Cattaneo time does not exist for the full 3D dynamics. This is the n = 2 case of L2.
4. Kinetic theory of fluids (phase 7): **M. H. Ernst, E. H. Hauge, J. M. J. van Leeuwen**, *Phys. Rev. A* **4**, 2055
   (1971), doi:10.1103/PhysRevA.4.2055 (**abstract**: t^(-d/2) tails of the velocity and of the kinetic parts of the
   shear-viscosity and heat-conductivity correlation functions); **J. J. Brey**, *J. Chem. Phys.* **79**, 4585 (1983),
   doi:10.1063/1.446374 (**abstract**: Burnett integrands have tails t^(1-d/2), so "the expressions for the Burnett
   transport coefficients are all divergent"); **N. G. Inozemtseva, B. I. Sadovnikov**, *Theor. Math. Phys.* **48**, 658
   (1981), doi:10.1007/BF01037993, and **51**, 409 (1982), doi:10.1007/BF01029268 (**metadata**: divergence of Burnett
   coefficients in 3D hard-sphere hydrodynamics); **M. H. Ernst, J. R. Dorfman**, *J. Stat. Phys.* **12**, 311 (1975),
   doi:10.1007/BF01012067 (**metadata**: nonanalytic dispersion relations); **Y. Pomeau, P. Résibois**, *Phys. Rep.*
   **19**, 63 (1975), doi:10.1016/0370-1573(75)90019-8 (**metadata**). **K. Saito, A. Dhar**, *Phys. Rev. Lett.* **104**,
   040601 (2010), doi:10.1103/PhysRevLett.104.040601 (**abstract**: finite kappa in a 3D anharmonic lattice; the tail
   exponent is not stated in the abstract).
5. Beyond-golden-rule damping of long-wavelength phonons (Akhiezer regime, omega tau_th << 1):
   **S. Simons**, *Proc. Phys. Soc.* **83**, 749 (1964), doi:10.1088/0370-1328/83/5/306 (**abstract**: phonon
   interaction theory extended to "one phonon having a wavelength larger than the mean free path of those with which it
   interacts", sound absorption at arbitrary temperature); **T. O. Woodruff, H. Ehrenreich**, *Phys. Rev.* **123**, 1553
   (1961), doi:10.1103/PhysRev.123.1553 (**metadata**); **A. Shvaika, M. Shpot, W. Schirmacher, T. Bryk, G. Ruocco**,
   arXiv:2104.13076 (2021) (**abstract**: "a sound attenuation proportional to omega^2 (Akhiezer-type damping)");
   **R. Brito-Orta, P. G. Klemens**, *J. Phys. Colloq.* **42**, C6-256 (1981), doi:10.1051/jphyscol:1981674
   (**abstract**: low-frequency longitudinal contributions to kappa "estimated using the Pomeranchuk mean free path and
   a cut-off given by the Akhieser theory"); **P. Jiang, L. Lindsay, Y. K. Koh** (item L1-8); **J. S. Jin**, *J. Appl.
   Phys.* **135**, 045103 (2024), doi:10.1063/5.0187032 (**abstract**: in that model Akhiezer damping reduces the Si conductivity
   by 12.6 % at 300 K and 16.8 % at 500 K); **V. Chiloyan, S. Huberman, A. A. Maznev, K. A. Nelson, G. Chen**,
   *Appl. Phys. Lett.* **116**, 163102 (2020), doi:10.1063/1.5139069 (**abstract**: nonthermal populations of long-MFP
   low-frequency phonons).
6. Low-frequency rate laws (historical): **C. Herring**, *Phys. Rev.* **95**, 954 (1954), doi:10.1103/PhysRev.95.954
   (**metadata**); **I. Pomeranchuk**, *Phys. Rev.* **60**, 820 (1941), doi:10.1103/PhysRev.60.820 (**metadata**); **S. Tamura,
   A. Sangu, H. J. Maris**, *Phys. Rev. B* **68**, 014302 (2003), doi:10.1103/PhysRevB.68.014302 (**abstract**: Herring
   rate Gamma = A omega^2 T^3 for L phonons well below the Debye temperature, with A = 0 along [001] of TeO2 —
   direction-dependent coefficients); **P. G. Klemens** (item L1-9); **A. Ward, D. A. Broido**, *Phys. Rev. B* **81**,
   085205 (2010), doi:10.1103/PhysRevB.81.085205 (**abstract**: first-principles forms of normal and umklapp relaxation
   times for Si and Ge).
7. Infrared divergences of hydrodynamic coefficients in other settings: **L. Sendra et al.**, *Phys. Rev. B* **106**,
   155301 (2022), doi:10.1103/PhysRevB.106.155301 (**read**, local copy: in 2D the integrals of the collective-limit
   coefficients diverge for flexural modes, "∝ integral 1/k dk which diverges for infinite samples", and finite
   discretisation "may not reach the region" concerned); **L. F. C. Pereira, D. Donadio**, *Phys. Rev. B* **87**, 125424
   (2013), doi:10.1103/PhysRevB.87.125424 (**metadata**: divergent kappa of strained graphene). Thermal viscosity and
   nonlocal length: Simoncelli, Marzari & Cepellotti, *Phys. Rev. X* **10**, 011019 (2020) (**text search**: computed
   on 27^3 grids, no infrared analysis found); Tur-Prats et al. 2025 (item L1-3); Beardo thesis Eq. 72 (item L1-2).
8. Frequency- and wavevector-dependent conductivity (**abstracts**): Y. K. Koh, D. G. Cahill, *Phys. Rev. B* **76**,
   075207 (2007), doi:10.1103/PhysRevB.76.075207; F. Yang, C. Dames, *Phys. Rev. B* **91**, 165311 (2015),
   doi:10.1103/PhysRevB.91.165311; S. G. Volz, *Phys. Rev. Lett.* **87**, 074301 (2001), doi:10.1103/PhysRevLett.87.074301;
   K. T. Regner et al., *Nat. Commun.* **4**, 1640 (2013), doi:10.1038/ncomms2630 (metadata); G. D. Mahan, F. Claro,
   *Phys. Rev. B* **38**, 1963 (1988), doi:10.1103/PhysRevB.38.1963; P. B. Allen, V. Perebeinos, *Phys. Rev. B* **98**,
   085427 (2018), doi:10.1103/PhysRevB.98.085427; P. B. Allen, *Phys. Rev. B* **97**, 134307 (2018),
   doi:10.1103/PhysRevB.97.134307; A. Kefayati, P. B. Allen, V. Perebeinos, *Phys. Rev. B* **105**, 235402 (2022),
   doi:10.1103/PhysRevB.105.235402; C. Hua, L. Lindsay, X. Chen, A. J. Minnich, *Phys. Rev. B* **100**, 085203 (2019),
   doi:10.1103/PhysRevB.100.085203. These treat kappa(omega) and kappa(q) through cut-off or suppression functions or
   full RTA/BTE solutions; none of the abstracts (nor the texts read: Hua & Lindsay 2020, Li, Jiang & Chen
   arXiv:2408.05654) addresses the existence of the low-frequency slope -kappa'(0) of a 3D crystal.

### What differs
* The criterion (divergent iff n >= 3/2 under constant velocity and classical capacity) is not stated in this form in
  any source found, but follows in one line from Vermeersch et al., Eq. 19 (same assumptions), from the
  Esfarjani-type power counting, and — for n = 2 — from the 3D long-time tail of item 3. It should be presented as a
  corollary with these attributions, not as a result. The project note `theory/aln/closure/reports/C-infrared.md`
  (16 September 2026) already derived the moment criteria (finite iff m·alpha < n in n dimensions); it is internal
  and does not change the attribution.
* Converged fixed-q low-frequency exponents for AlN were not found in the literature; they are new numerical data
  obtained with standard methods (single workflow, phono3py tetrahedron integration).
* **Contrary physics.** The fixed-q exponents of E25 are Fermi-golden-rule three-phonon linewidths, in which thermal
  phonons are treated as infinitely long-lived. For omega tau_th << 1 the damping of an acoustic wave enters the
  Akhiezer regime with attenuation ∝ omega^2 (item 5), and mode-coupling theory gives sound damping ∝ q^2 and a
  t^(-3/2) heat-current tail in 3D (item 3). Both give n = 2 >= 3/2 asymptotically. The bulk heat-flux memory time of
  the physical crystal is therefore expected to be infinite; a finite limit can hold at most for the three-phonon
  golden-rule operator. The crossover lies near 1/(2 pi tau_th) (OUR ESTIMATE: GHz to tens of GHz for AlN at 300 K;
  not computed), far below the 0.5–3 THz range of the campaign's meshes, which is why it cannot be seen there.
  Jiang, Lindsay & Koh (2016) already identified the omission of this damping as a possible source of overestimated
  long-MFP contributions in first-principles kappa; tau_mem weights those modes far more strongly than kappa
  (OUR INTERPRETATION).

**Verdict L2: criterion KNOWN (corollary of items 1–4); AlN exponents: new data, standard method; the inference
"finite bulk limit plausible" is CONTRADICTED for the physical crystal by items 3 and 5 and should be restricted to
the three-phonon golden-rule operator, where it remains UNRESOLVED.** The campaign statement S7 ("UNRESOLVED;
divergence PLAUSIBLE") is closer to the literature than the E25 status line ("divergence NOT SUPPORTED") once the
model level is stated.

---

## 5. L3 — film-thickness regulator and the cross-plane memory time

**Claim audited.** A film-thickness regulator (boundary scattering, Matthiessen 1/tau_L = 1/tau + |v_z|/L) makes the
cross-plane memory time mesh-converged for L up to ~1 um (AlN, 300 K: 5.3 ps at 30 nm to 35.6 ps at 1 um).

**Phases.** P1 "size-dependent thermal relaxation time thin film Cattaneo", "first principles hyperbolic heat
conduction thin film relaxation time", "second sound thin film first principles"; P2 Knudsen-number-dependent
effective coefficients, moment-ratio response times; P3 "effective relaxation time boundary scattering thin film
hydrodynamic", "cross-plane heat conduction thin films ab initio", FDTR/TDTR with hyperbolic or hydrodynamic models;
P4 Fuchs–Sondheimer theory (project note), Alvarez–Jou continued-fraction models; P5 references of Beardo 2020 and
Camacho de la Rosa 2025; P6 citing records of Beardo 2021 and Sendra 2021; P7 electron transport in metal films.

### Closest prior work

1. **F. X. Alvarez, D. Jou**, "Memory and nonlocal effects in heat transport: From diffusive to ballistic regimes",
   *Appl. Phys. Lett.* **90**, 083109 (2007), doi:10.1063/1.2645110, and "Size and frequency dependence of effective
   thermal conductivity in nanosystems", *J. Appl. Phys.* **103**, 094321 (2008), doi:10.1063/1.2913057 (**abstracts**).
   SOURCE SAYS: a generalised model with memory and nonlocal effects; effective conductivity depending on the Knudsen
   number and on frequency, for cross-plane layers and in-plane geometries. Phenomenological, not first-principles.
2. **A. Beardo, M. Calvo-Schwarzwälder, J. Camacho, T. G. Myers, P. Torres, L. Sendra, F. X. Alvarez, J. Bafaluy**,
   "Hydrodynamic Heat Transport in Compact and Holey Silicon Thin Films", *Phys. Rev. Appl.* **11**, 034003 (2019),
   doi:10.1103/PhysRevApplied.11.034003 (**abstract**: hydrodynamic equation with parameters calculated ab initio for
   Si films); **A. Beardo et al.**, PRB **101**, 075303 (2020) (item L1-7: FDTR, bulk KCM parameters with slip and
   temperature-jump boundary conditions); **P. Torres et al.** (item L1-5: Matthiessen boundary rate in the kinetic
   term, form factor with L_eff = 2.25 h for films).
3. **K. L. Krewer, W. Zhang, J. Arabski, G. Schmerber, E. Beaurepaire, M. Bonn, D. Turchinovich**,
   "Thickness-dependent electron momentum relaxation times in iron films", *Appl. Phys. Lett.* **116**, 102406 (2020),
   doi:10.1063/1.5142479; arXiv:2008.07852 (**abstract**). Neighbouring field. SOURCE SAYS: the macroscopic current
   response time (29 fs to 7 fs from thick to thin films) "depends on average and variation" of the distribution of
   microscopic relaxation times, and that variation "depends on film thickness because electrons with different
   relaxation times are affected differently by the confinement". This is the electronic analogue of a
   thickness-regularised tau_mem.
4. Ab initio cross-plane conduction in films: **B. Vermeersch, J. Carrete, N. Mingo**, *Appl. Phys. Lett.* **108**,
   193104 (2016), doi:10.1063/1.4948968 (**abstract**: suppression model matching Monte Carlo with ab initio phonons
   within ≤ 5 %; kappa_perp ~ ln L in single crystals); **C. Hua, A. J. Minnich**, *J. Appl. Phys.* **117**, 175306 (2015),
   doi:10.1063/1.4919432 (**abstract**: series solution of the frequency-dependent BTE for cross-plane transport);
   **J. Maassen, M. Lundstrom**, *J. Appl. Phys.* **117**, 135102 (2015), doi:10.1063/1.4916245, and **119**, 095102
   (2016), doi:10.1063/1.4942836 (**abstracts**: a simplified BTE recast exactly as Cattaneo/hyperbolic equations; "The
   correct physical boundary conditions ... are paramount").
5. Thermoreflectance with non-Fourier models and AlN films: **Z. Xiang, P. Jiang, R. Yang**, *J. Appl. Phys.* **132**,
   205104 (2022), doi:10.1063/5.0113257 (**abstract**: TDTR analysis with a phonon-hydrodynamic model, Si and Ge);
   **A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, *J. Appl. Phys.* **137**, 155103 (2025),
   doi:10.1063/5.0257299 (**abstract**: FDTR "high sensitivity to relaxation-time parameters"; a method to characterise
   the relaxation time "remains an open problem"); **M. S. B. Hoque et al.**, *Appl. Phys. Lett.* **125**, 262201
   (2024), doi:10.1063/5.0239769 (**abstract**: ballistic-to-diffusive transition in AlN films 1.6–2440 nm on sapphire;
   no memory time).

### What differs
* No source found computes an ab initio heat-flux memory (Cattaneo) time as a function of film thickness, for AlN or
  any material, nor shows that the film regulator restores Brillouin-zone convergence of that quantity.
* Size dependence of effective relaxation/memory parameters is known phenomenologically (item 1) and is built into
  ab initio hydrodynamic modelling of films through Matthiessen terms or boundary conditions (item 2); the moment-ratio
  mechanism is documented for electrons (item 3).

### Caveats raised by the literature
* The Matthiessen |v_z|/L regulator is heuristic. Ab initio cross-plane work uses Monte Carlo-calibrated suppression
  functions or full BTE solutions (item 4); Cattaneo-type film modelling depends on boundary conditions (item 4,
  Maassen & Lundstrom). The project note `theory/aln/closure/reports/C-infrared.md` (§4) records the double-counting
  risk of adding |v_z|/L before a boundary treatment. A film memory time is properly the first moment of the film's
  response function with stated boundary conditions.
* OUR INTERPRETATION: for L ≤ 1 um the regulator also removes the Akhiezer-regime problem of L2, since modes with
  omega tau_th << 1 have mean free paths far longer than the film thickness. This supports the choice of a
  thickness-regularised object.

**Verdict L3: PARTIALLY KNOWN (concept); the ab initio, mesh-converged film-regularised memory time for AlN is
POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE, incremental, and contingent on replacing or validating the Matthiessen
regulator against a cross-plane film solution.**

---

## 6. L4 — identifiability of the memory time in the force-constant-generated rate class

**Claim audited.** With three-phonon rates generated by symmetry-constrained third-order force constants perturbed
within the spread between two first-principles workflows, static data (all mode lifetimes within 10 %, kappa within
1 %) leave tau_mem undetermined by large factors on each discretised operator; the ambiguity resides in long-wavelength
acoustic modes; witnesses optimised on one mesh do not transfer to others.

**Phases.** P1 "uncertainty quantification lattice thermal conductivity force constants", "Bayesian uncertainty
interatomic force constants phonon", "sensitivity thermal conductivity third-order force constants"; P2 sloppy
models, practical identifiability, profile likelihood; P3 relaxation-time and MFP-spectrum reconstruction and its
uniqueness, exchange-correlation dependence, uncertainty of hydrodynamic/second-sound parameters, code benchmarks;
P4 Bayesian ensembles of interatomic potentials (2004), relaxation-spectrum inversion in rheology; P5 references of
Parks 2020 and Kurniawan 2022; P6 citing records of Parks 2020 (7) and of the Phonon Olympics (10); P7 systems
biology, interatomic-potential uncertainty, rheology.

### Closest prior work

1. **H. L. Parks, H.-Y. Kim, V. Viswanathan, A. J. H. McGaughey**, "Uncertainty quantification in first-principles
   predictions of phonon properties and lattice thermal conductivity", *Phys. Rev. Mater.* **4**, 083805 (2020),
   doi:10.1103/PhysRevMaterials.4.083805; arXiv:2005.08835. **Read** (abstract, introduction and results by text
   search). SOURCE SAYS: a 2000-member BEEF-vdW functional ensemble yields an ensemble of force constants and of kappa
   (Si: skewed distribution, mean 190 W m^-1 K^-1, standard deviation 36); kappa correlates with the TA sound speed
   (R^2 = 0.89) and mean Grüneisen parameter (0.85); "differences in ensemble predictions of thermal conductivity are
   correlated with the behavior of transverse-acoustic phonons with mean free paths between 100 and 300 nm". Static
   quantities only; harmonic and cubic constants vary together.
2. **A. J. H. McGaughey, L. Lindsay, H. Bao, T. Hamakawa, R. Juneja, S. Li, W. Li, R. Masuki, F. Meng, H. Meng,
   T. Pandey, C. Shao, J. Shiomi, T. Tadano, A. Togo, A. Wang, X. Zhang**, "Phonon Olympics: Phonon property and
   lattice thermal conductivity benchmarking from open-source packages", *J. Appl. Phys.* **138**, 135108 (2025),
   doi:10.1063/5.0289819 (**abstract**; data used by the campaign). SOURCE SAYS: teams' conductivities fall within
   ±15 % of their mean for each material; harmonic properties do not explain the variation; "It is the lifetimes
   associated with three-phonon scattering, obtained from the cubic force constants, that drive the variations".
3. Sensitivity of kappa to force constants and functionals: **Y. Srivastava, A. G. Gokhale, A. Jain**, *Phys. Rev. B*
   **112**, 195201 (2025), doi:10.1103/2n3p-hz62; arXiv:2509.14702 (**abstract**: solvers agree to ~1 % with identical
   constants; ~10 % between DFT packages; > 20 % between functionals); **A. Jain, A. J. H. McGaughey**, *Comput. Mater.
   Sci.* **110**, 115 (2015), doi:10.1016/j.commatsci.2015.08.014 (**metadata**); **H. Xie, X. Gu, H. Bao**, *Comput.
   Mater. Sci.* **138**, 368 (2017), doi:10.1016/j.commatsci.2017.06.047; arXiv:1706.08232 (**abstract**: errors
   introduced in third-order constants at typical first-principles force accuracy make kappa "quite different";
   imposing translational invariance does not always help; Grüneisen parameters are not a sufficient accuracy
   criterion); **H. Zhou, S. Zhou, Z. Hua, K. Bawane, T. Feng**, *Appl. Phys. Lett.* **123**, 192201 (2023),
   doi:10.1063/5.0173762 (**abstract**: fourth-order constants sensitive to energy-surface roughness of functionals,
   third-order smooth); **A. Krishnamoorthy et al.**, *Comput. Phys. Commun.* **254**, 107337 (2020),
   doi:10.1016/j.cpc.2020.107337 (**metadata**: Pareto-front uncertainty of force fields for kappa).
4. Sloppiness and identifiability: **R. N. Gutenkunst et al.**, *PLoS Comput. Biol.* **3**, e189 (2007),
   doi:10.1371/journal.pcbi.0030189 (**abstract**); **A. Raue et al.**, *Bioinformatics* **25**, 1923 (2009),
   doi:10.1093/bioinformatics/btp358 (**abstract**); **M. K. Transtrum et al.**, *J. Chem. Phys.* **143**, 010901 (2015),
   doi:10.1063/1.4923066 (**metadata**); **B. B. Machta, R. Chachra, M. K. Transtrum, J. P. Sethna**, *Science* **342**,
   604 (2013), doi:10.1126/science.1238723 (**metadata**); **Y. Kurniawan et al.**, *J. Chem. Phys.* **156**, 214103
   (2022), doi:10.1063/5.0084988 (**abstract**: interatomic potentials "are typically sloppy ... parameters are
   unidentifiable", profile-likelihood and Bayesian analysis); **S. L. Frederiksen, K. W. Jacobsen, K. S. Brown,
   J. P. Sethna**, *Phys. Rev. Lett.* **93**, 165501 (2004), doi:10.1103/PhysRevLett.93.165501 (**abstract**: Bayesian
   ensembles of fitted potentials to estimate errors of predicted properties).
5. Inverse problems for relaxation times: **A. J. Minnich**, *Phys. Rev. Lett.* **109**, 205901 (2012),
   doi:10.1103/PhysRevLett.109.205901 (**abstract**); **M. Forghani, N. G. Hadjiconstantinou, J.-P. M. Péraud**,
   *Phys. Rev. B* **94**, 155439 (2016), doi:10.1103/PhysRevB.94.155439 (**abstract**); **M. Forghani,
   N. G. Hadjiconstantinou**, *Phys. Rev. B* **97**, 195440 (2018), doi:10.1103/PhysRevB.97.195440, and *Appl. Phys.
   Lett.* **114**, 023106 (2019), doi:10.1063/1.5062846 (**metadata**), and "Comparison on two ways of phonon spectral
   property reconstruction from thermal spectroscopy experiment", *Int. J. Heat Mass Transf.* **175**, 121119 (2021),
   doi:10.1016/j.ijheatmasstransfer.2021.121119; arXiv:2103.01611 (**abstract**: "a given F(Lambda) does not correspond
   to a unique relaxation-times function ... different relaxation-time distributions will, in general, lead to
   different thermal responses"); **J. Cuffe et al.**, *Phys. Rev. B* **91**, 245423 (2015), doi:10.1103/PhysRevB.91.245423
   (**metadata**); P. Chen, I. M. Gamba, Q. Li, L. Wang, arXiv:2502.19533 (2025) (read in the 2 October campaign).
6. Dynamic and hydrodynamic predictions sensitive to scattering details: **Tur-Prats et al.** (item L1-3: l and tau
   "significantly reduced" with four-phonon scattering); **X. Li, Z. Han, X. Ruan, L. Shi**, *Phys. Rev. B* **112**,
   L100302 (2025), doi:10.1103/pp9w-jyb3 (**metadata**); **N. Malviya, N. K. Ravichandran**, *Phys. Rev. B* **108**,
   155201 (2023), doi:10.1103/PhysRevB.108.155201 (**abstract**: Callaway description fails for BAs and BSb);
   **J. Abou Haibeh, S. Huberman**, "Impact of exchange-correlation functionals on predictions of phonon hydrodynamics:
   a study of fluorides, chlorides, and hydrides", *Mater. Res. Express* **13**, 085901 (2026),
   doi:10.1088/2053-1591/ae53a0; arXiv:2311.05107 (**abstract**: "the selection of a functional impacts the prediction
   of thermal conductivity and the window for observation of phonon hydrodynamics").
7. Long-wavelength modes as the locus of first-principles disagreement: **P. Jiang, L. Lindsay, Y. K. Koh** (item
   L1-8: up to 30 % discrepancy between first-principles calculations for low-energy phonons in Si).
8. Rheology (phase 7, **metadata only**): **J. Honerkamp, J. Weese**, *Macromolecules* **22**, 4372 (1989),
   doi:10.1021/ma00201a036 (relaxation spectrum by regularisation); **A. R. Davies, R. S. Anderssen**, *J. Non-Newton.
   Fluid Mech.* **73**, 163 (1997), doi:10.1016/S0377-0257(97)00056-6 ("sampling localization" in relaxation-spectrum
   recovery). OUR INTERPRETATION: the weight-average relaxation time sum G_i lambda_i^2 / sum G_i lambda_i of linear
   viscoelasticity has the structure of tau_mem, and the long-time end of a relaxation spectrum is the part least
   determined by finite-window data — the same structure as L4's long-wavelength freedom. Not checked in the texts.

### What differs
* No source found constrains a symmetry-respecting cubic-force-constant class by static data (all lifetimes and kappa)
  and maps the attainable range of a dynamic moment (memory time, Cattaneo time, nonlocal length, viscosity).
* No source found locates the residual freedom in long-wavelength acoustic modes for a dynamic quantity, although
  items 1 and 7 show that first-principles uncertainty of kappa itself concentrates in long-MFP acoustic phonons.
* No source found tests whether extremal parameter sets transfer across Brillouin-zone discretisations.
* Relation to the earlier audit: L4 extends C2(a) (free-rate class) to the physically parametrised class; it is not
  covered by that campaign.

### Caveats raised by the literature and by L1–L2
* The audited quantity is infrared-controlled (L1) and possibly undefined for the physical crystal (L2). L4's large
  factors are properties of discretised operators; for publication the identifiability question should be posed for
  a mesh-converged, regularised quantity (film-regularised or frequency-cutoff memory time).
* The precision model is one workflow pair at one functional; items 2 and 3 document larger spreads (±15 % in kappa
  across packages; > 20 % across functionals).

**Verdict L4: PARTIALLY KNOWN in spirit (sloppiness, force-constant uncertainty of kappa, non-unique reconstruction);
the specific force-constant-class statement is POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE.**

---

## 7. Consequences for paper 1 (positioning)

1. **Attribution of the formula.** The RTA memory time is the heat-flux relaxation time of Sendra et al. (2021) and
   Beardo et al. (2021, Eq. 2). The conserving N/K is its full-operator generalisation (a conductivity-weighted mean
   relaxon lifetime in the language of Cepellotti & Marzari 2016; no source found defines that weighted mean
   explicitly). Present it as such.
2. **L1 as a convergence warning with numbers.** Cite Esfarjani et al. (Eq. 9) and Vermeersch et al. (Eqs. 19–20) for
   the power counting and Tur-Prats et al. for low-frequency dominance. State the implication for published single-mesh
   ab initio Cattaneo times and RTA nonlocal lengths as an interpretation, and test it on Ge (mesh sequence for Eq. 2 of
   Beardo et al.), which would turn it into a demonstrated statement.
3. **L2 as a corollary plus a model-level caveat.** State the n < 3/2 criterion as a corollary of Vermeersch et al.
   Eq. 19 under the same assumptions; cite the 3D long-time tail (Lepri, Livi & Politi Eqs. 134, 141; Ernst, Hauge &
   van Leeuwen; Brey) and Akhiezer-regime damping (Simons; Brito-Orta & Klemens; Jiang, Lindsay & Koh; Jin) for the
   physical crystal. Do not write that the bulk memory time is plausibly finite without "for the three-phonon
   golden-rule operator".
4. **L3 as the well-posed object.** The thickness-regularised memory time is the right target; replace or validate the
   Matthiessen regulator with a cross-plane film solution (Vermeersch et al. 2016 suppression function; Hua & Minnich
   2015 series) and state boundary conditions (Maassen & Lundstrom). Cite Alvarez & Jou and Krewer et al. as conceptual
   precedents.
5. **L4 in the identifiability literature.** Frame with Gutenkunst, Transtrum, Raue and Kurniawan; use Parks et al. and
   the Phonon Olympics as the precision context; restrict the claims to the stated discretised operators until the
   regularised quantity is analysed.

## 8. Unresolved prior-art uncertainty and limits of the search

1. Sendra et al. (2021) journal text and the Sendra thesis were not obtained (repository unreachable, Cloudflare);
   the RTA expressions are taken from the Beardo thesis, which reproduces them with attribution. A convergence remark
   in the journal paper cannot be excluded.
2. The published supplementary file of Beardo et al. (2021) was not inspected; arXiv v3 was used.
3. Read as abstract or metadata only: Herring 1954; Pomeranchuk 1941; Akhiezer 1939 (no DOI; cited by Klemens 1981);
   Woodruff & Ehrenreich 1961; Simons 1964; Brito-Orta & Klemens 1981; Klemens 1981; Pomeau & Résibois 1975; Ernst &
   Dorfman 1975; Inozemtseva & Sadovnikov 1981/82; Brey 1983; Jin 2024; Jiang, Lindsay & Koh 2016; Phonon Olympics
   2025 (beyond the campaign's own use of its data); the thermoreflectance and film papers of section 5; the
   identifiability, rheology and force-field uncertainty papers of section 6 except Parks et al.; Li et al. 2025.
4. Semantic Scholar was rate-limited (2/41 batch queries; no full-text snippet search). Full-text coverage of
   Elsevier, AIP and ACS journals depended on arXiv versions and author copies.
5. Not searched: Russian-language phonon-kinetics literature (Gurevich; Levinson; Gurzhi) and monographs (Ziman;
   Gurevich, *Transport in Phonon Systems*; Maris in *Physical Acoustics*; Struchtrup 2005). An older explicit statement
   that higher kinetic coefficients of the phonon gas diverge through soft long-wavelength phonons may exist there or
   in second-sound damping theory (Hardy 1970; Beck, Meier & Thellung 1974), not read.
6. Further single-mesh ab initio KCM parameter studies (e.g. GeTe, transition-metal dichalcogenides; titles found)
   were not inspected; they would add instances, not change the verdicts.
7. Alternative names were generated within this review (section 2). Additional terminology from the mathematical and
   theoretical-physics side (for example spectral-measure or Stieltjes-moment language for N/K, or "terminal relaxation
   time") should be searched if supplied.

## 9. Bibliography (all DOIs verified in Crossref; access level as stated in the text)

Abou Haibeh, J., Huberman, S. Mater. Res. Express 13, 085901 (2026), doi:10.1088/2053-1591/ae53a0; arXiv:2311.05107.
Allen, P. B. Phys. Rev. B 97, 134307 (2018), doi:10.1103/PhysRevB.97.134307.
Allen, P. B., Perebeinos, V. Phys. Rev. B 98, 085427 (2018), doi:10.1103/PhysRevB.98.085427.
Alvarez, F. X., Jou, D. Appl. Phys. Lett. 90, 083109 (2007), doi:10.1063/1.2645110.
Alvarez, F. X., Jou, D. J. Appl. Phys. 103, 094321 (2008), doi:10.1063/1.2913057.
Beardo, A. Generalized Hydrodynamic Heat Transport in Semiconductors, PhD thesis, Universitat Autònoma de Barcelona (2021), https://ddd.uab.cat/pub/tesis/2021/hdl_10803_673590/abr1de1.pdf.
Beardo, A. et al. Phys. Rev. Appl. 11, 034003 (2019), doi:10.1103/PhysRevApplied.11.034003.
Beardo, A., Hennessy, M. G., Sendra, L., Camacho, J., Myers, T. G., Bafaluy, J., Alvarez, F. X. Phys. Rev. B 101, 075303 (2020), doi:10.1103/PhysRevB.101.075303.
Beardo, A. et al. Sci. Adv. 7, eabg4677 (2021), doi:10.1126/sciadv.abg4677; arXiv:2007.05487.
Benenti, G., Donadio, D., Lepri, S., Livi, R. "Non-Fourier heat transport in nanosystems", arXiv:2212.09374 (text search; no relevant statement).
Brey, J. J. J. Chem. Phys. 79, 4585 (1983), doi:10.1063/1.446374.
Brito-Orta, R., Klemens, P. G. J. Phys. Colloq. 42, C6-256 (1981), doi:10.1051/jphyscol:1981674.
Camacho de la Rosa, A., Esquivel-Sirvent, R., Becerril, D. J. Appl. Phys. 137, 155103 (2025), doi:10.1063/5.0257299.
Cepellotti, A., Marzari, N. Phys. Rev. X 6, 041013 (2016), doi:10.1103/PhysRevX.6.041013; arXiv:1603.02608.
Cepellotti, A., Marzari, N. "Transport waves as crystal excitations", arXiv:1612.04317.
Chen, P., Gamba, I. M., Li, Q., Wang, L. arXiv:2502.19533 (2025).
Chiloyan, V., Huberman, S., Maznev, A. A., Nelson, K. A., Chen, G. Appl. Phys. Lett. 116, 163102 (2020), doi:10.1063/1.5139069.
Crawford, D. E., Zeng, Y., Vidal, J., Dong, J. Phys. Rev. B 112, 115201 (2025), doi:10.1103/p8wg-p1j3; arXiv:2501.00679.
Cuffe, J. et al. Phys. Rev. B 91, 245423 (2015), doi:10.1103/PhysRevB.91.245423.
Davies, A. R., Anderssen, R. S. J. Non-Newton. Fluid Mech. 73, 163 (1997), doi:10.1016/S0377-0257(97)00056-6.
Di Lucente, E., Simoncelli, M., Marzari, N. Adv. Phys. 75, 83–222 (2026), doi:10.1080/00018732.2026.2685376; arXiv:2608.13339.
Ernst, M. H., Dorfman, J. R. J. Stat. Phys. 12, 311 (1975), doi:10.1007/BF01012067.
Ernst, M. H., Hauge, E. H., van Leeuwen, J. M. J. Phys. Rev. A 4, 2055 (1971), doi:10.1103/PhysRevA.4.2055.
Esfarjani, K., Chen, G., Stokes, H. T. Phys. Rev. B 84, 085204 (2011), doi:10.1103/PhysRevB.84.085204; arXiv:1107.5288.
Forghani, M., Hadjiconstantinou, N. G., Péraud, J.-P. M. Phys. Rev. B 94, 155439 (2016), doi:10.1103/PhysRevB.94.155439.
Forghani, M., Hadjiconstantinou, N. G. Phys. Rev. B 97, 195440 (2018), doi:10.1103/PhysRevB.97.195440.
Forghani, M., Hadjiconstantinou, N. G. Appl. Phys. Lett. 114, 023106 (2019), doi:10.1063/1.5062846.
Forghani, M., Hadjiconstantinou, N. G. Int. J. Heat Mass Transf. 175, 121119 (2021), doi:10.1016/j.ijheatmasstransfer.2021.121119; arXiv:2103.01611.
Frederiksen, S. L., Jacobsen, K. W., Brown, K. S., Sethna, J. P. Phys. Rev. Lett. 93, 165501 (2004), doi:10.1103/PhysRevLett.93.165501.
Gutenkunst, R. N. et al. PLoS Comput. Biol. 3, e189 (2007), doi:10.1371/journal.pcbi.0030189.
Herring, C. Phys. Rev. 95, 954 (1954), doi:10.1103/PhysRev.95.954.
Honerkamp, J., Weese, J. Macromolecules 22, 4372 (1989), doi:10.1021/ma00201a036.
Hoque, M. S. B. et al. Appl. Phys. Lett. 125, 262201 (2024), doi:10.1063/5.0239769.
Hua, C., Lindsay, L. Phys. Rev. B 102, 104310 (2020), doi:10.1103/PhysRevB.102.104310; arXiv:2008.07596.
Hua, C., Lindsay, L., Chen, X., Minnich, A. J. Phys. Rev. B 100, 085203 (2019), doi:10.1103/PhysRevB.100.085203.
Hua, C., Minnich, A. J. Phys. Rev. B 89, 094302 (2014), doi:10.1103/PhysRevB.89.094302.
Hua, C., Minnich, A. J. J. Appl. Phys. 117, 175306 (2015), doi:10.1063/1.4919432.
Inozemtseva, N. G., Sadovnikov, B. I. Theor. Math. Phys. 48, 658 (1981), doi:10.1007/BF01037993; 51, 409 (1982), doi:10.1007/BF01029268.
Jain, A., McGaughey, A. J. H. Comput. Mater. Sci. 110, 115 (2015), doi:10.1016/j.commatsci.2015.08.014.
Jiang, P., Lindsay, L., Koh, Y. K. J. Appl. Phys. 119, 245705 (2016), doi:10.1063/1.4954674.
Jin, J. S. J. Appl. Phys. 135, 045103 (2024), doi:10.1063/5.0187032.
Johnson, J. A. et al. Phys. Rev. Lett. 110, 025901 (2013), doi:10.1103/PhysRevLett.110.025901 (median-MFP remark; context only).
Kefayati, A., Allen, P. B., Perebeinos, V. Phys. Rev. B 105, 235402 (2022), doi:10.1103/PhysRevB.105.235402.
Klemens, P. G. Int. J. Thermophys. 2, 55 (1981), doi:10.1007/BF00503574.
Koh, Y. K., Cahill, D. G. Phys. Rev. B 76, 075207 (2007), doi:10.1103/PhysRevB.76.075207.
Krewer, K. L. et al. Appl. Phys. Lett. 116, 102406 (2020), doi:10.1063/1.5142479; arXiv:2008.07852.
Krishnamoorthy, A. et al. Comput. Phys. Commun. 254, 107337 (2020), doi:10.1016/j.cpc.2020.107337.
Kurniawan, Y. et al. J. Chem. Phys. 156, 214103 (2022), doi:10.1063/5.0084988.
Lepri, S., Livi, R., Politi, A. Phys. Rep. 377, 1–80 (2003), doi:10.1016/S0370-1573(02)00558-6; arXiv:cond-mat/0112193.
Li, T., Jiang, B., Chen, Z. "AC thermal conductivity as a tool for solution mapping from diffusive to ballistic regime", arXiv:2408.05654 (2024).
Li, X., Han, Z., Ruan, X., Shi, L. Phys. Rev. B 112, L100302 (2025), doi:10.1103/pp9w-jyb3.
Lindsay, L., Broido, D. A., Reinecke, T. L. Phys. Rev. B 87, 165201 (2013), doi:10.1103/PhysRevB.87.165201 (ab initio kappa of AlN; context).
Maassen, J., Lundstrom, M. J. Appl. Phys. 117, 135102 (2015), doi:10.1063/1.4916245; 119, 095102 (2016), doi:10.1063/1.4942836.
Machta, B. B., Chachra, R., Transtrum, M. K., Sethna, J. P. Science 342, 604 (2013), doi:10.1126/science.1238723.
Mahan, G. D., Claro, F. Phys. Rev. B 38, 1963 (1988), doi:10.1103/PhysRevB.38.1963.
Malviya, N., Ravichandran, N. K. Phys. Rev. B 108, 155201 (2023), doi:10.1103/PhysRevB.108.155201; arXiv:2502.00337 (2025); arXiv:2606.17829 (2026).
McGaughey, A. J. H. et al. (Phonon Olympics) J. Appl. Phys. 138, 135108 (2025), doi:10.1063/5.0289819.
Minnich, A. J. Phys. Rev. Lett. 109, 205901 (2012), doi:10.1103/PhysRevLett.109.205901.
Parks, H. L., Kim, H.-Y., Viswanathan, V., McGaughey, A. J. H. Phys. Rev. Mater. 4, 083805 (2020), doi:10.1103/PhysRevMaterials.4.083805; arXiv:2005.08835.
Pereira, L. F. C., Donadio, D. Phys. Rev. B 87, 125424 (2013), doi:10.1103/PhysRevB.87.125424.
Pomeau, Y., Résibois, P. Phys. Rep. 19, 63 (1975), doi:10.1016/0370-1573(75)90019-8.
Pomeranchuk, I. Phys. Rev. 60, 820 (1941), doi:10.1103/PhysRev.60.820.
Raue, A. et al. Bioinformatics 25, 1923 (2009), doi:10.1093/bioinformatics/btp358.
Regner, K. T. et al. Nat. Commun. 4, 1640 (2013), doi:10.1038/ncomms2630 (kappa(omega) spectroscopy; context).
Saito, K., Dhar, A. Phys. Rev. Lett. 104, 040601 (2010), doi:10.1103/PhysRevLett.104.040601.
Sendra, L., Beardo, A., Torres, P., Bafaluy, J., Alvarez, F. X., Camacho, J. Phys. Rev. B 103, L140301 (2021), doi:10.1103/PhysRevB.103.L140301.
Sendra, L., Beardo, A., Bafaluy, J., Torres, P., Alvarez, F. X., Camacho, J. Phys. Rev. B 106, 155301 (2022), doi:10.1103/PhysRevB.106.155301.
Shvaika, A., Shpot, M., Schirmacher, W., Bryk, T., Ruocco, G. arXiv:2104.13076 (2021).
Simoncelli, M., Marzari, N., Cepellotti, A. Phys. Rev. X 10, 011019 (2020), doi:10.1103/PhysRevX.10.011019; arXiv:1906.09743.
Simons, S. Proc. Phys. Soc. 83, 749 (1964), doi:10.1088/0370-1328/83/5/306.
Srivastava, Y., Gokhale, A. G., Jain, A. Phys. Rev. B 112, 195201 (2025), doi:10.1103/2n3p-hz62; arXiv:2509.14702.
Tamura, S., Sangu, A., Maris, H. J. Phys. Rev. B 68, 014302 (2003), doi:10.1103/PhysRevB.68.014302.
Torres, P. et al. Phys. Rev. B 95, 165407 (2017), doi:10.1103/PhysRevB.95.165407; arXiv:1606.01149.
Transtrum, M. K. et al. J. Chem. Phys. 143, 010901 (2015), doi:10.1063/1.4923066.
Tur-Prats, J., Han, Z., Beardo, A., Ruan, X., Alvarez, F. X. Nano Lett. 25, 11203–11209 (2025), doi:10.1021/acs.nanolett.5c00855.
Vermeersch, B., Carrete, J., Mingo, N., Shakouri, A. Phys. Rev. B 91, 085202 (2015), doi:10.1103/PhysRevB.91.085202; arXiv:1406.7341.
Vermeersch, B., Mohammed, A. M. S., Pernot, G., Koh, Y. R., Shakouri, A. Phys. Rev. B 91, 085203 (2015), doi:10.1103/PhysRevB.91.085203.
Vermeersch, B., Carrete, J., Mingo, N. Appl. Phys. Lett. 108, 193104 (2016), doi:10.1063/1.4948968.
Volz, S. G. Phys. Rev. Lett. 87, 074301 (2001), doi:10.1103/PhysRevLett.87.074301 (frequency-dependent kappa by MD; context).
Ward, A., Broido, D. A. Phys. Rev. B 81, 085205 (2010), doi:10.1103/PhysRevB.81.085205.
Woodruff, T. O., Ehrenreich, H. Phys. Rev. 123, 1553 (1961), doi:10.1103/PhysRev.123.1553.
Xiang, Z., Jiang, P., Yang, R. J. Appl. Phys. 132, 205104 (2022), doi:10.1063/5.0113257.
Xie, H., Gu, X., Bao, H. Comput. Mater. Sci. 138, 368 (2017), doi:10.1016/j.commatsci.2017.06.047; arXiv:1706.08232.
Yang, F., Dames, C. Phys. Rev. B 91, 165311 (2015), doi:10.1103/PhysRevB.91.165311.
Zeng, Y., Dong, J. Phys. Rev. B 113, 235204 (2026), doi:10.1103/hx28-4ksc; arXiv:2507.05682.
Zhou, H., Zhou, S., Hua, Z., Bawane, K., Feng, T. Appl. Phys. Lett. 123, 192201 (2023), doi:10.1063/5.0173762.
