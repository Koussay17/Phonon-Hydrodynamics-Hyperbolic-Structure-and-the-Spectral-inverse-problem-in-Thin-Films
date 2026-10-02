# A2 — Physical and geometric admission before an AlN application

2026-09-28. Bounded second-generation assessment after reading B, C, D, L and PI-interface-map. No material calculation or code change. **The campaign supports conditional finite identities and regular-measure benchmarks; it has not established an AlN kinetic model or converged operator.** “True energy” below means the declared reference dispersion epsilon(q), not an experimentally exact dispersion.

## 1. Physical admission: separate three small-scale questions

**Quasiparticles.** Define gamma_s in ENERGY units as a spectral half-width. Require gamma_s/epsilon_s small on the retained domain for a sharp-mode approximation, and document shifts if bare harmonic energies are used. A small width compared with energy does not imply that nearby branches are dynamically distinguishable. Factors converting phono3py half-linewidths, angular frequencies and population decay rates must remain explicit. A numerical integration width sigma is not automatically gamma_s.

**Scalar populations.** For discarded same-q coherences, examine actual branch splittings delta E_ss' and the coherent/collisional couplings. A usual sufficient secular scale separation requires a coarse-graining time t_c with

    max(t_corr, hbar/epsilon, hbar/|delta E_ss'|) << t_c << t_coll

for the pairs one proposes to discard. This is a regime assumption, not a theorem supplied by the ratios. Exact degeneracy cannot satisfy the last splitting condition; scalar closure needs a separate invariant-diagonal/symmetry argument, or matrix occupations and phase-sensitive amplitudes. Approximate clustering does not authorize mixing unequal eigenfrequencies while retaining diagonal frequencies. Rotations must stay within the same q representation if q labels are kept. Passive basis covariance is not physical symmetry. The [Wigner transport derivation, discussion following Eq. (38)](https://arxiv.org/html/2112.06897) explicitly retains coherences and still assumes well-defined phonon excitations; invoking it alone does not derive our collision action.

**Delta-resonance approximation.** Even narrow modes can lie near a singular resonance. For a regular local mismatch, put g=|grad Delta| and H=||Hessian Delta||. A putative physical width w shifts the normal coordinate by order w/g; the local linear-root picture needs w H/g^2 small, as well as controlled variation of the interaction/Bose weight over that distance. Equivalently, inspect the weighted level density A(z)=integral_(Delta=z) KQ/|grad Delta| and its variation over the relevant spectral width. beta*w is a thermal diagnostic. These are local screens, not universal numerical thresholds. Summing individual widths to obtain a triad width itself assumes a spectral lineshape/convolution model. Lorentzian or finite-time tails do not inherit the Gaussian O(sigma^2) bias law.

A linewidth computed by the same unvalidated approximation is not independent validation of that approximation. No linewidth/splitting/curvature comparison for AlN was performed in this campaign.

## 2. Geometry and repeated sectors

For decay, eliminate momentum by b=p-a modulo G. In physical coordinates,

    grad_a Delta=hbar(v_b-v_a),
    grad_(p,a) Delta=hbar(v_p-v_b, v_b-v_a).

A vanishing fixed-parent derivative may be a chart failure or a singular slice of a regular joint surface. Use a different coordinate, including a parent coordinate, when the full gradient is nonzero. When it vanishes, classify the critical set and establish integrability of the WEIGHTED measure; do not cap a denominator. An empty sign-change scan is not an empty-set certificate: C explicitly missed pairs of regular near-critical roots. A/B/C/D also show that positivity, energy-null tests and tiny drift can all coexist with a wrong or missing measure.

Acoustic admission requires separate infrared estimates. At positive temperature, a linear acoustic branch has n approximately k_B T/(hbar c|q|) and n(1+n) approximately const/|q|^2. Check the capacity integral and the collision form independently, including interaction zeros, resonance Jacobians and test-function differences. Removing a ball around q=0 is a cutoff model until the omitted weighted contribution is bounded or shown to vanish. Exactly linear collinear resonances can be critical along whole sets. Branch crossings require smooth local projectors or an explicit matrix treatment, rather than assuming every sorted branch is differentiable.

**Identical daughters are not generically a finite continuum sector.** In the joint 2d-dimensional momentum domain, a=b (same branch and same q) lies on the d-dimensional diagonal p=2a modulo G. With a regular transverse resonance, its intersection has dimension d-1 inside a 2d-1 dimensional resonant surface, and hence zero coarea mass for a locally integrable density. Under smooth transverse conditions, a tube of radius ell around that intersection has mass O(ell^d). If transversality or integrability fails, this conclusion must be reexamined. At an identical-daughter point the fixed-parent gradient is necessarily zero, even if the joint surface remains regular.

The finite-oscillator repeated-mode Hamiltonian/Fock factors remain valid for their own discrete model. They do not justify adding a weighted atom to a continuum diagonal. Conversely, a surface quadrature node on that diagonal may represent a neighborhood of distinct continuum daughter states: attaching another repeated-mode factor solely because of that node's coordinates can change the quadrature. Declare whether an object is a finite oscillator event or a continuum weak-form evaluation. Same degenerate block does not mean identical mode; coincident legs use one shared basis and its symmetric tensor representation.

## 3. Energy and global quadrature admission

Require a single documented reference energy evaluation at all three off-grid legs, root solving against that same function, and a common entropy-variable reconstruction. B's finite proof applies when epsilon belongs to the test space and positive volume quadrature defines the capacity/energy functional. Solving with epsilon_h and reconstructing Bose factors with epsilon is a different scheme. Exact conservation of surrogate epsilon_h must be labeled as such; uniform energy error a only bounds physical triple mismatch by 3a. Physical force-constant/model error is additional to interpolation error.

Root tolerance, surface geometry, interaction interpolation, volume quadrature and basis resolution are separate error sources. Retain the coarea Jacobian, cover every connected component, merge periodic seams, and identify Normal versus Umklapp representatives without introducing a seam discontinuity in periodic energies. Test total measure and nonconstant smooth moments, not only invariants. Check the capacity matrix and collision action under refinement; a collection of disconnected quadrature triples can preserve energy while creating artificial invariants.

A global action needs one reversible-event enumeration convention and the complete updates of all participating modes. Daughter exchange, reciprocal partners, point-group orbits and repeated labels cannot each be blindly multiplied into a saved weight. An orbit factor suffices for an invariant scalar contraction; the full operator needs the correspondingly transformed event rows/test functions.

The [PI interface map](../PI-interface-map.md) establishes that phono3py's g[1] is a signed combination of channels; the tetrahedron path returns fixed-external-frequency integration weights, not a list of globally unique on-shell events. Nonzero weights at mesh labels do not certify resonance there. The existing equal-weight, distinct-index event API is narrower than the proposed weak-form representation. Increasing its resonance tolerance would not implement the latter.

## 4. Alternative formulations: what changes

* **Joint-surface coordinates and entropy moments:** useful reformulations of the declared scalar kinetic model. A joint atlas removes avoidable slice singularities; entropy moments make Bose stationarity and the energy test structural. They do not prove global admissibility of the finite moment ODE or dynamical convergence.
* **Symmetry quotient:** useful compression only after establishing a genuine symmetry of energies, amplitudes, occupations and represented observables. Preserve stabilizers, orbit incidence and coherent sewing. A passive rotation of a degenerate basis supplies no autonomous trace closure.
* **Renormalized dispersions:** potentially the appropriate quasiparticle model if derived self-consistently. They change the reference resonance surface. If the energy depends on occupations, conservation requires a compatible energy functional and its variation, not merely sum epsilon(n)n with an unchanged collision rule.
* **Finite-time or spectral kinetics:** potentially necessary physical extensions, not conserving replacements for the same on-shell population model. First-order transition amplitudes give the normalized energy kernel

      delta_t(Delta)=t/(2 pi hbar) sinc^2[Delta t/(2 hbar)].

  Its positive weight describes transition probability divided by observation time; it is not by itself a time-local Markov generator. The usual kinetic separation places t between correlation/oscillation and collision times. [Shi–Eyink, Eqs. (38)–(41)](https://arxiv.org/html/1507.08320v1) distinguish finite-time and Lorentzian regularizations from the limiting measure. At finite t, bare harmonic energy can exchange with interaction/correlation energy in a microscopic description. Merely inserting this positive off-shell kernel into the usual common-rate Bose flux still gives, at the physical Bose state,

      dot E=sum_e weight_e R_e Delta_e[1-exp(-beta Delta_e)] >=0,

  strictly positive for any active off-shell support. Spectral-function kinetics can instead impose frequency conservation on integrated spectral energies, but needs a consistent self-energy/vertex, equilibrium relation and conserved energy accounting. No such theory is supplied by replacing a Gaussian with a sinc or Lorentzian.

## 5. Decisive next tests and present boundary

1. **Geometry test:** manufacture a periodic joint surface containing a fixed-parent fold but no joint critical point. Integrate it using two independent coordinate choices; check normalization, component coverage, periodic seams and shifted meshes. Add an actually critical case and coupling zeros to test weighted, rather than unweighted, integrability.
2. **Representation test:** refine volume quadrature, surface quadrature and basis independently with exact-energy enrichment. Check nonlinear Bose drift, physical-energy defect, entropy production, capacity conditioning and action on smooth non-invariant test functions. Include a deliberately mismatched interpolant as a negative control.
3. **Repeated/coherent-sector test:** shrink a tube around a transverse daughter diagonal and verify vanishing measure. Separately rotate an exactly degenerate same-q sector and test the proposed scalar closure or retain its matrix action; a block norm alone is insufficient.
4. **Material admission evidence, still absent:** document reference dispersion/derivatives and their error, phase/counting conventions, relevant linewidth/splitting and geometric scales, and acoustic/critical exclusions with weighted remainder estimates before assigning a material interpretation. A first accepted patch would establish that patch only, not a full AlN operator.

**Result:** the regular scalar weak-form route remains a qualified candidate. The current evidence does not admit critical/acoustic/coherent sectors automatically, does not turn source mesh weights into physical events, and does not establish a material rate. No further calculations were undertaken for this report.