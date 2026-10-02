# A — A coarea measure with exact resonant quadrature

Independent first pass completed 2026-09-28, following the 2026-09-27 checkpoint. Read only this campaign question/assumptions and primary prior literature; no peer outputs. **Bounded benchmark arithmetic completed and saved.** No material inference, universal convergence theorem, or novelty claim.

## 1. Model, normalization, and the regular surface

Assume a scalar weak-coupling phonon kinetic approximation with dephased/resolved modes, Markov coarse-graining and the requisite factorized Bose statistics. This is not derived from a finite isolated Hamiltonian. The kinetic assumption includes correlation/coarse-graining times short compared with the collision time; mode coherences may be discarded only with a separately justified secular/dephasing limit. Work on a periodic d-dimensional Brillouin torus B, volume V_B, with C^2 branch energies epsilon_s(q)>=epsilon_min>0, fixed T>0, and a continuous nonnegative channel kernel K. Channel permutation/factorial counting is part of K and the channel enumeration contract, not another factor to infer from quadrature. For definiteness enumerate unordered daughter channels exactly once; an equivalent complete ordered daughter integral has its declared exchange factor.

Use normalized one-mode measure dq/V_B. The periodic momentum delta normalized for this measure is V_B sum_G delta^(d)(p-k-b-G). Eliminating b gives b=p-k modulo G and

    Delta(p,k)=epsilon_P(p)-epsilon_A(k)-epsilon_B(p-k),
    dmu_res = dp dk / V_B^2 * delta(Delta(p,k)).       (1)

The measure has units inverse energy. With physical wavevectors, grad_k Delta=hbar(v_B-v_A). For fixed p, if zero is a regular value,

    integral_B f(k) delta(Delta_p(k)) dk/V_B
      = integral_{Sigma_p} f(k) dH^(d-1)(k)
                              /[V_B |grad_k Delta|]. (2)

There is no free geometric normalization constant. Replacing energy by frequency requires the associated delta-function Jacobian. Identical-daughter or reciprocal counting must not be inserted again into (2).

For local coordinates k=(y,z), roots z_m(y) with |partial_z Delta|>=gamma>0 give weights 1/|partial_z Delta|. Use a finite collection of root charts with nonnegative partition functions chi_l summing to one ON the surface. Positive quadrature in y (and p for the weak form) gives positive weights alpha_e and triples lying exactly on Delta=0. Include every root; treat periodic seams once. A zero coordinate derivative is a chart failure, not necessarily a singular surface: change the solved coordinate if another derivative is nonzero.

For a compact regular surface, smooth integrands, uniformly nonsingular charts and consistently refined positive quadrature, this converges as a measure tested against fixed smooth functions. Constants depend on inverse powers of gamma and derivatives. This does not establish convergence of the full collision spectrum or nonlinear dynamics.

## 2. A finite entropy representation that retains the invariants

Surface quadrature alone does not specify how occupations at off-grid points are represented. Choose a finite common function space V_h=span{phi_j} over all branches that **contains the exact energy function epsilon** (enrich an interpolation space by epsilon if needed). Let epsilon=sum_j c_j phi_j. Interpolate the Bose entropy variable, not occupations independently:

    eta_h(q)=sum_j z_j phi_j(q)>0,
    n_h(q)=1/[exp(eta_h(q))-1].                       (3)

Choose positive one-mode quadrature masses m_l at points x_l. Define moments M_j=sum_l m_l phi_j(x_l)n_h(x_l), and capacity matrix C_ij=sum_l m_l phi_i phi_j n_h(1+n_h). Assume full evaluation rank, so C is positive definite. For each resonant triple e define

    r_e,j=phi_j(p_e)-phi_j(a_e)-phi_j(b_e),
    F_e=n_p(1+n_a)(1+n_b)-(1+n_p)n_a n_b,
    dot M=-sum_e alpha_e K_e r_e F_e,
    equivalently C dot z=sum_e alpha_e K_e r_e F_e.  (4)

All quantities at a collision triple use the SAME reconstruction (3). This is a proposed finite moment approximation, not an exact scalar closure of microscopic phonon dynamics.

Let s(n)=(1+n)log(1+n)-n log n, S_h=sum_l m_l s(n_h(x_l)), E_h=sum_l m_l epsilon(x_l)n_h(x_l)=c dot M. Then directly

    dot E_h=-sum_e alpha_e K_e Delta_e F_e=0,
    dot S_h=-sum_e alpha_e K_e (r_e dot z) F_e >=0.  (5)

The sign follows from F_e=R_e[exp(-r_e dot z)-1], R_e=(1+n_p)n_a n_b>0. For z=beta c, every F_e=0 separately: exact Bose detailed balance and nonlinear equilibrium stationarity. These claims concern the discrete quadrature functionals and hold while eta_h is positive; global preservation of the admissible coefficient domain is not proved here.

At equilibrium, set Q_e=n_p^0(1+n_a^0)(1+n_b^0). The linearized entropy-coordinate matrix is

    L_h=sum_e alpha_e K_e Q_e r_e r_e^T >=0,
    C_0 dot(delta z)=-L_h delta z,   L_h c=0.       (6)

Thus the metric as well as the collision form is declared. Replacing exact epsilon by an interpolant epsilon_h and solving Delta_h=0 conserves E_h built from that surrogate, not automatically the true energy. If ||epsilon_h-epsilon||_infinity<=a, physical triple mismatch is at most 3a.

For imperfect roots, c^T L_h c=sum_e alpha_e K_e Q_e Delta_e^2 is an explicit invariant-defect diagnostic; zero cannot arise through cancellations. Nonlinear energy drift is the first expression in (5) with nonzero Delta_e. Root residual rho translates into position error at most approximately rho/gamma in a regular chart. None of these bounds licenses calling a broadened tuple exactly resonant.

## 3. Solvable periodic benchmark

Take dimensionless wavevectors on [-pi,pi), an energy unit E_*>0 and three positive branches

    epsilon_P(p)=E_*(5+u),
    epsilon_A(k)=E_*(3+cos k),
    epsilon_B(b)=2 E_*.

For |u|<1, Delta=E_*(u-cos k), with roots k=+/-arccos u. The normalized fixed-parent and total resonance masses are both

    I(u)=1/[pi E_* sqrt(1-u^2)].                    (7)

Each root weighs 1/[2 pi E_* sqrt(1-u^2)]. More generally, integrate f by evaluating it at these two roots with these weights. The zero-detuning normalization is 1/(pi E_*), not one. This remains a benchmark in higher dimensions if the other torus coordinates are spectators.

A finite basis containing branchwise 1, cos q, sin q contains epsilon exactly. It tests (3)-(6), including Bose stationarity and capacity normalization, with analytic roots. Positive parent quadrature approximates any remaining p-dependence. Required checks: total mass (7), nonnegative dissipation, energy null vector, eventwise equilibrium flux, root residuals and convergence of a nonconstant smooth test function.

## 4. Decisive exclusions and a better manifold when slices fail

* At |u|>1 the resonance set is empty and the exact measure is zero. Positive-width Gaussian tails are a different finite measure.
* At u=1, Delta=E_*(1-cos k) is quadratic at k=0. Formula (7) diverges as u approaches 1 from below. A normalized Gaussian energy delta at u=1 gives mass proportional to 1/sqrt(E_* sigma); it has no finite sigma->0 measure here. Capping the coarea denominator invents a cutoff-dependent operator.
* Delta identically zero on an interval (for example collinear 1D exactly linear acoustic decay away from endpoints) is not an ordinary codimension-one resonance. Acoustic zero, singular Bose factors, crossings/degeneracy and non-smooth branch labels require separate analysis. A vanishing matrix element can sometimes regularize a weighted singularity, but that must be demonstrated quantitatively.
* Fixed-parent singularity does not imply that the global weak form is singular. In joint (p,k) coordinates,

      grad_(p,k) Delta = hbar(v_P-v_B, v_B-v_A).

  If daughters have equal velocity but the parent differs, solve a parent coordinate instead. The local model Delta=p-k^2 has singular fixed-p root density, while integral dp dk f delta(p-k^2)=integral dk f(k^2,k) is regular. A full (2d-1)-dimensional resonant atlas can therefore succeed when a fixed-parent atlas fails. This does not make an infinite pointwise parent rate finite.
* Treating every quadrature endpoint as an unrelated population unknown creates many disconnected reaction triples and artificial conserved quantities. Accurate quadrature of smooth samples alone does not certify dynamical convergence. The common field space in (3) addresses this representation problem, but its approximation/stability still requires audit.

For comparison only, resolving a Gaussian layer at a regular root requires h |grad Delta| much smaller than sigma before sigma tends to zero; the quadratic critical scale is instead h much smaller than sqrt(sigma/curvature). Neither is a universal joint scaling. The proposed regular-surface method has no physical broadening parameter: refine charts/quadrature and root tolerance after declaring the kinetic limit. Near a vanishing group-velocity difference the kinetic approximation itself may require reconsideration.

The regular coarea measure and velocity-degeneracy mechanism are established prior art: [Shi and Eyink, Resonance Van Hove Singularities in Wave Kinetics, Eq. (3) and Sec. 3](https://arxiv.org/html/1507.08320). The present benchmark and finite representation are conditional constructions; no novelty is claimed.

## 5. Executed bounded benchmark

Reproduction: [A-coarea-check.py](A-coarea-check.py), run with the pinned environment interpreter:

    D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe -B A-coarea-check.py

Units E_*=1, time_*=1; K=1, u=0.6, beta=0.7. The calculation uses 64 parent nodes, the two exact analytic roots, branchwise {1, cos q, sin q}, and 128 one-mode quadrature nodes per branch. A seeded small entropy-variable perturbation tests the nonlinear identities. All script assertions passed.

| Diagnostic | Result |
|---|---:|
| Analytic total resonance mass | 0.3978873577297384 |
| Quadrature mass relative error | 3.33e-16 |
| Maximum exact-energy mismatch / Bose forward-reverse relative error | 0 / 0 |
| Relative ||L_h c|| residual | 6.50e-17 |
| Minimum eigenvalue of L_h divided by ||L_h|| | -7.23e-17 (roundoff scale) |
| Minimum capacity eigenvalue | 0.0103263 |
| Nonlinear energy drift | 0 at reported precision |
| Entropy production from event sum | 1.0514127630978124e-5 |
| Relative agreement with one-mode entropy derivative | 1.78e-14 |
| Omit inverse group-gradient: relative mass error | 0.20 |

For f(p,k)=exp(2 cos p)(1+0.1 cos k), the exact integral is I(u)(1+0.1u) I_0(2), where I_0 is the modified Bessel function. Parent trapezoidal quadrature with 4, 8, 16 and 32 nodes gives relative errors 4.453e-2, 2.430e-5, 4.469e-14 and 0 at reported precision. This tests a smooth fixed integrand, not time-evolution convergence.

At the critical point u=1, use delta_sigma(E)=exp[-E^2/(2 sigma^2)]/(sqrt(2 pi) sigma). The substitution k=sqrt(sigma/E_*) x gives

    I_sigma ~ C/sqrt(E_* sigma),
    C=8^(1/4) Gamma(1/4)/(4 pi sqrt(2 pi))
     =0.19357700828395394.

The computed masses for sigma=0.1, 0.025, 0.00625, 0.0015625 are 0.619765, 1.227982, 2.450410, 4.898069. The corresponding relative errors in sqrt(sigma) I_sigma versus C fall from 1.245e-2 to 1.868e-4. Adaptive integration was performed in the rescaled x coordinate, which resolves the shrinking peak; reported rescaled quadrature errors were below 4.8e-14. This confirms the divergent scaling of this example rather than a finite resonant measure.

**Completed first-pass scope:** constructive regular-case algebra, analytic normalization, finite arithmetic, one fixed-integrand refinement sequence, and decisive singular counterexamples. General root topology certification, acoustic/critical weighted integrability, global moment-domain invariance, time discretization preserving the continuous-time identities, and dynamical/operator convergence remain unresolved. No additional calculations or peer comparisons were performed.