# Branch C: infrared moments and grazing-angle limits

Status: independent first pass, 2026-09-16. No other branch reports were read.
Method: radial asymptotics, causal scalar RTA response, angular boundary-layer
analysis, and independent quadrature benchmarks. These are not derivations of
a collision-matrix closure. Claims below are **derived under assumptions** and
await the campaign proof audit; numerical statements are explicitly finite-grid.

## 1. Scope and assumptions

Use `n` for spatial dimension and `L` for film thickness. At a fixed temperature
`T>0`, suppose an acoustic branch has, as `q=|q| -> 0`,

\[
 \omega(q,\Omega)=c(\Omega)q+o(q),\quad
 r(q,\Omega)=a(\Omega)q^\alpha[1+o(1)],\quad \alpha>0.
\]

Assume the relevant velocities have finite limits, `a` is bounded above and
away from zero, and angular weights are integrable. The directional heat
coupling and radial measure can then be written

\[
 C(q)v_i(q)^2\,d\mu(q)=b_i(\Omega)q^{n-1}[1+o(1)]dq\,d\Omega,
 \qquad b_i\ge0.
\]

The coefficient includes the Bose capacity limit `k_B`, normalization, and
directional velocity. It is positive on an angular set of nonzero measure.
Rates and velocities must be smooth enough for the stated uniform asymptotics.
There are no additional singularities elsewhere in the Brillouin zone. These
assumptions are a local model, **not established for the unresolved AlN modes**.
The 300 K total RTA rates are not identified with normal or resistive rates.

If an angular rate coefficient vanishes, the angular integrals below can
diverge independently of the radial criterion. If coupling behaves instead as
`q^beta`, replace `n` in the radial criteria by `n+beta`. For dispersion
`omega~q^s` with nonzero directional coupling, `beta=2s-2`. Temperature going
to zero is a separate, nonuniform limit; all radial results here hold at fixed
positive temperature.

## 2. Static and memory moments have different convergence criteria

Define

\[
 M_m=\int C(q)v_i(q)^2r(q)^{-m}\,d\mu(q),\quad
 \kappa(0)=M_1,\quad \tau_{\rm stat}=M_1/M_0,
 \quad\tau_{\rm mem}=M_2/M_1.
\]

Near zero,

\[
 M_m^{\rm IR}\sim\int b_i(\Omega)a(\Omega)^{-m}d\Omega
                    \int_0^{q_0}q^{n-1-m\alpha}\,dq.
\]

Therefore a nontrivial moment is finite iff `m alpha < n`; equality gives a
logarithmic divergence, and `m alpha > n` gives a power divergence.

| Quantity | Radial convergence condition |
|---|---|
| DC conductivity and static time | `alpha < n` |
| Memory slope/time, assuming finite DC | `2 alpha < n` |
| Second temporal derivative of conductivity | `3 alpha < n` |

An acoustic continuum can have finite conductivity without a finite memory
time or a positive collision spectral gap. In 3D, the example `r~q^2` gives
finite `M1`, divergent `M2`, and no analytic first-order expansion of the
continuum response at zero frequency. **This example is not a finding that
the actual AlN rates follow `q^2` below the mesh cutoff.**

The formal coefficient of `k^2` in a parity-paired scalar RTA constitutive
response `integral C v_i^2/(r+i k.v)` involves
`integral C v_i^2 (v.nhat)^2/r^3`; with linear acoustics and nondegenerate
couplings it requires `3 alpha<n`. This is a condition on that Taylor
expansion, not an identification of hydrodynamic viscosity.

### Mesh/cutoff scaling

For `q_min~1/N`, the omitted convergent moment scales as
`N^{-(n-m alpha)}`. A divergent cutoff moment instead scales as
`log N` at equality and `N^{m alpha-n}` above the threshold. Thus, if
`n/2<alpha<n`, `kappa_N` can approach a finite limit while
`tau_mem,N~N^{2 alpha-n}` diverges. At `alpha=n/2`, memory grows as `log N`.
If `alpha=n`, its ratio grows as `N^n/log N`; if `alpha>n`, it grows as
`N^alpha`. Constants depend on cutoff shape and quadrature.

This is why agreement with a published DC sum cannot certify memory
convergence, and why a mesh sequence must retain the same physical scattering
model. Deleting modes from one mesh is a sensitivity diagnostic, not refinement.

## 3. Continuum low-frequency response and long-time tail

The causal RTA conductivity for `Re p>0` is

\[
 \kappa(p)=\int\frac{C v_i^2}{r+p}\,d\mu
          =\int_0^\infty e^{-pt}K(t)\,dt,\qquad
 K(t)=\int C v_i^2e^{-rt}\,d\mu.
\]

Let `nu=n/alpha` and

\[
 B={1\over\alpha}\int b_i(\Omega)a(\Omega)^{-\nu}d\Omega.
\]

The change of variable `r=a q^alpha` gives transport-weighted rate density
`rho(r)~B r^{nu-1}` and

\[
 K(t)\sim B\Gamma(\nu)t^{-\nu}.
\]

For the leading power-law density and noninteger `nu`, the nonanalytic term
after all finite Taylor terms is

\[
 \kappa(p)=\sum_{\substack{k\ge0\\k<\nu-1}}(-p)^kM_{k+1}
            +{B\pi\over\sin(\pi\nu)}p^{\nu-1}
            +o(p^{\nu-1}).
\]

Here the powers use the causal branch, continued from positive real `p`;
harmonic forcing with `exp(-i Omega t)` uses `p=0^+-i Omega`. For `0<nu<1`
the sum is empty and the leading conductivity diverges. At integer `nu=j`,
the corresponding term is

\[
 (-1)^{j-1}B p^{j-1}\log(r_*/p),
\]

with finite Taylor terms of order less than `j-1`; the matching scale `r_*`
affects analytic terms, not the leading logarithmic coefficient. Subleading
rate-density corrections require their own asymptotic control; the displayed
leading singular term follows from the local power law after the indicated
finite subtractions.

In particular, for `1<nu<2`,

\[
 \kappa(0)-\kappa(p)\sim
 {B\pi\over\sin[\pi(\nu-1)]}p^{\nu-1}.
\]

In 3D with `alpha=2`, this is a square-root response, with a `t^{-3/2}`
kernel. Every finite positive-rate mesh instead produces a rational function
analytic around `p=0`, with nearest pole at `-r_min`. Mesh refinement can make
that apparent analytic window shrink to zero. To see continuum scaling at a
nonzero frequency one needs `r_min << |p| << r_asymptotic`, and convergence at
fixed `p` before the zero-frequency limit.

### Exact counterexample to using DC agreement as memory evidence

The dimensionless model `n=3`, `a=1`, `0<q<1`, `b=1` has

\[
 \kappa(p)=\int_0^1{q^2\over q^2+p}dq
 =1-\sqrt p\arctan(1/\sqrt p)
 =1-\frac\pi2\sqrt p+p+O(p^2).
\]

With a hard lower cutoff `1/N`, `M1=1-1/N` and `tau_mem=N` exactly.
Uniform radial midpoint and right-endpoint rules both produce `M1=1`
exactly for every `N`, while their second moments grow, respectively, as
`(pi^2/2)N` and `(pi^2/6)N`. Independent discretizations agree perfectly
on DC and both have divergent memory, with different finite-mesh coefficients.
The benchmark script checks these scalings and the finite-frequency exact
formula. This construction decisively separates moment convergence from
implementation agreement.

## 4. Grazing-angle boundary layer and the thin-film limit

Consider only the stated stationary in-plane RTA problem: identical reflecting
walls, constant in-plane temperature gradient, no transmission, fixed
specularity `0<=p_s<1`. For a shell with speed `v`, lifetime `tau`,
`Lambda=v tau`, set `epsilon=L/Lambda` and `mu=|v_z|/v`.
The implementation's suppression is

\[
 S(K,p_s)=1-{(1-p_s)K(1-e^{-1/K})\over1-p_s e^{-1/K}},
 \qquad K=\mu/\epsilon.
\]

For fixed `p_s<1` and `K>>1`,

\[
 S(K,p_s)\sim {1+p_s\over2(1-p_s)}K^{-1}.
\]

But `S(0,p_s)=1`: the grazing sector `mu=O(epsilon)` remains unsuppressed.
In an isotropic 3D shell the normalized in-plane conductivity is

\[
 R(\epsilon,p_s)={3\over2}\int_0^1(1-\mu^2)
                   S(\mu/\epsilon,p_s)d\mu.
\]

Matching the grazing sector to `epsilon<<mu<=1` gives

\[
 R(\epsilon,p_s)\sim
 {3(1+p_s)\over4(1-p_s)}\epsilon\log(1/\epsilon).
\]

This is a fixed-specularity limit; it is not uniform as `p_s->1`.
At exactly `p_s=1`, `R=1` for every thickness. The thin and specular limits
therefore do not commute. For diffuse walls,

\[
 R(\epsilon,0)=\frac34\epsilon
 [\log(1/\epsilon)+1-\gamma]+O(\epsilon^2).
\]

An independent exact expression is

\[
 R=1-{3\over2\epsilon}\left[{1\over4}-E_3(\epsilon)
                                         +E_5(\epsilon)\right],\qquad
 E_j(x)=\int_1^\infty e^{-xt}t^{-j}dt.
\]

It should be evaluated at high precision for small `epsilon`, because its
large terms cancel. The stable angular integral avoids that cancellation.
The leading diffuse asymptotic agrees with Sondheimer's original thin-film
formula, equation (21): [Sondheimer, Advances in Physics 1 (1952)]
(https://www.mit.edu/~levitov/8513/Sondheimer_AdvPhys1952.pdf).
This is a transport-equation analogy; it does not import electron material
parameters into phonon transport.

For smooth anisotropic dispersion, a nondegenerate codimension-one set
`v_z=0` with nonzero tangential coupling produces the same logarithmic
mechanism. Degenerate grazing sets, vanishing coupling, or a finite-measure
set of exactly two-dimensional modes can change the asymptotic. The supplied
AlN representatives do not alone establish all such continuum assumptions.

### Why a grazing node creates a false floor

For a fixed finite quadrature with total bulk transport fraction `w0` at
exactly `v_z=0`, every fixed `p_s<1` gives `R_mesh(L)->w0` as `L->0`.
For a regular continuum angular distribution, that set has measure zero and
`R(L)->0` (assuming integrability or the regularized radial behavior below).
Thus fixed-mesh thin-film convergence is nonuniform. Merely shifting a rule
to avoid an exactly grazing node replaces an artificial floor with possible
underestimation; it does not resolve the narrow angular layer.

Resolving each relevant shell requires angular spacing much smaller than
`L/(v tau)` for an ordinary uniform rule. A necessary weaker check for a
grazing-node rule is `w0 << epsilon log(1/epsilon)` if its floor is to be
negligible. Adaptive integration or local analytic cell integration is more
efficient than assigning a finite full-cell weight to the grazing velocity.

### Interaction with radial acoustic divergence

At fixed `L>0`, the low-q part of a diffuse-film shell behaves as
`q^{n-1} tau R(L/(v tau)) ~ const L q^{n-1} log(tau)`.
For `tau~q^{-alpha}`, this radial integral is finite for any finite
`alpha>0` and `n>0` under the angular assumptions above. Diffuse boundaries
can regularize a bulk static infrared divergence even though perfectly
grazing directions themselves never collide with a wall.

This regularization does not follow by adding `|v_z|/L` to the bulk rate and
then applying the same suppression: that would count the walls twice.
Nor does it justify applying the in-plane formula to cross-plane FDTR.

## 5. Quantitative diagnostics from the supplied AlN grid

Source: `theory/aln/rao_300K_modes.npz`, 793 irreducible q representatives,
12 branches, multiplicity sum `24^3=13824`, 300 K. Frequencies below
`0.4992881051 THz` have no positive-frequency samples. The lowest resolved
rate among nominally heat-carrying samples is `4.119943127e8 s^-1`, giving
`tau_max=2427.218 ps` and `r_min/(2 pi)=65.5709 MHz`. Tiny directional
velocity roundoff can formally make a nominally inactive mode active; the
minimum reported rate is therefore a conservative analytic-window bound.

| Finite-grid diagnostic | Basal | c axis |
|---|---:|---:|
| DC RTA conductivity (W m^-1 K^-1) | 309.775157 | 302.442616 |
| Memory time (ps) | 149.304006 | 231.655151 |
| Fraction of DC weight below 1 THz | 1.8243% | 3.6759% |
| Fraction of second-moment weight below 1 THz | 14.2846% | 32.0862% |
| Fraction of DC weight below 2 THz | 13.4839% | 16.1589% |
| Fraction of second-moment weight below 2 THz | 61.6075% | 74.4215% |
| Fraction of third-moment weight below 2 THz | 89.4419% | 96.0523% |
| Top ten irreducible modes' second-moment fraction | 56.5276% | 72.8560% |
| One largest irreducible mode's second-moment fraction | 13.0342% | 22.0109% |

Irreducible-mode counts include multiplicity in their weights and should not
be mistaken for independent full-zone observations. The two lowest-frequency
c-axis modes alone carry 2.7263% of DC but 28.5594% of the second moment.
Deleting all resolved modes below 2 THz changes the computed memory time to
66.2554 ps basal and 70.6741 ps c axis. This is a stress test of tail
dependence, **not** an estimate of the missing infrared contribution.

Descriptive fits of `log r` versus `log frequency`, using only the first
three branch indices and degeneracy weights, vary substantially with fitting
window. For branch 1 the fitted exponent changes from 1.286 (six points below
1.5 THz) to 1.630 (60 points below 4 THz); branch 3 changes from 1.337
(three points below 2 THz) to 2.084 (15 points below 4 THz). These windows
mix directions, may already depart from linear dispersion, and lack samples
approaching zero. They cannot determine which side of the `alpha=1.5`
memory-convergence threshold the unresolved 3D continuum occupies.

### Finite-grid dynamic response

Direct complex sums `kappa(-i Omega)=sum M1_mode/(1-i Omega tau_mode)` are
well defined at nonzero frequencies. A single pole fitted to the memory time
has the following relative error against that exact finite-grid RTA sum:

| Modulation frequency | Basal error | c-axis error |
|---|---:|---:|
| 1 MHz | 0.000240% | 0.000889% |
| 10 MHz | 0.0240% | 0.0882% |
| 100 MHz | 2.071% | 5.749% |
| 1 GHz | 42.42% | 61.39% |

These are model-reduction errors for the supplied discrete RTA spectrum,
not experimental errors and not evidence of continuum convergence. A small
relative error against the total response can hide a large relative error
in the much smaller dynamic correction; the synthetic benchmark records both.

### Actual grazing contamination

The exact-zero `v_z` nodes carry 12.9660024% of basal DC transport and
33.3808655% of its second moment. The threshold `|v_z|<1e-10 m/s` gives
12.9664860% and 33.3808765%, respectively. Thus the reported grazing weight
is not primarily an arbitrary floating-point threshold effect.

For diffuse boundaries on this fixed mesh:

| Thickness | Computed in-plane kappa (W m^-1 K^-1) | Near-zero-vz fraction of film result |
|---|---:|---:|
| 500 nm | 231.6431 | 17.3400% |
| 100 nm | 163.3198 | 24.5940% |
| 30 nm | 110.0318 | 36.5049% |
| 10 nm | 75.0790 | 53.4996% |
| 3 nm | 54.5262 | 73.6654% |
| 1 nm | 46.1615 | 87.0140% |

The exact finite-grid floor is approximately 40.16 W m^-1 K^-1. Values at
subatomic thicknesses in the JSON are only formal limit probes of the
discrete formula, not physical predictions. At 10 nm more than half the
computed film conductivity is assigned to the near-zero-normal-velocity
sector, making angular refinement indispensable. The fraction by itself
is not the error: nearby continuum grazing directions also conduct.

## 6. Independent numerical verification and reliability

Artifacts:

- `experiments/infrared_audit.py`: complete read-only data audit and synthetic checks.
- `experiments/infrared_diagnostics.json`: single-grid AlN diagnostics.
- `experiments/infrared_benchmarks.json`: independently solvable radial model and angular refinements.

Run with:

```powershell
python -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260916-210742-aln-spectral-closure\experiments\infrared_audit.py
```

The original repository is unchanged. Heat capacities were independently
evaluated through the hyperbolic-sine expression and agree with the repository
formula within `8.6e-16` relative. A separately written diffuse suppression
evaluation differs by at most `3.4e-16` absolute on the tested data. Positive
moment sums use `math.fsum`, with NumPy comparisons recorded, avoiding
subtractive cancellation. The Laplace secant is evaluated as
`sum A tau^2/(1+p tau)` rather than subtracting nearby conductivities.

For isotropic diffuse angular integrals, 70-digit exponential-integral values
were compared with a separately implemented, split adaptive angular integral
in double precision for `epsilon=1e-1,...,1e-8`. Maximum observed relative
disagreement was `2.74e-14`. The known thin-film two-term asymptotic reaches
relative error `4.68e-8` at `epsilon=1e-6` and `3.54e-10` at `1e-8`.

Fixed-rule refinement demonstrates a real under-resolution mechanism:

| epsilon | angular N | Gauss-Legendre relative error | Uniform trapezoid relative error |
|---|---:|---:|---:|
| 1e-3 | 32 | -4.317% | +374.36% |
| 1e-3 | 128 | -0.00718% | +72.875% |
| 1e-3 | 512 | +6.61e-8% | +9.295% |
| 1e-6 | 32 | -46.508% | +219403% |
| 1e-6 | 128 | -27.271% | +54804% |
| 1e-6 | 512 | -8.944% | +13662% |

The trapezoid includes the grazing endpoint, with exact artificial floor
`3/(4N)` at fixed N; Gauss avoids that node but misses the narrow sector when
its smallest angle is too large. Thus simply removing grazing samples is
not a defensible correction to the AlN data.

### Explicit numerical classification

- **Converged:** independent synthetic angular reference quadratures, within
  the tested range and displayed quantitative tolerance; algebraic moment and
  finite-frequency synthetic benchmarks are also independently checkable.
- **Stable finite sums:** AlN moments and finite-frequency discrete RTA sums;
  no sign of roundoff being the limiting uncertainty. Input-rate sum residual
  is about `9e-10` relative, consistent with finite source precision and
  irrelevant to the much larger convergence uncertainties.
- **Under-resolved:** fixed angular quadratures in sufficiently thin-film
  limits, demonstrated against an independent exact reference.
- **Inconclusive:** continuum convergence of actual AlN static/memory moments,
  finite-frequency response, and finite-thickness film results. Only one
  physical q mesh is supplied. High moments and thin films are especially
  exposed by the sensitivity diagnostics.
- **Not assessed:** collision-operator eigenvalue convergence, physical DFT
  convergence, temporal/spatial PDE discretization, interface boundary closure,
  actual surface specularity, and experimental adequacy. There is no time
  integrator or finite-domain PDE solver here, so CFL/temporal-truncation and
  artificial outer-boundary errors do not apply to these direct sums.

Consistency of the formulas with the stated RTA boundary problem, stability
of their floating-point evaluation, and convergence to the continuum are
distinct. The available evidence establishes the first two much more strongly
than the third. No physical AlN result should be called continuum converged
on this evidence alone.

## 7. Decisive next diagnostics and failed shortcuts

1. Recompute a physically consistent nested q-mesh sequence, e.g. 24, 36, 48,
   72 or finer, with scattering broadening/energy-conservation tolerances also
   converged. Record `M0,M1,M2,M3`, lowest active rates, cumulative frequency
   and lifetime fractions, and each directional response at fixed physical
   modulation frequency. More samples at the same grid do not substitute.
2. Plot mesh increments and fit an observed order only in an established
   asymptotic range. Test finite-limit power fits against logarithmic and
   divergent fits suggested above. Do not Richardson-extrapolate a divergent
   memory moment to a finite number.
3. Resolve acoustic rates by branch and direction on approach to Gamma;
   separate angular prefactors from radial exponents. Establish whether rate
   broadening introduces a false positive infrared floor. Use the dominant
   total-rate power, not an isotope-only exponent.
4. Converge the *dynamic correction* as well as total `kappa(p)`, and study
   the moving ratio `|p|/r_min`. Low-frequency curves on a single finite mesh
   can manufacture an apparently finite memory slope.
5. For each target film thickness, refine or integrate q cells crossing
   `v_z=0`; monitor the weight of the unresolved sector. Compare ordinary,
   shifted, and local angular/cell quadratures without relabelling them as
   independent physical meshes. State specularity and boundary model.

Failed shortcuts: inferring infrared convergence from reproduced published
DC values; assigning a physical continuum exponent from sparse fitted
frequency windows; treating low-frequency deletion as refinement; treating
absence of long mean free paths on one grid as absence in the continuum;
discarding grazing modes as a cure; adding a boundary rate before applying
the same walls' suppression; or presenting an in-plane reflecting-wall
calculation as cross-plane FDTR.

## 8. Strongest justified conclusion

The supplied data support accurate, reproducible **discrete RTA** moments and
boundary scenarios. They do not identify the continuum infrared exponent,
certify a finite continuum memory time, or resolve the thin-film grazing
limit. A finite static conductivity is compatible with nonanalytic dynamics
and divergent memory. For the present grid, the second moment is strongly
concentrated in the lowest resolved acoustic modes, and the thin-film response
develops a sizable finite-node grazing floor. These two obstructions require
separate radial and angular convergence tests.
