# Branch A: projection and resolvent elimination

Independent first-pass report, 2026-09-16. Read only the campaign question and assumptions, the designated AlN note, and the research protocol; no other branch report was consulted. These are conditional derivations, not identification of real-AlN hydrodynamic coefficients. The initially attempted write was blocked by a usage-limit rejection; this write follows the user's request to continue.

## 1. Variables and assumptions

Use entropy coordinates and the linearized kinetic equation

\[
 \partial_t y+V_i\partial_i y=-Cy+s,\qquad C=C^*\ge0,\quad V_i=V_i^*.
\]

Repeated spatial indices are summed. Work first with a finite full-Brillouin-zone discretization on a periodic spatial domain, at one uniform reference temperature. The velocity matrices are diagonal; the exact collision operator is generally not. Detailed balance and its exact invariants must be verified independently of modal lifetime data.

Let E have orthonormal columns spanning the proposed slow space; define Pi=EE*, Q=I-Pi, and y=Ea+w, Qw=w. Include every exact invariant that couples to the problem. In particular,

\[
 e_\mu=\epsilon_\mu\sqrt{W_\mu n^0_\mu(1+n^0_\mu)},\quad Ce=0,
 \qquad e^*e=k_BT^2 C_V.
\]

Here C_V denotes volumetric heat capacity. For real spatial wave vector k set

\[
 V_k=k_iV_i,\quad A=E^*CE,\quad K=QCE,\quad D=QCQ|_{Q\mathcal H},
\]
\[
 U=E^*V_kE,\quad B=QV_kE,\quad W=QV_kQ|_{Q\mathcal H},
 \quad X=K+iB,\quad Y=K^*+iB^*.
\]

All block coefficients have units of inverse time. Crucially, Y is not X* at fixed nonzero k; confusing them reverses the viscosity sign.

Two choices of E have different assumptions. A spectral subspace of the TOTAL operator has K=0 exactly. Its overlap with momentum determines its physical interpretation. The normal-collision kernel, span{e,P_x,P_y,P_z}, instead has K=QC_R E in general. A normal/resistive split helps identify mechanism; a complete total operator and its slow eigenvectors can identify dynamics even without that split.

## 2. Exact elimination

Let z be the causal Laplace variable and f=y(0)+s-hat(z). The block system is

\[
 \begin{pmatrix}z+A+iU&Y\\X&z+D+iW\end{pmatrix}
 \binom{a}{w}=\binom{f_P}{f_Q}.
\]

With R_Q=(z+D+iW)^{-1}, direct elimination gives

\[
 \boxed{S(z,k)=z+A+iU-YR_QX},\qquad
 Sa=f_P-YR_Qf_Q,\qquad w=R_Q(f_Q-Xa). \tag{A1}
\]

This is exact and closed for a, but frequency dependent. For zero source its time-domain form is

\[
 \dot a+(A+iU)a-\int_0^tY e^{-(D+iW)(t-r)}X a(r)\,dr
 =-Y e^{-(D+iW)t}w(0). \tag{A2}
\]

When K=0, the convolution term is PLUS the integral with kernel B* exp[-(D+iW)t] B. The right side records fast initial data. It cannot be dropped for arbitrary initial distributions.

## 3. Controlled local closure and its limits

ASSUMPTION: D >= gamma I with gamma>0. This is an operator coercivity condition, not a consequence of positive listed modal lifetimes.

For sigma=Re z>=0,

\[
 \|R_Q\|\le(\gamma+\sigma)^{-1},\qquad
 \|e^{-(D+iW)t}\|\le e^{-\gamma t}.
\]

The exact resolvent identity gives

\[
 R_Q-D^{-1}=-R_Q(z+iW)D^{-1},\qquad
 \|R_Q-D^{-1}\|\le
 \frac{|z|+\|W\|}{\gamma(\gamma+\sigma)}. \tag{A3}
\]

No convergent Neumann series is needed for this inequality. Replacing R_Q by D^{-1} yields

\[
 S_0=z+A_0+iU_0+B^*D^{-1}B,\tag{A4}
\]
\[
 A_0=A-K^*D^{-1}K,\qquad U_0=U-K^*D^{-1}B-B^*D^{-1}K.
\]

The collision Schur complement A_0 is positive semidefinite, U_0 is Hermitian, and B*D^{-1}B is positive semidefinite. Therefore this local Fourier generator is accretive: it does not introduce growing modes. Stability is distinct from approximation accuracy.

The quantitative error is

\[
 \boxed{\|S-S_0\|\le
 \frac{(\|K\|+\|B\|)^2(|z|+\|W\|)}
 {\gamma(\gamma+\sigma)}}. \tag{A5}
\]

For K=0 the squared prefactor is exactly ||B||^2. Since ||B||,||W|| <= |k| v_max, the retained k^2 term scales as k^2 v_max^2/gamma and its controlled correction has an additional factor (|z|+|k|v_max)/gamma.

For sigma>0 both inverse Schur operators have norm at most 1/sigma, so

\[
 \|S^{-1}-S_0^{-1}\|\le\|S-S_0\|/\sigma^2. \tag{A6}
\]

For the exact S, test the full accretive block operator on (a,-R_QXa) to obtain this inverse bound. Near a hydrodynamic pole, a better test is ||S_0^{-1}|| ||S-S_0||<1, using inverse perturbation. A small expansion parameter does not itself give a uniform relative response error near poles. Fast-sector forcing also requires the correction -Y(R_Q-D^{-1})f_Q.

### Slow coordinates need an alignment test

At k=0,

\[
 S(z,0)=A_0+z[I+K^*D^{-2}K]+O(z^2).
\]

Thus an order-one K/gamma changes the coefficient of the time derivative at leading order. Example:

\[
 C_{\rm pair}=\begin{pmatrix}1+\varepsilon&1\\1&1\end{pmatrix},
 \quad E=(1,0)^T,\quad D=K=1.
\]

Then S=z+1+epsilon-(z+1)^{-1}=epsilon+2z+O(z^2). The actual small decay rate is epsilon/2+O(epsilon^2), while static elimination predicts epsilon. Append an independent zero block if an exact energy invariant is desired. An arbitrary retained moment can fail by an order-one factor despite |z|/gamma tending to zero.

A sufficient controlled drift regime is C=C_N+C_R with C_N E=0, C_N restricted to Q >= gamma_N, C_R>=0, and rho=||C_R|| << gamma_N. Then

\[
 D\ge\gamma_N,\quad \|K\|\le\rho,\quad
 \|A-A_0\|\le\rho^2/\gamma_N,\quad
 \|K^*D^{-2}K\|\le(\rho/\gamma_N)^2.
\]

The global norm condition on C_R is sufficient and often stronger than needed. Slow-space residuals and the actual resolvent can supply sharper criteria.

## 4. Anisotropic coefficients and normalization

For K=0 put B_i=QV_iE. Then

\[
 \mathsf D_{ij}=B_i^*D^{-1}B_j,\qquad
 \partial_t a+A a+U_i\partial_i a
 -\mathsf D_{ij}\partial_i\partial_j a=0. \tag{A7}
\]

The minus sign becomes positive dissipation in Fourier space. The coefficient has units m^2/s and satisfies

\[
 k_ik_j a^*\mathsf D_{ij}a
 =\|D^{-1/2}k_i B_i a\|^2\ge0.
\]

Each tensor is computable through constrained cell problems D chi_{j beta}=QV_j E_beta. No dense inverse is needed, but collision action beyond the diagonal is needed.

For normal invariants define S_P=P*P and G_{i alpha}=e*V_i P_alpha. Reciprocal mode pairing q -> -q implies e even, P odd, V_i odd, hence e*P=e*V_i e=P*V_iP=0. Crystal inversion symmetry is not required; reciprocal pairing and collision parity must hold.

A displaced Bose distribution linearizes to

\[
 y_{\rm slow}=\frac{e\,\delta T}{k_BT^2}+\frac{Pu}{k_BT},
 \quad \delta U=C_V\delta T,
 \quad p=\frac{S_Pu}{k_BT},
 \quad j_i^{\rm drift}=\frac{G_{i\alpha}u_\alpha}{k_BT}.
\]

Heat current and momentum therefore need not have a common scalar proportionality. In the normal-only long-wavelength limit the two energy-coupled wave speeds are +/- c(n), with

\[
 c^2(n)=\frac{G_n S_P^{-1}G_n^*}{k_BT^2 C_V},
 \qquad G_n=n_iG_i. \tag{A8}
\]

This assumes independent momentum columns and nonzero coupling. Two transverse drift combinations have zero leading propagation speed. With resistance, underdamping must be checked from the full reduced dispersion matrix.

The leading resistance M=P*C_R P gives generalized rates of (M,S_P). These are not exact decay rates when QC_RP is nonzero; the previous rho/gamma_N bounds control that omission when applicable.

The distinct leading transport tensors are

\[
 \kappa^{\rm fast}_{ij}=\frac{e^*V_iQD^{-1}QV_je}{k_BT^2},\qquad
 \eta_{\alpha i,\beta j}=\frac{P_\alpha^*V_iQD^{-1}QV_jP_\beta}{k_BT}.
\]

They are exact static-closure coefficients when K=0, and leading normal-invariant coefficients under weak resistance. Nonzero K changes resistance and first-gradient coupling according to (A4). Kappa-fast excludes the retained drift channel and is not the total DC conductivity. Under reciprocal parity the thermal source is odd and the viscous source even. Conductivity consequently does not identify viscosity. Uniaxial symmetry reduces a second-rank conductivity to two components but does not reduce viscosity to those two scalars. Do not impose ordinary isotropic-fluid minor symmetries on crystal-momentum flux without deriving them.

## 5. Continuum and infrared obstruction

Finite-dimensional bounds require a gap uniform in any continuum limit claimed. One positive grid eigenvalue is insufficient. Acoustic collisions can become arbitrarily slow near zero wave number. Decoupled zero-frequency translations must be incorporated among invariants or removed with their vanishing continuum weight justified, never silently inverted.

For exact collisions the diagonal conductivity is a spectral Stieltjes integral,

\[
 \kappa_i(z)=\frac1{k_BT^2}\int_{(0,\infty)}
 \frac{d\mu_i(\lambda)}{\lambda+z},\qquad
 d\mu_i=\langle b_i,d\mathsf E_C(\lambda)b_i\rangle,
 \quad b_i=V_i e.
\]

A flux component in the invariant subspace adds a 1/z term and prevents finite DC conductivity. For m_j=int lambda^{-j} d mu, finite conductivity requires m_1 finite; a finite first memory derivative requires m_2 finite. An O(z^2) remainder is ensured by m_3 finite, using

\[
 \frac1{\lambda+z}=\frac1\lambda-\frac z{\lambda^2}
 +\frac{z^2}{\lambda^2(\lambda+z)}.
\]

A 3D diagonal acoustic illustration with density of states proportional to omega^2, finite low-frequency velocity and rate proportional to omega^p has m_j infrared convergent only for jp<3. At p=2 conductivity is finite but the memory slope diverges. This is a failure regime, not an inferred AlN rate law. Thus the reported finite-grid memory times need their own infrared and q-grid convergence test.

Absence of a global gap does not prove every effective description impossible: observable-weighted integrability can suffice. It removes the uniform exponential-memory estimate used here and requires another theorem and error norm.

## 6. Boundaries and characteristic scales

On a bounded sample,

\[
 \frac12\frac d{dt}\int_\Omega\|y\|^2+\int_\Omega y^*Cy
 =-\frac12\int_{\partial\Omega}y^*V_ny+
 \Re\int_\Omega y^*s.
\]

The outgoing-to-incoming reflection operator should be contractive in the |v_n|-weighted trace norm. Adiabatic reflection additionally requires zero total normal energy flux. Specular reflection, diffuse reflection with self-consistent wall temperature, and substrate transmission define different kinetic boundary problems.

Bulk projection does not determine temperature jumps or drift-slip coefficients. They require a kinetic half-space problem, reflection/transmission data, and the collision operator. Diffuse accommodation can cause viscous wall resistance; specular walls need not. No-slip is an asymptotic boundary result, not a universal phonon condition. Arbitrary fast initial populations also require the explicit term in (A2).

The main bulk scale ratios are |z|/gamma, |k|v_max/gamma, and slow-to-fast collision rates. A channel window schematically requires ell_N << d and d^2/nu << tau_R, with the actual anisotropic viscosity and wall slip. Cross-plane FDTR additionally involves penetration, beam, and interface scales; an in-plane reflecting-slab formula does not supply those conditions.

## 7. Non-identification from diagonal lifetimes

Frequencies and velocities determine equilibrium susceptibilities and streaming. Diagonal rates do not determine K,D, their eigenvectors, the fast gap, or momentum and viscous relaxation. For example,

\[
 C(a,b)=\begin{pmatrix}
 a+b&-a&0&-b\\-a&a+b&-b&0\\
 0&-b&a+b&-a\\-b&0&-a&a+b
 \end{pmatrix},\quad a,b>0,
\]

has invariant e=(1,1,1,1)^T and fixed diagonal d=a+b. The heat vector b_heat=(1,1,-1,-1)^T has eigenvalue 2b. Equal specified energies and velocities proportional to (1,1,-1,-1) therefore give conductivity proportional to 1/b although all diagonal lifetimes are unchanged. Reversal (1 3)(2 4) preserves C and reverses velocities. This is a counterexample within the operator assumptions, not a microscopic AlN scattering construction.

The conserving projected diagonal surrogate in the note also changes its actual diagonal:

\[
 (C_{\rm proj})_{\mu\mu}=r_\mu-r_\mu^2h_\mu^*(H^*DH)^+h_\mu.
\]

Its input r_mu is a bare model relaxation parameter, not generally its resulting collision diagonal. Preserving invariants and positivity does not validate the surrogate against the exact operator.

Next evidence needed: full-grid matrix-free collision action, invariants and slow-subspace residuals, a spectral/coercivity estimate or justified weighted alternative, converged cell problems and memory moments, and the physical boundary law. Normal/resistive actions help test the physical mechanism but do not replace the off-diagonal couplings.

## 8. Prior-art terminology and status

Equivalent search terms: Schur complement; Feshbach projection/partitioning; Mori-Zwanzig memory equation; adiabatic elimination; Chapman-Enskog/Hilbert expansion; kinetic cell problem; collision eigenmodes/relaxons; viscous heat equations; displaced Bose distribution; phonon drift; Peierls-Boltzmann hydrodynamics; Knudsen layer/Milne problem; Stieltjes transport response. No novelty claim is made.

Primary references checked during this branch:

- Mori, Transport, Collective Motion, and Brownian Motion (1965), [doi:10.1143/PTP.33.423](https://academic.oup.com/ptp/article/33/3/423/1925580): projection and memory equations for collective variables.
- Simoncelli, Marzari, Cepellotti, Generalization of Fourier's Law into Viscous Heat Equations (2020), [doi:10.1103/PhysRevX.10.011019](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011019): temperature/drift equations, thermal viscosity, and complementary parity sectors.
- Cepellotti and Marzari, Thermal transport in crystals as a kinetic theory of relaxons (2016), [arXiv:1603.02608](https://arxiv.org/abs/1603.02608): collective eigenvectors of the collision matrix.

DERIVED UNDER STATED ASSUMPTIONS: (A1)-(A8), including exact memory, resolvent error, stable static closure, conditional drift-wave speed, and infrared moment criteria.

FAILED APPROACHES: infer a spectral gap from diagonal lifetimes; infer viscosity from conductivity; assume normal projection makes QC_RE vanish; infer continuum memory convergence from one grid; discard boundary and initial layers without their own checks.

UNRESOLVED FOR AlN: collision action, slow spectrum, uniform gap or weighted alternative, continuum moment convergence, and interface law. No numerical real-AlN hydrodynamic closure or window is certified.

AUDIT TARGETS: Y versus X* sign; direct block inverse check; static Schur positivity; poor slow-coordinate example; diagonal-rate counterexample; and displaced-Bose normalization in (A8). No result here has yet received the campaign's independent proof audit.
