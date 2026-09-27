# D — pinned linewidth source, dimensions and event counting

**Independent first pass complete, 2026-09-26.** Read only question/assumptions and installed source; no peers, PI artifacts or material calculation. Source root: `D:\ResearchLab\envs\aln-phono3py-4.5.0\Lib\site-packages`. Versions phono3py = phonopy = 4.5.0. This report separates source identities from a conditional kinetic interpretation.

## 1. Source normalization: pp already contains 1/36 and 1/N_q

`reciprocal_to_normal.py` multiplies three eigenvectors without conjugation, applies inverse square-root masses, and divides by sqrt(f0 f1 f2). `interaction.py:176-189,1094` squares this contraction and converts it to **eV²**. With numerical masses in AMU, force constants in eV/angstrom³, frequencies f in THz, and N_q = product(mesh), its conversion is

\[
C_P=\frac{\hbar_{\rm SI}^{3}}
 {36\,8\,(2\pi\,10^{12})^3\,\mathrm{angstrom}^{6}\,m_u^3\,N_q}.
\]

Thus P = C_P |R|²/(f0 f1 f2), where R is the mass-scaled force-constant/eigenvector contraction. The 36 is (3!)²; 8 comes from three oscillator factors 1/2. N_q is the reciprocal mesh count, not the atom count or the force-constant supercell size. **Do not divide exported pp by N_q again.** For N_q = 27, C_P = 0.0001327959698104446.

An independent dimensional route uses oscillator lengths ell_i = sqrt(hbar/(2 m_i omega_i)) and the Hamiltonian coefficient Phi times product(ell_i)/(6 sqrt(N_q)). A scalar test with masses (27,14,27) AMU, frequencies (3,1,2) THz and Phi = 1.23 eV/angstrom³ gives squared coefficients 3.280864568656046e-9 and 3.2808645686560466e-9 eV² by the two routes: relative difference 2.22e-16.

## 2. Ordered gamma sum and cyclic/angular conversion

Define h_TH = h times 1e12 = 0.00413566733 eV/THz and hbar_ps = 0.0006582118985531608 eV ps. `imag_self_energy.py:201-208,768-795` implements

\[
\gamma_0=C_\gamma\sum_t w_t\sum_{jk}P^t_{0jk}
\left[(n_j+n_k+1)D(f_0-f_j-f_k)
 +(n_j-n_k)\{D(f_0+f_j-f_k)-D(f_0-f_j+f_k)\}\right],
\]

\[
C_\gamma=18\pi/h_{\rm TH}^2
 =3.306215697005045\times10^6\;\mathrm{THz}^2/\mathrm{eV}^2.
\]

Here D has units 1/THz. `triplets.py:246-250` forms exactly these Gaussian arguments when broadening is selected; the exact-resonance reference instead uses a spectral delta density. A finite isolated resonant tuple does not supply a finite value of delta(0) or imply irreversible kinetics.

Gamma is a cyclic-frequency half-linewidth in THz. The implemented lifetime convention is `conductivity/utils.py:90`: inverse tau = **4 pi gamma**, numerically in ps^-1. The factors are the amplitude-versus-population factor 2 and the cyclic-to-angular factor 2 pi.

`triplets.py:85,130,392-395` permits daughter-wavevector swapping and counts reciprocal-grid representatives, with sum(w_t) = N_q. Its weights reconstruct the ordered fixed-external-wavevector sum. They are not boson factorial moments or automatically counts of globally unique reversible events. Daughter exchange can be represented partly by these weights and partly by the ordered band loop. An importer cannot blindly add another exchange factor or count a reverse reaction separately.

## 3. Conditional distinct-daughter event correspondence

Assume Born/golden-rule Markov kinetics, diagonal populations with factorized mode distributions, exact resonance, and correctly sewn mode bases. Also assume the relevant cubic coefficients are permutation consistent. A finite Hamiltonian alone does not establish these approximations.

In the convention X_q = a_q + a†_-q, the monomial a_p a†_a a†_b uses all-incoming labels (p,-a,-b). Its coefficient A_pab is the **sum of the six ordered coefficients** for distinct full mode labels. Only when those coefficients agree is A_pab = 6V and |A_pab|² = 36P. Its Fock-state matrix element includes sqrt(N_p(N_a+1)(N_b+1)). Hermitian conjugation supplies the reverse process.

The required reversal/sewing maps are not supplied by pp. At a degenerate wavevector, eigenvectors at -q can differ from the conjugates at q by a unitary matrix. Nor does the source's numerical conversion prove permutation consistency of supplied force constants. Consequently 6V is conditional; absent these checks, coherently summing the relevant ordered coefficients requires amplitude data.

Changing the golden-rule density from energy to ordinary frequency gives

\[
\frac{2\pi}{\hbar_{\rm ps}}|A|^2\delta(E)
 =\frac{|A|^2}{\hbar_{\rm ps}^2}D_f.
\]

For one unique reversible event p <-> a+b with distinct daughters,

\[
\kappa=\frac{36P}{\hbar_{\rm ps}^2}D_f
       =8\pi C_\gamma P D_f,\qquad
F=\kappa[n_p(1+n_a)(1+n_b)-(1+n_p)n_an_b].
\]

The extra factor two relative to 4 pi C_gamma is the two daughter orders in the linewidth sum. At equilibrium, let Q = n_p(1+n_a)(1+n_b), W_i = n_i(1+n_i), and nu = (-1,1,1). The declared scalar linearization is L = kappa Q nu nu^T W^-1, with dot(delta n) = -L delta n. For an isolated distinct-mode event, its three diagonal coefficients agree with the corresponding 4 pi gamma contributions, including the absorption channels. Calling the source's Python Bose-sum kernel on synthetic frequencies (3,1,2) THz at 300 K gives ratios (1.0000000000000002, 0.9999999999999997, 1.0).

## 4. Repeated daughters change both counting and the correspondence

For p <-> a+a, there are three distinct ordered tensor placements, so permutation consistency gives A_paa = 3V. The forward/reverse Fock factors are N_p(N_a+1)(N_a+2) and (N_p+1)N_a(N_a-1). They are not obtained by replacing two independent daughter means by the same mean.

Under the additional **single-mode geometric-moment closure**, the second factorial moments are 2 n_a² and 2(1+n_a)². Hence

\[
\kappa_{aa}=\frac{18P}{\hbar_{\rm ps}^2}D_f
 =4\pi C_\gamma P D_f,\quad
F_{aa}=\kappa_{aa}[n_p(1+n_a)^2-(1+n_p)n_a^2],\quad
\nu=(-1,2).
\]

The one repeated daughter order in the parent linewidth agrees with this parent coefficient. However, the full daughter population derivative includes both occurrences: its diagonal coefficient is 4 kappa_aa (n_a-n_p). The ordered linewidth kernel gives 4 pi gamma_a = 2 kappa_aa (n_a-n_p). Thus the blanket identification L_ii = 4 pi gamma_i fails by a factor two for the repeated daughter in this declared closure; coincident occurrence contributions must remain explicit. This is a model/counting distinction, not a demonstrated defect in phono3py's self-energy calculation.

The direct source-kernel test at (2,1) THz gives event-diagonal / inverse-lifetime ratios **(1,2)** within roundoff. Summing the thermal Fock distribution through occupation 1000 verifies the two factorial moments within 6.67e-16 relative error; omitted probability is 2.85e-70. Geometric closure is an additional ansatz, not an established invariant family of the exact cubic dynamics.

## Evidence and limits

`D-source-factor-check.py` and `.json` persist constants, source SHA-256 hashes, the independent oscillator-length calculation and source-kernel arithmetic. The synthetic checks use a formal common spectral density of 1/THz to compare coefficients; they do not assign a physical finite delta(0). Run the script with the pinned interpreter `D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe`.

Classification: roundoff-limited finite arithmetic; no material, continuum or broadening convergence test. The ordered linewidth is source-defined. Its conversion to a unique reversible event requires verified sewing/permutation conventions, occurrence counting and an explicitly chosen kinetic closure. Broadening must not silently replace those assumptions or exact detailed balance.
