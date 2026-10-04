# 00 — Question (workstream B: AlN transport baseline)

## Objective
Build and validate a first-principles lattice thermal-conductivity pipeline for wurtzite AlN
(phono3py/phonopy 4.5.0) against the Phonon Olympics benchmark before any new physics, and
export the data needed by the downstream workstreams (collision operator, film BTE).

## Claims to test
1. With the phono3py team's settings, phono3py 4.5.0 reproduces their 300 K values
   (RTA 253/232, LBTE 285/271 W/(m K), in-plane/cross-plane; J. Appl. Phys. 138, 135108 (2025)).
2. The 300 K conductivity is converged with respect to the q mesh and the Brillouin-zone
   integration (tetrahedron vs Gaussian), for RTA and LBTE; the residual error is quantified.
3. kappa(T), 100–1000 K, agrees with the Olympics phono3py table and is placed relative to the
   ShengBTE table.
4. Natural-abundance isotope scattering changes kappa(300 K) only weakly; size quantified.
5. The exported collision matrix is documented exactly (representation, weights, averaging,
   symmetrization, units); its symmetry and energy (even-sector) null-vector behaviour are
   measured, given the repository finding that phono3py's matrix is conductivity-equivalent
   (odd sector) but not the physical operator on even populations.
6. At least one number is reproduced by an independent route.

## Fixed inputs
Phonon Olympics author repository, commit 0640f07735059be9717a7565c2a0f22dc0da7a17,
`Aluminum Nitride/phono3py/AlN_kappa_input_files`: POSCAR, BORN, fc2.hdf5 (5x5x3 supercell),
fc3.hdf5 (3x3x2 supercell). Local copy:
`D:\ResearchLab\datasets\phonon-olympics\0640f077...\AlN_phono3py` (SHA-256 checked by every run).

## Phonon Olympics phono3py-team settings (from `AlN summary final.xlsx` and the kappa sheets)
phonopy 2.12.0 / phono3py 2.1.0; q mesh 31x31x17 (Gamma-centred); linear tetrahedron method;
isotopically pure; NAC by the Gonze method (BORN: eps = diag(4.4805, 4.4805, 4.7120),
Z*(Al) = diag(2.5139, 2.5139, 2.6781)); direct LBTE solution (Chaput); group velocities from the
dynamical-matrix derivative; supercells 5x5x3 (fc2) and 3x3x2 (fc3, all triplets, translational
invariance imposed within phono3py); 300 K values 252.992/231.947 (RTA), 285.048/271.262 (LBTE).
