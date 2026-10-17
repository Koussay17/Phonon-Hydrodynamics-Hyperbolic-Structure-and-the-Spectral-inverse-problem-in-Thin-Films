"""Write phono3py displacement dataset and pw.x inputs for one functional.
Pair cutoff 3.9 A: keeps the Al-N shell at 3.65 A (PBEsol) / 3.67 A (PBE) for both functionals; next shell 4.38 A.
usage: gen_disp.py <relax.out or reference> <functional> <outdir> <ecutwfc> <ecutrho>"""
import sys, re, os, numpy as np
from phonopy.structure.atoms import PhonopyAtoms
from phono3py import Phono3py
BOHR = 0.529177210903
src, func, out, ecw, ecr = sys.argv[1:6]
if src == "reference":  # PBEsol structure of the reference ShengBTE workflow (3.harmonic__input__scf.in)
    lat = np.array([[5.888098425987400, 0, 0], [-2.944049212993700, 5.099242816888617, 0],
                    [0, 0, 9.424851438595599]]) * BOHR
    frac = np.array([[1/3, 2/3, 0.9995896991], [2/3, 1/3, 0.4995896991],
                     [1/3, 2/3, 0.3815261244], [2/3, 1/3, 0.8815261244]])
else:  # final structure of a vc-relax
    t = open(src).read()
    assert "End final coordinates" in t, "relaxation not finished"
    blk = t[t.rindex("Begin final coordinates"):t.rindex("End final coordinates")]
    alat = float(re.search(r"CELL_PARAMETERS \(alat=\s*([\d.]+)\)", blk).group(1))
    rows = blk.split("CELL_PARAMETERS")[1].splitlines()[1:4]
    lat = np.array([[float(x) for x in r.split()] for r in rows]) * alat * BOHR
    pos = blk.split("ATOMIC_POSITIONS (crystal)")[1].splitlines()[1:5]
    frac = np.array([[float(x) for x in r.split()[1:4]] for r in pos])
cell = PhonopyAtoms(symbols=["Al", "Al", "N", "N"], cell=lat, scaled_positions=frac)
p = Phono3py(cell, supercell_matrix=np.diag([3, 3, 2]), primitive_matrix="auto")
p.generate_displacements(distance=0.03, cutoff_pair_distance=3.9)
os.makedirs(out, exist_ok=True)
p.save(os.path.join(out, "phono3py_disp.yaml"))
head = open(sys.argv[0].replace("gen_disp.py", "scf_head.in")).read().replace("@FUNC@", func).replace("@ECUTW@", ecw).replace("@ECUTR@", ecr)
def write(name, sc):
    s = head.replace("@NAT@", str(len(sc)))
    s += "CELL_PARAMETERS angstrom\n" + "".join("  %.10f %.10f %.10f\n" % tuple(v) for v in sc.cell)
    s += "ATOMIC_POSITIONS angstrom\n" + "".join("%s %.10f %.10f %.10f\n" % (e, *x) for e, x in zip(sc.symbols, sc.positions))
    s += "K_POINTS gamma\n"
    open(os.path.join(out, name), "w").write(s)
write("disp-00000.in", p.supercell)
n = 0
for i, sc in enumerate(p.supercells_with_displacements):
    if sc is not None:
        write("disp-%05d.in" % (i + 1), sc); n += 1
np.savetxt(os.path.join(out, "lattice_angstrom.txt"), lat)
print(func, "supercells:", n, "+ perfect")
