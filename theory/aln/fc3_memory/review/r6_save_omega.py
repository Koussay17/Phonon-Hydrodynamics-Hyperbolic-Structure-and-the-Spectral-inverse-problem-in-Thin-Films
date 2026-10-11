"""Review R6 (part 1, phono3py env): save phono3py's symmetrised reducible collision matrix Omega (5x5x3, sigma 0.1,
300 K, all triplets) and X vectors for given fc3 tensors, for an entry-by-entry comparison with the event operator.
usage: python -B r6_save_omega.py <outdir> name=path ...
"""
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("RAYON_NUM_THREADS", "2")
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402

import r4_p3_lbte as R  # noqa: E402  (installs the capture hooks)
import aln_common as common  # noqa: E402
from p3_events import load_ph3  # noqa: E402

outdir = Path(sys.argv[1])
outdir.mkdir(parents=True, exist_ok=True)
common.verify_inputs()
ph3 = load_ph3(True, [5, 5, 3], 0.1)
fc3_ref = np.array(ph3.fc3)
for spec in sys.argv[2:]:
    name, path = spec.split("=", 1)
    fc3 = fc3_ref if path == "REF" else np.load(path)
    ph3.fc3 = np.ascontiguousarray(fc3, dtype="double")
    ph3.init_phph_interaction()
    R.STASH.clear()
    ph3.run_thermal_conductivity(is_LBTE=True, temperatures=[300.0], is_reducible_collision_matrix=True,
                                 log_level=0)
    Om = R.STASH["Omega"][0, 0]
    n = int(np.sqrt(Om.size))
    bz = ph3.phph_interaction.bz_grid
    np.savez_compressed(outdir / f"omega_{name}.npz", Omega=Om.reshape(n, n), X=R.STASH["X"].reshape(n, 3),
                        kappa=np.array(ph3.thermal_conductivity.kappa)[0, 0], grg2bzg=np.array(bz.grg2bzg),
                        freqs_bz=R.STASH["freqs"])
    print(name, "saved", flush=True)
