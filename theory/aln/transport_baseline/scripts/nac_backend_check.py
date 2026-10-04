"""Which backend's polar (Gonze-Lee) phonons are closer to the cutoff-converged values?

Frequencies on the 9x9x5 BZ grid with lang C / Rust at G_cutoff x1 (default) and x2
(Lambda fixed); reports max |f - f_ref| with f_ref = Rust x2 and C x2.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common
from phono3py.cui.load import load

def freqs(lang, factor):
    kw = dict(unitcell_filename=str(common.DATASET/"POSCAR"), supercell_matrix=[3,3,2], phonon_supercell_matrix=[5,5,3],
              fc3_filename=str(common.DATASET/"fc3.hdf5"), fc2_filename=str(common.DATASET/"fc2.hdf5"),
              born_filename=str(common.DATASET/"BORN"), is_nac=True, symmetrize_fc=False, log_level=0, lang=lang)
    ph = load(**kw)
    if factor != 1.0:
        ph.mesh_numbers = [1,1,1]; ph.init_phph_interaction()
        dm = ph.dynamical_matrix; dm.run(np.array([1.0,0,0])/3)
        _, _, gc, _, lam = dm.Gonze_nac_dataset
        ph = load(**kw); ph.nac_params = dict(ph.nac_params, G_cutoff=gc*factor, Lambda=lam)
    ph.mesh_numbers = [9,9,5]; ph.init_phph_interaction(); ph.run_phonon_solver()
    f, _, _ = ph.get_phonon_data()
    dm = ph.dynamical_matrix; dm.run(np.array([1.0,0,0])/3)
    _, _, gc, gl, lam = dm.Gonze_nac_dataset
    return np.array(f), {"G_cutoff": float(gc), "Lambda": float(lam), "G_count": int(len(gl))}

res = {}
F = {}
for lang in ("C", "Rust"):
    for fac in (1.0, 2.0):
        F[(lang, fac)], res[f"{lang}_x{fac}"] = freqs(lang, fac)
ref = F[("Rust", 2.0)]
for k, v in F.items():
    res[f"maxdiff_{k[0]}_x{k[1]}_vs_Rust_x2"] = float(np.abs(v - ref).max())
    res[f"maxdiff_{k[0]}_x{k[1]}_vs_C_x2"] = float(np.abs(v - F[("C", 2.0)]).max())
print(json.dumps(res, indent=1))
common.write_json(common.CAMPAIGN/"results"/"nac_backend_check_m995.json", res)
