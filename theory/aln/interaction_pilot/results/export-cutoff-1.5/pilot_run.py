"""One-grid-point interaction-export pilot, not converged conductivity."""
from pathlib import Path
import os, sys, json, platform, hashlib, time, importlib.metadata
import numpy as np
import h5py
import phono3py
from phono3py.file_IO import write_phonon_to_hdf5

ROOT = Path(__file__).resolve().parent
RUN = ROOT/"runs"/"gp1-m333"
os.chdir(RUN)
started=time.time()
assert phono3py.__version__ == "4.5.0"
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
input_record={}
for name in ("POSCAR","BORN","fc2.hdf5","fc3.hdf5"):
    input_record[name]={"sha256":digest(name),"bytes":Path(name).stat().st_size}
for name in ("fc2.hdf5","fc3.hdf5"):
    with h5py.File(name) as f:
        print(name,{k:list(f[k].shape) for k in f},flush=True)
print("Loading pinned inputs without force-constant symmetrization",flush=True)
ph=phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
    phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",
    fc3_filename="fc3.hdf5",born_filename="BORN",is_nac=True,
    symmetrize_fc=False,log_level=1,lang="Rust")
ph.nac_params=dict(ph.nac_params, G_cutoff=1.7947217950209668, Lambda=0.26615773936469767)
ph.mesh_numbers=[3,3,3]
ph.init_phph_interaction()
ph.run_phonon_solver()
freq,eig,address=ph.get_phonon_data()
write_phonon_to_hdf5(freq,eig,address,ph.mesh_numbers,bz_grid=ph.grid)
print("Phonons",freq.shape,"min",freq.min(),"max",freq.max(),flush=True)
ph.sigmas=[0.1]
ph.run_thermal_conductivity(temperatures=[300],grid_points=[1],
    is_full_pp=True,write_pp=True,write_gamma=True,write_gamma_detail=True,
    is_N_U=True,write_kappa=False,log_level=1)
record={"scope":"single-point coarse-mesh interaction/self-energy export; no conductivity or hydrodynamic conclusion",
    "python":platform.python_version(),"packages":{k:importlib.metadata.version(k) for k in
    ["phono3py","phonopy","phonors","numpy","scipy","h5py","spglib"]},
    "inputs":input_record,"settings":{"mesh":[3,3,3],"grid_points":[1],"temperature_K":300,
    "sigma_THz":0.1,"backend":"Rust","symmetrize_fc":False,"NAC":True},
    "frequency_THz":{"shape":list(freq.shape),"min":float(freq.min()),"max":float(freq.max())},
    "elapsed_seconds":time.time()-started,"script_sha256":digest(__file__)}
record["outputs"]={p.name:{"bytes":p.stat().st_size,"sha256":digest(p)}
    for p in RUN.glob("*.hdf5") if p.name not in ("fc2.hdf5","fc3.hdf5")}
(ROOT/"pilot-run.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
print(json.dumps(record,indent=2),flush=True)
