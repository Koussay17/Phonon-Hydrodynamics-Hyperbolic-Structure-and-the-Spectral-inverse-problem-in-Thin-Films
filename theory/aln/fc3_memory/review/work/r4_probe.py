import os, sys, time
os.environ["OMP_NUM_THREADS"]="2"; os.environ["RAYON_NUM_THREADS"]="2"
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261005-fc3-memory\review")
import numpy as np
import r4_p3_lbte as R
from p3_events import load_ph3
ph3 = load_ph3(True, [5,5,3], 0.1)
t0=time.time()
ph3.init_phph_interaction()
ph3.run_thermal_conductivity(is_LBTE=True, temperatures=[300.0], is_reducible_collision_matrix=True, log_level=0)
tc = ph3.thermal_conductivity
print("time", time.time()-t0)
print("kappa", np.array(tc.kappa)[0,0])
print("gamma shape", np.array(tc.gamma).shape, "grid_points", np.array(tc.grid_points).shape, np.array(tc.grid_points)[:10])
print("Omega", R.STASH["Omega"].shape, "X", R.STASH["X"].shape, "freqs", R.STASH["freqs"].shape)
itr = ph3.phph_interaction; bz = itr.bz_grid
print("grg2bzg", np.array(bz.grg2bzg)[:10], "sigma_cutoff", ph3.sigma_cutoff)
print([a for a in dir(tc) if not a.startswith('__')][:80])
