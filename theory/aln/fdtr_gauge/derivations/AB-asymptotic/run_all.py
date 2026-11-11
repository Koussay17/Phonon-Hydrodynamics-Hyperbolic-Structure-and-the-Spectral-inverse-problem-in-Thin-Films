"""Rerun every AB-A check in order; log exit codes, runtimes and output to results/run_all.log."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["b2_lowfreq_series.py", "b3_highfreq_reflections.py", "b4_rank_strata.py",
           "b5_xi1_sensitivity_and_bs_inversion.py", "a1_phase_sensitivity.py", "a2_nuisance_collinearity.py",
           "a3_delay_variance_and_lumped.py", "x1_twin_crosscheck.py", "x2_beam_vs_1d.py"]
env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2", PYTHONWARNINGS="ignore")
log = []
for s in SCRIPTS:
    t0 = time.time()
    r = subprocess.run([sys.executable, s], cwd=HERE, env=env, capture_output=True, text=True)
    dt = time.time() - t0
    log.append("=== %s  exit %d  %.1f s\n%s" % (s, r.returncode, dt, r.stdout[-6000:]))
    print("%-42s exit %d  %.1f s" % (s, r.returncode, dt))
    if r.returncode:
        print(r.stderr[-3000:])
with open(os.path.join(HERE, "results", "run_all.log"), "w") as fh:
    fh.write("\n".join(log))
