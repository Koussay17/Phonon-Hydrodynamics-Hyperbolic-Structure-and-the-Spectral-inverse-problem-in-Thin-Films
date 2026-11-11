#!/bin/bash
# E30: film-regularised event-operator memory time on 7x7x5, 9x9x6, 12x12x8 (REF, workflows, witnesses)
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
run() { out=$1; shift; [ -f "$out" ] || $PY -B scripts_c/film_event.py "$@" > "${out%.json}.log" 2>&1; }
run results_c/film/film_m775.json results_c/m775x/events_m775_s0.1.npz results_c/m775x results_c/film/film_m775.json \
  alamode shengbte orb0.1_c0_max orb0.1_c2_max cal1_c0_max cal1_c2_max cal1_c0_min cal1_c2_min cal2_c0_max cal2_c2_max \
  orb0.1_J_c0_max orb0.1_J_c2_max cal1nz_J_c0_max cal1nz_J_c2_max
run results_c/film/film_m996w.json results_c/m996/events_m996_s0.1.npz review/na_scratch/m996w results_c/film/film_m996w.json \
  orb0.1_c0_max orb0.1_c2_max c1_x_cont02 c1_z_cont02
run results_c/film/film_m996.json results_c/m996/events_m996_s0.1.npz results_c/m996 results_c/film/film_m996.json \
  alamode shengbte cal1_c0_max
run results_c/film/film_m996j.json results_c/m996/events_m996_s0.1.npz results_c/m996/joint results_c/film/film_m996j.json \
  orb0.1_J_c0_max orb0.1_J_c2_max cal1nz_J_c0_max cal1nz_J_c2_max
run results_c/film/film_m12128.json review/na_scratch/events_m12128_s0.1.npz review/na_scratch/m12128w results_c/film/film_m12128.json \
  orb0.1_c2_max c1_z_cont02
echo done > results_c/film/e30.done
