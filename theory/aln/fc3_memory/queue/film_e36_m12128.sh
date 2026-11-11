#!/bin/bash
# E36: 12x12x8 film metrics of the two-mesh optima, rerun with the low-memory film_event.py after the first attempt ran
# out of memory (three-mesh optimisation resident). Starts when film_transfer2.sh has finished.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
while [ ! -f results_c/film/e36.done ]; do sleep 120; done
NAMES="orb0.1_m2L30max orb0.1_m2L30min orb0.1_m2L100max orb0.1_m2L100min orb0.1_m2L300max orb0.1_m2L300min cal1nz_m2L100max cal1nz_m2L100min"
[ -f results_c/film/film2_m12128.json ] || $PY -B scripts_c/film_event.py review/na_scratch/events_m12128_s0.1.npz review/na_scratch/film2_pp results_c/film/film2_m12128.json $NAMES > results_c/film/film2_m12128.log 2>&1
echo done > results_c/film/e36_m12128.done
