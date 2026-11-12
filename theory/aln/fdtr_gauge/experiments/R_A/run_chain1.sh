#!/bin/bash
# R-A chain 1 (sequential, one heavy process at a time): full 7x7x5 evaluation, light RTA routes, widths.
cd /d/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/experiments/R_A
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
$PY -u -B ra1_film_event_ext.py m775 > out/ra1_m775.log 2>&1; echo "ra1 m775 EXIT $?" >> out/chain1.log
OMP_NUM_THREADS=1 $PY -u -B ra1c_ref_rta.py > out/ra1c.log 2>&1; echo "ra1c EXIT $?" >> out/chain1.log
OMP_NUM_THREADS=1 $PY -u -B ra2_rta_film_c.py > out/ra2.log 2>&1; echo "ra2 EXIT $?" >> out/chain1.log
OMP_NUM_THREADS=1 $PY -u -B ra2b_rta_taul_meshes.py > out/ra2b.log 2>&1; echo "ra2b EXIT $?" >> out/chain1.log
$PY -u -B ra3_widths.py > out/ra3.log 2>&1; echo "ra3 EXIT $?" >> out/chain1.log
echo CHAIN1 DONE >> out/chain1.log
