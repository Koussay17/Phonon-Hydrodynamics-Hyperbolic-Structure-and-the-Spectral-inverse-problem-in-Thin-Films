#!/bin/bash
# E33: three-mesh robust film runs (5x5x3 + 7x7x5 + 9x9x6 data and objective), one at a time (9x9x6 model about 10 GB).
# Resumable: finished runs skipped.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=results_c/film_opt
one() {  # name tag eps col sense L
  [ -f $O/$1.json ] && return
  $PY -B scripts_c/robust_mesh.py --tag $2 --eps $3 --col $4 --sense $5 --meshes m553,m775,m996 --objective m553,m775,m996 --film $6 --out $O/$1 > $O/$1.out 2>&1
}
for L in 100 30; do for s in max min; do one m3_orb0.1_z_L${L}_$s orb 0.1 2 $s $L; done; done
for t in cal1nz cal3wnz cal3w; do for s in max min; do one m3_${t}_z_L100_$s $t 1 2 $s 100; done; done
for s in max min; do one m3_orb0.1_z_L300_$s orb 0.1 2 $s 300; done
echo done > $O/film_opt3.done
