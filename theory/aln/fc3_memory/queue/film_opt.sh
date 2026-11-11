#!/bin/bash
# E31: attainable range of the film memory time tau_z(L) (robust_mesh.py --film). Resumable: finished runs skipped.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
export PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
O=results_c/film_opt
one() {  # name tag eps col sense meshes objective L
  [ -f $O/$1.json ] && return
  $PY -B scripts_c/robust_mesh.py --tag $2 --eps $3 --col $4 --sense $5 --meshes $6 --objective $7 --film $8 --out $O/$1 > $O/$1.out 2>&1
}
export -f one; export O
laneA() {  # 5x5x3 single-mesh runs, 3 at a time
  {
  for L in 30 300 1000; do for s in max min; do echo "m553_orb0.1_z_L${L}_$s orb 0.1 2 $s m553 m553 $L"; done; done
  for t in cal1 cal3w cal3wnz; do for s in max min; do echo "m553_${t}_z_L100_$s $t 1 2 $s m553 m553 100"; done; done
  for s in max min; do echo "m553_orb0.2_z_L100_$s orb 0.2 2 $s m553 m553 100"; done
  } | xargs -P 3 -L 1 bash -c 'one "$@"' _
}
laneB() {  # two-mesh robust runs (5x5x3 + 7x7x5 data and objective), one at a time (about 4 GB each)
  for L in 100 30 300; do for s in max min; do one m553m775_orb0.1_z_L${L}_$s orb 0.1 2 $s m553,m775 m553,m775 $L; done; done
  for s in max min; do one m553m775_cal1nz_z_L100_$s cal1nz 1 2 $s m553,m775 m553,m775 100; done
}
laneA & laneB & wait
echo done > $O/film_opt.done
