#!/bin/bash
# E39: independent multistart of the three-mesh film optima (E33), per-cluster class orb 0.1, 100 nm, max and min.
# Starts rand:1 and rand:2 (uniform in the box, robust_mesh.py --start rand:SEED) against the zero/warm starts of E33.
# Each optimum is then checked on 7x7x5, 9x9x6 and the held-out 12x12x8 exactly as in E37 (m3_validate2.sh).
# Names: m3_orb0.1rs<SEED>_z_L100_<sense> (tag prefix "orb" selects the scaled coordinates in witness_verify.py).
# One run at a time (the three-mesh model holds about 10-15 GB). Resumable: finished runs and checks are skipped.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
P3=/c/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
O=results_c/film_opt; W=results_c/m553/witness; R=results_c/film3
validate() {  # name (e.g. orb0.1rs1_m3L100max) u.npy
  local n=$1 u=$2
  [ -f $R/$n.ok ] && return
  export OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
  [ -f $W/fc3_$n.npy ] || $PY -B scripts_c/witness_verify.py build x results_c/m553 $u $n >> $R/build.log 2>&1
  for spec in "7 7 5 results_c/m775x events_m775_s0.1.npz m775" "9 9 6 results_c/m996 events_m996_s0.1.npz m996" \
              "12 12 8 review/na_scratch events_m12128_s0.1.npz m12128"; do
    set -- $spec; d=$4/film2_pp
    [ -f $d/pp_r01_$n.npz ] || $P3 -B scripts_c/p3_events.py --mesh $1 $2 $3 --events $4/$5 --outdir $d --fc3 $n=$W/fc3_$n.npy >> $d/p3.log 2>&1
    [ -f $R/${n}_$6.json ] || $PY -B scripts_c/film_event.py $4/$5 $d $R/${n}_$6.json $n > $R/${n}_$6.log 2>&1
    [ -f $R/${n}_$6_transfer.json ] || $PY -B scripts_c/mesh_transfer_check.py $4/$5 $d $R/${n}_$6_transfer.json $n > $R/${n}_$6_transfer.log 2>&1
  done
  [ -f $R/${n}_m12128.json ] && [ -f $R/${n}_m12128_transfer.json ] && echo done > $R/$n.ok
}
for seed in 1 2; do
  for s in max min; do
    b=m3_orb0.1rs${seed}_z_L100_$s
    if [ ! -f $O/$b.json ]; then
      export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
      $PY -B scripts_c/robust_mesh.py --tag orb --eps 0.1 --col 2 --sense $s --meshes m553,m775,m996 \
          --objective m553,m775,m996 --film 100 --start rand:$seed --out $O/$b > $O/$b.out 2>&1
    fi
    [ -f $O/${b}_u.npy ] && validate orb0.1rs${seed}_m3L100$s $O/${b}_u.npy
  done
done
echo done > $R/e39.ok
