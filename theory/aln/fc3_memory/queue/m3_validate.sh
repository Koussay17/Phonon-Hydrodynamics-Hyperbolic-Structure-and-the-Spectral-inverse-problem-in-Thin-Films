#!/bin/bash
# E37: held-out 12x12x8 validation of the three-mesh film optima (E33), each as soon as its run has written _u.npy.
# Per tensor: witness build -> phono3py amplitudes on the 12x12x8 events -> film_event (low-memory) -> mesh_transfer_check.
# Also evaluated on the fitted meshes 7x7x5 and 9x9x6 for the reproduction check. Resumable; ends after film_opt3.done.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
P3=/c/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
export OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
O=results_c/film_opt; W=results_c/m553/witness; R=results_c/film3
mkdir -p $R
one() {  # u.npy
  b=$(basename $1 _u.npy)                    # m3_orb0.1_z_L100_max
  t=${b#m3_}; t=${t%%_z_*}; L=${b##*_L}; L=${L%_*}; s=${b##*_}
  n="${t}_m3L${L}${s}"                       # orb0.1_m3L100max
  [ -f $R/$n.done ] && return
  [ -f $W/fc3_$n.npy ] || $PY -B scripts_c/witness_verify.py build x results_c/m553 $1 $n >> $R/build.log 2>&1
  for spec in "7 7 5 results_c/m775x events_m775_s0.1.npz m775" "9 9 6 results_c/m996 events_m996_s0.1.npz m996" \
              "12 12 8 review/na_scratch events_m12128_s0.1.npz m12128"; do
    set -- $spec; d=$4/film2_pp              # REF amplitudes already there (E36)
    [ -f $d/pp_r01_$n.npz ] || $P3 -B scripts_c/p3_events.py --mesh $1 $2 $3 --events $4/$5 --outdir $d --fc3 $n=$W/fc3_$n.npy >> $d/p3.log 2>&1
    [ -f $R/${n}_$7.json ] || $PY -B scripts_c/film_event.py $4/$5 $d $R/${n}_$7.json $n > $R/${n}_$7.log 2>&1
    [ -f $R/${n}_$7_transfer.json ] || $PY -B scripts_c/mesh_transfer_check.py $4/$5 $d $R/${n}_$7_transfer.json $n > $R/${n}_$7_transfer.log 2>&1
  done
  [ -f $R/${n}_m12128.json ] && [ -f $R/${n}_m12128_transfer.json ] && echo done > $R/$n.done
}
while :; do
  for u in $O/m3_*_u.npy; do [ -f "$u" ] && one $u; done
  [ -f $O/film_opt3.done ] && { for u in $O/m3_*_u.npy; do one $u; done; break; }
  sleep 600
done
echo done > $R/e37.done
