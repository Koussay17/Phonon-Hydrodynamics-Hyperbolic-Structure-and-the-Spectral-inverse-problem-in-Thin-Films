#!/bin/bash
# E36: transfer of the two-mesh (5x5x3 + 7x7x5) film optima (E31) to 7x7x5 (reproduction check), 9x9x6 and 12x12x8.
# Resumable: existing tensors, amplitudes and results are skipped.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
P3=/c/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
export OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
O=results_c/film_opt; W=results_c/m553/witness
NAMES=""
for spec in "orb0.1 orb0.1_z_L30" "orb0.1 orb0.1_z_L100" "orb0.1 orb0.1_z_L300" "cal1nz cal1nz_z_L100"; do
  set -- $spec
  for s in max min; do
    L=${2##*_L}; n="$1_m2L${L}${s}"; NAMES="$NAMES $n"
    [ -f $W/fc3_$n.npy ] || $PY -B scripts_c/witness_verify.py build x results_c/m553 $O/m553m775_$2_${s}_u.npy $n >> $O/transfer2_build.log 2>&1
  done
done
for spec in "7 7 5 results_c/m775x events_m775_s0.1.npz static/pp_r01_REF.npz m775" \
            "9 9 6 results_c/m996 events_m996_s0.1.npz static/pp_r01_REF.npz m996" \
            "12 12 8 review/na_scratch events_m12128_s0.1.npz m12128w/pp_r01_REF.npz m12128"; do
  set -- $spec; d=$4/film2_pp; mkdir -p $d
  [ -f $d/pp_r01_REF.npz ] || cp $4/$6 $d/pp_r01_REF.npz
  todo=""; for n in $NAMES; do [ -f $d/pp_r01_$n.npz ] || todo="$todo --fc3 $n=$W/fc3_$n.npy"; done
  [ -z "$todo" ] || $P3 -B scripts_c/p3_events.py --mesh $1 $2 $3 --events $4/$5 --outdir $d $todo >> $d/p3.log 2>&1
  [ -f results_c/film/film2_$7.json ] || $PY -B scripts_c/film_event.py $4/$5 $d results_c/film/film2_$7.json $NAMES > results_c/film/film2_$7.log 2>&1
  [ -f results_c/film/film2_$7_transfer.json ] || $PY -B scripts_c/mesh_transfer_check.py $4/$5 $d results_c/film/film2_$7_transfer.json $NAMES > results_c/film/film2_$7_transfer.log 2>&1
done
echo done > results_c/film/e36.done
