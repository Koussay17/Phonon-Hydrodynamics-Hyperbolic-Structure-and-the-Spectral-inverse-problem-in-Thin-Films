#!/bin/bash
# E32: transfer of the 5x5x3 film-optimised tensors (E31) to 5x5x3 (check), 7x7x5 and 9x9x6; film_event evaluation.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
P3=/c/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
export OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
W=results_c/m553/witness
NAMES="orb0.1_fL30_max orb0.1_fL30_min orb0.1_fL100_max orb0.1_fL100_min orb0.1_fL300_max orb0.1_fL300_min orb0.1_fL1000_max orb0.1_fL1000_min cal1_fL100_max cal1_fL100_min"
args=""; for n in $NAMES; do args="$args --fc3 $n=$W/fc3_$n.npy"; done
for spec in "5 5 3 results_c/m553 events_m553_s0.1.npz pp_r01_REF.npz m553" "7 7 5 results_c/m775x events_m775_s0.1.npz static/pp_r01_REF.npz m775" "9 9 6 results_c/m996 events_m996_s0.1.npz static/pp_r01_REF.npz m996"; do
  set -- $spec; d=$4/filmw_pp; mkdir -p $d
  [ -f $d/pp_r01_REF.npz ] || cp $4/$6 $d/pp_r01_REF.npz
  todo=""; for n in $NAMES; do [ -f $d/pp_r01_$n.npz ] || todo="$todo --fc3 $n=$W/fc3_$n.npy"; done
  [ -z "$todo" ] || $P3 -B scripts_c/p3_events.py --mesh $1 $2 $3 --events $4/$5 --outdir $d $todo >> $d/p3.log 2>&1
  [ -f results_c/film/filmw_$7.json ] || $PY -B scripts_c/film_event.py $4/$5 $d results_c/film/filmw_$7.json $NAMES > results_c/film/filmw_$7.log 2>&1
done
echo done > results_c/film/e32.done
