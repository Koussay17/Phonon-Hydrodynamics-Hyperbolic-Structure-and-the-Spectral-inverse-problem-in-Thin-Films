#!/bin/sh
# Master queue for blocks A and B (resumable: every step skips work already on disk). Two lanes:
#   lane L (linearisations, 6 single-thread workers) and lane R (chunked RTA mesh sequence, 6 threads).
# Start:  sh /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory/queue/run_all.sh
C=/d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
S=$C/scripts_c
R=D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c
PY=C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
P3=C:/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
LOG=$C/queue/logs
mkdir -p $LOG
export ALN_DATASET='D:\ResearchLab\datasets\phonon-olympics\0640f07735059be9717a7565c2a0f22dc0da7a17\AlN_phono3py'
cd $S

lane_L() {
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  echo "$(date) L: 9x9x6 parts" >> $LOG/run_all.log
  seq 0 15 | xargs -P 6 -I{} sh -c "[ -f $R/m996/lin_G_{}.npy ] || $PY -B linearize_parts.py part $R/m996/events_m996_s0.1.npz $R/m996 {} 16 > $R/m996/lin_part_{}.log 2>&1"
  if [ ! -f $R/m996/merge.done ]; then
    echo "$(date) L: 9x9x6 merge" >> $LOG/run_all.log
    OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 LIN_DTYPE=float32 $PY -B linearize_parts.py merge $R/m996 16 $R/m553 orb cal1nz cal3w cal3wnz > $LOG/merge_m996.log 2>&1 && echo done > $R/m996/merge.done
  fi
  echo "$(date) L: 7x7x5 parts" >> $LOG/run_all.log
  seq 0 5 | xargs -P 6 -I{} sh -c "[ -f $R/m775x/lin_G_{}.npy ] || $PY -B linearize_parts.py part $R/m775x/events_m775_s0.1.npz $R/m775x {} 6 > $R/m775x/lin_part_{}.log 2>&1"
  if [ ! -f $R/m775x/merge3w.done ]; then
    echo "$(date) L: 7x7x5 merge" >> $LOG/run_all.log
    OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 $PY -B linearize_parts.py merge $R/m775x 6 $R/m553 cal3w cal3wnz > $LOG/merge_m775.log 2>&1 && echo done > $R/m775x/merge3w.done
  fi
  echo "$(date) L: 5x5x3 three-workflow classes are built by orbit_scaling.py (already on disk)" >> $LOG/run_all.log
  if [ ! -f $R/m996/temperature_witnesses.json ]; then
    echo "$(date) L: kappa(T) 9x9x6" >> $LOG/run_all.log
    W=$C/review/na_scratch/m996w
    OMP_NUM_THREADS=4 $PY -B $C/review/r10_temperature.py $R/m996/events_m996_s0.1.npz $R/m996/pp_r01_REF.npz $R/m996/temperature_witnesses.json 100,200,300,500,800 cal1nz_J_x=$R/m996/joint/pp_r01_cal1nz_J_c0_max.npz cal1nz_J_z=$R/m996/joint/pp_r01_cal1nz_J_c2_max.npz orb_J_x=$R/m996/joint/pp_r01_orb0.1_J_c0_max.npz orb_J_z=$R/m996/joint/pp_r01_orb0.1_J_c2_max.npz orb_E15_z=$W/pp_r01_orb0.1_c2_max.npz c1_z_cont02=$W/pp_r01_c1_z_cont02.npz alamode=$R/m996/pp_r01_alamode.npz shengbte=$R/m996/pp_r01_shengbte.npz > $LOG/temperature_m996.log 2>&1
  fi
  echo "$(date) L: done" >> $LOG/run_all.log
}

lane_R() {
  for m in "48 48 26" "56 56 30" "64 64 35"; do
    d=$R/rta/m$(echo $m | tr -d ' ')
    [ -f $d/assembled.json ] && continue
    echo "$(date) R: RTA $m" >> $LOG/run_all.log
    $P3 -B rta_chunks.py run --mesh $m --outdir $d --chunk 64 --threads 6 >> $LOG/rta_$(echo $m | tr -d ' ').log 2>&1 && \
      $PY -B rta_chunks.py assemble --outdir $d > $d/assemble.log 2>&1
  done
  echo "$(date) R: done" >> $LOG/run_all.log
}

echo "$(date) run_all start (pid $$)" >> $LOG/run_all.log
lane_L &
lane_R &
wait
echo "$(date) run_all finished" >> $LOG/run_all.log
