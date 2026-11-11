#!/bin/sh
# Block A job queue (independent of any interactive session). Sequential; each job logs to queue/logs.
# 1) wait for the running 48x48x26 RTA run; 2) RTA 56x56x30 (tetrahedron); 3) RTA 64x64x35.
set -u
Q=/d/ResearchLab/orchestration/campaigns/20261005-fc3-memory/queue
mkdir -p $Q/logs
export ALN_CAMPAIGN='D:\ResearchLab\orchestration\campaigns\20261005-fc3-memory\results\rta'
export OMP_NUM_THREADS=12 RAYON_NUM_THREADS=12
P=C:/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
RK=D:/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline/scripts/run_kappa.py
RUNS=/d/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results/rta/runs
echo "$(date) queue A start" >> $Q/logs/queue_A.log
until grep -q '"status": "finished"' $RUNS/rta-m484826-tetra/run.json 2>/dev/null; do sleep 120; done
echo "$(date) 48x48x26 finished" >> $Q/logs/queue_A.log
for m in "56 56 30" "64 64 35"; do
  id="rta-m$(echo $m | tr -d ' ')-tetra"
  [ -f $RUNS/$id/run.json ] && grep -q '"status": "finished"' $RUNS/$id/run.json && continue
  echo "$(date) start $id" >> $Q/logs/queue_A.log
  $P -B $RK --run-id $id --mesh $m --method rta --temps 300 --note "E27: tau_mem(RTA) mesh sequence" > $Q/logs/$id.log 2>&1
  echo "$(date) end $id rc=$?" >> $Q/logs/queue_A.log
done
echo "$(date) queue A done" >> $Q/logs/queue_A.log
