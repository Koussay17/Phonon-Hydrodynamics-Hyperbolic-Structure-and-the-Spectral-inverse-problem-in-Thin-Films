#!/bin/bash
# Even-sector operator checks (small meshes), 2 threads, run beside the main queue.
D=/d/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline
PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
export OMP_NUM_THREADS=2 RAYON_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONUNBUFFERED=1
cd $D
for spec in "5 5 3:tetra:" "7 7 5:tetra:" "9 9 5:tetra:" "5 5 3:s0.1:--sigma 0.1" "7 7 5:s0.1:--sigma 0.1" "9 9 5:s0.1:--sigma 0.1" "9 9 5:s0.2:--sigma 0.2"; do
  m=${spec%%:*}; rest=${spec#*:}; tag=${rest%%:*}; extra=${rest#*:}
  name="physop_m${m// /}_${tag}"
  [ -f results/$name.json ] && continue
  $PY -X utf8 -B scripts/physical_operator_check.py --mesh $m $extra --json results/$name.json > runs/logs/$name.log 2>&1
  echo "$(date -Iseconds) $name rc=$?" >> runs/logs/queue_physop.log
done
