#!/bin/bash
# Parallel low-thread queue (3 threads) run beside queue_main: reproduction mesh series and
# small production meshes. Uses the same done/running markers as queue_side.
set -u
D=/d/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline
PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
export OMP_NUM_THREADS=3 RAYON_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 PYTHONUNBUFFERED=1
cd "$D"
QLOG=runs/logs/queue_par.log
step() {
  local name="$1"; shift
  if [ -f "runs/logs/$name.done" ] || [ -f "runs/logs/$name.running" ]; then echo "$(date -Iseconds) skip $name" >> $QLOG; return; fi
  touch "runs/logs/$name.running"
  echo "$(date -Iseconds) start $name" >> $QLOG
  "$@" > "runs/logs/$name.log" 2>&1
  local rc=$?
  echo "$(date -Iseconds) end $name rc=$rc" >> $QLOG
  [ $rc -eq 0 ] && touch "runs/logs/$name.done"; rm -f "runs/logs/$name.running"
}
OOC="$PY -X utf8 -B scripts/lbte_ooc.py"
RK="bash scripts/launch.sh"
for m in "19 19 11" "23 23 13" "27 27 15"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-v2C $OOC all --out runs/ooc-m${tag}-v2C --mesh $m --temps 300 --lang C --no-r0-average \
       --gv-delta-q 1e-5 --block 16 --rows-per-chunk 1024
done
for m in "11 11 7" "15 15 9" "19 19 11" "23 23 13"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-prod $OOC all --out runs/ooc-m${tag}-prod --mesh $m --temps 300 --lang Rust \
       --block 16 --rows-per-chunk 1024
done
step prod-rta-m15159-gauss $RK prod-rta-m15159-gauss 3 --mesh 15 15 9 --method rta --temps 300 \
     --sigmas 0.05 0.1 0.2 --note "production setting, Gaussian smearing widths in THz"
echo "$(date -Iseconds) queue_par finished" >> $QLOG
