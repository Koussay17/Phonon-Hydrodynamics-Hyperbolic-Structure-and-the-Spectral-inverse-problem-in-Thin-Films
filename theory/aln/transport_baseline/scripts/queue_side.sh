#!/bin/bash
# Second queue: mesh series (reproduction + production), Gaussian smearing, larger mesh.
# Starts when queue_main has finished (avoids CPU contention).
set -u
D=/d/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline
PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
export OMP_NUM_THREADS=6 RAYON_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 PYTHONUNBUFFERED=1
cd "$D"
QLOG=runs/logs/queue_side.log
until grep -q "queue_main finished" runs/logs/queue_main.log 2>/dev/null; do sleep 60; done
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
# reproduction-by-mesh (phono3py-2.1.0-like settings), LBTE + RTA via rows+PCG
for m in "19 19 11" "23 23 13" "27 27 15"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-v2C $OOC all --out runs/ooc-m${tag}-v2C --mesh $m --temps 300 --lang C --no-r0-average \
       --gv-delta-q 1e-5 --block 16 --rows-per-chunk 1024
done
# production mesh series (tetrahedron), 300 K
for m in "11 11 7" "15 15 9" "19 19 11" "23 23 13" "27 27 15"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-prod $OOC all --out runs/ooc-m${tag}-prod --mesh $m --temps 300 --lang Rust \
       --block 16 --rows-per-chunk 1024
done
# Gaussian smearing (production), RTA, three widths per run (pp reused across widths)
for m in "15 15 9" "19 19 11" "23 23 13"; do
  tag=$(echo $m | tr -d ' ')
  step prod-rta-m${tag}-gauss $RK prod-rta-m${tag}-gauss 6 --mesh $m --method rta --temps 300 \
       --sigmas 0.05 0.1 0.2 --note "production setting, Gaussian smearing widths in THz"
done
# Gaussian LBTE (production) at 19x19x11 for two widths
for s in 0.1 0.2; do
  step ooc-m191911-prod-s$s $OOC all --out runs/ooc-m191911-prod-s$s --mesh 19 19 11 --temps 300 \
       --lang Rust --sigma $s --block 16 --rows-per-chunk 1024
done
# beyond the Olympics mesh: production RTA and LBTE at 35x35x19
step prod-rta-m353519 $RK prod-rta-m353519 6 --mesh 35 35 19 --method rta --temps 300 \
     --note "production setting, mesh beyond the Olympics choice"
step ooc-m353519-prod $OOC all --out runs/ooc-m353519-prod --mesh 35 35 19 --temps 300 --lang Rust \
     --block 8 --rows-per-chunk 512
echo "$(date -Iseconds) queue_side finished" >> $QLOG
