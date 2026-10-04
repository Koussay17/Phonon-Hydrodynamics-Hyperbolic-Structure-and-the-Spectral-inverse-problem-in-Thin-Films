#!/bin/bash
# Final sequential queue (6 threads). Launched detached (PowerShell Start-Process) after the
# first queues were terminated by memory-pressure process reaping at ~12:40 on 3 Oct.
# Steps already marked .done are skipped; out-of-core runs resume from progress.json.
set -u
D=/d/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline
PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
export OMP_NUM_THREADS=6 RAYON_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 PYTHONUNBUFFERED=1
cd "$D"
QLOG=runs/logs/queue_final.log
step() {
  local name="$1"; shift
  if [ -f "runs/logs/$name.done" ]; then echo "$(date -Iseconds) skip $name" >> $QLOG; return; fi
  echo "$(date -Iseconds) start $name" >> $QLOG
  "$@" >> "runs/logs/$name.log" 2>&1
  local rc=$?
  echo "$(date -Iseconds) end $name rc=$rc" >> $QLOG
  [ $rc -eq 0 ] && touch "runs/logs/$name.done"
}
OOC="$PY -X utf8 -B scripts/lbte_ooc.py"
RK="bash scripts/launch.sh"
ALLT="100 125 150 175 200 225 250 275 300 400 500 600 700 800 900 1000"
step ooc-m313117-v2C  $OOC all --out runs/ooc-m313117-v2C --mesh 31 31 17 --temps 300 --lang C \
     --no-r0-average --gv-delta-q 1e-5 --block 16 --rows-per-chunk 1024
step ooc-m313117-prod $OOC all --out runs/ooc-m313117-prod --mesh 31 31 17 --temps 100 200 300 500 1000 \
     --lang Rust --block 16 --rows-per-chunk 1024 --keep-raw 2
step prod-rta-m313117-iso $RK prod-rta-m313117-iso 6 --mesh 31 31 17 --method rta --temps $ALLT --isotope \
     --note "production setting, natural-abundance isotope scattering"
step ooc-m313117-prod-iso $OOC solve --out runs/ooc-m313117-prod --temp-index 2 --rows-per-chunk 1024 \
     --gamma-iso-from runs/prod-rta-m313117-iso/kappa-m313117.hdf5
for m in "11 11 7" "15 15 9" "19 19 11" "23 23 13" "27 27 15"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-prod $OOC all --out runs/ooc-m${tag}-prod --mesh $m --temps 300 --lang Rust --block 16 --rows-per-chunk 1024
done
for m in "23 23 13" "27 27 15"; do
  tag=$(echo $m | tr -d ' ')
  step ooc-m${tag}-v2C $OOC all --out runs/ooc-m${tag}-v2C --mesh $m --temps 300 --lang C --no-r0-average \
       --gv-delta-q 1e-5 --block 16 --rows-per-chunk 1024
done
for m in "15 15 9" "19 19 11"; do
  tag=$(echo $m | tr -d ' ')
  step prod-rta-m${tag}-gauss $RK prod-rta-m${tag}-gauss 6 --mesh $m --method rta --temps 300 \
       --sigmas 0.05 0.1 0.2 --note "production setting, Gaussian smearing widths in THz"
done
step ooc-m191911-prod-s0.1 $OOC all --out runs/ooc-m191911-prod-s0.1 --mesh 19 19 11 --temps 300 --lang Rust \
     --sigma 0.1 --block 16 --rows-per-chunk 1024
echo "$(date -Iseconds) queue_final finished" >> $QLOG
