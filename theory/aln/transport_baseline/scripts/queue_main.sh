#!/bin/bash
# Sequential production queue for workstream B (run in background).
# Each job writes its own log under runs/logs; a job is skipped if its marker exists.
set -u
D=/d/ResearchLab/orchestration/campaigns/20261002-aln-transport-baseline
PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
export OMP_NUM_THREADS=6 RAYON_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 PYTHONUNBUFFERED=1
cd "$D"
QLOG=runs/logs/queue_main.log
step() {  # step <name> <command...>
  local name="$1"; shift
  if [ -f "runs/logs/$name.done" ]; then echo "$(date -Iseconds) skip $name" >> $QLOG; return; fi
  echo "$(date -Iseconds) start $name" >> $QLOG
  "$@" > "runs/logs/$name.log" 2>&1
  local rc=$?
  echo "$(date -Iseconds) end $name rc=$rc" >> $QLOG
  [ $rc -eq 0 ] && touch "runs/logs/$name.done"
}
OOC="$PY -X utf8 -B scripts/lbte_ooc.py"
RK="bash scripts/launch.sh"
ALLT="100 125 150 175 200 225 250 275 300 400 500 600 700 800 900 1000"

# Q2/Q3: production RTA at the Olympics mesh, all temperatures (N/U split; isotope)
step prod-rta-m313117-NU    $RK prod-rta-m313117-NU 6 --mesh 31 31 17 --method rta --temps $ALLT --is-N-U \
     --note "production setting (phono3py 4.5.0 defaults: Rust, make_r0_average=True, analytic gv), tetrahedron, no isotope, N/U split"
step prod-rta-m313117-iso   $RK prod-rta-m313117-iso 6 --mesh 31 31 17 --method rta --temps $ALLT --isotope \
     --note "production setting, natural-abundance isotope scattering"
# Q4: reproduction LBTE at 31x31x17 (phono3py-2.1.0-like settings) via out-of-core route
step ooc-m313117-v2C        $OOC all --out runs/ooc-m313117-v2C --mesh 31 31 17 --temps 300 --lang C \
     --no-r0-average --gv-delta-q 1e-5 --block 16 --rows-per-chunk 1024
# Q5: production LBTE at 31x31x17, temperature sweep; keep pre-assembly matrix at 300 K (index 2)
step ooc-m313117-prod       $OOC all --out runs/ooc-m313117-prod --mesh 31 31 17 \
     --temps 100 200 300 400 500 600 800 1000 --lang Rust --block 16 --rows-per-chunk 1024 --keep-raw 2
# Q9: isotope LBTE at 31x31x17, 300 K (diagonal isotope term added in the solve)
step ooc-m313117-prod-iso   $OOC solve --out runs/ooc-m313117-prod --temp-index 2 --rows-per-chunk 1024 \
     --gamma-iso-from runs/prod-rta-m313117-iso/kappa-m313117.hdf5
echo "$(date -Iseconds) queue_main finished" >> $QLOG
