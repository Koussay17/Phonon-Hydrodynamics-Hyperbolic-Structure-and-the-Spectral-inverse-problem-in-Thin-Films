#!/bin/bash
# Block C: second functional (PBE) against PBEsol, 3x3x2 supercell, phono3py displacements (0.03 A, pair cutoff 3.9 A).
# Resumable: every step skips work already finished. Launch from PowerShell (see campaign E38), not with nohup.
set -u
B=/home/koussay/blockc; cd $B
export OMP_NUM_THREADS=1
PW="mpirun --use-hwthread-cpus -np 4 pw.x -nk 1"
PY=$B/venv/bin/python
log() { echo "$(date -u +%FT%TZ) $*" >> $B/progress.log; }
done_ok() { [ -f "$1" ] && grep -q "JOB DONE" "$1"; }
runpw() {  # dir input output
  ( cd $1 && mkdir -p tmp && /usr/bin/time -v $PW -in $2 > $3.part 2> $3.time && mv $3.part $3; rm -rf tmp )
}
# 1. PBE relaxation of the primitive cell
mkdir -p relax
done_ok relax/relax_pbe.out || { log "relax pbe start"; cp relax_pbe.in relax/; runpw relax relax_pbe.in relax_pbe.out; log "relax pbe end"; }
done_ok relax/relax_pbe.out || { log "relax pbe FAILED"; exit 1; }
# 2. cutoff test (pre-registered rule in ecut_check.py); 120/600 reference = pbe/disp-0000{0,1}.out
if [ ! -f ecut.choice ]; then
  for e in "80 480" "60 480"; do
    t=ecut_test/${e/ /_}; mkdir -p $t
    for f in disp-00000 disp-00001; do
      sed "s/ecutwfc .*/ecutwfc         = ${e% *}/; s/ecutrho .*/ecutrho         = ${e#* }/" pbe/$f.in > $t/$f.in
      done_ok $t/$f.out || { t0=$(date +%s); runpw $t $f.in $f.out; log "ecut test $e $f $(( $(date +%s)-t0 )) s"; }
    done
  done
  $PY ecut_check.py >> $B/progress.log 2>&1 || { log "ecut check FAILED"; exit 1; }
fi
read EW ER < ecut.choice
if [ ! -f inputs_${EW}_${ER}.done ]; then
  if [ "$EW" != 120 ]; then mkdir -p pbe_e120_invalid && mv pbe/*.out pbe/*.time pbe_e120_invalid/ 2>/dev/null; log "cutoff $EW/$ER chosen; invalid 120/600 outputs moved to pbe_e120_invalid/"; fi
  $PY gen_disp.py relax/relax_pbe.out pbe pbe $EW $ER >> $B/progress.log 2>&1 && \
  $PY gen_disp.py reference pbesol pbesol $EW $ER >> $B/progress.log 2>&1 && touch inputs_${EW}_${ER}.done
fi
# 2b. code-path check at the production cutoff (kcheck.py): Gamma-trick vs complex k-point path for disp-00001
if [ ! -f kcheck.ok ]; then
  mkdir -p kcheck
  sed 's/^K_POINTS gamma/K_POINTS automatic\n1 1 1 0 0 0/' pbe/disp-00001.in > kcheck/disp-00001.in
  done_ok pbe/disp-00001.out    || { t0=$(date +%s); runpw pbe disp-00001.in disp-00001.out; log "pbe disp-00001.in $(( $(date +%s)-t0 )) s"; }
  done_ok kcheck/disp-00001.out || { t0=$(date +%s); runpw kcheck disp-00001.in disp-00001.out; log "kcheck complex path $(( $(date +%s)-t0 )) s"; }
  $PY kcheck.py >> $B/progress.log 2>&1 || { log "kcheck FAILED: stop"; exit 1; }
fi
# 3. supercell SCFs: PBE branch first (a complete PBE fc3 as early as possible), then the PBEsol control
for d in pbe pbesol; do
  for i in $d/disp-*.in; do
    f=$(basename $i); o=${f%.in}.out
    done_ok $d/$o && continue
    t0=$(date +%s); runpw $d $f $o; log "$d $f rc=$? $(( $(date +%s)-t0 )) s"
  done
  log "$d branch complete"
done
# 4. PBEsol relaxation with the same settings (consistency check of the reference structure)
done_ok relax/relax_pbesol.out || { cp relax_pbesol.in relax/; runpw relax relax_pbesol.in relax_pbesol.out; log "relax pbesol end"; }
log "pipeline finished"
