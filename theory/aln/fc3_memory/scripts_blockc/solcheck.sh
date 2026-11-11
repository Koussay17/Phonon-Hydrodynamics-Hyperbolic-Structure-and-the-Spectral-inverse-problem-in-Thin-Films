#!/bin/bash
# PBEsol code-path checks (added 11 Nov): folding test and Gamma-vs-complex path, before any further PBEsol production.
B=/home/koussay/blockc; cd $B; export OMP_NUM_THREADS=1
PW="mpirun --use-hwthread-cpus -np 4 pw.x -nk 1"
log() { echo "$(date -u +%FT%TZ) $*" >> $B/progress.log; }
mkdir -p solcheck && cd solcheck
# 1. primitive PBEsol (reference structure) at 60/480 on Gamma-centred 3x3x2
sed -e 's/f60_480/sol60/g' -e 's/pbe-n-kjpaw/pbesol-n-kjpaw/g' -e 's/celldm(1)=5.88809843/celldm(1)=5.888098425987400/' $B/fold_check/f60_480.in \
 | awk 'BEGIN{c=0} /CELL_PARAMETERS/{print "CELL_PARAMETERS alat"; print "   1.000000000   0.000000000   0.000000000"; print "  -0.500000000   0.866025404   0.000000000"; print "   0.000000000   0.000000000   1.600661069"; c=3; next} c>0{c--; next} {print}' \
 | sed -e 's/0.9998140361/0.9995896991/; s/0.4998140361/0.4995896991/; s/0.3813017874/0.3815261244/; s/0.8813017874/0.8815261244/' > prim_sol60.in
grep -q 'pbesol' prim_sol60.in || { log "solcheck: input build failed"; exit 1; }
[ -f prim_sol60.out ] || { $PW -in prim_sol60.in > prim_sol60.out 2>&1; log "solcheck primitive PBEsol 60/480 done"; }
# 2. perfect supercell PBEsol, complex path, default mixing
sed 's/^K_POINTS gamma/K_POINTS automatic\n1 1 1 0 0 0/' $B/pbesol/disp-00000.in > sc0_complex.in
[ -f sc0_complex.out ] || { mkdir -p tmp; $PW -in sc0_complex.in > sc0_complex.out 2>&1; rm -rf tmp; log "solcheck disp-00000 complex path done"; }
# 3. disp-00001 PBEsol complex path (Gamma result exists in pbesol/disp-00001.out)
sed 's/^K_POINTS gamma/K_POINTS automatic\n1 1 1 0 0 0/' $B/pbesol/disp-00001.in > d1_complex.in
[ -f d1_complex.out ] || { mkdir -p tmp; $PW -in d1_complex.in > d1_complex.out 2>&1; rm -rf tmp; log "solcheck disp-00001 complex path done"; }
log "solcheck finished"
