#!/bin/bash
# R-S chain 3 (after chain 2): 7x7x5 slabs for REF and the delta-extreme film witnesses (regularised full-operator
# delta on 7x7x5, both thicknesses, from rs4_filmw_static.json), black walls; fits; mesh-trend summary.
cd "$(dirname "$0")"
PY=C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
FS=${1:-sub4}
EXT=$($PY -c "import json; d=json.load(open('out/rs4_filmw_static.json'))['m775']['extremes']; print(' '.join(sorted({d[L][k] for L in d for k in ('max','min')})))")
echo "extremes: $EXT" > out/rs2_m775_filmw_${FS}.log
$PY -B -W ignore rs2_slab.py m775 filmw full,rta ${FS}-nowalls REF $EXT >> out/rs2_m775_filmw_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs2b_fits.py m775 filmw $FS REF $EXT > out/rs2b_m775_filmw_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs5_trend.py $FS > out/rs5_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs6_proxy.py $FS > out/rs6_${FS}.log 2>&1 || exit 1
echo CHAIN3_DONE >> out/rs5_${FS}.log
