#!/bin/bash
# R-S chain 2: REF + 12 witnesses on 7x7x5 (full operator with wall variants, then RTA), fits, margins;
# then the mesh-trend inputs (film witnesses: static times on 5x5x3 and 7x7x5, slabs on 5x5x3).
# usage: ./run_rs_chain2.sh <freqset, e.g. sub3>      One heavy process at a time, two BLAS threads.
cd "$(dirname "$0")"
PY=C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
FS=$1
W12="orb0.1_m3L30max orb0.1_m3L30min orb0.1_m3L100max orb0.1_m3L100min orb0.1_m3L300max orb0.1_m3L300min cal1nz_m3L100max cal1nz_m3L100min cal3wnz_m3L100max cal3wnz_m3L100min cal3w_m3L100max cal3w_m3L100min"
WW="cal1_fL100_max cal1_fL100_min orb0.1_fL30_max orb0.1_fL30_min orb0.1_fL100_max orb0.1_fL100_min orb0.1_fL300_max orb0.1_fL300_min orb0.1_fL1000_max orb0.1_fL1000_min"
$PY -B -W ignore rs2_slab.py m775 film2 full $FS REF $W12 > out/rs2_${FS}_full.log 2>&1 || exit 1
$PY -B -W ignore rs2_slab.py m775 film2 rta $FS REF $W12 > out/rs2_${FS}_rta.log 2>&1 || exit 1
$PY -B -W ignore rs2b_fits.py m775 film2 $FS REF $W12 > out/rs2b_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs3_margins.py $FS > out/rs3_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs4_filmw_static.py > out/rs4_filmw_static.log 2>&1 || exit 1
$PY -B -W ignore rs2_slab.py m553 filmw full,rta ${FS}-nowalls REF $WW > out/rs2_m553_${FS}.log 2>&1 || exit 1
$PY -B -W ignore rs2b_fits.py m553 filmw $FS REF $WW > out/rs2b_m553_${FS}.log 2>&1 || exit 1
echo CHAIN2_DONE >> out/rs2b_m553_${FS}.log
