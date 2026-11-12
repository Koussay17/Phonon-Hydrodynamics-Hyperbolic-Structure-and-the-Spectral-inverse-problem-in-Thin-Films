#!/bin/bash
# R-S chain 1: REF on the full FDTR frequency set (exact), fits on exact data, AAA interpolation check.
# One heavy process at a time, two BLAS threads (set in rs_core).
cd "$(dirname "$0")"
PY=C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
$PY -B -W ignore rs2_slab.py m775 film2 full,rta all REF > out/rs2_all_REF.log 2>&1 || exit 1
$PY -B -W ignore rs2b_fits.py m775 film2 all REF > out/rs2b_all_REF.log 2>&1 || exit 1
$PY -B -W ignore rs2c_interp_check.py > out/rs2c.log 2>&1 || exit 1
echo CHAIN1_DONE >> out/rs2c.log
