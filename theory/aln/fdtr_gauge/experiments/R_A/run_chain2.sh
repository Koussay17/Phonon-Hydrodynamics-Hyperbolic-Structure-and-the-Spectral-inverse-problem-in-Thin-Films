#!/bin/bash
# R-A chain 2 (sequential): sensitivity box (JAX), band check of the kinetic kernel on 7x7x5, delta analyses, tables.
cd /d/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/experiments/R_A
PY=/c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe
$PY -u -B ra4_box.py > out/ra4.log 2>&1; echo "ra4 EXIT $?" >> out/chain2.log
$PY -u -B ra1d_freq_kernel.py m775 300,500 > out/ra1d_m775.log 2>&1; echo "ra1d m775 EXIT $?" >> out/chain2.log
$PY -u -B ra5_delta.py > out/ra5.log 2>&1; echo "ra5 EXIT $?" >> out/chain2.log
$PY -u -B ra5b_observable.py > out/ra5b.log 2>&1; echo "ra5b EXIT $?" >> out/chain2.log
$PY -u -B ra9_tables.py > out/ra9_tables.txt 2> out/ra9.err; echo "ra9 EXIT $?" >> out/chain2.log
echo CHAIN2 DONE >> out/chain2.log
