#!/bin/sh
# E25: long-wavelength acoustic linewidths, convergence in the integration mesh
P=C:/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
S=D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/scripts_c/ir_linewidths.py
O=D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results/ir_lw
Q="--q 1/48,0,0 --q 2/48,0,0 --q 3/48,0,0 --q 4/48,0,0 --q 6/48,0,0 --q 8/48,0,0 --q 0,0,1/24 --q 0,0,2/24 --q 0,0,3/24 --q 0,0,4/24 --q 0,0,6/24"
$P -B $S --mesh 48 48 24 $Q --out $O/lw_m48.json --threads 12
$P -B $S --mesh 96 96 48 $Q --q 1/96,0,0 --q 0,0,1/48 --out $O/lw_m96.json --threads 12
$P -B $S --mesh 144 144 72 $Q --q 1/144,0,0 --q 2/144,0,0 --q 0,0,1/72 --q 0,0,2/72 --out $O/lw_m144.json --threads 12
