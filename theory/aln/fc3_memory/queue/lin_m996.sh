#!/bin/sh
# 9x9x6 exact linearisation: 16 row chunks, 6 concurrent single-thread workers
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory/scripts_c
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
seq 0 15 | xargs -P 6 -I{} sh -c '[ -f D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/m996/lin_G_{}.npy ] || C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe -B linearize_parts.py part D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/m996/events_m996_s0.1.npz D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/m996 {} 16 > D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/m996/lin_part_{}.log 2>&1'
echo done > D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/m996/lin_parts.done
