#!/bin/bash
# E39 with a memory guard (12 November): the first E39 attempt lost all four runs to allocation failures while other
# work held memory. Each run now starts only when at least 16 GB of physical memory is free, and a run that ends without
# its .json (failure) is archived as .failed<k>.out and retried, up to 4 attempts. Validation as in e39_multistart.sh.
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
O=results_c/film_opt
free_gb() { powershell -NoProfile -c "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB)"; }
for seed in 1 2; do
  for s in max min; do
    b=m3_orb0.1rs${seed}_z_L100_$s
    for k in 1 2 3 4; do
      [ -f $O/$b.json ] && break
      until [ "$(free_gb | tr -d '\r')" -ge 16 ]; do sleep 300; done
      echo "$(date) start $b attempt $k" >> queue/e39_guarded.log
      OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        /c/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe -B scripts_c/robust_mesh.py --tag orb --eps 0.1 \
        --col 2 --sense $s --meshes m553,m775,m996 --objective m553,m775,m996 --film 100 --start rand:$seed --out $O/$b \
        > $O/$b.out 2>&1
      [ -f $O/$b.json ] || { mv $O/$b.out $O/$b.failed$k.out; echo "$(date) $b attempt $k failed" >> queue/e39_guarded.log; }
    done
  done
done
# validation of every finished run (same chain as E37), then the end marker
bash queue/e39_multistart.sh >> queue/e39_guarded.log 2>&1
