#!/bin/bash
# Relaunch the resumable queues after a reboot (each skips work already on disk). Run from Git Bash.
# run_all.sh, queue_A.sh, film_opt3.sh and m3_validate2.sh are finished (film_opt3.done, e37.ok); only E39 remains here.
# Block C (E38) runs in WSL and is relaunched from PowerShell (see E38).
cd /d/ResearchLab/orchestration/campaigns/20261005-fc3-memory
[ -f results_c/film3/e39.ok ] || nohup bash queue/e39_guarded.sh > queue/e39_guarded.nohup 2>&1 &
echo "$(date) relaunched after reboot" >> queue/logs/restarts.log
