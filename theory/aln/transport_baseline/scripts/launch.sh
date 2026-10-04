#!/bin/bash
# Launch one phono3py run with fixed thread settings and a log file.
# usage: launch.sh <run-id> <threads> [run_kappa.py args except --run-id]
# Interpreter: the isolated D: environment when mounted, otherwise the C: replica
# built from the same pinned requirements (phono3py/phonopy 4.5.0, phonors 0.3.0,
# numpy 2.5.3, scipy 1.18.1, h5py 3.16.0, spglib 2.7.0, Python 3.14.7).
set -u
RUN_ID="$1"; shift
NT="$1"; shift
SCRIPTS="$(cd "$(dirname "$0")" && pwd -W)"
if [ -x /d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe ]; then
  PY=/d/ResearchLab/envs/aln-phono3py-4.5.0/Scripts/python.exe
else
  PY=/c/Users/Koussay/ResearchLab/envs/aln-phono3py-4.5.0-replica/Scripts/python.exe
fi
CAMP="$("$PY" -c "import sys; sys.path.insert(0, r'$SCRIPTS'); import common; print(common.CAMPAIGN.as_posix())")"
mkdir -p "$CAMP/runs/logs"
LOG="$CAMP/runs/logs/${RUN_ID}.log"
export OMP_NUM_THREADS="$NT" RAYON_NUM_THREADS="$NT" OPENBLAS_NUM_THREADS="$NT" MKL_NUM_THREADS="$NT"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
echo "# launch $(date -Iseconds) threads=$NT python=$PY campaign=$CAMP args: $*" > "$LOG"
"$PY" -X utf8 -B "$SCRIPTS/run_kappa.py" --run-id "$RUN_ID" "$@" >> "$LOG" 2>&1
RC=$?
echo "# exit code $RC at $(date -Iseconds)" >> "$LOG"
exit $RC
