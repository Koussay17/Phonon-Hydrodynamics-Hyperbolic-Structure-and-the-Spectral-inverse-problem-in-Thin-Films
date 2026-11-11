"""Task 3 complement: per-cluster-orbit bounded minima at eps = 0.05 and 0.2 (same settings as t3_batch.py).
usage: python -B t3_batch_min.py <events.npz> <resultdir> <outdir> <nworkers> [maxouter maxiter]
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import t3_batch  # noqa: E402

t3_batch.TASKS = [(0, "min", 0.05), (2, "min", 0.05), (0, "min", 0.2), (2, "min", 0.2)]

if __name__ == "__main__":
    t3_batch.main()
