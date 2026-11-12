"""Second job queue, safe to run beside run_queue.py: a job is started only if at most one other cx_* job runs; 15 s
after the start the count is re-checked and, if more than two jobs are running, this queue terminates its own child
(a subprocess it created) and retries later.  Nothing outside this queue is ever touched.
Usage: python run_queue2.py <queue name>
"""
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_queue import QUEUES, running_jobs, DIR, PY

env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
for job in QUEUES[sys.argv[1]]:
    tag = "_".join(a.replace(".py", "").replace(".json", "") for a in job)
    while True:
        while running_jobs() > 1:
            time.sleep(23)
        log = open(os.path.join(DIR, f"queue_{tag}.log"), "w")
        t0 = time.time()
        p = subprocess.Popen([PY] + job, cwd=DIR, stdout=log, stderr=subprocess.STDOUT, env=env)
        time.sleep(15)
        if running_jobs() > 2 and p.poll() is None:
            p.terminate()
            p.wait()
            log.close()
            print(time.strftime("%H:%M:%S"), "oversubscribed, retrying later:", " ".join(job), flush=True)
            time.sleep(60)
            continue
        print(time.strftime("%H:%M:%S"), "start", " ".join(job), flush=True)
        rc = p.wait()
        log.close()
        print(time.strftime("%H:%M:%S"), "done", " ".join(job), "rc", rc, "%.0fs" % (time.time() - t0), flush=True)
        break
    time.sleep(5)
