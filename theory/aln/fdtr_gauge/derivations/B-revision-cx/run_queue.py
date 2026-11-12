"""Sequential job queue for the counterexample scripts, throttled so that at most two compute jobs of this branch run
at the same time (each job is single-threaded).  A job is launched only when at most one other cx_* interpreter is
running.  Read-only process listing via PowerShell; nothing is ever stopped by this script.
Usage: python run_queue.py <queue name>
"""
import os
import subprocess
import sys
import time

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx"
PY = r"C:/Users/Koussay/ResearchLab/envs/ultimate-research/Scripts/python.exe"
QUEUES = {
    "A": [
        ["cx_06_joint_pairs.py", "entropic", "4", "10"],
        ["cx_06_joint_pairs.py", "entropic", "5", "10"],
        ["cx_06_joint_pairs.py", "local", "5", "10"],
        ["cx_05_multistart.py", "aln", "20000"],
        ["cx_12_homotopy_to_aln.py", "entropic", "band1", "20"],
        ["cx_12_homotopy_to_aln.py", "entropic", "band2", "20"],
        ["cx_13_four_freq_random_targets.py", "entropic", "30", "2000"],
        ["cx_13_four_freq_random_targets.py", "local", "30", "2000"],
        ["cx_11_practical_profiles.py"],
        ["cx_02b_bare_free.py"],
        ["cx_05b_random_triples.py", "generic", "30", "4000", "entropic"],
        ["cx_12_homotopy_to_aln.py", "local", "band1", "10"],
        ["cx_07b_kernel_solve.py", "mixed", "300"],
        ["cx_07b_kernel_solve.py", "free", "300"],
    ],
    "B": [
        ["cx_09b_certify_pairs.py", "cx_06_joint_local_N4.json", "5", "cx_09b_certified_local_N4.json"],
        ["cx_14_two_freq_entropic.py"],
        ["cx_06_joint_pairs.py", "entropic", "6", "10", "moderate"],
        ["cx_06_joint_pairs.py", "local", "6", "10", "moderate"],
        ["cx_16_aln_conditioning.py"],
        ["cx_15_aln_region_pairs.py", "entropic", "3", "6"],
        ["cx_15_aln_region_pairs.py", "local", "3", "6"],
        ["cx_15_aln_region_pairs.py", "entropic", "4", "6"],
    ],
    "C": [
        ["cx_07b_kernel_solve.py", "tied2", "300"],
    ],
    "E": [
        ["cx_09_certify.py", "cases_aln_N3.json", "cx_09_certified_aln_N3.json"],
        ["cx_05_multistart.py", "aln2", "20000"],
    ],
    "D": [
        ["cx_09b_certify_pairs.py", "cx_06_joint_entropic_N4.json", "5", "cx_09b_certified_entropic_N4.json"],
        ["cx_09b_certify_pairs.py", "cx_06_joint_entropic_N5.json", "5", "cx_09b_certified_entropic_N5.json"],
    ],
}


def running_jobs():
    cmd = ("Get-CimInstance Win32_Process -Filter \"Name = 'python.exe'\" | Where-Object { $_.CommandLine -match 'cx_' "
           "-and $_.WorkingSetSize -gt 20MB } | ForEach-Object { $_.CommandLine }")
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", cmd], capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return 99
    return len([l for l in out.splitlines() if l.strip() and "run_queue.py" not in l])


if __name__ == "__main__":
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for job in QUEUES[sys.argv[1]]:
        while running_jobs() > 1:
            time.sleep(20)
        tag = "_".join(a.replace(".py", "") for a in job)
        log = open(os.path.join(DIR, f"queue_{tag}.log"), "w")
        t0 = time.time()
        print(time.strftime("%H:%M:%S"), "start", " ".join(job), flush=True)
        rc = subprocess.run([PY] + job, cwd=DIR, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
        log.close()
        print(time.strftime("%H:%M:%S"), "done", " ".join(job), "rc", rc, "%.0fs" % (time.time() - t0), flush=True)
        time.sleep(5)
