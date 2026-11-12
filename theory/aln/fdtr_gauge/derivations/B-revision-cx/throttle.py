"""Cooperative throttle: a job waits (light state, no numerical imports) while at least `max_others` other cx_* compute
jobs of this branch are running.  Read-only process listing; never stops anything."""
import os
import subprocess
import time


def others_running():
    me = os.getpid()
    parent = os.getppid()
    cmd = ("Get-CimInstance Win32_Process -Filter \"Name = 'python.exe'\" | Where-Object { $_.CommandLine -match 'cx_' "
           "-and $_.WorkingSetSize -gt 20MB } | ForEach-Object { \"$($_.ProcessId) $($_.ParentProcessId)\" }")
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", cmd], capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return 99
    n = 0
    for line in out.splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        pid, ppid = int(parts[0]), int(parts[1])
        if pid in (me, parent) or ppid == me:
            continue
        n += 1
    return n


def wait_slot(max_others=2, period=30):
    while others_running() >= max_others:
        time.sleep(period)
