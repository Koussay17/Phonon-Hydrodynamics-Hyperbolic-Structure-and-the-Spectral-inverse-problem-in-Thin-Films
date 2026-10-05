"""Shared helpers for the AlN transport-baseline campaign (workstream B).

Inputs: hash-pinned Phonon Olympics phono3py files for AlN (author repository
McGaughey-Lab/Phonon-Olympics, commit 0640f07735059be9717a7565c2a0f22dc0da7a17,
folder "Aluminum Nitride/phono3py/AlN_kappa_input_files").  Every run refuses to
start if any SHA-256 differs from the download manifest.

Locations are resolved in this order:
  1. environment variables ALN_CAMPAIGN / ALN_DATASET,
  2. the D: campaign and dataset directories (primary),
  3. the C: staging copy (used only while D: is unavailable; same hashes).
"""
from __future__ import annotations

import ctypes
import ctypes.wintypes as wt
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import time
from pathlib import Path

D_CAMPAIGN = Path(r"D:\ResearchLab\orchestration\campaigns\20261002-aln-transport-baseline")
D_DATASET = Path(r"D:\ResearchLab\datasets\phonon-olympics"
                 r"\0640f07735059be9717a7565c2a0f22dc0da7a17\AlN_phono3py")
C_STAGING = Path(r"C:\Users\Koussay\ResearchLab\experiments\20261002-aln-transport-baseline-staging")


def _resolve(env, primary, fallback):
    if os.environ.get(env):
        return Path(os.environ[env])
    if primary.exists():
        return primary
    return fallback


CAMPAIGN = _resolve("ALN_CAMPAIGN", D_CAMPAIGN, C_STAGING)
DATASET = _resolve("ALN_DATASET", D_DATASET, C_STAGING / "inputs")

PINNED_SHA256 = {
    "POSCAR": "f4ffc5437f2d7c00457744c62bb88e9de91f21b0cda12cdb0db6e25407e909e4",
    "BORN": "e3783b344d87c80cd1615bd83e42531a28b849aeebb21ddc2054c63b3ca63003",
    "fc2.hdf5": "3fee31949075f3facf6f8510e4529bbf275a427f447c4c6ba29eec57e780763b",
    "fc3.hdf5": "d0237fe50df7e61d0db28976835bc86bc1f1bdedaf71a96bc6f6bf55e72be900",
}
SUPERCELL_FC3 = [3, 3, 2]
SUPERCELL_FC2 = [5, 5, 3]


def sha256(path: os.PathLike | str, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def verify_inputs() -> dict:
    rec = {}
    for name, ref in PINNED_SHA256.items():
        p = DATASET / name
        d = sha256(p)
        if d != ref:
            raise RuntimeError(f"Input hash mismatch for {p}: {d} != {ref}")
        rec[name] = {"path": str(p), "sha256": d, "bytes": p.stat().st_size}
    return rec


def software_versions() -> dict:
    pk = {}
    for k in ["phono3py", "phonopy", "phonors", "numpy", "scipy", "h5py", "spglib", "symfc"]:
        try:
            pk[k] = importlib.metadata.version(k)
        except importlib.metadata.PackageNotFoundError:
            pk[k] = None
    return {"python": sys.version, "executable": sys.executable,
            "platform": platform.platform(), "processor": platform.processor(),
            "packages": pk,
            "env_threads": {k: os.environ.get(k) for k in
                            ["OMP_NUM_THREADS", "RAYON_NUM_THREADS",
                             "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]}}


class _PMC(ctypes.Structure):
    _fields_ = [("cb", wt.DWORD), ("PageFaultCount", wt.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


def memory_info() -> dict:
    """Current and peak working set / private commit of this process (GiB)."""
    try:
        pmc = _PMC()
        pmc.cb = ctypes.sizeof(_PMC)
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.GetCurrentProcess.restype = wt.HANDLE
        k32.K32GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(_PMC), wt.DWORD]
        k32.K32GetProcessMemoryInfo.restype = wt.BOOL
        h = k32.GetCurrentProcess()
        ok = k32.K32GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb)
        if not ok:
            return {}
        return {"working_set_GB": pmc.WorkingSetSize / 2**30,
                "peak_working_set_GB": pmc.PeakWorkingSetSize / 2**30,
                "private_GB": pmc.PagefileUsage / 2**30,
                "peak_private_GB": pmc.PeakPagefileUsage / 2**30}
    except Exception:  # pragma: no cover
        return {}


class Timer:
    def __init__(self):
        self.t0 = time.perf_counter()
        self.marks = {}

    def mark(self, name):
        self.marks[name] = time.perf_counter() - self.t0
        return self.marks[name]


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default) + "\n", encoding="utf-8")


def _json_default(o):
    import numpy as np
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))
