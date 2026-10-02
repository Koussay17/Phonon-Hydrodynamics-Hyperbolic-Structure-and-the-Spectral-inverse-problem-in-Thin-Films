"""Make both import conventions used by the suite work under bare `pytest`.

Most tests insert `src/` and import modules directly; a few import `src.<module>`.
Adding both the repository root and `src/` keeps `pytest` and `python -m pytest`
equivalent.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT, ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
