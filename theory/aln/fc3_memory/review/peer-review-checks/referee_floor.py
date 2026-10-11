"""Referee check: tau / tau_CS and the rigorous floor (1 - dK) tau_CS / tau (dK = 1 %) per mesh, reference fc3."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from referee_check import Mesh, ROOT, FOURPI  # noqa: E402

out = {}
for name, ev, pp in (("m553", ROOT / "m553" / "events_m553_s0.1.npz", ROOT / "m553" / "pp_r01_REF.npz"),
                     ("m775", ROOT / "m775x" / "events_m775_s0.1.npz", ROOT / "m775x" / "pp_r01_REF.npz"),
                     ("m996", ROOT / "m996" / "events_m996_s0.1.npz", ROOT / "m996" / "pp_r01_REF.npz")):
    M = Mesh(ev)
    g, _ = M.rates(pp)
    m = M.metrics(g, want_x=True)
    rec = {}
    for ax, c in (("x", 0), ("z", 2)):
        b = M.b[:, c]
        tcs = m[f"K_{ax}"] / float(b @ b) / FOURPI
        rec[ax] = {"tau": m[f"tau_{ax}"], "tauCS": tcs, "tau_over_tauCS": m[f"tau_{ax}"] / tcs,
                   "floor_ratio": 0.99 * tcs / m[f"tau_{ax}"]}
    out[name] = rec
    print(name, json.dumps(rec), flush=True)
Path(Path(__file__).resolve().parent / "referee_floor.json").write_text(json.dumps(out, indent=1))
