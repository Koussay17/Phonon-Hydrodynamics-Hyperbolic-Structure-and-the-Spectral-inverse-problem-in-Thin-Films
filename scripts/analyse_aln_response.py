"""Uniform bulk diagonal-RTA conductivity, not a simulated FDTR signal."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.spectral import mode_heat_capacity, rta_moments


def pole_errors(s, weights, tau, tau_static):
    """Cancellation-free errors; s != 0 and a nonzero current are required.

    weights are modal DC contributions, all nonnegative.
    Returns exact response and errors relative to the total response and
    to its dynamic correction. Matching the memory pole cancels the linear
    term, so its difference is evaluated as a variance-weighted expression.
    """
    weights, tau = np.asarray(weights), np.asarray(tau)
    k0 = float(np.sum(weights))
    tm = float(np.sum(weights * tau) / k0)
    den = 1 + s * tau
    response = np.sum(weights / den)
    dynamic = -s * np.sum(weights * tau / den)
    static_diff = s / (1 + s * tau_static) * np.sum(
        weights * (tau_static - tau) / den)
    memory_diff = s**2 / (1 + s * tm)**2 * np.sum(
        weights * (tau - tm)**2 / den)
    return response, {
        "static_pole_relative_complex_error": float(abs(static_diff) / abs(response)),
        "memory_pole_relative_complex_error": float(abs(memory_diff) / abs(response)),
        "static_pole_relative_dynamic_correction_error": float(abs(static_diff) / abs(dynamic)),
        "memory_pole_relative_dynamic_correction_error": float(abs(memory_diff) / abs(dynamic)),
    }


def calculate(repo):
    source = repo / "theory/aln/rao_300K_modes.npz"
    with np.load(source, allow_pickle=False) as z:
        c = (mode_heat_capacity(z["omega_rad_s"], 300)
             * z["degeneracy"][:, None] / z["degeneracy"].sum()
             / z["cell_volume_m3"])
        v, r = z["velocity_m_s"], z["total_rate_s"]
        tau = np.divide(1., r, out=np.zeros_like(r), where=r > 0)
    frequencies = [1e3, 1e5, 1e6, 1e7, 5e7, 1e8, 2e8, 1e9]
    result = {
        "scope": ("Uniform bulk diagonal RTA at 300 K on one 24^3 mesh; "
                  "no spatial streaming, FDTR boundary response or N/U inference. "
                  "Continuum convergence is unestablished."),
        "harmonic_convention": "exp(+i omega t)",
        "error_method": "Algebraic differences; memory-pole variance identity avoids low-frequency cancellation",
        "provenance": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "input_relative_path": source.relative_to(repo).as_posix(),
            "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "directions": {},
    }
    for name, v2 in [("basal", (v[:, :, 0]**2 + v[:, :, 1]**2) / 2),
                     ("c_axis", v[:, :, 2]**2)]:
        moments = rta_moments(c, v2, r)
        weights = moments.pop("contributions")
        rows = []
        for f in frequencies:
            response, errors = pole_errors(2j * np.pi * f, weights, tau,
                                          moments["tau_static_s"])
            rows.append({
                "frequency_Hz": f,
                "kappa_real_W_mK": float(response.real),
                "kappa_imag_W_mK": float(response.imag),
                "phase_degree": float(np.angle(response, deg=True)),
                **errors,
            })
        result["directions"][name] = {"moments": moments, "response": rows}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    args = parser.parse_args()
    result = calculate(args.repo)
    (args.repo / "theory/aln/bulk_response.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for name, d in result["directions"].items():
        print(name, d["response"][-2])


if __name__ == "__main__":
    main()
