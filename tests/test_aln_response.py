"""Independent high-precision checks of cancellation-free bulk diagnostics."""
import importlib.util
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/analyse_aln_response.py"
SPEC = importlib.util.spec_from_file_location("aln_response", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_identical_lifetimes_have_exact_single_pole():
    for frequency in (1e-12, 1., 1e12):
        response, errors = MOD.pole_errors(
            1j * frequency, np.array([1., 2., 3.]), np.array([2., 2., 2.]), 2.)
        np.testing.assert_allclose(response, 6 / (1 + 2j * frequency), rtol=1e-15)
        assert all(value == 0 for value in errors.values())


@pytest.mark.parametrize("frequency", [1e-12, 1e-6, 1., 100.])
def test_reduction_errors_match_independent_high_precision_direct_subtraction(frequency):
    weights = np.array([1., 2., 7.])
    tau = np.array([.25, 3., 100.])
    ts = float(3000 / 1421)
    response, errors = MOD.pole_errors(1j * frequency, weights, tau, ts)
    with mp.workdps(100):
        s = mp.j * mp.mpf(frequency)
        w, t = list(map(mp.mpf, weights)), list(map(mp.mpf, tau))
        k0 = sum(w)
        exact = sum(a / (1 + s * b) for a, b in zip(w, t))
        tm = sum(a * b for a, b in zip(w, t)) / k0
        static = k0 / (1 + s * mp.mpf(ts))
        memory = k0 / (1 + s * tm)
        expected = {
            "static_pole_relative_complex_error": abs(exact - static) / abs(exact),
            "memory_pole_relative_complex_error": abs(exact - memory) / abs(exact),
            "static_pole_relative_dynamic_correction_error": abs(exact - static) / abs(exact - k0),
            "memory_pole_relative_dynamic_correction_error": abs(exact - memory) / abs(exact - k0),
        }
        np.testing.assert_allclose(response, complex(exact), rtol=2e-15)
        for name, value in expected.items():
            np.testing.assert_allclose(errors[name], float(value), rtol=3e-14, atol=0)
