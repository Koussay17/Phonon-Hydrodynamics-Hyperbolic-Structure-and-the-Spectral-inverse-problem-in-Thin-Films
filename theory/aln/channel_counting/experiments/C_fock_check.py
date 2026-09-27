"""Independent finite Fock-space channel-counting checks (campaign approach C).

Convention: H3 = (1/6) sum_ordered T[l,m,n] X_l X_m X_n,
X_l = a_l + a^dagger_{bar(l)}; T symmetric and reciprocal-conjugate.
This script never imports material code or converts amplitudes to irreversible rates.
Requires Python >=3.10, NumPy, SciPy. Run with python -B C_fock_check.py.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import platform
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy import sparse

T = 1.0 + 2.0j


def sparse_norm(operator):
    return float(np.sqrt(np.sum(np.abs(operator.data) ** 2)))


def ladder_operators(n_modes, local_dim):
    local = sparse.diags(np.sqrt(np.arange(1, local_dim)), 1,
                         shape=(local_dim, local_dim), format="csr")
    identity = sparse.eye(local_dim, format="csr")
    result = []
    for mode in range(n_modes):
        operator = sparse.csr_matrix([[1.0]])
        for position in range(n_modes):
            operator = sparse.kron(
                operator, local if position == mode else identity, format="csr"
            )
        result.append(operator)
    return result


def tensor_hamiltonian(fields, incoming, reverse, unique):
    dimension = fields[0].shape[0]
    operator = sparse.csr_matrix((dimension, dimension), dtype=complex)
    reciprocal = tuple(reverse[index] for index in incoming)
    for indices, coefficient in [(incoming, T), (reciprocal, T.conjugate())]:
        permutations = list(itertools.permutations(indices))
        if unique:
            permutations = sorted(set(permutations))
        for ordering in permutations:
            operator = operator + (coefficient / 6.0) * (
                fields[ordering[0]] @ fields[ordering[1]] @ fields[ordering[2]]
            )
    operator.eliminate_zeros()
    return operator


def energy_projection(operator, harmonic_energies):
    coo = operator.tocoo()
    keep = harmonic_energies[coo.row] == harmonic_energies[coo.col]
    projected = sparse.coo_matrix(
        (coo.data[keep], (coo.row[keep], coo.col[keep])), shape=coo.shape
    ).tocsr()
    projected.eliminate_zeros()
    return projected


def complex_pair(value):
    return [float(value.real), float(value.imag)]


def run_model(repeated):
    if repeated:
        names = ["p", "d", "bar_p", "bar_d"]
        reverse = [2, 3, 0, 1]
        incoming = (0, 3, 3)
        energies = np.array([2, 1, 2, 1], dtype=int)
        local_dim = 4
        daughter_indices = (1, 1)
        coefficient = T / 2.0
    else:
        names = ["p", "a", "b", "bar_p", "bar_a", "bar_b"]
        reverse = [3, 4, 5, 0, 1, 2]
        incoming = (0, 4, 5)
        energies = np.array([3, 1, 2, 3, 1, 2], dtype=int)
        local_dim = 3
        daughter_indices = (1, 2)
        coefficient = T
    n_modes = len(names)
    dims = (local_dim,) * n_modes
    annihilators = ladder_operators(n_modes, local_dim)
    creators = [operator.getH().tocsr() for operator in annihilators]
    fields = [annihilators[j] + creators[reverse[j]] for j in range(n_modes)]
    raw = tensor_hamiltonian(fields, incoming, reverse, unique=True)
    occupations = np.array(list(itertools.product(range(local_dim), repeat=n_modes)))
    harmonic_energies = occupations @ energies
    projected = energy_projection(raw, harmonic_energies)

    a, b = daughter_indices
    decay = creators[a] @ creators[b] @ annihilators[0]
    reciprocal_decay = (
        creators[reverse[a]] @ creators[reverse[b]] @ annihilators[reverse[0]]
    )
    collected = (
        coefficient * decay + coefficient.conjugate() * decay.getH()
        + coefficient.conjugate() * reciprocal_decay
        + coefficient * reciprocal_decay.getH()
    ).tocsr()
    difference = projected - collected
    samples = []
    counts = itertools.product([1, 2], [0, 1]) if repeated else itertools.product([1, 2], [0, 1], [0, 1])
    for count in counts:
        initial = [0] * n_modes
        initial[0] = count[0]
        initial[1] = count[1]
        if not repeated:
            initial[2] = count[2]
        final = initial.copy()
        final[0] -= 1
        final[a] += 1
        final[b] += 1
        row = np.ravel_multi_index(tuple(final), dims)
        col = np.ravel_multi_index(tuple(initial), dims)
        if repeated:
            ladder = math.sqrt(count[0] * (count[1] + 1) * (count[1] + 2))
        else:
            ladder = math.sqrt(count[0] * (count[1] + 1) * (count[2] + 1))
        predicted = coefficient * ladder
        measured = complex(raw[row, col])
        reverse_measured = complex(raw[col, row])
        reciprocal_initial = [initial[reverse[j]] for j in range(n_modes)]
        reciprocal_final = [final[reverse[j]] for j in range(n_modes)]
        reciprocal_measured = complex(raw[
            np.ravel_multi_index(tuple(reciprocal_final), dims),
            np.ravel_multi_index(tuple(reciprocal_initial), dims)
        ])
        samples.append({
            "initial": initial, "final": final,
            "matrix_amplitude": complex_pair(measured),
            "ladder_prediction": complex_pair(predicted),
            "absolute_error": abs(measured - predicted),
            "reverse_conjugation_error": abs(reverse_measured - predicted.conjugate()),
            "reciprocal_conjugation_error": abs(reciprocal_measured - predicted.conjugate()),
            "energy_change": int(np.dot(energies, np.subtract(final, initial))),
            "total_number_change": int(sum(final) - sum(initial)),
        })
    vacuum_sample = samples[0]
    vacuum_amplitude = complex(*vacuum_sample["matrix_amplitude"])
    coo = projected.tocoo()
    commutator_residual = float(np.sqrt(np.sum(np.abs(
        coo.data * (harmonic_energies[coo.row] - harmonic_energies[coo.col])
    ) ** 2)))
    checks = {
        "names": names, "local_dimension": local_dim, "matrix_dimension": local_dim**n_modes,
        "harmonic_energies": energies.tolist(), "incoming_tuple": list(incoming),
        "unique_ordered_tuple_count": len(set(itertools.permutations(incoming))),
        "collected_coefficient": complex_pair(coefficient),
        "raw_nonzero_entries": raw.nnz, "projected_nonzero_entries": projected.nnz,
        "hermiticity_absolute_frobenius_error": sparse_norm(raw - raw.getH()),
        "projected_vs_collected_absolute_frobenius_error": sparse_norm(difference),
        "projected_vs_collected_relative_frobenius_error": sparse_norm(difference) / sparse_norm(collected),
        "energy_commutator_frobenius_error": commutator_residual,
        "max_sample_amplitude_error": max(item["absolute_error"] for item in samples),
        "max_reverse_conjugation_error": max(item["reverse_conjugation_error"] for item in samples),
        "max_reciprocal_conjugation_error": max(item["reciprocal_conjugation_error"] for item in samples),
        "parent_one_daughter_vacuum_squared_amplitude": abs(vacuum_amplitude) ** 2,
        "samples": samples,
    }
    if repeated:
        duplicated = tensor_hamiltonian(fields, incoming, reverse, unique=False)
        row = np.ravel_multi_index(tuple(vacuum_sample["final"]), dims)
        col = np.ravel_multi_index(tuple(vacuum_sample["initial"]), dims)
        bad_amplitude = complex(duplicated[row, col])
        checks["duplicated_permutation_negative_control"] = {
            "amplitude_ratio": complex_pair(bad_amplitude / vacuum_amplitude),
            "squared_amplitude_ratio": abs(bad_amplitude / vacuum_amplitude) ** 2,
            "matrix_vs_twice_correct_absolute_frobenius_error": sparse_norm(duplicated - 2.0 * raw),
        }
    assert checks["projected_vs_collected_relative_frobenius_error"] < 1e-14
    assert checks["max_sample_amplitude_error"] < 1e-13
    assert checks["max_reverse_conjugation_error"] < 1e-13
    assert checks["max_reciprocal_conjugation_error"] < 1e-13
    assert checks["hermiticity_absolute_frobenius_error"] < 1e-12
    assert commutator_residual == 0.0
    assert all(item["energy_change"] == 0 and item["total_number_change"] == 1 for item in samples)
    return checks


def rational_moment_checks():
    result = []
    for ratio, cutoffs in [
        (Fraction(1, 2), [4, 8, 16, 32, 64]),
        (Fraction(3, 4), [8, 16, 32, 64, 128]),
    ]:
        mean = ratio / (1 - ratio)
        expected_annihilation = 2 * mean**2
        expected_creation = 2 * (1 + mean)**2
        convergence = []
        for cutoff in cutoffs:
            mass = sum((1 - ratio) * ratio**n for n in range(cutoff + 1))
            annihilation = sum((1 - ratio) * ratio**n * n * (n - 1)
                               for n in range(cutoff + 1))
            creation = sum((1 - ratio) * ratio**n * (n + 1) * (n + 2)
                          for n in range(cutoff + 1))
            assert 1 - mass == ratio**(cutoff + 1)
            assert 0 <= annihilation <= expected_annihilation
            assert 0 <= creation <= expected_creation
            convergence.append({
                "cutoff": cutoff, "omitted_probability": float(1 - mass),
                "annihilation_partial_moment": float(annihilation),
                "creation_partial_moment": float(creation),
                "annihilation_relative_error": float((expected_annihilation - annihilation) / expected_annihilation),
                "creation_relative_error": float((expected_creation - creation) / expected_creation),
            })
        assert convergence[-1]["annihilation_relative_error"] < 1e-11
        assert convergence[-1]["creation_relative_error"] < 1e-11
        result.append({
            "geometric_ratio_exact": str(ratio), "mean_exact": str(mean),
            "annihilation_moment_exact": str(expected_annihilation),
            "creation_moment_exact": str(expected_creation),
            "convergence": convergence,
        })
    # Independent mean-one distributions disprove a closure based on the mean alone.
    state_control = {
        "common_mean": 1,
        "fock_N1_annihilation_factorial_moment": 0,
        "fock_N1_creation_factorial_moment": 6,
        "geometric_mean1_annihilation_factorial_moment": 2,
        "geometric_mean1_creation_factorial_moment": 8,
    }
    # Exact detailed-balance examples: q_p=q_a*q_b or q_p=q_d**2.
    np_dist, na, nb = Fraction(1, 5), Fraction(1), Fraction(1, 2)
    np_rep, nd = Fraction(1, 3), Fraction(1)
    dist_forward = np_dist * (1 + na) * (1 + nb)
    dist_reverse = (1 + np_dist) * na * nb
    rep_forward = np_rep * 2 * (1 + nd)**2
    rep_reverse = (1 + np_rep) * 2 * nd**2
    assert dist_forward == dist_reverse
    assert rep_forward == rep_reverse
    return {
        "geometric_partial_sums": result, "mean_only_negative_control": state_control,
        "exact_detailed_balance": {
            "distinct_forward_and_reverse": str(dist_forward),
            "repeated_factorial_forward_and_reverse": str(rep_forward),
            "residuals_exact": [0, 0],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = {
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpreter": sys.executable, "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "tensor_amplitude": complex_pair(T),
        "hamiltonian_convention": "(1/6) sum_ordered T[l,m,n] X_l X_m X_n",
        "distinct": run_model(False), "repeated": run_model(True),
        "thermal": rational_moment_checks(),
        "limitations": [
            "No material/source normalization or exported pp mapping is tested.",
            "Energy projection is algebraic; finite matrices imply no irreversible rate.",
            "Selected ladder paths stay below truncation boundaries.",
            "Both tested incoming orbits are disjoint from their reciprocal orbit.",
            "Thermal sums are unnormalized partial sums, with the omitted mass reported.",
            "Geometric independent occupations are an extra kinetic closure assumption.",
        ],
        "all_assertions_passed": True,
    }
    results["vacuum_squared_amplitude_repeated_over_distinct"] = (
        results["repeated"]["parent_one_daughter_vacuum_squared_amplitude"]
        / results["distinct"]["parent_one_daughter_vacuum_squared_amplitude"]
    )
    destination = args.output_dir / "C_fock_check.json"
    destination.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(destination), "script_sha256": results["script_sha256"],
        "all_assertions_passed": True,
        "distinct_relative_matrix_error": results["distinct"]["projected_vs_collected_relative_frobenius_error"],
        "repeated_relative_matrix_error": results["repeated"]["projected_vs_collected_relative_frobenius_error"],
        "vacuum_squared_amplitude_ratio": results["vacuum_squared_amplitude_repeated_over_distinct"],
        "max_thermal_final_relative_error": max(
            max(case["convergence"][-1]["annihilation_relative_error"],
                case["convergence"][-1]["creation_relative_error"])
            for case in results["thermal"]["geometric_partial_sums"]
        ),
    }, indent=2))


if __name__ == "__main__":
    main()

