"""D independent phase-information toys: NumPy only; no material inputs.

Run with --output ABSOLUTE_JSON_PATH to persist diagnostics.
K rows label one daughter block; columns label parent modes. The third leg is a
singleton. Hbar = k_B T = 1; parent frequencies 2, daughter frequencies 1.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import numpy as np


def floats(x):
    return np.asarray(x, dtype=float).tolist()


def gram_data(k):
    g = k.conj().T @ k
    return {
        "P": floats(np.abs(k)**2),
        "block_norm_squared": float(np.sum(np.abs(k)**2)),
        "G_real": floats(g.real),
        "G_imaginary": floats(g.imag),
        "G_eigenvalues": floats(np.linalg.eigvalsh(g)),
        "G_diagonal": floats(np.diag(g).real),
    }


def thermal_population_operator(p):
    # Exact resonant reactions a_i <-> b_j + c. Mode order a0,a1,b0,b1,c.
    # The positive common equilibrium Bose flux factor is absorbed into time.
    frequencies = np.array([2., 2., 1., 1., 1.])
    n = 1 / np.expm1(frequencies)
    susceptibility = n * (1 + n)
    onsager = np.zeros((5, 5))
    for j in range(2):
        for i in range(2):
            s = np.zeros(5)
            s[i], s[2+j], s[4] = -1., 1., 1.
            onsager += p[j, i] * np.outer(s, s)
    generator = -onsager / susceptibility[None, :]
    root = np.sqrt(susceptibility)
    entropy_operator = onsager / np.outer(root, root)
    # Block totals and uniform occupation embedding; B @ R is identity.
    block = np.array([[1., 1., 0., 0., 0.],
                      [0., 0., 1., 1., 0.],
                      [0., 0., 0., 0., 1.]])
    embed = block.T @ np.diag([.5, .5, 1.])
    defect = block @ generator @ (np.eye(5) - embed @ block)
    hidden = np.array([1., -1., 0., 0., 0.])
    record = {
        "mode_frequencies": floats(frequencies),
        "equilibrium_occupations": floats(n),
        "susceptibilities": floats(susceptibility),
        "P": floats(p),
        "L": floats(generator),
        "B": floats(block),
        "R": floats(embed),
        "projected_block_generator": floats(block @ generator @ embed),
        "closure_defect": floats(defect),
        "closure_defect_spectral_norm": float(np.linalg.norm(defect, 2)),
        "hidden_perturbation": floats(hidden),
        "hidden_block_totals": floats(block @ hidden),
        "hidden_block_derivative": floats(block @ generator @ hidden),
        "energy_conservation_residual": float(np.linalg.norm(frequencies @ generator)),
        "entropy_operator_symmetry_residual": float(np.linalg.norm(entropy_operator - entropy_operator.T)),
        "entropy_operator_eigenvalues": floats(np.linalg.eigvalsh(entropy_operator)),
    }
    return record


def covariance_evolution(k):
    # Explicit additional model: d(delta N)/dt = -{G,delta N}/2, G=K^*K.
    # Thermal Markov reservoir; N_eq=n_B(2) I. Delta N(0)=I/2; scale it small
    # for a perturbation. This is not a derivation of a full phonon generator.
    g = k.conj().T @ k
    lam, vec = np.linalg.eigh(g)
    rows = []
    for t in (0., .1, .5, 1., 2.):
        propagator = (vec * np.exp(-.5 * t * lam)) @ vec.conj().T
        perturbation = .5 * propagator @ propagator.conj().T
        rows.append({"t": t, "total_perturbation": float(np.trace(perturbation).real),
                     "offdiagonal_absolute": float(abs(perturbation[0, 1]))})
    return {"initial_total": 1., "initial_total_derivative": float(-.5*np.trace(g).real),
            "initial_total_second_derivative": float(.5*np.trace(g@g).real), "trajectory": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    h = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2.)
    k_plus = np.ones((2, 2))
    k_minus = np.array([[1., 1.], [1., -1.]])
    family = []
    for theta in (0., np.pi/3, np.pi/2, np.pi):
        k = np.array([[1., 1.], [1., np.exp(1j*theta)]])
        actual = np.linalg.eigvalsh(k.conj().T @ k)
        expected = np.array([2-2*abs(np.cos(theta/2)), 2+2*abs(np.cos(theta/2))])
        family.append({"theta": float(theta), "Gram_eigenvalues": floats(actual),
                       "expected_eigenvalues": floats(expected),
                       "maximum_eigenvalue_error": float(np.max(np.abs(actual-expected))),
                       "modulus_squared_error_from_ones": float(np.max(np.abs(np.abs(k)**2-1)))})
    cases = {}
    for label, k in (("plus", k_plus), ("minus", k_minus)):
        rotated = k @ h
        g = k.T @ k
        prediction_without_interference = np.abs(k)**2 @ np.abs(h)**2
        cases[label] = {
            "base": gram_data(k), "parent_Hadamard_rotated": gram_data(rotated),
            "Gram_covariance_residual": float(np.linalg.norm(rotated.T@rotated-h.T@g@h)),
            "phase_blind_rotated_P_prediction": floats(prediction_without_interference),
            "phase_blind_prediction_error_frobenius": float(np.linalg.norm(np.abs(rotated)**2-prediction_without_interference)),
            "base_population_operator": thermal_population_operator(np.abs(k)**2),
            "rotated_population_operator": thermal_population_operator(np.abs(rotated)**2),
            "explicit_thermal_covariance_model": covariance_evolution(k),
        }
    triple_cases = {}
    for label in ("all_plus", "one_minus"):
        tensor = np.ones((2, 2, 2))
        if label == "one_minus":
            tensor[1, 1, 1] = -1.
        grams = []
        for axis in range(3):
            unfolded = np.moveaxis(tensor, axis, 0).reshape(2, -1)
            grams.append(floats(np.linalg.eigvalsh(unfolded@unfolded.T)))
        transformed = np.einsum("ia,jb,kc,abc->ijk", h, h, h, tensor)
        triple_cases[label] = {
            "tensor": floats(tensor), "elementwise_squared_modulus": floats(np.abs(tensor)**2),
            "block_norm_squared": float(np.sum(np.abs(tensor)**2)),
            "Gram_eigenvalues_by_leg": grams,
            "all_leg_Hadamard_rotated_P": floats(np.abs(transformed)**2),
            "rotated_block_norm_squared": float(np.sum(np.abs(transformed)**2)),
            "permutation_symmetry_residual": float(max(np.linalg.norm(tensor-np.transpose(tensor,p))
                for p in ((0,2,1), (1,0,2), (2,1,0)))),
        }
    block_plus = np.array(cases['plus']['rotated_population_operator']['projected_block_generator'])
    block_minus = np.array(cases['minus']['rotated_population_operator']['projected_block_generator'])
    result = {
        "scope": "Independent finite synthetic constructions; exact degeneracy assumed; no AlN data or peer results read.",
        "environment": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__,
                        "platform": platform.platform(), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "continuous_phase_family": family,
        "two_by_two_cases": cases,
        "three_leg_cases": triple_cases,
        "rotated_projected_block_generator_difference": float(np.linalg.norm(block_plus-block_minus)),
        "conventions": {
            "K_indices": "K[daughter,parent], with singleton third leg; tensor entries are V[parent,daughter,0]=K[daughter,parent].",
            "G": "K conjugate-transpose times K",
            "population_generator": "L = -sum(P[j,i] s_ij s_ij^T) diag(n*(1+n))^-1; common Bose equilibrium factor absorbed into time.",
            "exact_closure_test": "B L (I - R B) = 0, where B sums block populations and R uniformly embeds block totals.",
            "covariance_model": "Additional specified thermal Markov model only, not a derivation from all cubic phonon dynamics.",
        },
    }
    encoded = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        if not args.output.is_absolute():
            raise ValueError('--output must be an absolute path')
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded)


if __name__ == '__main__':
    main()
