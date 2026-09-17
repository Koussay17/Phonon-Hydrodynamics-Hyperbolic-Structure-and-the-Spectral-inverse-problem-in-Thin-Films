"""Fixed-diagonal collision counterexamples; exact, dense, spectral, and 80-dps checks.

Run: py constructive_collision.py
Only NumPy, SymPy, and mpmath are required. All collision rates are dimensionless.
The current-response convention omits the common positive k_B T^2 prefactor.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent


def walsh_signs(n: int) -> np.ndarray:
    assert n > 0 and n & (n - 1) == 0
    return np.array(
        [[(-1) ** ((i & j).bit_count()) for j in range(n)] for i in range(n)],
        dtype=float,
    )


def graph4(a: float, b: float, c: float) -> np.ndarray:
    # Opposite edges carry the same weight. Every degree is a+b+c.
    return np.array(
        [[a+b+c, -a, -b, -c], [-a, a+b+c, -c, -b],
         [-b, -c, a+b+c, -a], [-c, -b, -a, a+b+c]], dtype=float,
    )


def response_dense(collision: np.ndarray, current: np.ndarray) -> dict:
    n = len(current)
    energy = np.ones(n) / np.sqrt(n)
    lifted = collision + np.outer(energy, energy)
    x = np.linalg.solve(lifted, current)
    dc = float(current @ x)
    m2 = float(x @ x)
    return {
        "dc": dc, "m2": m2, "memory": m2 / dc,
        "energy_residual": float(np.linalg.norm(collision @ energy)),
        "solve_relative_residual": float(np.linalg.norm(collision @ x - current)
                                         / np.linalg.norm(current)),
        "diagonal_error": float(np.max(np.abs(np.diag(collision) - 1.0))),
        "energy_gauge_error": float(abs(energy @ x)),
    }


def family8_matrix(t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # A is the operator in the orthonormal Walsh basis; mode zero is energy.
    a = np.diag([0.0, t*t, 2.0, 2.0-t*t, 1.0, 1.0, 1.0, 1.0])
    a[1, 2] = a[2, 1] = t
    a[4, 7] = a[7, 4] = -t
    signs = walsh_signs(8)
    collision = signs @ a @ signs.T / 8
    return collision, signs[:, 2] / np.sqrt(8), signs / np.sqrt(8)


def exact_symbolic_checks() -> dict:
    t, p = sp.symbols("t p", positive=True)
    signs = sp.Matrix(walsh_signs(8).astype(int))
    a = sp.diag(0, t*t, 2, 2-t*t, 1, 1, 1, 1)
    a[1, 2] = a[2, 1] = t
    a[4, 7] = a[7, 4] = -t
    c = sp.simplify(signs * a * signs.T / 8)
    assert all(sp.simplify(c[i, i] - 1) == 0 for i in range(8))
    assert sp.simplify(c * sp.ones(8, 1)) == sp.zeros(8, 1)
    assert c == c.T
    # Inversion permutes physical modes i -> i xor 3.
    assert all(sp.simplify(c[i, j] - c[i ^ 3, j ^ 3]) == 0
               for i in range(8) for j in range(8))
    block = sp.Matrix([[t*t, t], [t, 2]])
    x = sp.simplify(block.inv() * sp.Matrix([0, 1]))
    assert x == sp.Matrix([-1/t, 1])
    assert sp.simplify((x.T*x)[0] - (1 + 1/t**2)) == 0
    kp = sp.factor((block + p*sp.eye(2)).inv()[1, 1])
    kp_expected = (p+t*t)/(p*p+(2+t*t)*p+t*t)
    assert sp.simplify(kp-kp_expected) == 0
    # Independent exact solve in original physical coordinates at a rational t.
    current = signs[:, 2] / sp.sqrt(8)
    energy = sp.ones(8, 1) / sp.sqrt(8)
    rational_c = c.subs(t, sp.Rational(1, 2))
    physical_x = (rational_c + energy*energy.T).inv() * current
    exact_dc = sp.simplify((current.T * physical_x)[0])
    exact_m2 = sp.simplify((physical_x.T * physical_x)[0])
    assert exact_dc == 1 and exact_m2 == 5
    # The complete three-current DC tensor is unchanged as well.
    currents_basis = sp.zeros(8, 3)
    for j, index in enumerate((2, 5, 6)):
        currents_basis[index, j] = 1
    lifted_a = a + sp.diag(1, 0, 0, 0, 0, 0, 0, 0)
    tensor_solution = sp.simplify(lifted_a.inv()*currents_basis)
    dc_tensor = sp.simplify(currents_basis.T*tensor_solution)
    m2_tensor = sp.simplify(tensor_solution.T*tensor_solution)
    assert dc_tensor == sp.eye(3)
    assert m2_tensor == sp.diag(1+1/t**2, 1, 1)
    physical_currents = signs.extract(range(8), (2, 5, 6))/sp.sqrt(8)
    physical_tensor_solution = (rational_c+energy*energy.T).inv()*physical_currents
    assert sp.simplify(physical_currents.T*physical_tensor_solution) == sp.eye(3)
    assert sp.simplify(physical_tensor_solution.T*physical_tensor_solution) == sp.diag(5, 1, 1)
    # Four-mode, equal-DC pair: exact rational physical matrices.
    ca = sp.Rational(4, 3)*sp.eye(4) - sp.ones(4, 4)/3
    cb = sp.Matrix([[1, -sp.Rational(1, 2), 0, -sp.Rational(1, 2)],
                    [-sp.Rational(1, 2), 1, -sp.Rational(1, 2), 0],
                    [0, -sp.Rational(1, 2), 1, -sp.Rational(1, 2)],
                    [-sp.Rational(1, 2), 0, -sp.Rational(1, 2), 1]])
    b4 = sp.Matrix([1, 0, 0, -1])
    ea = sp.ones(4, 4)/4
    graph_results = []
    for operator in (ca, cb):
        solution = (operator+ea).inv()*b4
        dc = sp.simplify((b4.T*solution)[0])
        memory = sp.simplify((solution.T*solution)[0]/dc)
        graph_results.append({"dc": str(dc), "memory": str(memory)})
    assert graph_results == [{"dc": "3/2", "memory": "3/4"},
                             {"dc": "3/2", "memory": "5/6"}]
    return {"fixed_diagonal": True, "energy_conservation": True,
            "inversion_symmetry": True, "response": str(kp),
            "dc": "1", "memory": "1 + t**(-2)", "dc_tensor": str(dc_tensor), "m2_tensor": str(m2_tensor),
            "exact_physical_solve_t_half": {"dc": str(exact_dc), "m2": str(exact_m2)},
            "graph_pair": graph_results}


def high_precision(t_text: str) -> dict:
    mp.mp.dps = 80
    t = mp.mpf(t_text)
    a = mp.diag([0, t*t, 2, 2-t*t, 1, 1, 1, 1])
    a[1, 2] = a[2, 1] = t
    a[4, 7] = a[7, 4] = -t
    signs = mp.matrix(walsh_signs(8).astype(int).tolist())
    c = signs*a*signs.T/8
    e = mp.matrix([1/mp.sqrt(8)]*8)
    b = signs[:, 2]/mp.sqrt(8)
    x = mp.lu_solve(c+e*e.T, b)
    dc = (b.T*x)[0]
    tau = (x.T*x)[0]/dc
    exact_tau = 1+1/(t*t)
    residual = mp.norm(c*x-b)/mp.norm(b)
    return {"t": t_text, "decimal_precision": 80,
            "dc": mp.nstr(dc, 30), "memory": mp.nstr(tau, 30),
            "dc_abs_error": mp.nstr(abs(dc-1), 12),
            "memory_relative_error": mp.nstr(abs(tau/exact_tau-1), 12),
            "relative_residual": mp.nstr(residual, 12)}


def main() -> None:
    result = {"environment": {"python": platform.python_version(),
                               "numpy": np.__version__, "sympy": sp.__version__,
                               "mpmath": mp.__version__},
              "symbolic": exact_symbolic_checks()}
    b = np.array([1.0, 0.0, 0.0, -1.0])
    graph_pair = []
    for name, operator in (("complete", graph4(1/3, 1/3, 1/3)),
                            ("cycle", graph4(0.5, 0, 0.5))):
        record = {"name": name, **response_dense(operator, b)}
        record["spectrum"] = np.linalg.eigvalsh(operator).tolist()
        record["response_at_p_1"] = float(b @ np.linalg.solve(operator+np.eye(4), b))
        graph_pair.append(record)
    result["same_dc_graph_pair"] = graph_pair
    # Fixed degree graphs with a current that observes the arbitrarily slow mode.
    slow_graph = []
    bslow = np.array([1, -1, 1, -1], dtype=float)/2
    for delta in (1.0, 0.1, 0.01, 0.001):
        operator = graph4(0, 1-delta/2, delta/2)
        slow_graph.append({"delta": delta, **response_dense(operator, bslow),
                           "spectrum": np.linalg.eigvalsh(operator).tolist()})
    result["graph_varying_dc_and_slow_rate"] = slow_graph
    # Same graph family, orthogonal current: the hidden slow rate is unobserved.
    hidden_graph = []
    bfast = np.array([1, 1, -1, -1], dtype=float)/2
    for delta in (1.0, 0.1, 0.01, 0.001):
        operator = graph4(0, 1-delta/2, delta/2)
        hidden_graph.append({"delta": delta, **response_dense(operator, bfast),
                             "spectrum": np.linalg.eigvalsh(operator).tolist()})
    result["graph_identical_response_hidden_slow_rate"] = hidden_graph
    # Eigensolver and dense inverse checks from the original physical matrix.
    doubles = []
    for t in (0.75, 0.5, 0.2, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
        operator, current, h = family8_matrix(t)
        record = {"t": t}
        exact_tau = 1+1/t**2
        slow_exact = 2*t*t/(2+t*t+np.sqrt(4+t**4))
        record["slow_eigenvalue_exact"] = slow_exact
        record["memory_exact"] = exact_tau
        record["max_positive_offdiagonal"] = float(np.max(operator-np.diag(np.diag(operator))))
        try:
            record["dense"] = response_dense(operator, current)
            record["dense"]["dc_abs_error"] = abs(record["dense"]["dc"]-1)
            record["dense"]["memory_relative_error"] = abs(record["dense"]["memory"]/exact_tau-1)
        except np.linalg.LinAlgError as error:
            record["dense"] = {"failure": str(error)}
        q = h[:, 1:]
        eigenvalues, eigenvectors = np.linalg.eigh(q.T @ operator @ q)
        weights = (eigenvectors.T @ (q.T @ current))**2
        dc = float(np.sum(weights/eigenvalues))
        memory = float(np.sum(weights/eigenvalues**2)/dc)
        record["spectral"] = {"min_eigenvalue": float(eigenvalues.min()), "dc": dc,
                               "memory": memory, "dc_abs_error": abs(dc-1),
                               "memory_relative_error": abs(memory/exact_tau-1)}
        # Direct resolvent at a fixed frequency stays benign when t is small.
        p = 0.3
        actual = float(current @ np.linalg.solve(operator+p*np.eye(8), current))
        expected = (p+t*t)/(p*p+(2+t*t)*p+t*t)
        record["resolvent_p_point3_abs_error"] = abs(actual-expected)
        doubles.append(record)
    result["family8_double_precision"] = doubles
    result["family8_high_precision"] = [high_precision(t) for t in ("0.01", "0.000001", "1e-10", "1e-20")]
    # Report errors, do not silently accept ill-conditioned double-precision data.
    trusted = [record for record in doubles if record["t"] >= 1e-3]
    assert max(record["dense"]["dc_abs_error"] for record in trusted) < 1e-8
    assert max(record["dense"]["memory_relative_error"] for record in trusted) < 1e-8
    assert max(record["resolvent_p_point3_abs_error"] for record in doubles) < 1e-12
    target = ROOT / "constructive_collision_results.json"
    target.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(target), "graph_pair": graph_pair,
                      "family8_doubles": doubles,
                      "high_precision": result["family8_high_precision"]}, indent=2))


if __name__ == "__main__":
    main()
