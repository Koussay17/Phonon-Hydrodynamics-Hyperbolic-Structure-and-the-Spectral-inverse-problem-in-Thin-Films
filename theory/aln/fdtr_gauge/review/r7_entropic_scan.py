"""Random scan of the 4-parameter rank under the 'entropic' interface temperature (theta = (1 + taul p) T on the GK side),
with the 'local' convention as control. Nondimensional (bs = 1, xi = 1 at the point, frequencies omega xi^2 in 1e-3..1e3).
Double precision; derivatives by central differences in log-parameters (relative step 1e-6); column-normalised SVD.
An exact rank drop would show as s4 ~ 1e-9 or below; generic points give s4 ~ 1e-3..1e-1."""
import json
import numpy as np

rng = np.random.default_rng(20261112)
w = np.logspace(-3, 3, 61)
p = 1j * w


def lnZ(lq, conv):
    xi, b, tR, tl = np.exp(lq)
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(p)
    L = 1 / S
    if conv == "entropic":
        L = L / (1 + tl * p)
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    return np.log((1 + r * E) / ((1 - r * E) * B))


def s4(lq, conv, h=1e-6):
    cols = []
    for j in range(4):
        e = np.zeros(4); e[j] = h
        d = (lnZ(lq + e, conv) - lnZ(lq - e, conv)) / (2 * h)
        cols.append(np.concatenate([d.real, d.imag]))
    M = np.array(cols).T
    nr = np.linalg.norm(M, axis=0)
    nr[nr == 0] = 1.0
    M = M / nr
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[-1] / s[0])


out = {}
for conv in ("local", "entropic"):
    vals = {"interior": [], "resonance": [], "resonance_b=bs": [], "b=bs": []}
    for _ in range(300):
        lb = rng.uniform(-1.5, 1.5)
        lt = rng.uniform(np.log(1e-2), np.log(1e1), size=2)
        vals["interior"].append(s4(np.array([0.0, lb, lt[0], lt[1]]), conv))
        vals["resonance"].append(s4(np.array([0.0, lb, lt[0], lt[0]]), conv))
        vals["resonance_b=bs"].append(s4(np.array([0.0, 0.0, lt[0], lt[0]]), conv))
        vals["b=bs"].append(s4(np.array([0.0, 0.0, lt[0], lt[1]]), conv))
    out[conv] = {k: {"min": float(np.min(v)), "median": float(np.median(v)), "n_below_1e-8": int(np.sum(np.array(v) < 1e-8))}
                 for k, v in vals.items()}
print(json.dumps(out, indent=1))
with open(__file__.replace(".py", ".json"), "w") as fh:
    json.dump(out, fh, indent=1)
