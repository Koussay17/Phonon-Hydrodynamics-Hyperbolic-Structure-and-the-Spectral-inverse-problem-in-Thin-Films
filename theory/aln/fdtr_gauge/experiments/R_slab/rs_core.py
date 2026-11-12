"""Core of the full-operator kinetic slab (campaign 20261111-fdtr-gauge-memory, branch R-S).

Model (ASSUMPTIONS). Linearised phonon BTE in a c-oriented film 0 < z < d, frequency domain e^{pt}, symmetrised
variables f (f = theta e at equilibrium at temperature theta; e_i = sqrt(c_i), c_i mode capacity in J/(m^3 K)):
    V f' = -(p + Chat) f,   V = diag(v_z),   Chat = 4 pi 1e12 C(g) [1/s],   Chat = Chat^T >= 0,   Chat e = 0,
C(g) the orbit-tied conserving three-phonon event operator of fc3-memory (aln_geometry / fc3class conventions,
re-derived from the raw event list as in R_A_audit/au3_event775.py). No Matthiessen regulator.
Walls diffuse and black: f_i(0) = T_A e_i for v_z,i > 0, f_i(d) = T_B e_i for v_z,i < 0 (T_A, T_B the wall emission
temperatures; same convention as the RTA slab au6_slab.py).
Heat flux q = b.f (b = V e), energy density U = e.f, energy temperature theta = U / C_tot.

Parity. P: (q_x, q_y, q_z, j) -> (q_x, q_y, -q_z, j) (time reversal o C2z; stored grid map 15): P V P = -V,
P Chat P = Chat, P e = e (orbit tying makes the last two exact). Moving modes come in pairs (+, -); the plane q_z = 0
(v_z = 0 exactly) holds the local modes. With s = (f+ + f-)/sqrt2, o = (f+ - f-)/sqrt2:
    v o' = -[(p + C_e) s + sqrt2 C_p0 f_0],   0 = sqrt2 C_0p s + (p + C_00) f_0,   v s' = -(p + C_o) o,
C_e = C_pp + C_pm, C_o = C_pp - C_pm. Eliminating f_0:  v o' = -A_e s,  A_e = p + C_e - 2 C_p0 (p + C_00)^{-1} C_0p.
Modes s = phi e^{-lam z}, o = psi e^{-lam z}:  A_e phi = mu B phi,  B = v A_o^{-1} v,  mu = lam^2,
psi = lam A_o^{-1} v phi  (A_o = p + C_o). Exact energy structure: A_e e_p = p w_E, w_E = e_p - C_p0 (p + C_00)^{-1} e_0.

Solvers (all exact for the discrete mode set, no z discretisation except the doubling check):
  slab_paired   : eig of M = v^-1 A_o v^-1 A_e (LAPACK geev), diffusive eigenpair polished by Rayleigh-quotient
                  inverse iteration on the pencil (A_e, B); production route.
  slab_unpaired : first-order eigenproblem V^-1 A_S on all moving modes (no parity used, eigenvalues lam not lam^2).
  slab_doubling : adding-doubling (invariant imbedding) z-march: Taylor propagator of a layer of thickness d/2^K,
                  K doublings of the reflection/transmission pair; no eigendecomposition.
  slab_paired_real_svd : real p only; lam = 1/sigma of Z = R_O^-T ... (well-conditioned check of the eigenvalues).
"""
from __future__ import annotations

import os
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import scipy.sparse as sps  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c")
EVENTS = {"m775": FC3 / "m775x/events_m775_s0.1.npz", "m553": FC3 / "m553/events_m553_s0.1.npz"}
PPDIR = {("m775", "film2"): FC3 / "m775x/film2_pp", ("m775", "filmw"): FC3 / "m775x/filmw_pp",
         ("m553", "filmw"): FC3 / "m553/filmw_pp"}
WIT12 = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
        [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
WITW = ["cal1_fL100_max", "cal1_fL100_min"] + [f"orb0.1_fL{L}_{s}" for L in (30, 100, 300, 1000) for s in ("max", "min")]
KB_SI = 1.380649e-23
FOURPI = 4.0 * np.pi
PARITY_MAP = 15


# ----------------------------------------------------------------------------------------------- geometry
def geometry(mesh):
    Z = np.load(EVENTS[mesh])
    f = Z["freqs"]
    Nq, nb = f.shape
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(valid.sum())
    eps = f.ravel()[valid]
    kT = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    vz = Z["gv"].reshape(Nq * nb, 3)[valid][:, 2] * 100.0                 # THz*A -> m/s
    x = eps / kT
    nocc = 1.0 / np.expm1(x)
    dsc = np.sqrt(nocc * (1 + nocc))
    cap = KB_SI * x * x * nocc * (1 + nocc) / (Nq * float(Z["volume"]) * 1e-30)      # J/(m^3 K) per mode
    pf, af, bf = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    P, A, B = newidx[pf], newidx[af], newidx[bf]
    assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
    rep = A == B
    lp = np.sqrt((eps[A] + eps[B]) / eps[P])
    ld = 1.0 / lp
    vP = -lp / dsc[P]
    vA = np.where(rep, 2 * ld, ld) / dsc[A]
    vB = np.where(rep, 0.0, ld / dsc[B])
    s = np.sinh(x / 2)
    kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (2 * Nq * float(Z["KB"]) * float(Z["T"]) ** 2)
    # orbits: connected components of the image graph under the 24 stored grid maps
    NF = Nq * nb
    key = (pf * NF + np.minimum(af, bf)) * NF + np.maximum(af, bf)
    order = np.argsort(key)
    skey = key[order]
    m = key.size
    rows, cols = [np.arange(m)], [np.arange(m)]
    n_missing = 0
    for R in Z["rot_maps"]:
        def img(full):
            return R[full // nb] * nb + full % nb
        ip_, ia_, ib_ = img(pf), img(af), img(bf)
        k2 = (ip_ * NF + np.minimum(ia_, ib_)) * NF + np.maximum(ia_, ib_)
        pos = np.clip(np.searchsorted(skey, k2), 0, m - 1)
        found = skey[pos] == k2
        n_missing += int((~found).sum())
        rows.append(np.arange(m)[found])
        cols.append(order[pos[found]])
    Gr = sps.coo_matrix((np.ones(sum(r.size for r in rows)), (np.concatenate(rows), np.concatenate(cols))), shape=(m, m))
    n_orb, orb = connected_components(Gr, directed=False)
    # parity permutation of the valid modes
    R = Z["rot_maps"][PARITY_MAP]
    addr, msh = Z["addr"], Z["mesh"]
    tgt = np.stack([addr[:, 0] % msh[0], addr[:, 1] % msh[1], (-addr[:, 2]) % msh[2]], 1)
    assert np.array_equal(addr[R] % msh, tgt), "map 15 is not (qx, qy, -qz)"
    full = np.arange(Nq * nb)
    perm_full = R[full // nb] * nb + full % nb
    perm = newidx[perm_full[valid]]
    assert (perm >= 0).all()
    zero = np.abs(vz) < 1e-12
    ip = np.where((vz > 0) & ~zero)[0]
    im = perm[ip]
    i0 = np.where(zero)[0]
    assert np.all(vz[im] < 0) and np.allclose(vz[im], -vz[ip], rtol=1e-12, atol=1e-9)
    assert np.allclose(eps[im], eps[ip], rtol=1e-13) and np.all(perm[i0] == i0)
    return {"mesh": mesh, "eps": eps, "vz": vz, "cap": cap, "e": np.sqrt(cap), "P": P, "A": A, "B": B, "rep": rep,
            "vP": vP, "vA": vA, "vB": vB, "s": s, "kf": kf, "orb": orb, "n_orb": int(n_orb),
            "n_missing_images": n_missing, "n": int(eps.size), "m": int(m), "perm": perm, "ip": ip, "im": im, "i0": i0,
            "C_tot": float(cap.sum()), "Nq": Nq}


def load_pp(mesh, group, name):
    Pz = np.load(PPDIR[(mesh, group)] / f"pp_r01_{name}.npz")
    return Pz["pp"], float(Pz["conv"]), Pz["g0w"]


def rates(G, pp, conv, g0w):
    P, A, B, rep, s = G["P"], G["A"], G["B"], G["rep"], G["s"]
    W = np.where(rep, conv * g0w / (8 * s[P] * s[A] * s[A]), conv * g0w / (4 * s[P] * s[A] * s[np.where(rep, A, B)]))
    g = W * pp
    sums = np.bincount(G["orb"], weights=g, minlength=G["n_orb"])
    cnt = np.bincount(G["orb"], minlength=G["n_orb"])
    return (sums / cnt)[G["orb"]]


def assemble(G, g):
    """Dense symmetrised operator in 1/s (physical rate = 4 pi x operator in THz)."""
    n = G["n"]
    idx = [G["P"], G["A"], G["B"]]
    val = [G["vP"], G["vA"], G["vB"]]
    r_, c_, w_ = [], [], []
    for i in range(3):
        for j in range(3):
            r_.append(idx[i]); c_.append(idx[j]); w_.append(g * val[i] * val[j])
    C = sps.coo_matrix((np.concatenate(w_), (np.concatenate(r_), np.concatenate(c_))), shape=(n, n)).toarray()
    C = 0.5 * (C + C.T)
    return C * (FOURPI * 1e12)


def operator(G, mesh_group, name):
    pp, conv, g0w = load_pp(G["mesh"], mesh_group, name)
    return assemble(G, rates(G, pp, conv, g0w))


def rta_operator(G, C):
    """Energy-conserving RTA on the diagonal of C: relaxation of each mode to one local temperature."""
    r = np.diag(C).copy()
    e = G["e"]
    re = r * e
    return np.diag(r) - np.outer(re, re) / (e @ re)


def regulated(G, C, L, variant="QDQ"):
    """C + regulator |v_z|/L (1/s) in its energy-conserving completions (au3: E-QDQ, E-phys)."""
    d = np.abs(G["vz"]) / L
    e = G["e"]
    eh = e / np.linalg.norm(e)
    if variant == "QDQ":
        h = d * eh
        return C + np.diag(d) - np.outer(eh, h) - np.outer(h, eh) + (eh @ h) * np.outer(eh, eh)
    g = d * e
    return C + np.diag(d) - np.outer(g, g) / (e @ g)


def kinetic_times(G, Ch):
    """K (W/mK), tau_R, tau_l (s) of a conserving operator with kernel span{e} (KKT bordered solves)."""
    e, vz = G["e"], G["vz"]
    n = e.size
    K_ = np.zeros((n + 1, n + 1))
    K_[:n, :n] = Ch
    K_[:n, n] = e
    K_[n, :n] = e
    lu = sla.lu_factor(K_, check_finite=False)
    b = vz * e
    x = sla.lu_solve(lu, np.append(b, 0.0), check_finite=False)[:n]
    K = float(b @ x)
    eh = e / np.linalg.norm(e)
    u = vz * x
    u = u - eh * (eh @ u)
    y = sla.lu_solve(lu, np.append(u, 0.0), check_finite=False)[:n]
    M = float(u @ y)
    return {"K": K, "tauR": float(x @ x) / K, "taul": float(e @ e) * M / K ** 2,
            "delta": float(x @ x) / K - float(e @ e) * M / K ** 2}


# ----------------------------------------------------------------------------------------------- paired blocks
def blocks(G, C):
    ip, im, i0 = G["ip"], G["im"], G["i0"]
    Cpp = C[np.ix_(ip, ip)]
    Cpm = C[np.ix_(ip, im)]
    sym = {"pp_vs_mm": float(np.max(np.abs(C[np.ix_(im, im)] - Cpp)) / np.max(np.abs(Cpp))),
           "pm_vs_mp": float(np.max(np.abs(C[np.ix_(im, ip)] - Cpm)) / np.max(np.abs(Cpp))),
           "p0_vs_m0": float(np.max(np.abs(C[np.ix_(im, i0)] - C[np.ix_(ip, i0)])) / np.max(np.abs(Cpp))) if i0.size else 0.0,
           "energy_residual": float(np.linalg.norm(C @ G["e"]) / np.linalg.norm(np.diag(C) * G["e"]))}
    e = G["e"]
    return {"Ce": Cpp + Cpm, "Co": Cpp - Cpm, "Cpp": Cpp, "Cpm": Cpm, "Cp0": C[np.ix_(ip, i0)], "C00": C[np.ix_(i0, i0)],
            "ep": e[ip], "e0": e[i0], "v": G["vz"][ip], "C_tot": float(e @ e), "sym": sym}


def blocks_from_parts(ep, e0, v, Cpp, Cpm, Cp0, C00):
    """Paired blocks for synthetic tests (Cpp, Cpm symmetric; Cp0 (n x n0))."""
    return {"Ce": Cpp + Cpm, "Co": Cpp - Cpm, "Cpp": Cpp, "Cpm": Cpm, "Cp0": Cp0, "C00": C00, "ep": ep, "e0": e0,
            "v": v, "C_tot": float(2 * ep @ ep + e0 @ e0), "sym": {}}


def _schur_even(Bk, p):
    """A_e (energy-structure corrected), w_E, X = (p + C00)^-1 C0p."""
    n = Bk["ep"].size
    n0 = Bk["e0"].size
    Ae = Bk["Ce"].astype(complex) + p * np.eye(n)
    ep, e0 = Bk["ep"], Bk["e0"]
    if n0:
        lu00 = sla.lu_factor(Bk["C00"] + p * np.eye(n0), check_finite=False)
        X = sla.lu_solve(lu00, Bk["Cp0"].T.astype(complex), check_finite=False)
        Ae = Ae - 2.0 * (Bk["Cp0"] @ X)
        wE = ep - Bk["Cp0"] @ sla.lu_solve(lu00, e0.astype(complex), check_finite=False)
    else:
        X = np.zeros((0, n), complex)
        wE = ep.astype(complex)
    Ae = 0.5 * (Ae + Ae.T)
    # enforce A_e e_p = p w_E exactly (symmetric rank-two correction; removes rounding of the Schur complement in the
    # conserved direction, ~eps * cond(C00) * |C_p0| against p)
    r = p * wE - Ae @ ep
    ne = np.linalg.norm(ep)
    eh = ep / ne
    corr = float(np.linalg.norm(r) / max(np.linalg.norm(p * wE), 1e-300))
    Ae = Ae + (np.outer(r, eh) + np.outer(eh, r) - (eh @ r) * np.outer(eh, eh)) / ne
    return Ae, wE, X, corr


def wall_operator(wall, ep, bp):
    """Wall reflection operator W (incident -> re-emitted, paired labelling) and emission weight tau_w.
    ("black",): W = 0, tau_w = 1 (diffuse black wall, the audit's convention);
    ("specular", r): fraction r specularly reflected (P is the specular map of a c-plane wall), rest absorbed and
                     re-emitted at the reservoir temperature: W = r I, tau_w = 1 - r;
    ("diffuse", t): gray transmission t; the reflected part re-emitted diffusely at the temperature that conserves its
                    energy flux: W = (1 - t) e_p b_p^T / (b_p.e_p), tau_w = t."""
    kind = wall[0]
    if kind == "black":
        return None, 1.0
    if kind == "specular":
        r = float(wall[1])
        return (lambda X: r * X), 1.0 - r
    if kind == "diffuse":
        t = float(wall[1])
        c = (1.0 - t) / float(bp @ ep)
        return (lambda X: c * (np.outer(ep, bp @ X) if X.ndim == 2 else ep * (bp @ X))), t
    raise KeyError(kind)


def _bc_solve(Phi, Psi, lam, d, ep, bp, uS, C_tot, p, wall=("black",)):
    """Two-port from the modal basis: returns Y, Th, diagnostics. Walls at z = 0 and d (mirror images):
    (I - W) s + (I + W) o = sqrt2 tau_w T_A e_p at z = 0;  (I - W) s - (I + W) o = sqrt2 tau_w T_B e_p at z = d."""
    W, tw = wall_operator(wall, ep, bp)
    if W is None:
        Phs, Pss = Phi, Psi
    else:
        Phs, Pss = Phi - W(Phi), Psi + W(Psi)
    E = np.exp(-lam * d)
    Pp = (Phs + Pss) + (Phs - Pss) * E[None, :]
    Pm = (Phs + Pss) - (Phs - Pss) * E[None, :]
    luP = sla.lu_factor(Pp, check_finite=False)
    luM = sla.lu_factor(Pm, check_finite=False)
    rt2 = np.sqrt(2.0)
    out = []
    for TA, TB in ((1.0, 0.0), (0.0, 1.0)):
        sg = sla.lu_solve(luP, rt2 * tw * (TA + TB) * ep.astype(complex), check_finite=False)
        dl = sla.lu_solve(luM, rt2 * tw * (TA - TB) * ep.astype(complex), check_finite=False)
        a, b = 0.5 * (sg + dl), 0.5 * (sg - dl)
        o0 = Psi @ (a - E * b)
        od = Psi @ (E * a - b)
        s0 = Phi @ (a + E * b)
        sd = Phi @ (E * a + b)
        q0 = rt2 * (bp @ o0)
        qd = rt2 * (bp @ od)
        # integral of s over the film: (1 - E)/lam, with the small-lam limit d
        fac = np.where(np.abs(lam * d) < 1e-8, d * (1 - lam * d / 2), (1 - E) / np.where(lam == 0, 1, lam))
        intS = Phi @ ((a + b) * fac)
        U = uS @ intS
        if W is None:
            r0 = s0 + o0 - rt2 * TA * ep
            r1 = sd - od - rt2 * TB * ep
        else:
            r0 = (s0 - W(s0)) + (o0 + W(o0)) - rt2 * tw * TA * ep
            r1 = (sd - W(sd)) - (od + W(od)) - rt2 * tw * TB * ep
        out.append((q0, qd, U, (uS @ s0) / C_tot, (uS @ sd) / C_tot,
                    float(np.max(np.abs(r0))), float(np.max(np.abs(r1)))))
    (q0A, qdA, UA, t0A, tdA, bcA0, bcA1), (q0B, qdB, UB, t0B, tdB, bcB0, bcB1) = out
    Y = np.array([[q0A, q0B], [qdA, qdB]])
    Th = np.array([[t0A, t0B], [tdA, tdB]])
    ebal = max(abs((q0A - qdA) - p * UA), abs((q0B - qdB) - p * UB)) / max(abs(q0A), 1e-300)
    diag = {"energy_balance_rel": float(ebal), "reciprocity_rel": float(abs(Y[0, 1] + Y[1, 0]) / abs(Y[0, 1])),
            "symmetry_rel": float(abs(Y[0, 0] + Y[1, 1]) / abs(Y[0, 0])),
            "bc_residual_rel": float(max(bcA0, bcA1, bcB0, bcB1) / (np.sqrt(2) * np.max(ep))),
            "pivot_ratio_Pminus": float(np.abs(np.diag(luM[0])).max() / np.abs(np.diag(luM[0])).min())}
    return Y, Th, diag


def slab_paired(Bk, p, ds, polish=True, return_modes=False, walls=None):
    """Production solver. ds: iterable of thicknesses (m). Returns list of (Y, Th, diag) over ds (walls None: black)
    or, with walls = [wall specs], a list over walls of such lists; and the eigen diagnostics."""
    ep, v = Bk["ep"], Bk["v"]
    n = ep.size
    Ao = Bk["Co"].astype(complex) + p * np.eye(n)
    luo = sla.lu_factor(Ao, check_finite=False)
    Ae, wE, X, corr = _schur_even(Bk, p)
    M = (Ao / v[:, None]) @ (Ae / v[:, None])
    mu, Phi = sla.eig(M, check_finite=False, overwrite_a=True)
    del M
    bp = v * ep
    xo = sla.lu_solve(luo, bp.astype(complex), check_finite=False)
    Kp = 2 * (bp @ xo)
    mu_est = p * (2 * (ep @ wE)) / Kp
    j = int(np.argmin(np.abs(mu - mu_est)))
    info = {"schur_energy_correction_rel": corr, "mu_d_geev": complex(mu[j]), "mu_d_estimate": complex(mu_est)}
    Bm = v[:, None] * sla.lu_solve(luo, np.diag(v).astype(complex), check_finite=False)
    Bm = 0.5 * (Bm + Bm.T)
    if polish:
        ph = Phi[:, j].copy()
        mu_j = mu[j]
        hist = []
        for _ in range(6):
            with warnings.catch_warnings():          # the shifted pencil is singular by design (inverse iteration)
                warnings.simplefilter("ignore", sla.LinAlgWarning)
                y = sla.solve(Ae - mu_j * Bm, Bm @ ph, check_finite=False)
            ph = y / np.linalg.norm(y)
            mu_new = (ph @ (Ae @ ph)) / (ph @ (Bm @ ph))
            hist.append(abs(mu_new - mu_j) / abs(mu_new))
            mu_j = mu_new
            if hist[-1] < 1e-15:
                break
        mu[j] = mu_j
        Phi[:, j] = ph
        info["polish_history"] = hist
        info["mu_d_polished"] = complex(mu_j)
        info["polish_shift_rel"] = float(abs(mu_j - info["mu_d_geev"]) / abs(mu_j))
    lam = np.sqrt(mu.astype(complex))
    lam = np.where(lam.real < 0, -lam, lam)
    Psi = lam[None, :] * sla.lu_solve(luo, v[:, None] * Phi, check_finite=False)
    nrm = np.sqrt(np.sum(np.abs(Phi) ** 2, 0) + np.sum(np.abs(Psi) ** 2, 0))
    Phi, Psi = Phi / nrm, Psi / nrm
    # eigen-residuals on the pencil (A_e, B) (backward-error form)
    nA, nB = np.linalg.norm(Ae, 1), np.linalg.norm(Bm, 1)
    R_ = Ae @ Phi - (Bm @ Phi) * mu[None, :]
    res = np.linalg.norm(R_, axis=0) / ((nA + np.abs(mu) * nB) * np.linalg.norm(Phi, axis=0))
    info.update({"eig_backward_max": float(res.max()), "eig_backward_diffusive": float(res[j]),
                 "lam_min_abs": float(np.abs(lam).min()), "lam_max_abs": float(np.abs(lam).max()), "j_diff": j})
    uS = np.sqrt(2.0) * wE
    if walls is None:
        out = [_bc_solve(Phi, Psi, lam, d, ep, bp, uS, Bk["C_tot"], p) for d in ds]
    else:
        out = [[_bc_solve(Phi, Psi, lam, d, ep, bp, uS, Bk["C_tot"], p, wall=w) for d in ds] for w in walls]
    if return_modes:
        return out, info, {"lam": lam, "Phi": Phi, "Psi": Psi}
    return out, info


# ----------------------------------------------------------------------------------------------- independent checks
def slab_unpaired(G, C, p, ds):
    """First-order eigenproblem on all moving modes (no parity), A_S g = lam V g; two-port by direct BC solve."""
    vz, e = G["vz"], G["e"]
    mov = np.abs(vz) >= 1e-12
    i0 = np.where(~mov)[0]
    im_ = np.where(mov)[0]
    A = C[np.ix_(im_, im_)].astype(complex) + p * np.eye(im_.size)
    if i0.size:
        lu00 = sla.lu_factor(C[np.ix_(i0, i0)] + p * np.eye(i0.size), check_finite=False)
        X = sla.lu_solve(lu00, C[np.ix_(i0, im_)].astype(complex), check_finite=False)
        A = A - C[np.ix_(im_, i0)] @ X
        wE = e[im_] - X.T @ e[i0]                    # U = wE . f_mov  (f_0 = -X f_mov)
    else:
        wE = e[im_].astype(complex)
    V = vz[im_]
    lam, Gm = sla.eig(A / V[:, None], check_finite=False)
    pos = lam.real > 0
    assert pos.sum() == (V > 0).sum(), (pos.sum(), (V > 0).sum())
    lp_, ln_ = lam[pos], lam[~pos]
    Gp, Gn = Gm[:, pos], Gm[:, ~pos]
    Gp = Gp / np.linalg.norm(Gp, axis=0)
    Gn = Gn / np.linalg.norm(Gn, axis=0)
    up, dn = V > 0, V < 0
    ev = e[im_]
    b = V * ev
    out = []
    for d in ds:
        Ep, En = np.exp(-lp_ * d), np.exp(ln_ * d)          # e^{-lam d} (Re>0), e^{+lam d} (Re<0)
        # f(z) = Gp a e^{-lp z} + Gn c e^{-ln (z - d)}
        M = np.block([[Gp[up], Gn[up] * En[None, :]], [Gp[dn] * Ep[None, :], Gn[dn]]])
        lu = sla.lu_factor(M, check_finite=False)
        res = []
        for TA, TB in ((1.0, 0.0), (0.0, 1.0)):
            rhs = np.concatenate([TA * ev[up], TB * ev[dn]]).astype(complex)
            sol = sla.lu_solve(lu, rhs, check_finite=False)
            a, c = sol[:lp_.size], sol[lp_.size:]
            f0 = Gp @ a + Gn @ (c * En)
            fd = Gp @ (a * Ep) + Gn @ c
            res.append((b @ f0, b @ fd, (wE @ f0) / G["C_tot"], (wE @ fd) / G["C_tot"]))
        Y = np.array([[res[0][0], res[1][0]], [res[0][1], res[1][1]]])
        Th = np.array([[res[0][2], res[1][2]], [res[0][3], res[1][3]]])
        out.append((Y, Th))
    return out


def slab_doubling(Bk, p, d, K=23, order=10, return_RT=False, wall=("black",)):
    """Adding-doubling z-march (no eigendecomposition). Layer h = d / 2^K, propagator exp(h G) by Taylor series of
    the given order on [f+; f-] (G = [[-v^-1 At, -v^-1 Bt], [v^-1 Bt, v^-1 At]], At = p + C_pp - S, Bt = C_pm - S,
    S = C_p0 (p + C00)^-1 C_0p); reflection R and the defect Tt = I - T are doubled K times:
        X = (I - R^2)^-1,  R' = R + T R X T,  T' = T X T,  Tt' = 2 Tt - Tt^2 - T R^2 X T  (T = I - Tt)."""
    ep, v = Bk["ep"], Bk["v"]
    n = ep.size
    n0 = Bk["e0"].size
    if n0:
        lu00 = sla.lu_factor(Bk["C00"] + p * np.eye(n0), check_finite=False)
        X0 = sla.lu_solve(lu00, Bk["Cp0"].T.astype(complex), check_finite=False)
        S = Bk["Cp0"] @ X0
        wfac = Bk["Cp0"] @ sla.lu_solve(lu00, Bk["e0"].astype(complex), check_finite=False)
    else:
        S = np.zeros((n, n), complex)
        wfac = np.zeros(n, complex)
    At = Bk["Cpp"] + p * np.eye(n) - S
    Bt = Bk["Cpm"] - S
    Gm = np.block([[-At / v[:, None], -Bt / v[:, None]], [Bt / v[:, None], At / v[:, None]]])
    if K is None:                              # layer with |h G|_1 <= 0.05
        K = int(np.ceil(np.log2(max(d * np.linalg.norm(Gm, 1) / 0.05, 1.0))))
    h = d / 2 ** K
    Gm *= h
    gnorm = float(np.linalg.norm(Gm, 1))
    # Taylor series of exp(Gm) - I (keep the defect from the identity for precision)
    term = Gm.copy()
    Ph = Gm.copy()
    for k in range(2, order + 1):
        term = term @ Gm / k
        Ph += term
    tail = float(np.linalg.norm(term, 1))
    del Gm, term
    # Ph = exp(hG) - I ; blocks of Phi = I + Ph
    P11, P12, P21, P22 = Ph[:n, :n], Ph[:n, n:], Ph[n:, :n], Ph[n:, n:]
    I = np.eye(n)
    lu22 = sla.lu_factor(I + P22, check_finite=False)
    R = -sla.lu_solve(lu22, P21, check_finite=False)                      # f-(0) = R f+(0) + Tu f-(h)
    Tu = sla.lu_solve(lu22, I.astype(complex), check_finite=False)
    Ttd = -(P11 + P12 @ R)                                                # I - T_down
    Ttu = sla.lu_solve(lu22, P22, check_finite=False)                     # I - T_up = (I + P22)^-1 P22
    Rb = P12 @ Tu
    sym = {"R_top_vs_bottom": float(np.max(np.abs(R - Rb)) / max(np.max(np.abs(R)), 1e-300)),
           "T_down_vs_up": float(np.max(np.abs(Ttd - Ttu)) / max(np.max(np.abs(Ttu)), 1e-300))}
    Tt = 0.5 * (Ttd + Ttu)
    T = I - Tt
    R = 0.5 * (R + Rb)
    del Ph, P11, P12, P21, P22
    for _ in range(K):
        R2 = R @ R
        luX = sla.lu_factor(I - R2, check_finite=False)
        XT = sla.lu_solve(luX, T, check_finite=False)                     # X T
        RXT = R @ XT
        Rn = R + T @ RXT
        Ttn = 2 * Tt - Tt @ Tt - T @ (R2 @ XT)
        R, Tt = Rn, Ttn
        T = I - Tt
    # walls: f+(0) = tau_w TA ep + W f-(0), f-(d) = tau_w TB ep + W f+(d);  f-(0) = R f+(0) + T f-(d), f+(d) = T f+(0) + R f-(d)
    bp = v * ep
    wE = ep - wfac
    Wf, tw = wall_operator(wall, ep, bp)
    W = np.zeros((n, n)) if Wf is None else Wf(np.eye(n))
    Msys = np.block([[I - W @ R, -W @ T], [-W @ T, I - W @ R]])
    lus = sla.lu_factor(Msys, check_finite=False)
    out = []
    for TA, TB in ((1.0, 0.0), (0.0, 1.0)):
        sol = sla.lu_solve(lus, np.concatenate([tw * TA * ep, tw * TB * ep]).astype(complex), check_finite=False)
        fp0, fmd = sol[:n], sol[n:]
        fm0 = R @ fp0 + T @ fmd
        fpd = T @ fp0 + R @ fmd
        q0 = bp @ (fp0 - fm0)
        qd = bp @ (fpd - fmd)
        U0 = wE @ (fp0 + fm0)
        Ud = wE @ (fpd + fmd)
        out.append((q0, qd, U0 / Bk["C_tot"], Ud / Bk["C_tot"]))
    Y = np.array([[out[0][0], out[1][0]], [out[0][1], out[1][1]]])
    Th = np.array([[out[0][2], out[1][2]], [out[0][3], out[1][3]]])
    info = {"h": h, "K": K, "order": order, "norm_hG": gnorm, "taylor_tail": tail, **sym}
    if return_RT:
        return Y, Th, info, R, T
    return Y, Th, info


def slab_eigs_real_svd(Bk, p):
    """Real p > 0: lam = singular values of Y = R_e v^-1 R_o^T (A_e = R_e^T R_e, A_o = R_o^T R_o, Cholesky)."""
    assert np.isreal(p) and p > 0
    n = Bk["ep"].size
    Ae, wE, X, corr = _schur_even(Bk, p)
    Ae = Ae.real
    Ao = Bk["Co"] + p * np.eye(n)
    Re_ = sla.cholesky(Ae, lower=False)
    Ro = sla.cholesky(Ao, lower=False)
    Ymat = Re_ @ (Ro.T / Bk["v"][:, None])
    return np.sort(sla.svdvals(Ymat))


def slab_eigs_real_paired(Bk, p):
    n = Bk["ep"].size
    v = Bk["v"]
    Ao = Bk["Co"] + p * np.eye(n)
    Ae, wE, X, corr = _schur_even(Bk, p)
    M = (Ao / v[:, None]) @ (Ae.real / v[:, None])
    mu = sla.eigvals(M)
    return np.sort(np.sqrt(mu.astype(complex)).real)
