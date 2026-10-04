"""Inventory (paths, sizes, SHA-256, shapes, units) of the data exported for downstream work.

Production setting (phono3py 4.5.0 defaults: Rust kernels, make_r0_average=True, analytic
group velocities incl. Gonze-Lee NAC derivative), tetrahedron, 31x31x17, no isotope unless
stated.  Writes results/saved_data_inventory.json.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

C = common.CAMPAIGN
R = C / "runs"


def entry(p, describe):
    p = Path(p)
    if not p.exists():
        return {"path": str(p), "missing": True}
    e = {"path": str(p), "bytes": p.stat().st_size, "sha256": common.sha256(p), "content": describe}
    if p.suffix == ".npy":
        m = np.load(p, mmap_mode="r")
        e["shape"], e["dtype"] = list(m.shape), str(m.dtype)
    elif p.suffix == ".npz":
        with np.load(p) as z:
            e["arrays"] = {k: list(z[k].shape) for k in z.files}
    elif p.suffix == ".hdf5":
        with h5py.File(p, "r") as h:
            e["datasets"] = {k: list(h[k].shape) for k in h.keys()}
    return e


def main():
    ooc = R / "ooc-m313117-prod"
    inv = {"setting": "production: phono3py 4.5.0 defaults (lang=Rust, make_r0_average=True, analytic gv), "
                      "linear tetrahedron, mesh 31x31x17 (N=16337, 864 irreducible points), NAC Gonze-Lee, no isotope",
           "conventions": {
               "frequency": "THz (ordinary frequency f = omega/2pi)",
               "gamma": "phono3py imaginary self-energy = half linewidth, THz; tau = 1/(4 pi gamma) ps",
               "group_velocity": "THz*Angstrom (x100 -> m/s)",
               "heat_capacity": "eV/K per mode (phonopy constants)",
               "weights": "integer star size of each irreducible point; sum = 16337",
               "assembled_matrix": "irreducible kappa-star representation, index (i,b,a): ir point i (864), band b (12), "
                                   "Cartesian a (3); n=31104. A = sym(P A_w P), A_w[(iba),(jcd)] = w_i w_j "
                                   "[sum_R Omega_{(q_i b),(R q_j c)} R_ad + delta_ij delta_bc sum_{R in stab(q_i)} Gamma_ib R_ad], "
                                   "w_i = sqrt(star_i/24), P = degeneracy averaging, sym = (M+M^T)/2. Units THz. "
                                   "Off-diagonal part is Chaput's conductivity-equivalent operator (all three triplet "
                                   "channels with + sign, no momentum reversal): valid for odd (current-carrying) "
                                   "populations only; it does NOT annihilate the energy vector (see results/physop_*.json).",
               "kappa_from_matrix": "X_(iba) = w_i v_iba f THzToEv/(4 kB T^2 sinh(f THzToEv/2kBT)); Y = A^+ X; "
                                    "kappa = conv kB T^2/N sum_ib sum_R [(R x)(R y)^T + transpose], N = 16337, "
                                    "conv = phono3py get_unit_to_WmK()/V = 1602.17733/(2 pi V[A^3]) (phonopy e = 1.60217733e-19 C), "
                                    "kB = 8.617338256808316e-05 eV/K, THzToEv = 4.13566733e-3; sum over the 24 rotations",
               "pre_assembly_matrix": "A_w before degeneracy averaging and symmetrization (kept for asymmetry studies)",
               "gamma_N_U": "phono3py is_N_U split of the RTA gamma (Normal / Umklapp by G=0 vs G!=0 of the BZ-folded "
                            "triplet), same mesh and setting, all temperatures in the file",
           },
           "files": {}}
    inv["files"]["meta"] = entry(ooc / "meta.npz", "ir grid points, star weights, w=sqrt(star/24), rotated grid points, "
                                 "Cartesian rotations (24), frequencies, group velocities (ir points), temperatures, mesh, volume")
    inv["files"]["gamma_rows_stage"] = entry(ooc / "gamma.npy", "gamma (T, ir, band) from the collision-row computation, THz")
    temps = list(np.load(ooc / "meta.npz")["temps"]) if (ooc / "meta.npz").exists() else []
    for k, T in enumerate(temps):
        inv["files"][f"assembled_matrix_T{T:g}"] = entry(ooc / f"A_T{k}.npy", f"assembled symmetrized collision matrix at {T:g} K")
        inv["files"][f"solution_T{T:g}"] = entry(ooc / f"Y_T{k}.npy", "LBTE solution Y = A^+ X (ir, band, 3)")
        inv["files"][f"solve_record_T{T:g}"] = entry(ooc / f"solve_T{k}.json", "PCG record, kappa, residuals")
    inv["files"]["pre_assembly_matrix_T300"] = entry(ooc / "A_T2_prephono3pyavg_presym.npy",
                                                      "A_w at 300 K before degeneracy averaging/symmetrization")
    inv["files"]["solve_record_T300_isotope"] = entry(ooc / "solve_T2_iso.json", "PCG record with isotope diagonal")
    inv["files"]["rta_NU_all_T"] = entry(R / "prod-rta-m313117-NU" / "kappa-m313117.hdf5",
                                         "phono3py RTA output: frequency, group_velocity, gv_by_gv, heat_capacity, gamma, "
                                         "gamma_N, gamma_U, mode_kappa, kappa, weight, qpoint, grid_point, temperature")
    inv["files"]["rta_isotope_all_T"] = entry(R / "prod-rta-m313117-iso" / "kappa-m313117.hdf5",
                                              "phono3py RTA output with natural-abundance isotope scattering (gamma_isotope)")
    common.write_json(C / "results" / "saved_data_inventory.json", inv)
    print(json.dumps({k: {kk: v.get(kk) for kk in ("path", "bytes", "sha256", "shape", "missing")}
                      for k, v in inv["files"].items()}, indent=1))


if __name__ == "__main__":
    main()
