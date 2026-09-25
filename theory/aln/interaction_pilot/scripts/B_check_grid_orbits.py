"""Independent finite-group partition check; does not read or compute pp values."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import inspect
import json
from collections import Counter
from pathlib import Path
import h5py
import numpy as np
import spglib
from phonopy.interface.vasp import read_vasp
from phonopy.phonon.grid import BZGrid
from phono3py.phonon3.triplets import get_triplets_at_q

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--run-dir",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    base=args.run_dir
    paths={name:base/name for name in
        ["POSCAR","phonon-m333.hdf5","pp-m333-g1-s0.1.hdf5",
         "gamma_detail-m333-g1-s0.1.hdf5"]}
    with h5py.File(paths["phonon-m333.hdf5"],"r") as f:
        mesh=f["mesh"][:].astype(np.int64)
        address=f["grid_address"][:].astype(np.int64)
        frequency=f["frequency"][:]
        phonon_keys=sorted(f.keys())
    with h5py.File(paths["pp-m333-g1-s0.1.hdf5"],"r") as f:
        triplet=f["triplet"][:].astype(np.int64)
        full=f["triplet_all"][:].astype(np.int64)
        weights=f["weight"][:].astype(np.int64)
        pp_shape=list(f["pp"].shape)
        pp_keys=sorted(f.keys())
    with h5py.File(paths["gamma_detail-m333-g1-s0.1.hdf5"],"r") as f:
        assert np.array_equal(f["triplet"][:],triplet)
        assert np.array_equal(f["triplet_all"][:],full)
        assert np.array_equal(f["weight"][:],weights)
    cell=read_vasp(paths["POSCAR"])
    symprec=1e-5
    dataset=spglib.get_symmetry_dataset(
        (cell.cell,cell.scaled_positions,cell.numbers),symprec=symprec)
    assert dataset is not None
    grid=BZGrid(mesh,lattice=cell.cell,symmetry_dataset=dataset,
                is_time_reversal=True,store_dense_gp_map=True,lang="C")
    assert np.array_equal(grid.D_diag,mesh)
    assert np.array_equal(grid.PS,[0,0,0])
    assert np.array_equal(grid.Q,np.eye(3,dtype=np.int64))
    assert np.array_equal(grid.addresses,address)
    n=int(np.prod(mesh))
    def idx(a):
        b=np.asarray(a,dtype=np.int64)%mesh
        return int(b@np.array([1,mesh[0],mesh[0]*mesh[1]],dtype=np.int64))
    regular=np.array([[i%mesh[0],(i//mesh[0])%mesh[1],
                       i//(mesh[0]*mesh[1])] for i in range(n)],dtype=np.int64)
    b2g=np.array([idx(a) for a in address],dtype=np.int64)
    assert np.array_equal(b2g,grid.bzg2grg)
    assert np.array_equal(b2g[grid.grg2bzg],np.arange(n))
    fibers={i:np.flatnonzero(b2g==i).tolist() for i in range(n)}
    assert all(fibers.values())
    assert int(grid.gp_map[0])==0 and int(grid.gp_map[-1])==len(address)
    assert np.array_equal(np.diff(grid.gp_map),[len(fibers[i]) for i in range(n)])
    alias_frequency_error=max(
        float(np.max(np.abs(frequency[fibers[i]]-frequency[fibers[i][0]])))
        for i in range(n))
    identity=tuple(range(n))
    def compose(p,q):
        return tuple(p[q[i]] for i in range(n))
    def closure(generators):
        out={identity}
        queue=[identity]
        while queue:
            p=queue.pop()
            for gen in generators:
                candidate=compose(gen,p)
                if candidate not in out:
                    out.add(candidate)
                    queue.append(candidate)
        return out
    def orbits(generators):
        pending=set(range(n))
        out=[]
        while pending:
            orbit={min(pending)}
            queue=list(orbit)
            while queue:
                i=queue.pop()
                for gen in generators:
                    j=gen[i]
                    if j not in orbit:
                        orbit.add(j)
                        queue.append(j)
            pending-=orbit
            out.append(tuple(sorted(orbit)))
        return out
    # Construct reciprocal address rotations independently from direct-cell
    # spglib operations; equal diagonal mesh makes Q=D transforms scalar.
    reciprocal=[]
    for direct in dataset.rotations:
        rfloat=np.linalg.inv(direct).T
        r=np.rint(rfloat).astype(np.int64)
        assert np.allclose(rfloat,r,rtol=0,atol=1e-12)
        reciprocal.extend([r,-r])
    rotations={tuple(r.reshape(-1).tolist()):r for r in reciprocal}
    perms={tuple(idx(r@a) for a in regular) for r in rotations.values()}
    assert all(sorted(p)==list(range(n)) for p in perms)
    assert closure(perms)==perms
    library_perms={tuple(idx(r@a) for a in regular) for r in grid.rotations}
    assert perms==library_perms
    anchor_bz=int(triplet[0,0])
    assert anchor_bz==1
    anchor=int(b2g[anchor_bz])
    ar=regular[anchor]
    little={p for p in perms if p[anchor]==anchor}
    assert closure(little)==little
    swap=tuple(idx(-ar-a) for a in regular)
    assert compose(swap,swap)==identity
    group=closure(little|{swap})
    classes=orbits(little|{swap})
    hclasses=orbits(little)
    class_of={r:c for c in classes for r in c}
    canon=np.array([min(class_of[r]) for r in range(n)],dtype=np.int64)
    hcanon=np.array([min(c) for r in range(n) for c in hclasses if r in c],
                    dtype=np.int64)
    assert len(hcanon)==n
    rows=[]
    represented=[]
    for j,(t,w) in enumerate(zip(triplet,weights)):
        assert np.all((t>=0)&(t<len(address)))
        assert int(b2g[t[0]])==anchor
        summed=address[t].sum(axis=0)
        assert np.all(summed%mesh==0)
        r=int(b2g[t[1]])
        assert int(b2g[t[2]])==swap[r]
        orbit=class_of[r]
        assert int(w)==len(orbit)
        represented.append(orbit)
        stabilizer=sum(p[r]==r for p in group)
        assert len(group)==len(orbit)*stabilizer
        rows.append({"row":j,"triplet_bz":t.tolist(),
                     "second_regular":r,"orbit_regular":list(orbit),
                     "export_weight":int(w),"enumerated_weight":len(orbit),
                     "stabilizer_size":int(stabilizer),
                     "reciprocal_G":(summed//mesh).tolist()})
    assert len(set(represented))==len(classes)==len(triplet)==8
    assert set(represented)==set(classes)
    assert int(weights.sum())==n==27
    full_seconds=[]
    for t in full:
        assert int(b2g[t[0]])==anchor
        assert np.all(address[t].sum(axis=0)%mesh==0)
        r=int(b2g[t[1]])
        assert int(b2g[t[2]])==swap[r]
        full_seconds.append(r)
    assert sorted(full_seconds)==list(range(n))
    assert len(full)==n
    # Native reconstruction is a separate comparison, not input to the
    # independent orbit enumeration above.
    native_t,native_w,native_map,native_q=get_triplets_at_q(
        anchor_bz,grid,is_time_reversal=True,swappable=True,lang="C")
    assert np.array_equal(native_t,triplet)
    assert np.array_equal(native_w,weights)
    def partition_from_map(mapping):
        return {tuple(np.flatnonzero(mapping==rep).tolist())
                for rep in np.unique(mapping)}
    assert partition_from_map(native_map)==set(classes)
    assert partition_from_map(native_q)==set(hclasses)
    assert np.array_equal(native_map[native_map],native_map)
    assert np.array_equal(native_q[native_q],native_q)
    tr=tuple(idx(-a) for a in regular)
    assert compose(tr,tr)==identity
    tr_fixed=sum(tr[i]==i for i in range(n))
    assert tr_fixed==int(np.prod(np.gcd(2,mesh)))==1
    neg_anchor=tr[anchor]
    assert neg_anchor!=anchor
    negative_little={p for p in perms if p[neg_anchor]==neg_anchor}
    negative_swap=tuple(idx(-regular[neg_anchor]-a) for a in regular)
    negative_classes=set(orbits(negative_little|{negative_swap}))
    assert {tuple(sorted(tr[i] for i in c)) for c in classes}==negative_classes
    burnside_fixed=[sum(p[i]==i for i in range(n)) for p in group]
    assert sum(burnside_fixed)==len(group)*len(classes)
    package_versions={p:importlib.metadata.version(p)
                      for p in ["phono3py","phonopy","spglib","numpy","h5py"]}
    result={
        "status":"PASS",
        "method":"spglib direct rotations -> reciprocal permutations; independent integer closure and BFS; native triplet maps used only afterward",
        "run_directory":str(base.resolve()),
        "versions":package_versions,"symprec":symprec,
        "spacegroup_number":int(dataset.number),
        "spacegroup_international":dataset.international,
        "backend_for_grid_reconstruction":"C",
        "grid":{"mesh":mesh.tolist(),"D_diag":grid.D_diag.tolist(),
                "P":grid.P.tolist(),"Q":grid.Q.tolist(),"PS":grid.PS.tolist(),
                "bz_rows":len(address),"regular_points":n,
                "fiber_size_histogram":dict(Counter(len(v) for v in fibers.values())),
                "max_alias_frequency_abs_error_THz":alias_frequency_error,
                "bzg2grg":b2g.tolist(),"grg2bzg":grid.grg2bzg.tolist(),
                "dense_gp_map":grid.gp_map.tolist()},
        "anchor":{"bz":anchor_bz,"regular":anchor,"address":address[anchor_bz].tolist(),
                  "time_reverse_regular":neg_anchor,
                  "canonical_time_reverse_bz":int(grid.grg2bzg[neg_anchor])},
        "group_counts":{"direct_operations":len(dataset.rotations),
                        "reciprocal_with_time_reversal":len(perms),
                        "little_group":len(little),"little_plus_swap":len(group),
                        "little_orbits":len(hclasses),"triplet_orbits":len(classes),
                        "burnside_fixed_counts_sorted":sorted(burnside_fixed)},
        "rows":rows,"all_27_full_rows_bijective":True,
        "independent_map_triplets":canon.tolist(),"independent_map_q":hcanon.tolist(),
        "native_map_triplets":native_map.tolist(),"native_map_q":native_q.tolist(),
        "native_maps_match_partitions":True,
        "negative_anchor_grid_partition_matches":True,
        "time_reversal_fixed_points":tr_fixed,
        "triplet_and_weights_equal_in_pp_and_gamma_detail":True,
        "pp_shape_not_values":pp_shape,
        "exported_dataset_names":{"phonon":phonon_keys,"pp":pp_keys},
        "file_sha256":{name:sha256(path) for name,path in paths.items()},
        "source_sha256":{str(Path(inspect.getsourcefile(fn)).resolve()):
                         sha256(inspect.getsourcefile(fn))
                         for fn in [BZGrid,get_triplets_at_q]},
        "limitations":[
            "Q, P, PS, both grid maps and rotation lists were reconstructed, not present as dedicated export datasets.",
            "The reproduction uses symprec=1e-5, time reversal enabled, swappable=True and the C backend; does not establish Rust equivalence.",
            "The -q anchor interaction tensor was not exported or compared; only its grid orbit partition was derived.",
            "No band permutation, degenerate subspace sewing matrix, eigenvector gauge, or interaction reciprocity was tested.",
            "Momentum closure is not energy resonance or proof of a physical collision operator."
        ]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"group_counts":result["group_counts"],
                     "anchor":result["anchor"],"rows":rows,
                     "alias_frequency_error_THz":alias_frequency_error,
                     "output":str(args.output)},indent=2))

if __name__=="__main__":
    main()
