"""Bounded algebraic sewing tests; no physical rate or scalar closure claim."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np

def relative(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))

def channels(v,vr,sewing):
    s0,s1,s2=sewing
    forward=6*np.einsum("ijk,jb,kc->ibc",v,s1.conj(),s2.conj(),optimize=True)
    inverse=6*np.einsum("abc,ia->ibc",vr,s0.conj(),optimize=True)
    return forward,inverse

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--material-array",type=Path)
    args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(20260927)
    def unitary(n):
        z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        return np.linalg.qr(z)[0]
    sizes=(3,4,2)
    v=rng.normal(size=sizes)+1j*rng.normal(size=sizes)
    s=[unitary(n) for n in sizes]
    vr=np.einsum("ijk,ia,jb,kc->abc",v.conj(),*s,optimize=True)
    f,g=channels(v,vr,s)
    error=relative(g,f.conj())
    # Deliberate wrong conjugation must be detectable, not just a real-gauge test.
    wrong=6*np.einsum("ijk,jb,kc->ibc",v,s[1],s[2],optimize=True)
    negative=relative(g,wrong.conj())
    checks={"complex_unitary_channel_identity":error<1e-12,
            "wrong_conjugation_detected":negative>0.1}
    result={"scope":"finite tensor algebra; exact frequency-preserving sewing assumed for physical channels",
            "synthetic":{"inverse_vs_forward_conjugate":error,
                         "wrong_conjugation_relative_error":negative},"checks":checks}
    if args.material_array:
        with np.load(args.material_array,allow_pickle=False) as z:
            v=z["aligned_0"];vr=z["reversed_amplitude"]
            sewing=[z[f"block_sewing{i}"] for i in range(3)]
            f,g=channels(v,vr,sewing)
            error=relative(g,f.conj())
            # No claim of general unequal-frequency unitary mode freedom.
            ff=z["original_frequencies"];fr=z["reversed_frequencies"]
            intertwining=[float(np.linalg.norm(ff[i,:,None]*u-u*fr[i,None,:]))
                          for i,u in enumerate(sewing)]
        result["material"]={"array_sha256":hashlib.sha256(args.material_array.read_bytes()).hexdigest(),
            "inverse_vs_forward_conjugate":error,
            "frequency_intertwining_residual_THz":intertwining,
            "physical_rows":{"parent":1,"daughters":[26,14]},
            "scope":"one selected tensor; approximate numerical frequency blocks; no exact resonant channel admitted"}
        checks["selected_material_channel_identity"]=error<1e-10
        checks["frequency_intertwining"]=max(intertwining)<1e-10
    result["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output_dir/"channel-sewing.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    if not all(checks.values()):raise RuntimeError("Sewing check failed")
if __name__=="__main__":
    main()
