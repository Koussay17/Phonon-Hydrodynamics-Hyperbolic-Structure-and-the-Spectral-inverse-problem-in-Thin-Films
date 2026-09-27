"""Reproduce finite channel checks; optional pinned-source and selected material audits."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/"theory/aln/channel_counting"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--source-python",type=Path,help="Pinned phono3py 4.5.0 interpreter")
    parser.add_argument("--pilot-root",type=Path,help="Optional retained corrected-cutoff pilot")
    args=parser.parse_args()
    if args.pilot_root and not args.source_python:
        parser.error("--pilot-root requires --source-python")
    output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    shutil.copytree(ARCHIVE/"experiments",output/"experiments")
    shutil.copy2(ARCHIVE/"D-source-factor-check.py",output/"D-source-factor-check.py")
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="2")
    records=[]
    checks={}
    def run(script,extra=(),interpreter=None):
        completed=subprocess.run([str(interpreter or sys.executable),"-B",str(output/script),
                                  *map(str,extra)],cwd=output,env=env,capture_output=True,timeout=240)
        (output/(Path(script).stem+".log")).write_bytes(completed.stdout+b"\nSTDERR\n"+completed.stderr)
        records.append({"script":script,"returncode":completed.returncode})
        if completed.returncode:raise RuntimeError(f"{script} failed: see its saved log")
        print(script,"completed",flush=True)
    run("experiments/C_fock_check.py",["--output-dir",output/"fock"])
    run("experiments/C2_closure_check.py",["--output-dir",output/"closure"])
    run("experiments/channel_sewing_check.py",["--output-dir",output/"synthetic-sewing"])
    c=json.loads((output/"fock/C_fock_check.json").read_text())
    c2=json.loads((output/"closure/C2_closure_check.json").read_text())
    s=json.loads((output/"synthetic-sewing/channel-sewing.json").read_text())
    checks["fock_matrices"]=max(c[k]["projected_vs_collected_relative_frobenius_error"]
                               for k in ("distinct","repeated"))<1e-12
    checks["vacuum_ratio"]=abs(c["vacuum_squared_amplitude_repeated_over_distinct"]-0.5)<1e-12
    checks["closure_exact_assertions"]=c2["all_assertions_passed"]
    checks["complex_sewing_and_negative_control"]=all(s["checks"].values())
    if args.source_python:
        run("D-source-factor-check.py",interpreter=args.source_python)
        d=json.loads((output/"D-source-factor-check.json").read_text())
        baseline=json.loads((ARCHIVE/"D-source-factor-check.json").read_text())
        checks["pinned_source_hashes"]=d["source_hashes"]==baseline["source_hashes"]
        checks["pinned_versions"]=all(d["versions"][k]=="4.5.0" for k in ("phono3py","phonopy"))
        checks["independent_oscillator_units"]=d["unit_checks"]["pp_independent_relative_difference"]<1e-12
        for kind,target in (("distinct",[1,1,1]),("repeated",[1,2])):
            ratios=d["conditional_channel_checks"][kind]["event_diagonal_over_inverse_lifetime"]
            checks[kind+"_conditional_linewidth_comparison"]=len(ratios)==len(target) and max(
                abs(x-y) for x,y in zip(ratios,target))<1e-12
    if args.pilot_root:
        run("experiments/material_permutations.py",["--pilot-root",args.pilot_root.resolve(),
            "--output-dir",output/"material"],args.source_python)
        m=json.loads((output/"material/permutation-audit.json").read_text())
        checks["selected_material_permutations"]=max(
            row["amplitude_relative_frobenius"] for row in m["permutations"])<1e-10
        checks["selected_pp_reconstruction"]=max(
            row["pp_global_scaled_difference"] for row in m["permutations"])<1e-10
        checks["selected_reverse_contraction"]=m["time_reversal"]["block_only_relative_frobenius"]<1e-10
        run("experiments/channel_sewing_check.py",["--output-dir",output/"material-sewing",
            "--material-array",output/"material/permutation-audit.npz"])
        ms=json.loads((output/"material-sewing/channel-sewing.json").read_text())
        checks["selected_material_channel_sewing"]=all(ms["checks"].values())
    status={"scope":"bounded algebra and selected tensor; not converged material kinetics",
            "checks":checks,"commands":records,"source_requested":bool(args.source_python),
            "material_requested":bool(args.pilot_root)}
    (output/"reproduction-status.json").write_text(json.dumps(status,indent=2)+"\n")
    if not all(checks.values()):raise RuntimeError("Failed scientific check: see reproduction-status.json")
    print(json.dumps(status,indent=2))
if __name__=="__main__":
    main()
