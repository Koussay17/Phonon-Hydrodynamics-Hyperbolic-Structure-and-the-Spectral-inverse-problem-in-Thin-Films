"""Reproduce bounded finite checks and optionally the selected material contraction."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/"theory/aln/degenerate_action"

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--pilot-root",type=Path)
    p.add_argument("--material-python",type=Path)
    args=p.parse_args()
    output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    shutil.copytree(ARCHIVE/"experiments",output/"experiments")
    shutil.copy2(ARCHIVE/"D-phase-information.py",output/"D-phase-information.py")
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="2")
    records=[]
    def run(script,extra=(),interpreter=None):
        completed=subprocess.run([str(interpreter or sys.executable),"-B",str(output/script),*map(str,extra)],
            env=env,cwd=output,capture_output=True,timeout=180)
        (output/(Path(script).stem+".log")).write_bytes(completed.stdout+b"\nSTDERR\n"+completed.stderr)
        records.append({"script":script,"returncode":completed.returncode})
        if completed.returncode:
            raise RuntimeError(f"{script} failed; see saved log")
        print(script,"completed",flush=True)
        return completed.stdout
    run("experiments/C_lumpability_check.py")
    run("D-phase-information.py",["--output",output/"D-phase-information.json"])
    run("experiments/second_memory_check.py",["--repo-root",ROOT,"--output-dir",output/"experiments"])
    regression=json.loads(run("experiments/check_frequency_blocks.py"))
    c=json.loads((output/"experiments/C_lumpability_results.json").read_text())
    m=json.loads((output/"experiments/second_memory_check.json").read_text())
    d=json.loads((output/"D-phase-information.json").read_text())
    checks={
        "C_checks":all(c["checks"].values()),
        "memory_checks":all(m["checks"].values()),
        "pinned_repository_action_compared":m["repo"]["status"]=="compared",
        "clustering_regressions":all(regression["checks"].values()),
        "phase_family_exact_reference":max(row["maximum_eigenvalue_error"] for row in d["continuous_phase_family"])<1e-12,
        "phase_same_projected_generator":d["rotated_projected_block_generator_difference"]<1e-12,
    }
    if args.pilot_root:
        run("experiments/reconstruct_complex_triplet.py",
            ["--pilot-root",args.pilot_root.resolve(),"--output-dir",output/"material"],
            args.material_python)
        v=json.loads((output/"material/complex-triplet.json").read_text())
        checks["selected_pp_reconstruction"]=max(row["pp_max_scaled_error"] for row in v["comparisons"])<1e-10
        checks["selected_block_transport"]=v["block_only_transport_relative_frobenius"]<1e-10
    status={"scope":"bounded finite and selected contraction checks; not material collision validation",
            "checks":checks,"commands":records,"material_requested":bool(args.pilot_root)}
    (output/"reproduction-status.json").write_text(json.dumps(status,indent=2)+"\n")
    if not all(checks.values()):raise RuntimeError("scientific check failed; see reproduction-status.json")
    print(json.dumps(status,indent=2))

if __name__=="__main__":
    main()
