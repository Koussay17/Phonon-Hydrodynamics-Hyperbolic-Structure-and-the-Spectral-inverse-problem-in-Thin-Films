"""Replay bounded resonance/entropy benchmarks without material inputs."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/"theory/aln/resonance_measure"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    # Archived checks use assert: an optimized interpreter would strip them and
    # still report success (28 Sep numerical audit). Refuse it explicitly.
    if sys.flags.optimize:
        raise SystemExit("reproduce_resonance_measure: do not run under python -O")
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    # Copy executable inputs only; no stale JSON can pass as a fresh result.
    for source in ARCHIVE.rglob("*.py"):
        target=out/source.relative_to(ARCHIVE)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="2")
    env.pop("PYTHONOPTIMIZE",None)
    commands=[]
    def run(name,extra=()):
        result=subprocess.run([sys.executable,"-B",str(out/name),*map(str,extra)],
            cwd=out,env=env,capture_output=True,timeout=180)
        (out/(Path(name).stem+".log")).write_bytes(result.stdout+b"\nSTDERR\n"+result.stderr)
        commands.append({"script":name,"returncode":result.returncode})
        if result.returncode:raise RuntimeError(f"{name} failed; see saved log")
        print(name,"completed",flush=True)
        return result.stdout
    a=json.loads(run("branches/A-coarea-check.py"))
    (out/"A-results.json").write_text(json.dumps(a,indent=2)+"\n")
    run("scripts/B_weak_form_check.py")
    b=json.loads((out/"runs/B_weak_form_check.json").read_text())
    run("experiments/C_resonance_benchmark.py",["--output-dir",out/"C-results"])
    c=json.loads((out/"C-results/C_resonance_benchmark.json").read_text())
    run("D-limits-check.py")
    d=json.loads((out/"D-limits-check.json").read_text())
    checks={
        "A_mass":a["mass_relative_error"]<1e-12,
        "A_energy_kernel":a["L_energy_relative_residual"]<1e-12,
        "A_entropy_identity":a["entropy_identity_relative_error"]<1e-10,
        "B_exact_scalar_calculus":all(v=="0" for v in b["exact_scalar_residuals"].values()),
        "B_interpolation_negative_control":b["interpolation"]["occupation_interpolation_equilibrium_drift"]=="5/7",
        "B_energy_kernel":b["interpolation"]["energy_residual"]==["0","0"],
        "C_assertions":c["all_assertions_passed"],
        "D_independent_poisson_sums":max(r["absolute_discrepancy"] for r in d["linear_gaussian_comb"])<1e-12,
        "D_positive_heating":all(r["energy_drift"]>0 for r in d["equilibrium_heating"]["rows"]),
        "D_regular_heating_order":abs(d["equilibrium_heating"]["rows"][-1]["observed_energy_drift_order"]-2)<0.02,
        "D_small_drift_wrong_measure":any(r["energy_drift"]<1e-12 and abs(r["unit_weight_measure"]-2)>1
            for r in d["small_drift_is_not_measure_convergence"]),
    }
    run("experiments/C2_weak_form.py",["--output-dir",out/"C2-results"])
    c2=json.loads((out/"C2-results/C2_weak_form.json").read_text())
    checks["C2_finite_weak_form_assertions"]=c2["all_assertions_passed"]
    checks["C2_separated_variables_finite_and_accurate"]=all(
        r["mobility_relative_error"]<1e-12 and r["flux_relative_error"]<=r["flux_tolerance"]
        for r in c2["separated_entropy_variables"])
    checks["C2_thermal_identity_negative_control"]=(
        c2["thermal_identity_negative_control"]["distributivity_check_passes"]=={"1.0":True,"1.3":False})
    report={"scope":"bounded analytic/synthetic checks; no material or general dynamic convergence claim",
        "checks":checks,"commands":commands}
    (out/"reproduction-status.json").write_text(json.dumps(report,indent=2)+"\n")
    if not all(checks.values()):raise RuntimeError("Scientific check failed")
    print(json.dumps(report,indent=2))
if __name__=="__main__":
    main()
