"""Reproduce a bounded export pilot on D or another user-selected workspace.

Requires theory/aln/interaction_pilot/requirements-pilot.txt in an isolated
Python environment. A completed run retains a KNOWN FAILED validation check;
successful process completion is not certification of a physical operator.
"""
from pathlib import Path
import argparse, hashlib, importlib.metadata, json, os, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/"theory/aln/interaction_pilot"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="New directory; existing paths are never overwritten")
    args = parser.parse_args()
    for package, version in (("phono3py", "4.5.0"), ("phonopy", "4.5.0"), ("phonors", "0.3.0")):
        if importlib.metadata.version(package) != version:
            raise RuntimeError(f"{package}=={version} required for this pinned pilot")
    expected = json.loads((ARCHIVE/"pilot-run.json").read_text(encoding="utf-8"))["inputs"]
    for name, entry in expected.items():
        data = (args.input_dir/name).read_bytes()
        if len(data) != entry["bytes"] or hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError(f"Input does not match pinned source: {name}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for p in ARCHIVE.glob("*.py"):
        shutil.copy2(p, output/p.name)
    shutil.copytree(ARCHIVE/"scripts", output/"scripts")
    for name in ("gp1-m333", "pair-m333"):
        run = output/"runs"/name
        run.mkdir(parents=True)
        for file in expected:
            shutil.copy2(args.input_dir/file, run/file)
    env = dict(os.environ, OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2",
               RAYON_NUM_THREADS="2", PYTHONDONTWRITEBYTECODE="1")
    commands = [
        ["pilot_run.py"], ["pilot_pair.py"], ["inspect_pilot.py"],
        ["reconstruct_gamma.py"], ["pair_checks.py"], ["harmonic_controls.py"],
        ["matrix_control.py"], ["cutoff_control.py"], ["review_numerical_arrays.py"],
        ["scripts/B_check_grid_orbits.py", "--run-dir", str(output/"runs/gp1-m333"),
         "--output", str(output/"B_grid_orbits.json")],
    ]
    executed = []
    for command in commands:
        script = command[0]
        result = subprocess.run([sys.executable, "-B", str(output/script), *command[1:]],
                                cwd=output, env=env, capture_output=True, timeout=180)
        log = output/(Path(script).stem+".log")
        log.write_bytes(result.stdout+b"\nSTDERR\n"+result.stderr)
        expected_failure = False
        if script == "inspect_pilot.py" and result.returncode == 1:
            report = json.loads((output/"pilot-inspection.json").read_text())
            failed = [k for k,v in report["checks"].items() if not v]
            expected_failure = failed == ["BZ_duplicate_frequencies"]
        executed.append({"script":script, "returncode":result.returncode,
                         "expected_failed_validation":expected_failure})
        print(script, result.returncode,
              "(known validation failure)" if expected_failure else "", flush=True)
        if result.returncode and not expected_failure:
            raise RuntimeError(f"{script} failed; inspect {log}")
    status = {
        "status":"pilot_executed_with_unresolved_validation_failure",
        "failed_check":"BZ_duplicate_frequencies",
        "scope":"No grid convergence, canonical event operator, or material transport validation",
        "commands":executed}
    (output/"reproduction-status.json").write_text(json.dumps(status,indent=2)+"\n")
    print(json.dumps(status,indent=2))


if __name__ == "__main__":
    main()
