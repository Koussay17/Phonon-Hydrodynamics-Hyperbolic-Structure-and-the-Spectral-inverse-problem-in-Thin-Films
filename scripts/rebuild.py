"""Rebuild authored figures and PDFs; leave third-party reference PDFs untouched."""
from pathlib import Path
import argparse
import contextlib
import io
import json
import os
import re
import runpy
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def figures():
    os.environ["MPLBACKEND"] = "Agg"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.show = lambda: None
    prior = Path.cwd()
    try:
        os.chdir(ROOT/"notebooks")
        nbpath = ROOT/"notebooks/01_cattaneo_energy_trajectory.ipynb"
        nb = json.loads(nbpath.read_text(encoding="utf-8"))
        namespace = {"__name__":"__main__"}
        count=0
        for cell in nb["cells"]:
            if cell["cell_type"] != "code":
                continue
            count+=1
            captured=io.StringIO()
            with contextlib.redirect_stdout(captured):
                exec(compile("".join(cell["source"]),str(nbpath),"exec"),namespace)
            cell["execution_count"]=count
            cell["outputs"]=[{"output_type":"stream","name":"stdout",
                              "text":captured.getvalue().splitlines(keepends=True)}]
        nbpath.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+"\n",encoding="utf-8")
        plt.close("all")
        print("FIGURE 01",flush=True)
        for script in sorted((ROOT/"notebooks").glob("*.py")):
            output=io.StringIO()
            with contextlib.redirect_stdout(output):
                runpy.run_path(str(script),run_name="__main__")
            (ROOT/".build"/(script.stem+".txt")).write_text(output.getvalue(),encoding="utf-8")
            plt.close("all")
            print("FIGURE",script.stem,flush=True)
        result=subprocess.run([sys.executable,"-X","utf8","-B",
                               str(ROOT/"scripts/analyse_experiment.py")],
                              cwd=ROOT,capture_output=True,check=True)
        (ROOT/".build/experimental_analysis.txt").write_bytes(result.stdout)
        print("FIGURE 06 and synthetic experimental analysis",flush=True)
        result=subprocess.run([sys.executable,"-X","utf8","-B",
                               str(ROOT/"scripts/analyse_aln_spectrum.py")],
                              cwd=ROOT,capture_output=True,check=True)
        (ROOT/".build/aln_spectral_analysis.txt").write_bytes(result.stdout)
        print("FIGURES 07-08 and spectral AlN analysis",flush=True)
        result=subprocess.run([sys.executable,"-X","utf8","-B",
                               str(ROOT/"scripts/analyse_aln_response.py")],
                              cwd=ROOT,capture_output=True,check=True)
        (ROOT/".build/aln_bulk_response.txt").write_bytes(result.stdout)
        print("Bulk RTA response diagnostic",flush=True)
    finally:
        os.chdir(prior)

def pdfs(selected=None, install_missing=False):
    sources=sorted((ROOT/"notes").glob("*.tex"))
    if selected:
        sources=[s for s in sources if s.stem in selected]
        if len(sources)!=len(set(selected)):
            raise ValueError("unknown note name")
    # Keep the existing published PDFs until every requested source compiles.
    for source in sources:
        engine="pdflatex" if source.stem=="camacho_2025" else "xelatex"
        executable=shutil.which(engine)
        if executable is None:
            raise RuntimeError(f"{engine} is required")
        target=ROOT/".build"/source.stem
        target.mkdir(exist_ok=True)
        options=[]
        if "miktex" in executable.lower():
            options=["--enable-installer" if install_missing else "--disable-installer","--miktex-disable-maintenance",
                     "--miktex-disable-diagnose"]
        command=[executable,*options,"-interaction=nonstopmode","-halt-on-error",
                 "-file-line-error",f"-output-directory={target}",source.name]
        for pass_number in (1,2):
            result=subprocess.run(command,cwd=source.parent,capture_output=True,timeout=180)
            log=result.stdout.decode("utf-8",errors="replace")
            (target/f"pass{pass_number}.txt").write_text(log,encoding="utf-8")
            if result.returncode:
                raise RuntimeError(f"{source.name}: build failed\n"+log[-4500:])
        log=(target/(source.stem+".log")).read_text(encoding="utf-8",errors="replace")
        warnings=[line for line in log.splitlines() if
                  any(x in line for x in ("Overfull","Missing character","Warning:"))]
        (target/"warnings.txt").write_text("\n".join(warnings),encoding="utf-8")
        print("PDF",source.stem, f"({len(warnings)} warnings)",flush=True)
    for source in sources:
        shutil.copy2(ROOT/".build"/source.stem/(source.stem+".pdf"),source.with_suffix(".pdf"))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--figures-only",action="store_true")
    group.add_argument("--pdfs-only",action="store_true")
    parser.add_argument("--note",action="append",help="TeX stem to rebuild, with --pdfs-only")
    parser.add_argument("--install-missing",action="store_true",help="Allow MiKTeX to download missing packages")
    args=parser.parse_args()
    (ROOT/".build").mkdir(exist_ok=True)
    if not args.pdfs_only:
        figures()
    if not args.figures_only:
        pdfs(args.note,args.install_missing)

if __name__=="__main__":
    main()
