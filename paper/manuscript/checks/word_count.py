"""Approximate word counts of the manuscript sections.

Counted: running text words; each inline formula counts as one word. Excluded: comments, displayed
equations (counted separately), labels, citation keys and cross-reference keys.

Usage (from this directory):  python -B word_count.py
"""
import re
from pathlib import Path

DISPLAY = r"(equation|equation\*|multline|gather|align|split|aligned|gathered)"


def count(text):
    text = re.sub(r"(?<!\\)%.*", "", text)                       # comments
    n_disp = len(re.findall(r"\\begin\{(equation|multline|gather|align)\*?\}", text))
    text = re.sub(r"\\begin\{" + DISPLAY + r"\}.*?\\end\{\1\}", " ", text, flags=re.S)
    text = re.sub(r"\\\[.*?\\\]", " ", text, flags=re.S)
    text = re.sub(r"\$[^$]+\$", " MATH ", text)                  # inline math = one word
    text = re.sub(r"\\(cite|ref|eqref|label)\{[^}]*\}", " ", text)
    text = re.sub(r"\\begin\{[^}]*\}(\[[^]]*\])?|\\end\{[^}]*\}", " ", text)
    text = re.sub(r"\\item\[[^]]*\]", " ", text)
    text = re.sub(r"\\[a-zA-Z@]+\*?", " ", text)                 # command names (arguments kept)
    words = re.findall(r"[A-Za-z0-9][A-Za-z0-9'\-]*", text)
    return len(words), n_disp


root = Path(__file__).resolve().parent.parent
rows = []
main = (root / "main.tex").read_text(encoding="utf-8")
abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", main, flags=re.S).group(1)
rows.append(("abstract (main.tex)",) + count(abstract))
for name in ["introduction", "kinetics", "identifiability", "frontface", "aln", "fdtr", "conclusions", "notation"]:
    rows.append((f"sections/{name}.tex",) + count((root / "sections" / f"{name}.tex").read_text(encoding="utf-8")))
total = sum(r[1] for r in rows)
for r in rows:
    print(f"{r[0]:32s} {r[1]:6d} words   {r[2]:3d} displayed equations")
print(f"{'total':32s} {total:6d} words")
