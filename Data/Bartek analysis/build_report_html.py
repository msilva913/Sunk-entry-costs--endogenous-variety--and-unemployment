"""
build_report_html.py — Convert report.md to a self-contained HTML file
=======================================================================
Uses pandoc for Markdown→HTML conversion (handles LaTeX math correctly)
and Python post-processing to inline all referenced images as base64 data
URIs, producing a fully portable single-file HTML report.

Why pandoc instead of the Python `markdown` library:
  The Python markdown library does not understand LaTeX math delimiters.
  It converts `_{subscript}` inside `$...$` to `<em>` italic tags, breaking
  every equation. Pandoc parses math blocks first and passes them verbatim
  to MathJax, avoiding this completely.

Image path resolution (for each ![...](path) reference):
  1. Relative to the report.md file itself
  2. data/results/           (IRF plots, persistence plots)
  3. data/instruments/       (permanence_ratios_by_supersector.png)
  4. Current working directory

Math rendering: MathJax 3 loaded from CDN (requires internet on first open;
browser-cached thereafter).

Usage
-----
  cd "Data/Bartek analysis"
  python build_report_html.py

  Optional arguments:
    --input   path to report.md         (default: report.md)
    --output  path to output HTML file  (default: report.html)
    --results path to results directory (default: data/results)

Dependencies
------------
  pandoc  (system install — https://pandoc.org/installing.html)
    Windows: winget install JohnMacFarlane.Pandoc
    macOS:   brew install pandoc
    Linux:   apt install pandoc
  No extra Python packages needed (standard library only).

Run
---
  python build_report_html.py
"""

import argparse
import base64
import mimetypes
import os
import re
import subprocess
import sys
from pathlib import Path

# ── working directory ─────────────────────────────────────────────────────────
try:
    os.chdir(Path(__file__).resolve().parent)
except NameError:
    os.chdir(
        Path.home()
        / "Documents/GitHub/Sunk_entry_costs_endogenous_variety_unemployment"
        / "Data/Bartek analysis"
    )

MATHJAX_URL = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Convert report.md to self-contained HTML (pandoc + base64 images)"
)
parser.add_argument("--input",   default="report.md",   help="Input Markdown file")
parser.add_argument("--output",  default="report.html",  help="Output HTML file")
parser.add_argument("--results", default="data/results", help="Results directory for image lookup")
args = parser.parse_args()

INPUT_PATH  = Path(args.input)
OUTPUT_PATH = Path(args.output)
RESULTS_DIR = Path(args.results)
INSTR_DIR   = Path("data/instruments")

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
if not INPUT_PATH.exists():
    sys.exit(f"ERROR: {INPUT_PATH} not found.\nRun from the 'Data/Bartek analysis' directory.")

try:
    result = subprocess.run(["pandoc", "--version"],
                            capture_output=True, text=True, check=True)
    pandoc_ver = result.stdout.splitlines()[0]
    print(f"Using {pandoc_ver}")
except (FileNotFoundError, subprocess.CalledProcessError):
    sys.exit(
        "ERROR: pandoc not found.\n"
        "Install from https://pandoc.org/installing.html\n"
        "  Windows: winget install JohnMacFarlane.Pandoc\n"
        "  macOS:   brew install pandoc\n"
        "  Linux:   apt install pandoc"
    )

# ---------------------------------------------------------------------------
# Step 1: pandoc Markdown -> HTML
# ---------------------------------------------------------------------------
print(f"\nConverting {INPUT_PATH} with pandoc...")

pandoc_cmd = [
    "pandoc",
    "--standalone",
    "--from", "markdown+tex_math_dollars+tex_math_single_backslash+smart",
    "--to", "html5",
    f"--mathjax={MATHJAX_URL}",
    "--toc",
    "--toc-depth=3",
    str(INPUT_PATH),
]

result = subprocess.run(pandoc_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("pandoc stderr:")
    print(result.stderr)
    sys.exit("ERROR: pandoc conversion failed.")

for line in result.stderr.splitlines():
    if "WARNING" in line and "fetch" not in line:
        print(f"  pandoc: {line}")

html = result.stdout

# ---------------------------------------------------------------------------
# Step 2: inject custom CSS before </head>
# ---------------------------------------------------------------------------
CUSTOM_CSS = """<style>
*,*::before,*::after{box-sizing:border-box;}
body{font-family:"Georgia","Times New Roman",serif;font-size:16px;line-height:1.8;
     color:#1a1a1a;background:#f5f5f5;margin:0;}
/* layout */
header[role=banner],body>nav#TOC,body>section,body>main,article,
.page-wrap{max-width:940px;margin:0 auto;background:#fff;
           box-shadow:0 2px 16px rgba(0,0,0,.09);}
header[role=banner]{padding:36px 52px 24px;border-bottom:2px solid #2c3e50;}
h1.title{font-size:1.8rem;color:#2c3e50;margin-bottom:4px;}
.subtitle{font-size:1.05rem;color:#555;margin-bottom:6px;}
.author,.date{font-size:0.9rem;color:#888;font-style:italic;}
body>nav#TOC,section,main,article{padding:0 52px 60px;}
/* TOC */
#TOC{background:#f4f7fb;border-left:4px solid #2c7be5;padding:20px 28px;
     margin:28px 0 36px;border-radius:0 4px 4px 0;}
#TOC>ul{padding-left:16px;}
#TOC li{margin:4px 0;}
#TOC a{color:#2c7be5;text-decoration:none;}
#TOC a:hover{text-decoration:underline;}
nav#TOC::before{content:"Contents";font-weight:700;font-size:1rem;
                color:#2c3e50;display:block;margin-bottom:10px;}
/* headings */
h1{font-size:1.55rem;color:#2c3e50;margin:2.2rem 0 0.8rem;
   border-bottom:1px solid #ddd;padding-bottom:5px;}
h2{font-size:1.25rem;color:#2c3e50;margin:1.9rem 0 0.6rem;}
h3{font-size:1.08rem;color:#34495e;margin:1.5rem 0 0.5rem;}
h4{font-size:0.98rem;color:#555;margin:1.2rem 0 0.4rem;}
/* text */
p{margin:0.7rem 0;}
ul,ol{margin:0.4rem 0 0.4rem 1.8rem;}
li{margin:0.2rem 0;}
hr{border:none;border-top:1px solid #ddd;margin:1.8rem 0;}
/* code */
code{font-family:"Consolas","Monaco",monospace;font-size:0.84em;
     background:#f0f3f7;padding:2px 5px;border-radius:3px;color:#c0392b;}
pre{background:#f0f3f7;border-left:3px solid #2c7be5;padding:14px 18px;
    overflow-x:auto;border-radius:0 4px 4px 0;margin:1rem 0;}
pre code{background:none;color:#1a1a1a;padding:0;font-size:0.9em;}
/* tables */
table{border-collapse:collapse;width:100%;margin:1.2rem 0;font-size:0.9rem;}
th{background:#2c3e50;color:#fff;padding:8px 12px;text-align:left;font-weight:600;}
td{padding:6px 12px;border-bottom:1px solid #e0e0e0;}
tr:nth-child(even) td{background:#f7f9fc;}
tr:hover td{background:#eef3fb;}
/* images */
img{max-width:100%;height:auto;display:block;margin:1.8rem auto;
    border:1px solid #ddd;border-radius:4px;box-shadow:0 2px 10px rgba(0,0,0,.12);}
/* figure captions */
figure{margin:0;}
figcaption{text-align:center;font-size:0.87rem;color:#555;
           margin-top:0.4rem;margin-bottom:1.6rem;line-height:1.5;font-style:italic;}
/* math */
.math.display{overflow-x:auto;margin:1rem 0;}
mjx-container{overflow-x:auto;max-width:100%;}
/* references */
#references p,section#references p{padding-left:2.2em;text-indent:-2.2em;
  margin:0.35rem 0;font-size:0.91rem;line-height:1.6;}
/* print */
@media print{
  body{background:#fff;font-size:11pt;}
  body>*{box-shadow:none;max-width:100%;padding:0;}
  nav#TOC{border:1px solid #ccc;}
  pre{border:1px solid #ccc;}
  img{box-shadow:none;border:1px solid #ccc;}
}
</style>
"""
html = html.replace("</head>", CUSTOM_CSS + "\n</head>", 1)

# ---------------------------------------------------------------------------
# Step 3: embed all local images as base64 data URIs
# ---------------------------------------------------------------------------
SEARCH_DIRS = [
    INPUT_PATH.parent,
    RESULTS_DIR,
    INSTR_DIR,
    Path("."),
]

def _find_image(src: str):
    p = Path(src)
    if p.is_absolute() and p.exists():
        return p
    for d in SEARCH_DIRS:
        c1 = d / p.name   # filename only
        if c1.exists():
            return c1
        c2 = d / p        # path as written
        if c2.exists():
            return c2
    return None

def _make_data_uri(found: Path) -> str:
    mime, _ = mimetypes.guess_type(str(found))
    if mime is None:
        mime = "image/png"
    data = base64.b64encode(found.read_bytes()).decode("ascii")
    kb   = found.stat().st_size / 1024
    print(f"  Embedded: {found.name}  ({kb:.0f} KB)")
    return f"data:{mime};base64,{data}"

def _replace_src(m: re.Match) -> str:
    full = m.group(0)
    src  = m.group(1)
    if src.startswith(("data:", "http://", "https://")):
        return full
    found = _find_image(src)
    if found is None:
        print(f"  WARNING image not found: {src}")
        return full
    return full.replace(src, _make_data_uri(found), 1)

print("\nEmbedding images...")
html = re.sub(r'<img\b[^>]*\bsrc="([^"]+)"', _replace_src, html)
html = re.sub(r"<img\b[^>]*\bsrc='([^']+)'", _replace_src, html)

# ---------------------------------------------------------------------------
# Step 4: write output
# ---------------------------------------------------------------------------
OUTPUT_PATH.write_text(html, encoding="utf-8")
size_kb = OUTPUT_PATH.stat().st_size / 1024
print(f"\nSaved: {OUTPUT_PATH}  ({size_kb:.0f} KB)")
print(f"Open:  file:///{OUTPUT_PATH.resolve().as_posix()}")
print()
print("Note: MathJax loads from CDN (needs internet on first open, then cached).")