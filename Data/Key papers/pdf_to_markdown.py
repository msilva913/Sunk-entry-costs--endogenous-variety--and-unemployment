"""
pdf_to_markdown.py -- Mirror the PDF paper library as searchable Markdown
=========================================================================
Converts every PDF in this folder to a Markdown file under ``markdown/``.
The PDFs are never modified; the Markdown is a derived, disposable mirror.

Why this exists
---------------
The non-negotiable rules in ``context/principles.md`` (N11-N14) forbid citing a
paper we have not read and require a verbatim quotation for any claim about the
literature.  Verifying a quote by scrolling a PDF is slow and error prone.  A
plain-text mirror makes the library greppable, so a claim can be checked in one
command:

    grep -rn "product destruction" "Key papers/markdown/"

Page markers are preserved for exactly this reason.  A quote is only useful in
the draft if it can be attributed to a page.

Fidelity over prettiness
------------------------
The output is deliberately close to the source text rather than well-formed
Markdown, because the default copy is the one quotes are taken from.  Line
breaks, hyphenation, and column layout are preserved.

Ligatures are the one exception, and they are always expanded.  A PDF storing
the single glyph U+FB01 means "fi", so leaving it makes ``grep finding`` miss
the word.  Expanding restores the characters the author wrote.

Two opt-in flags trade fidelity for searchability:

``--dehyphenate``  joins words split by an end-of-line hyphen.
``--reflow``       joins lines within a paragraph, so a grep for a phrase
                   spanning a line break matches.  Table-like blocks are left
                   alone.

Either flag stamps a warning in the file header.  Do not quote from output
generated with them; regenerate without flags first.

Usage
-----
    python pdf_to_markdown.py              # convert new or changed PDFs
    python pdf_to_markdown.py --force      # reconvert everything
    python pdf_to_markdown.py --reflow     # search-friendly, not quote-safe
    python pdf_to_markdown.py --outdir ../Notes/paper_text

Requires the ``pdftotext`` executable (xpdf or poppler).  On this machine it
ships with Git for Windows and is already on PATH.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import shutil
import subprocess
import sys
from pathlib import Path

BASE_DIR   = Path(__file__).resolve().parent
DEFAULT_OUT = BASE_DIR / "markdown"

# pdftotext emits a form feed between pages. That is the only reliable page
# delimiter it gives us, so it becomes the page marker.
_FORM_FEED = "\f"


def find_pdftotext() -> str:
    """Locate the pdftotext executable, or explain how to get one."""
    exe = shutil.which("pdftotext")
    if exe:
        return exe
    sys.exit(
        "pdftotext not found on PATH.\n"
        "It ships with Git for Windows (mingw64/bin) and with poppler-utils.\n"
        "Install poppler, or add the Git mingw64 bin directory to PATH."
    )


def slugify(name: str) -> str:
    """Filesystem-friendly stem, so long journal filenames stay greppable."""
    s = re.sub(r"[^\w\s-]", "", name).strip()
    s = re.sub(r"[\s_]+", "_", s)
    # Journal-supplied filenames run long. Cap the stem so the full path stays
    # under the Windows 260-character limit even in a deep output directory.
    return s[:80].rstrip("_")


def extract_text(pdf: Path, exe: str) -> str:
    """Run pdftotext with layout preserved. Returns raw text including \\f."""
    # -layout keeps column structure, which matters for tables in the papers.
    # -enc UTF-8 avoids mojibake in author names and math.
    proc = subprocess.run(
        [exe, "-layout", "-enc", "UTF-8", str(pdf), "-"],
        capture_output=True,
    )
    if proc.returncode != 0:
        err = proc.stderr.decode("utf-8", "replace").strip()
        raise RuntimeError(f"pdftotext failed ({proc.returncode}): {err}")
    return proc.stdout.decode("utf-8", "replace")


# Typographic ligatures are a rendering artifact, not content: the PDF means
# "finding" even though it stores one glyph. Expanding them restores the real
# characters, so grep for "finding" matches. Always applied.
_LIGATURES = {
    "ﬀ": "ff", "ﬁ": "fi", "ﬂ": "fl",
    "ﬃ": "ffi", "ﬄ": "ffl", "ﬅ": "st", "ﬆ": "st",
}


def normalize_ligatures(text: str) -> str:
    for glyph, plain in _LIGATURES.items():
        text = text.replace(glyph, plain)
    return text


def dehyphenate(text: str) -> str:
    """Join words broken across a line by a trailing hyphen.

    Only applied on request. It improves search but silently edits the text,
    so the result must not be treated as quotable source.
    """
    return re.sub(r"(\w)-\n\s*(\w)", r"\1\2", text)


def reflow(text: str) -> str:
    """Join lines within a paragraph so sentences occupy one line.

    ``pdftotext -layout`` preserves the source line breaks, which keeps tables
    readable but means a grep for a phrase spanning a line break finds nothing.
    Reflowing fixes phrase search at the cost of the column layout, so it is
    opt-in. Blank lines and indented or short lines (headings, table rows) are
    left alone.
    """
    out = []
    for para in text.split("\n\n"):
        lines = para.split("\n")
        # Table-like blocks have lots of internal whitespace; do not touch them.
        if sum("   " in ln for ln in lines) > len(lines) / 2:
            out.append(para)
            continue
        out.append(re.sub(r"\s*\n\s*", " ", para).strip())
    return "\n\n".join(out)


def to_markdown(raw: str, source: Path, exe: str,
                dehyph: bool, do_reflow: bool) -> str:
    """Wrap the extracted text in page markers and a provenance header."""
    pages = raw.split(_FORM_FEED)
    # A trailing form feed leaves an empty final element.
    if pages and not pages[-1].strip():
        pages.pop()

    body_parts = []
    for i, page in enumerate(pages, start=1):
        page = normalize_ligatures(page.rstrip())
        if dehyph:
            page = dehyphenate(page)
        if do_reflow:
            page = reflow(page)
        # Squeeze runs of blank lines; pdftotext is generous with them.
        page = re.sub(r"\n{3,}", "\n\n", page)
        body_parts.append(f"<!-- page {i} -->\n\n{page}\n")

    header = (
        f"# {source.stem}\n\n"
        f"> Auto-generated from `{source.name}` by `pdf_to_markdown.py`.\n"
        f"> Do not edit. Regenerate instead. The PDF is the source of record.\n"
        f">\n"
        f"> - pages: {len(pages)}\n"
        f"> - converted: {_dt.date.today().isoformat()}\n"
        f"> - tool: {Path(exe).name}"
        f"{', de-hyphenated' if dehyph else ''}"
        f"{', reflowed' if do_reflow else ''}"
        f"{' (edited for search: not safe to quote verbatim)' if (dehyph or do_reflow) else ''}\n\n"
        f"---\n\n"
    )
    return header + "\n".join(body_parts)


def needs_update(pdf: Path, md: Path, force: bool) -> bool:
    if force or not md.exists():
        return True
    return pdf.stat().st_mtime > md.stat().st_mtime


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT,
                    help="destination directory (default: ./markdown)")
    ap.add_argument("--force", action="store_true",
                    help="reconvert even if the Markdown is up to date")
    ap.add_argument("--dehyphenate", action="store_true",
                    help="join hyphenated line breaks; improves search, "
                         "makes the output unsafe to quote from")
    ap.add_argument("--reflow", action="store_true",
                    help="join lines within paragraphs so phrase search works "
                         "across line breaks; drops the column layout")
    args = ap.parse_args()

    # Journal filenames carry typographic dashes and accents that the default
    # Windows console codepage cannot encode. Print should never be what fails.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass

    exe = find_pdftotext()
    pdfs = sorted(BASE_DIR.glob("*.pdf"))
    if not pdfs:
        print(f"No PDFs found in {BASE_DIR}")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)

    converted = skipped = failed = 0
    for pdf in pdfs:
        md = args.outdir / f"{slugify(pdf.stem)}.md"
        if not needs_update(pdf, md, args.force):
            skipped += 1
            continue
        try:
            raw = extract_text(pdf, exe)
        except RuntimeError as e:
            print(f"  FAILED  {pdf.name}: {e}")
            failed += 1
            continue

        if not raw.strip():
            # Scanned papers with no text layer need OCR, which is out of scope.
            print(f"  EMPTY   {pdf.name}: no text layer, needs OCR")
            failed += 1
            continue

        md.write_text(to_markdown(raw, pdf, exe, args.dehyphenate, args.reflow),
                      encoding="utf-8")
        # --outdir may point outside this folder, so relative_to can fail.
        try:
            shown = md.relative_to(BASE_DIR)
        except ValueError:
            shown = md
        print(f"  ok      {pdf.name} -> {shown}")
        converted += 1

    print(f"\n{converted} converted, {skipped} up to date, {failed} failed "
          f"({len(pdfs)} PDFs in {BASE_DIR.name})")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
