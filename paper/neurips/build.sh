#!/bin/bash
# Build the NeurIPS-style preprint (paper/neurips/paper.md -> paper/neurips/paper_neurips.pdf) and its
# LaTeX source bundle for arXiv (paper/neurips/arxiv_source.zip). NeurIPS 2026 style, [preprint] option.
# Fails on any LaTeX error or missing glyph; reports where the references start so the 9-page
# main-text limit can be checked.
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
B=build; rm -rf "$B"; mkdir -p "$B/figs"
python3 "$ROOT/paper/prep.py" paper.md "$B/paper_prep.md"
pandoc "$B/paper_prep.md" -o "$B/main.tex" -s --template neurips.latex
cp kit/neurips_2026.sty "$B/"
pend=$(grep -c "\[\[PENDING" paper.md || true); [ "$pend" -eq 0 ] || echo "note: $pend PENDING placeholders remain in paper.md" >&2
python3 -c 'import re,sys;[print(f) for f in sorted(set(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", open(sys.argv[1]).read())))]' "$B/main.tex" |
while read -r f; do [ -f "$ROOT/$f" ] || { echo "missing figure $f" >&2; exit 1; }; mkdir -p "$B/$(dirname "$f")"; cp "$ROOT/$f" "$B/$f"; done
(cd "$B" && for i in 1 2; do xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null || { grep -m5 '^!' main.log >&2; exit 1; }; done)
miss=$(grep -c 'Missing character' "$B/main.log" || true)
[ "$miss" -eq 0 ] || { grep 'Missing character' "$B/main.log" | sort | uniq -c >&2; echo "missing glyphs: $miss" >&2; exit 1; }
cp "$B/main.pdf" paper_neurips.pdf
pages=$(pdfinfo paper_neurips.pdf | awk '/^Pages/{print $2}')
refp=""; for i in $(seq 1 "$pages"); do pdftotext -f "$i" -l "$i" paper_neurips.pdf - | grep -qx 'References' && { refp=$i; break; }; done
echo "paper_neurips.pdf: $pages pages; references start on page ${refp:-?} (main text must end by page 9)"
(cd "$B" && rm -f ../arxiv_source.zip && zip -qr ../arxiv_source.zip main.tex neurips_2026.sty figs)
echo "arxiv_source.zip: $(du -h arxiv_source.zip | cut -f1)"
