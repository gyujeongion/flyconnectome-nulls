#!/bin/bash
# build bioRxiv long-form PDF and supplement PDF from the markdown sources (pandoc + xelatex)
set -e
cd "$(dirname "$0")/.."
COMMON=(--pdf-engine=xelatex -V mainfont="STIXGeneral" -V mathfont="STIX Two Math" -V geometry:margin=2.4cm -V fontsize=11pt -V linestretch=1.15 -V colorlinks=true -V linkcolor=blue -V urlcolor=blue --resource-path=.)
python3 paper/prep.py PAPER.md paper/_paper.md; python3 paper/prep.py SUPPLEMENTARY.md paper/_supp.md
pandoc paper/_paper.md -o paper/biorxiv.pdf "${COMMON[@]}" -V documentclass=article -H paper/header.tex
pandoc paper/_supp.md -o paper/supplementary.pdf "${COMMON[@]}" -V documentclass=article -H paper/header.tex
python3 paper/prep.py paper/alife.md paper/_alife.md
pandoc paper/_alife.md -o paper/alife.pdf "${COMMON[@]}" -V documentclass=article -V classoption=twocolumn -V fontsize=10pt -V geometry:margin=1.9cm -H paper/header_twocol.tex
echo built
