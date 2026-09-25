#!/usr/bin/env python3
"""Post hoc: how the corrected-grid equivalence verdicts depend on the bound delta.

delta = 0.10 was fixed after the first grid (PAPER.md 2.6). This reports, for each connectome-minus-
control comparison of the corrected grid, the smallest symmetric bound its registered 90 % interval
satisfies (max of |lower|, |upper|), and the verdict at several bounds. It reads the intervals from
results/equivalence.txt, the file the paper's Table 2 is taken from, so the numbers cannot drift from
the paper; it computes nothing new from the runs.
"""
import json, os, re, sys

DELTAS = [0.05, 0.075, 0.10, 0.15, 0.20]
SRC = "results/equivalence.txt"
ROW = re.compile(r"^\s+A - (N\d)\s+median ([+-]\d\.\d+).*?90% CI \[([+-]\d\.\d+),([+-]\d\.\d+)\]")

cell, out = None, {}
in_grid = False
for line in open(SRC):
    if line.startswith("=== corrected grid"):
        in_grid = True; continue
    if line.startswith("===") and in_grid:
        break
    if not in_grid:
        continue
    head = line.strip()
    if head in ("P0T0", "P1T0", "P0T1", "P1T1", "pooled (seed means)"):
        cell = "seed mean" if head.startswith("pooled") else head
        continue
    m = ROW.match(line)
    if m and cell:
        n, med, lo, hi = m.group(1), *map(float, m.groups()[1:])
        need = max(abs(lo), abs(hi))
        out.setdefault(cell, {})[n] = {"median": med, "ci90": [lo, hi], "smallest_bound": round(need, 3),
                                       "equivalent_at": {str(d): bool(lo > -d and hi < d) for d in DELTAS}}

n_rows = sum(len(v) for v in out.values())
if n_rows != 20 or set(out) != {"P0T0", "P1T0", "P0T1", "P1T1", "seed mean"}:
    sys.exit(f"parse check failed: {n_rows} rows, cells {sorted(out)} (expected 20 rows in 5 cells)")

lines = ["Corrected grid, last probe: smallest symmetric bound each 90 % interval satisfies (post hoc)",
         f"{'cell':10s} {'ctrl':4s} {'median':>7s} {'90% CI':>17s} {'bound':>6s}  " + "  ".join(f"d={d:<5}" for d in DELTAS)]
for c in ("seed mean", "P0T0", "P1T0", "P0T1", "P1T1"):
    for n in ("N4", "N5", "N1", "N2"):
        r = out[c][n]
        lines.append(f"{c:10s} {n:4s} {r['median']:+7.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}] {r['smallest_bound']:6.3f}  "
                     + "  ".join(f"{'eq' if r['equivalent_at'][str(d)] else '--':7s}" for d in DELTAS))
print("\n".join(lines))
os.makedirs("results", exist_ok=True)
json.dump(out, open("results/delta_sensitivity.json", "w"), indent=1)
open("results/delta_sensitivity.txt", "w").write("\n".join(lines) + "\n")
print("\nwrote results/delta_sensitivity.json and .txt")
