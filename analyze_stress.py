"""Ceiling-effect probe: evolved elites re-evaluated in harder worlds (no further evolution)."""
import json, glob, numpy as np
from scipy.stats import wilcoxon
NULLS = ["N1", "N2", "N3"]
rows = [json.loads(l) for f in glob.glob("runs/paper_*/assay.jsonl") for l in open(f)]
TAGS = ["base", "stress_pred0.12", "stress_pred0.15", "stress_pred0.18", "stress_toxic0.7", "stress_toxic1.0"]
for scheme in ["global", "pergroup", "distmatch"]:
    D = {}
    for r in rows:
        if r["homeo"] == scheme: D.setdefault(r["seed"], {})[r["assay"]] = r["res"]
    seeds = sorted(s for s, v in D.items() if all(t in v for t in TAGS))
    if not seeds: continue
    print(f"--- ceiling-effect stress grid, {scheme} (n={len(seeds)}) ---")
    for tag in TAGS:
        A = np.array([D[s][tag]["A"]["fit"] for s in seeds])
        surv = np.median([D[s][tag]["A"]["surv"] for s in seeds])
        line = f"  {tag:18s} A {np.median(A):+.3f} (surv {surv:.2f})"
        for n in NULLS:
            nv = np.array([D[s][tag][n]["fit"] for s in seeds]); d = A - nv
            p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
            line += f" | {n} {np.median(nv):+.3f} diff {np.median(d):+.3f} {int((d > 0).sum())}/{len(d)} p={p:.3f}"
        print(line)
