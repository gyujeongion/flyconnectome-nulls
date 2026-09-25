"""Is the recurrent interior doing work? Generation-500 samples re-evaluated with the interior zeroed or scrambled."""
import json, glob, os, numpy as np
from seeds import registered  # registered grids are seeds 0-9; see seeds.py
from scipy.stats import wilcoxon
CONDS = ["A", "N1", "N2", "N4", "N5"]; TAGS = ["intact", "interior_zeroed", "interior_scrambled"]
D = {}
for r in sorted(registered(glob.glob("runs/ecofix_P*_s*"))):
    f = f"{r}/assay.jsonl"
    if not os.path.exists(f): continue
    for l in open(f):
        d = json.loads(l)
        if d["assay"] in TAGS and d["gen"] == 500: D[(r, d["assay"])] = d["res"]
runs = sorted({k[0] for k in D if all((k[0], t) in D for t in TAGS)})
print(f"runs with all three blocks: {len(runs)}/40")
seeds = sorted({int(r.split('_s')[-1]) for r in runs})
def seedmean(fn): 
    out = {}
    for r in runs: out.setdefault(int(r.split('_s')[-1]), []).append(fn(r))
    return np.array([np.mean(v) for _, v in sorted(out.items())])
print("\nfitness (median over runs) and loss when the interior is removed:")
for t in TAGS:
    print(f"  {t:20s} " + " ".join(f"{c} {np.median([D[(r, t)][c]['fit'] for r in runs]):.2f}" for c in CONDS))
for t in TAGS[1:]:
    print(f"\n  relative loss from intact, {t} (seed means over all ecologies, n={len(seeds)}):")
    rel = {c: seedmean(lambda r, c=c, t=t: (D[(r, 'intact')][c]['fit'] - D[(r, t)][c]['fit']) / max(D[(r, 'intact')][c]['fit'], 1e-6)) for c in CONDS}
    print("    " + " | ".join(f"{c} {np.median(v):.0%}" for c, v in rel.items()))
    for n in CONDS[1:]:
        d = rel["A"] - rel[n]
        print(f"      A - {n}: {np.median(d):+.3f} p={wilcoxon(d).pvalue:.4f}")
