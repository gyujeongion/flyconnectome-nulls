"""Mechanism of the generation-0 advantage: wiring vs engineered innate readout vs vision."""
import json, numpy as np
from scipy.stats import wilcoxon
rng = np.random.default_rng(0); NULLS = ["N1", "N2", "N3"]
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for r, i in enumerate(o): run = max(run, (m - r) * ps[i]); adj[i] = min(1, run)
    return adj
def ci(d, B=10000):
    d = np.asarray(d); return np.percentile(np.median(rng.choice(d, (B, len(d))), 1), [2.5, 97.5])
R = [json.loads(l) for l in open("runs/gen0/gen0assay.jsonl")]
for scheme in ["global", "pergroup", "distmatch"]:
    D = {}
    for r in R:
        if r["homeo"] == scheme: D.setdefault(r["seed"], {})[r["assay"]] = r["res"]
    seeds = sorted(s for s, v in D.items() if len(v) >= 6)
    if not seeds: continue
    print(f"\n========== generation 0, {scheme} (seeds {seeds}) ==========")
    for tag in ["normal", "no_innate_readout", "blind", "anosmic", "blind_and_no_innate", "random_walk"]:
        A = np.array([D[s][tag]["A"]["fit"] for s in seeds]); ps, parts = [], []
        for n in NULLS:
            d = np.array([D[s][tag]["A"]["fit"] - D[s][tag][n]["fit"] for s in seeds])
            ps.append(wilcoxon(d).pvalue if np.any(d != 0) else 1.0); parts.append((n, np.array([D[s][tag][n]["fit"] for s in seeds]), d))
        adj = holm(ps)
        print(f"  {tag:22s} A {np.median(A):+.3f} | " + " | ".join(
            f"{n} {np.median(nv):+.3f} diff {np.median(d):+.3f} CI[{ci(d)[0]:+.3f},{ci(d)[1]:+.3f}] {int((d>0).sum())}/{len(d)} p={pa:.3f}"
            for (n, nv, d), pa in zip(parts, adj)))
