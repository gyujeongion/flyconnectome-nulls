"""Pre-registered analysis for PROTOCOL.md (paper-v1)."""
import json, glob, os, sys, numpy as np
from scipy.stats import wilcoxon
rng = np.random.default_rng(0)
NULLS = ["N1", "N2", "N3"]; ABL = ["no_plasticity", "anosmic", "blind"]
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for r, i in enumerate(o): run = max(run, (m - r) * ps[i]); adj[i] = min(1, run)
    return adj
def boot_ci(d, B=10000):
    d = np.asarray(d); bs = np.median(rng.choice(d, (B, len(d))), 1); return np.percentile(bs, [2.5, 97.5])
def load(scheme):
    out = {}
    for f in sorted(glob.glob(f"runs/paper_{scheme}_s*/probe.jsonl")):
        run = os.path.dirname(f)
        if not os.path.exists(f"{run}/done.flag"): continue
        s = int(run.split("_s")[-1]); P = [json.loads(l) for l in open(f)]
        out[s] = P
    return out
res_all = {}
for scheme in ["global", "pergroup", "distmatch"]:
    D = load(scheme)
    if not D: continue
    print(f"\n==================== scheme: {scheme}  (completed seeds: {sorted(D)}) ====================")
    metrics = {}
    for s, P in D.items():
        gens = np.array([p["gen"] for p in P])
        for c in ["A"] + NULLS:
            y = np.array([p[f"{c}.normal"]["fit"] for p in P])
            last = P[-1]
            metrics.setdefault(c, {}).setdefault(s, dict(
                AUC=float(np.trapezoid(y, gens) / (gens[-1] - gens[0])) if len(gens) > 1 else float(y[0]),
                gen0=float(y[0]), final=float(y[-1]),
                **{f"abl_{a}": float(last[f"{c}.normal"]["fit"] - last[f"{c}.{a}"]["fit"]) for a in ABL}))
    rows = []
    for mname in ["AUC", "gen0", "final"] + [f"abl_{a}" for a in ABL]:
        seeds = sorted(set.intersection(*[set(metrics[c]) for c in ["A"] + NULLS]))
        A = np.array([metrics["A"][s][mname] for s in seeds])
        ps, parts = [], []
        for n in NULLS:
            d = np.array([metrics["A"][s][mname] - metrics[n][s][mname] for s in seeds])
            p = wilcoxon(d).pvalue if len(d) >= 5 and np.any(d != 0) else float("nan")
            ps.append(p if np.isfinite(p) else 1.0)
            lo, hi = boot_ci(d) if len(d) >= 3 else (np.nan, np.nan)
            parts.append((n, np.array([metrics[n][s][mname] for s in seeds]), d, lo, hi))
        adj = holm(ps)
        tag = "PRIMARY " if mname == "AUC" else "        "
        print(f"{tag}{mname:18s} n={len(seeds)} A median {np.median(A):.3f}")
        for (n, nv, d, lo, hi), p, pa in zip(parts, ps, adj):
            print(f"            vs {n}: null median {np.median(nv):.3f} | diff median {np.median(d):+.3f} 95%CI [{lo:+.3f},{hi:+.3f}] | A>null {int((d > 0).sum())}/{len(d)} | p={p:.4f} p_holm={pa:.4f}")
            rows.append(dict(scheme=scheme, metric=mname, null=n, n=len(d), A_median=float(np.median(A)), null_median=float(np.median(nv)),
                             diff_median=float(np.median(d)), ci=[float(lo), float(hi)], A_gt=int((d > 0).sum()), p=float(p), p_holm=float(pa)))
    res_all[scheme] = rows
json.dump(res_all, open("runs/paper_stats.json", "w"), indent=1)
print("\nsaved runs/paper_stats.json")
