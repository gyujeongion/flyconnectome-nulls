"""Calibration-sensitivity analysis of pre-evolution circuit metrics.
Paired by seed: diff = A - null. Wilcoxon signed-rank (two-sided), Holm correction over the 3 nulls per (scheme, metric)."""
import json, sys, numpy as np
from scipy.stats import wilcoxon
path = sys.argv[1] if len(sys.argv) > 1 else "runs/sens/circuit.jsonl"
R = [json.loads(l) for l in open(path)]
spec = lambda v: float(v.split()[0])
METRICS = {
    "cos_KC (lower=better odor separation)": (lambda r: r["cos_KC"], "lower"),
    "cos_PN": (lambda r: r["cos_PN"], "lower"),
    "DAN_sugar response": (lambda r: r["DAN_sugar"], "higher"),
    "DAN_malaise response": (lambda r: r["DAN_malaise"], "higher"),
    "learning specificity eta=-0.05": (lambda r: spec(r["spec_eta-0.05"]), "higher"),
    "loom escape (L-R turn; negative=correct)": (lambda r: r["turn_loomL_minus_loomR"], "lower"),
}
schemes = sorted({r["homeo"] for r in R}); nulls = ["N1", "N2", "N3"]
def holm(ps):
    order = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * ps[i]); adj[i] = min(1.0, run)
    return adj
out = []
for sch in schemes:
    print(f"\n===== calibration scheme: {sch} =====")
    cal = {c: [r["calib"] for r in R if r["homeo"] == sch and r["cond"] == c] for c in ["A"] + nulls}
    for c, L in cal.items():
        if L:
            keys = [k for k in L[0] if isinstance(L[0][k], (int, float))]
            print(f"  calib {c}: " + " ".join(f"{k}={np.mean([x[k] for x in L]):.3g}" for k in keys))
    for mname, (f, better) in METRICS.items():
        vals = {c: {r["seed"]: f(r["res"]) for r in R if r["homeo"] == sch and r["cond"] == c} for c in ["A"] + nulls}
        seeds = sorted(set(vals["A"]) & set.intersection(*[set(vals[c]) for c in nulls])) if all(vals[c] for c in nulls) else []
        if len(seeds) < 5:
            print(f"  {mname}: n={len(seeds)} (skip)"); continue
        ps, rows = [], []
        for c in nulls:
            d = np.array([vals["A"][s] - vals[c][s] for s in seeds])
            try: p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
            except Exception: p = 1.0
            ps.append(p); rows.append((c, d))
        adj = holm(np.array(ps))
        A = np.array([vals["A"][s] for s in seeds])
        line = f"  {mname:42s} A={np.median(A):+.3f} [{A.min():+.3f},{A.max():+.3f}]"
        for (c, d), p, pa in zip(rows, ps, adj):
            nv = np.array([vals[c][s] for s in seeds])
            wins = int(np.sum(d < 0)) if better == "lower" else int(np.sum(d > 0))
            line += f" | {c} {np.median(nv):+.3f} A-better {wins}/{len(d)} p_holm={pa:.3f}"
            out.append(dict(scheme=sch, metric=mname, null=c, n=len(d), A_median=float(np.median(A)), null_median=float(np.median(nv)),
                            median_diff=float(np.median(d)), A_better=wins, p=float(p), p_holm=float(pa)))
        print(line)
json.dump(out, open(path.replace(".jsonl", "_stats.json"), "w"), indent=1)
