"""Post-hoc assays: is the nulls' olfaction dependence associative learning? does the innate advantage transfer?"""
import json, glob, numpy as np
from scipy.stats import wilcoxon
rng = np.random.default_rng(0); NULLS = ["N1", "N2", "N3"]
def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for r, i in enumerate(o): run = max(run, (m - r) * ps[i]); adj[i] = min(1, run)
    return adj
def ci(d, B=10000):
    d = np.asarray(d); return np.percentile(np.median(rng.choice(d, (B, len(d))), 1), [2.5, 97.5])
rows = [json.loads(l) for f in glob.glob("runs/paper_*/assay.jsonl") for l in open(f)]
for scheme in ["global", "pergroup", "distmatch"]:
    D = {}
    for r in rows:
        if r["homeo"] != scheme: continue
        D.setdefault(r["seed"], {})[r["assay"]] = r["res"]
    seeds = sorted(s for s, v in D.items() if len(v) >= 6)
    if not seeds: continue
    print(f"\n========== {scheme} (seeds {seeds}) ==========")
    def val(s, a, c, k="fit"): return D[s][a][c][k]
    METRICS = {
        "plasticity benefit (base - base_noplast)": lambda s, c: val(s, "base", c) - val(s, "base_noplast", c),
        "reversal cost (base - reversal)": lambda s, c: val(s, "base", c) - val(s, "reversal", c),
        "reversal x plasticity interaction": lambda s, c: (val(s, "base", c) - val(s, "reversal", c)) - (val(s, "base_noplast", c) - val(s, "reversal_noplast", c)),
        "extra toxic bites after reversal (h2)": lambda s, c: val(s, "reversal", c, "bad_h2") - val(s, "base", c, "bad_h2"),
        "transfer: toxicity 0.9": lambda s, c: val(s, "transfer_toxic0.9", c),
        "transfer: wider plume 0.12": lambda s, c: val(s, "transfer_plume0.12", c),
        "base fitness (elites, gen 250)": lambda s, c: val(s, "base", c),
    }
    for mname, f in METRICS.items():
        A = np.array([f(s, "A") for s in seeds]); ps, parts = [], []
        for n in NULLS:
            d = np.array([f(s, "A") - f(s, n) for s in seeds])
            p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
            ps.append(p); parts.append((n, np.array([f(s, n) for s in seeds]), d))
        adj = holm(ps)
        print(f"  {mname:42s} A median {np.median(A):+.3f}")
        for (n, nv, d), p, pa in zip(parts, ps, adj):
            lo, hi = ci(d)
            print(f"      vs {n}: null {np.median(nv):+.3f} | diff {np.median(d):+.3f} CI[{lo:+.3f},{hi:+.3f}] A>null {int((d>0).sum())}/{len(d)} p_holm={pa:.4f}")
