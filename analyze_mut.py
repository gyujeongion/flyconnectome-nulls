"""PROTOCOL_MUT.md analysis. Unit = run (40); per run, condition value = mean over its 8 parent genotypes."""
import glob, os, numpy as np
from scipy.stats import wilcoxon, spearmanr
import json

COND = ["A", "N1", "N2", "N3"]; GENS = [0, 250, 500]
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj

def per_parent(clone, mut):
    """clone [MM], mut [MM] -> metrics for one parent"""
    q95 = np.percentile(clone, 95); base = clone.mean(); d = mut - base
    return dict(p_ben=float((mut > q95).mean()), up=float(np.percentile(d, 90)), low=float(np.percentile(d, 10)),
                med=float(np.median(d)), up_rel=float(np.percentile(d, 90) / (abs(base) + 1e-3)),
                low_rel=float(np.percentile(d, 10) / (abs(base) + 1e-3)), base=float(base))

runs = sorted(r for r in glob.glob("runs/eco_*_s*") if os.path.isdir(r) and not r.endswith(".log"))
M = {}   # (run, gen, scale_idx, cond) -> metric dict (mean over parents)
for r in runs:
    for g in GENS:
        f = f"{r}/mut_g{g:06d}.npz"
        if not os.path.exists(f): continue
        Z = np.load(f); cl, mu, cond = Z["clone"], Z["mut"], Z["cond"]
        for si in range(mu.shape[0]):
            for c in COND:
                ix = np.where(cond == c)[0]
                ms = [per_parent(cl[i], mu[si, i]) for i in ix]
                M[(r, g, si, c)] = {k: float(np.mean([m[k] for m in ms])) for k in ms[0]}
SC = [0.5, 1.0, 2.0]
done = sorted({k[0] for k in M if k[1] == 500})
print(f"runs with all checkpoints: {len(done)}/40")

def test(g, si, key):
    rs = [r for r in done if all((r, g, si, c) in M for c in COND)]
    A = np.array([M[(r, g, si, "A")][key] for r in rs]); ps, txt = [], []
    for n in COND[1:]:
        d = A - np.array([M[(r, g, si, n)][key] for r in rs])
        p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0; ps.append(p)
        txt.append((n, np.median(d), int((d > 0).sum()), len(d)))
    adj = holm(ps)
    vals = " ".join(f"{c} {np.median([M[(r, g, si, c)][key] for r in rs]):.3f}" for c in COND)
    return f"   {key:7s} | {vals} || " + " | ".join(f"A-{n} {m:+.3f} {k_}/{k} p_holm {a:.4f}" for (n, m, k_, k), a in zip(txt, adj))

for g in GENS:
    for si, s in enumerate(SC):
        print(f"\n=== generation {g}, mutation magnitude x{s} (medians over runs; A-null paired over runs) ===")
        for key in ["p_ben", "up", "up_rel", "med", "low", "low_rel", "base"]:
            print(test(g, si, key))

print("\n=== P4 (exploratory): gen-0 (A-N1) P(beneficial) at x1 vs gen-600 (A-N1) final fitness gap, Spearman over runs ===")
xs, ys = [], []
for r in done:
    P = [json.loads(l) for l in open(f"{r}/probe.jsonl")]
    fin = P[-1]["A.normal"]["fit"] - P[-1]["N1.normal"]["fit"]
    xs.append(M[(r, 0, 1, "A")]["p_ben"] - M[(r, 0, 1, "N1")]["p_ben"]); ys.append(fin)
rho, p = spearmanr(xs, ys); print(f"   rho {rho:+.3f} p {p:.4f} n {len(xs)}")

print("\n=== CI and equivalence (post hoc): run-level A-null difference in P(beneficial), bootstrap 95% CI of the median; equivalence margin +-0.02 ===")
rng_ = np.random.default_rng(0)
for g in GENS:
    for si, s in enumerate(SC):
        rs = [r for r in done if all((r, g, si, c) in M for c in COND)]
        txt = []
        for n in ["N1", "N2"]:
            d = np.array([M[(r, g, si, "A")]["p_ben"] - M[(r, g, si, n)]["p_ben"] for r in rs])
            b = np.median(rng_.choice(d, (10000, len(d))), 1); lo, hi = np.percentile(b, [2.5, 97.5])
            txt.append(f"A-{n} {np.median(d):+.4f} CI[{lo:+.4f},{hi:+.4f}] {'within' if lo > -0.02 and hi < 0.02 else 'NOT within'} +-0.02")
        print(f"  gen {g} x{s}: " + " | ".join(txt))
print("\n=== median mutation effect (offspring cost), relative to parent, x1 ===")
for g in GENS:
    rs = [r for r in done if all((r, g, 1, c) in M for c in COND)]
    print(f"  gen {g}: " + " ".join(f"{c} {np.median([M[(r, g, 1, c)]['med'] / (abs(M[(r, g, 1, c)]['base']) + 1e-3) for r in rs]):+.3f}" for c in COND))
