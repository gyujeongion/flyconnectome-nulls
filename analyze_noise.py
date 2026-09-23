"""Defect-28 robustness (post hoc): original (biased) vs zero-mean turning noise, gen-0 and gen-500 population samples, home ecology.
Seed-level inference: each seed's 4 ecologies averaged (n = 10); per-ecology medians descriptive."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon
COND = ["A", "N1", "N2", "N3"]
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj
D = {}   # (run, gen, z) -> res
for r in sorted(glob.glob("runs/eco_P*_s*")):
    f = f"{r}/assay.jsonl"
    if not os.path.isdir(r) or not os.path.exists(f): continue
    for l in open(f):
        d = json.loads(l)
        if d["assay"].startswith("noise_zm"): D[(r, d["gen"], int(d["assay"][-1]))] = d["res"]
runs = sorted({k[0] for k in D}); print("runs:", len(runs))
lab = lambda r: {"P0.0_T0.0": "P0T0", "P0.2_T0.0": "P1T0", "P0.0_T1.0": "P0T1", "P0.2_T1.0": "P1T1"}[r.split("eco_P")[1].rsplit("_s", 1)[0].replace("T", "T")[:9]]
seed = lambda r: int(r.rsplit("_s", 1)[1])
def seedmean(fn):
    out = {}
    for r in runs: out.setdefault(seed(r), []).append(fn(r))
    return np.array([np.mean(v) for _, v in sorted(out.items())])
for g in [0, 500]:
    print(f"\n=== generation {g} ===")
    for z in [0, 1]:
        print(f" noise {'biased (as evolved)' if z == 0 else 'zero-mean'}: " + " ".join(
            f"{c} {np.median([D[(r, g, z)][c]['fit'] for r in runs]):.3f}" for c in COND))
    # drop caused by removing the drift, per condition (seed means)
    drop = {c: seedmean(lambda r, c=c: D[(r, g, 0)][c]["fit"] - D[(r, g, 1)][c]["fit"]) for c in COND}
    print(" fitness lost when drift removed (seed means): " + " ".join(f"{c} {np.median(v):+.3f}" for c, v in drop.items()))
    ps, txt = [], []
    for n in COND[1:]:
        d = drop["A"] - drop[n]; ps.append(wilcoxon(d).pvalue); txt.append(f"A-{n} {np.median(d):+.3f} ({int((d > 0).sum())}/10)")
    print("   extra loss of connectome vs control: " + " | ".join(f"{t} p_holm {a:.4f}" for t, a in zip(txt, holm(ps))))
    # does the ordering survive zero-mean noise?
    ps, txt = [], []
    for n in COND[1:]:
        d = seedmean(lambda r, n=n: D[(r, g, 1)]["A"]["fit"] - D[(r, g, 1)][n]["fit"]); ps.append(wilcoxon(d).pvalue)
        txt.append(f"A-{n} {np.median(d):+.3f} A ahead {int((d > 0).sum())}/10")
    print("   connectome minus control under zero-mean noise: " + " | ".join(f"{t} p_holm {a:.4f}" for t, a in zip(txt, holm(ps))))
    print("   per ecology (zero-mean, median A-N1 / A-N2 over seeds): " + " ".join(
        f"{e}: {np.median([D[(r, g, 1)]['A']['fit'] - D[(r, g, 1)]['N1']['fit'] for r in runs if e in r]):+.2f}/"
        f"{np.median([D[(r, g, 1)]['A']['fit'] - D[(r, g, 1)]['N2']['fit'] for r in runs if e in r]):+.2f}"
        for e in ["P0.0_T0.0", "P0.2_T0.0", "P0.0_T1.0", "P0.2_T1.0"]))
