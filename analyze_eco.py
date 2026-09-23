"""paper-v2 ecology grid analysis (pre-registered in PROTOCOL_ECO.md).

Primary : per-ecology paired A - null on common-garden probe AUC, and the
          interaction (P1 gap - P0 gap) tested per seed.
Scale    : gaps are reported both absolute and relative (gap / cell mean fitness),
           because the toxic cells have lower absolute fitness overall.
Secondary: generation-0 and final fitness, ablation effects.
"""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon

COND = ["A", "N1", "N2", "N3"]
CELLS = [("0.0", "0.0"), ("0.2", "0.0"), ("0.0", "1.0"), ("0.2", "1.0")]
FIT = {"AUC", "gen0", "final"}   # relative gap only meaningful for fitness levels, not ablation deltas near 0
LBL = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}


def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj


def boot(d, n=10000, seed=0):
    r = np.random.default_rng(seed)
    b = r.choice(d, (n, len(d)), replace=True).mean(1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def load(pv, pt):
    """{seed: {cond: dict of metrics}}"""
    out = {}
    for f in sorted(glob.glob(f"runs/eco_P{pv}_T{pt}_s*/probe.jsonl")):
        run = os.path.dirname(f)
        if not os.path.exists(f"{run}/done.flag"):
            continue
        P = [json.loads(l) for l in open(f)]
        gens = np.array([p["gen"] for p in P], float)
        d = {}
        for c in COND:
            y = np.array([p[f"{c}.normal"]["fit"] for p in P], float)
            last = P[-1]
            d[c] = dict(AUC=float(np.trapezoid(y, gens) / (gens[-1] - gens[0])), gen0=float(y[0]), final=float(y[-1]),
                        abl_anosmic=float(last[f"{c}.normal"]["fit"] - last[f"{c}.anosmic"]["fit"]),
                        abl_blind=float(last[f"{c}.normal"]["fit"] - last[f"{c}.blind"]["fit"]),
                        abl_noplast=float(last[f"{c}.normal"]["fit"] - last[f"{c}.no_plasticity"]["fit"]))
        out[int(run.split("_s")[-1])] = d
    return out


def gaps(D, metric):
    """per seed: {null: absolute gap}, plus the cell's mean level for normalisation"""
    seeds = sorted(D)
    lvl = np.array([np.mean([D[s][c][metric] for c in COND]) for s in seeds])
    g = {n: np.array([D[s]["A"][metric] - D[s][n][metric] for s in seeds]) for n in COND[1:]}
    return seeds, g, lvl


def test(tag, d):
    p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
    lo, hi = boot(d)
    return p, f"{tag} median {np.median(d):+.4f} CI[{lo:+.3f},{hi:+.3f}] {int((d > 0).sum())}/{len(d)}"


DATA = {cell: load(*cell) for cell in CELLS}
print("runs found:", {LBL[c]: len(DATA[c]) for c in CELLS})

for metric in ["AUC", "gen0", "final", "abl_anosmic", "abl_blind", "abl_noplast"]:
    print(f"\n================ {metric} ================")
    for cell in CELLS:
        D = DATA[cell]
        if not D:
            continue
        seeds, g, lvl = gaps(D, metric)
        ps, txt = zip(*[test(n, g[n]) for n in COND[1:]])
        adj = holm(ps)
        print(f" {LBL[cell]} (n={len(seeds)}) level {lvl.mean():.3f} | " +
              " | ".join(f"{t} p_holm {a:.4f}" + (f" rel {np.median(g[n]) / lvl.mean():+.3f}" if metric in FIT else "")
                         for t, a, n in zip(txt, adj, COND[1:])))
    # pre-registered interaction: predator effect on the gap, within matched toxicity
    print(" -- interaction (P1 gap - P0 gap), per seed, matched toxicity --")
    for pt in ["0.0", "1.0"]:
        D0, D1 = DATA[("0.0", pt)], DATA[("0.2", pt)]
        common = sorted(set(D0) & set(D1))
        if len(common) < 3:
            continue
        ps, txt = [], []
        for n in COND[1:]:
            d = np.array([(D1[s]["A"][metric] - D1[s][n][metric]) - (D0[s]["A"][metric] - D0[s][n][metric]) for s in common])
            p, t = test(n, d); ps.append(p); txt.append(t)
        adj = holm(ps)
        print(f"    toxicity {pt} (n={len(common)}) | " + " | ".join(f"{t} p_holm {a:.4f}" for t, a in zip(txt, adj)))
    # H2: toxicity effect on the gap, matched predator, absolute and relative
    print(" -- toxicity effect (T1 gap - T0 gap), per seed, matched predator --")
    for pv in ["0.0", "0.2"]:
        D0, D1 = DATA[(pv, "0.0")], DATA[(pv, "1.0")]
        common = sorted(set(D0) & set(D1))
        if len(common) < 3:
            continue
        l0 = np.mean([[D0[s][c][metric] for c in COND] for s in common])
        l1 = np.mean([[D1[s][c][metric] for c in COND] for s in common])
        ps, txt, rel = [], [], []
        for n in COND[1:]:
            a0 = np.array([D0[s]["A"][metric] - D0[s][n][metric] for s in common])
            a1 = np.array([D1[s]["A"][metric] - D1[s][n][metric] for s in common])
            p, t = test(n, a1 - a0); ps.append(p); txt.append(t)
            rel.append(np.median(a1) / (l1 + 1e-9) - np.median(a0) / (l0 + 1e-9))
        adj = holm(ps)
        print(f"    predator {pv} (n={len(common)}) | " +
              " | ".join(f"{t} p_holm {a:.4f}" + (f" rel_delta {r:+.3f}" if metric in FIT else "") for t, a, r in zip(txt, adj, rel)))
