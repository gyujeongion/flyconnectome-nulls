#!/usr/bin/env python3
"""Post hoc, not registered (research log 59). Two questions about the n = 20 extension of the corrected grid.

1. Were the nulls of seeds 10-19 built like those of seeds 0-9? Olfactory-to-descending direct share of the calibrated
   weights each run actually used (generation-0 genomes), and whether that share predicts a null's lead over the
   connectome across seeds (Spearman, per null and pooled over N1 and N2).
2. How many seeds would a per-ecology equivalence verdict need? Plug-in simulation: experiments of n seeds drawn from
   the 20 observed paired differences, 90 % percentile-bootstrap interval of the median, share inside +-0.10. This
   treats the observed distribution as the truth, so it is a rough guide, not a power analysis for a new design.
"""
import json
import numpy as np
from scipy.stats import spearmanr

B = np.load("data/brain.npz", allow_pickle=True)
orn = np.array(sorted(set(int(i) for i in B["in_orn"]))); dn = np.array(sorted(set(int(i) for i in B["dn"])))
ECO = {"0.0": "P0T0", "0.2": "P1T0"}


def share(W):
    W = np.abs(W.astype(np.float64)); return float(W[np.ix_(dn, orn)].sum() / W[:, orn].sum())


def last(r):
    p = None
    for l in open(f"runs/{r}/probe.jsonl"):
        p = l
    return json.loads(p)


out = {"shortcut_share": {}, "share_vs_lead": {}, "seeds_needed": {}}
rng = np.random.default_rng(1)
for pv, lab in ECO.items():
    runs = [f"ecofix_P{pv}_T0.0_s{s}" for s in range(20)]
    P = [last(r) for r in runs]
    Z = [np.load(f"runs/{r}/elite_g000000.npz") for r in runs]
    pooled = ([], [])
    for c in ("N1", "N2", "N4", "N5"):
        sh = np.array([share(z[f"{c}_W"][0]) for z in Z])
        lead = np.array([p[f"{c}.normal"]["fit"] - p["A.normal"]["fit"] for p in P])
        out["shortcut_share"][f"{lab}/{c}"] = {"s0-9_median_pct": 100 * float(np.median(sh[:10])), "s10-19_median_pct": 100 * float(np.median(sh[10:])),
                                              "range_pct": [100 * float(sh.min()), 100 * float(sh.max())]}
        if c in ("N1", "N2"):
            rho, p = spearmanr(sh, lead); out["share_vs_lead"][f"{lab}/{c}"] = {"rho": float(rho), "p": float(p), "n": 20}
            pooled[0].extend(sh); pooled[1].extend(lead)
        d = -lead                                    # connectome minus null
        need = {}
        for n in (20, 40, 60, 100):
            ok = 0
            for _ in range(1000):
                x = rng.choice(d, n); b = np.median(x[rng.integers(0, n, (400, n))], 1)
                lo, hi = np.percentile(b, [5, 95]); ok += (lo > -0.1) and (hi < 0.1)
            need[n] = ok / 1000
        out["seeds_needed"][f"{lab}/{c}"] = {"median": float(np.median(d)), "sd": float(d.std(ddof=1)), "p_inside_by_n": need}
    rho, p = spearmanr(*pooled); out["share_vs_lead"][f"{lab}/N1+N2"] = {"rho": float(rho), "p": float(p), "n": 40}
json.dump(out, open("results/posthoc_n20.json", "w"), indent=1)
print(json.dumps(out, indent=1))
