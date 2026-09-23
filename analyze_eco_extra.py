"""paper-v2 follow-up (exploratory, labelled as such in the paper):
(1) crossover generation: first probe generation at which the connectome's common-garden fitness
    falls to or below each null's, per seed and ecology;
(2) cross-ecology common garden on generation-500 elites (secondary endpoint, PROTOCOL_ECO.md);
(3) key numbers for manuscript provenance (results/eco_key.json)."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon

COND = ["A", "N1", "N2", "N3"]
CELLS = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}

def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj

def runs(pv, pt):
    return sorted(r for r in glob.glob(f"runs/eco_P{pv}_T{pt}_s*") if os.path.exists(f"{r}/done.flag"))

print("=== (1) crossover generation (median over seeds; 'never' = connectome stays ahead through gen 550) ===")
cross = {}
for (pv, pt), lab in CELLS.items():
    out = []
    for n in COND[1:]:
        gs = []
        for r in runs(pv, pt):
            P = [json.loads(l) for l in open(f"{r}/probe.jsonl")]
            g = [p["gen"] for p in P if p[f"A.normal"]["fit"] <= p[f"{n}.normal"]["fit"] and p["gen"] > 0]
            gs.append(g[0] if g else np.inf)
        gs = np.array(gs, float); fin = np.isfinite(gs)
        cross[(lab, n)] = gs
        med = np.median(gs)
        out.append(f"{n} median {('never' if not np.isfinite(med) else int(med))} crossed {fin.sum()}/{len(gs)}")
    print(f" {lab}: " + " | ".join(out))

print("\n=== (2) cross-ecology common garden, generation-500 elites (fitness A - null, per seed) ===")
rows = []
for (pv, pt), lab in CELLS.items():
    for r in runs(pv, pt):
        f = f"{r}/assay.jsonl"
        if not os.path.exists(f): continue
        last = {}
        for l in open(f):
            d = json.loads(l)
            if d["gen"] == 500 and d["assay"].startswith("eco_P"): last[d["assay"]] = d["res"]   # keep last write
        if len(last) == 4: rows.append((lab, int(r.split("_s")[-1]), last))
n_runs = len(rows)
print(f" runs with complete cross-ecology assay: {n_runs}/40")
TEST = {"eco_P0.0_T0.0": "P0T0", "eco_P0.2_T0.0": "P1T0", "eco_P0.0_T1.0": "P0T1", "eco_P0.2_T1.0": "P1T1"}
cross_tab = {}
for home in CELLS.values():
    for tkey, tlab in TEST.items():
        R = [x for x in rows if x[0] == home]
        if len(R) < 3: continue
        ps, txt = [], []
        for n in COND[1:]:
            d = np.array([x[2][tkey]["A"]["fit"] - x[2][tkey][n]["fit"] for x in R])
            p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0; ps.append(p)
            txt.append((n, np.median(d), int((d > 0).sum()), len(d)))
            cross_tab[(home, tlab, n)] = float(np.median(d))
        adj = holm(ps)
        print(f" evolved in {home} -> tested in {tlab}: " + " | ".join(
            f"{n} {m:+.3f} {g}/{k} p_holm {a:.4f}" for (n, m, g, k), a in zip(txt, adj)))

# H3: does the connectome's relative standing depend on the test ecology matching home? (home vs away gap)
print("\n=== (3) H3 home-vs-away: A-null gap at home minus mean gap away (per seed) ===")
for home in CELLS.values():
    R = [x for x in rows if x[0] == home]
    if len(R) < 3: continue
    hk = [k for k, v in TEST.items() if v == home][0]
    ps, txt = [], []
    for n in COND[1:]:
        d = np.array([(x[2][hk]["A"]["fit"] - x[2][hk][n]["fit"]) -
                      np.mean([x[2][k]["A"]["fit"] - x[2][k][n]["fit"] for k in TEST if k != hk]) for x in R])
        p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0; ps.append(p); txt.append((n, np.median(d), int((d > 0).sum()), len(d)))
    adj = holm(ps)
    print(f" {home}: " + " | ".join(f"{n} {m:+.3f} {g}/{k} p_holm {a:.4f}" for (n, m, g, k), a in zip(txt, adj)))
