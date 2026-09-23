"""Analyses requested in council manuscript review round 1 (all post hoc, labelled so in the paper)."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon
COND = ["A", "N1", "N2", "N3"]
CELLS = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}
rng = np.random.default_rng(0)
def boot(d, n=10000):
    b = rng.choice(d, (n, len(d)), replace=True); m = np.median(b, 1)
    return np.percentile(m, 2.5), np.percentile(m, 97.5)
def runs(pv, pt): return sorted(r for r in glob.glob(f"runs/eco_P{pv}_T{pt}_s*") if os.path.exists(f"{r}/done.flag"))
def probes(r): return [json.loads(l) for l in open(f"{r}/probe.jsonl")]

print("=== (a) sustained crossing: first probe after which the control stays >= connectome at every later probe ===")
for (pv, pt), lab in CELLS.items():
    out = []
    for n in COND[1:]:
        gs = []
        for r in runs(pv, pt):
            P = probes(r); d = np.array([p[f"{n}.normal"]["fit"] - p["A.normal"]["fit"] for p in P]); g = [p["gen"] for p in P]
            k = next((i for i in range(1, len(d)) if np.all(d[i:] >= 0)), None)
            gs.append(g[k] if k is not None else np.inf)
        gs = np.array(gs, float); fin = np.isfinite(gs)
        out.append(f"{n} sustained {fin.sum()}/{len(gs)}, median gen {int(np.median(gs[fin])) if fin.any() else 'never'}")
    print(f" {lab}: " + " | ".join(out))

print("\n=== (b) fitness components at the final probe (medians over seeds; top-16 common garden, intact) ===")
for (pv, pt), lab in CELLS.items():
    row = []
    for c in COND:
        L = [probes(r)[-1][f"{c}.normal"] for r in runs(pv, pt)]
        row.append(f"{c}: fit {np.median([x['fit'] for x in L]):.2f} surv {np.median([x['surv'] for x in L]):.2f} "
                   f"good {np.median([x['good'] for x in L]):.1f} bad {np.median([x['bad'] for x in L]):.1f}")
    print(f" {lab} | " + " | ".join(row))

print("\n=== (c) absolute ablation costs at the final probe (intact - ablated, median over seeds) ===")
for (pv, pt), lab in CELLS.items():
    row = []
    for c in COND:
        L = [probes(r)[-1] for r in runs(pv, pt)]
        f = lambda k: np.median([x[f"{c}.normal"]["fit"] - x[f"{c}.{k}"]["fit"] for x in L])
        row.append(f"{c}: smell {f('anosmic'):+.2f} vision {f('blind'):+.2f} plast {f('no_plasticity'):+.2f} readouts {f('random_walk'):+.2f}")
    print(f" {lab} | " + " | ".join(row))

print("\n=== (d) population samples vs top-16: mean clone fitness of stored samples / top-16 probe fitness, gen 500 ===")
rat = {c: [] for c in COND}
for (pv, pt), lab in CELLS.items():
    for r in runs(pv, pt):
        f = f"{r}/mut_g000500.npz"
        if not os.path.exists(f): continue
        Z = np.load(f); P = {p["gen"]: p for p in probes(r)}
        for c in COND:
            rat[c].append(Z["clone"][Z["cond"] == c].mean() / P[500][f"{c}.normal"]["fit"])
print(" " + " | ".join(f"{c} median ratio {np.median(v):.3f} [IQR {np.percentile(v,25):.3f}-{np.percentile(v,75):.3f}]" for c, v in rat.items()))
for n in COND[1:]:
    d = np.array(rat["A"]) - np.array(rat[n]); print(f"   A-{n} ratio diff median {np.median(d):+.3f} p {wilcoxon(d).pvalue:.4f}")

