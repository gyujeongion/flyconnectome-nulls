"""PROTOCOL_FIX analysis: corrected turning noise + interface-preserving controls N4/N5.
Registered inference at seed level (n=10, per ecology and averaged over ecologies)."""
import json, glob, os, numpy as np
from seeds import registered  # registered grids are seeds 0-9; see seeds.py
from scipy.stats import wilcoxon
CONDS = ["A", "N1", "N2", "N4", "N5"]
CELLS = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}
rng = np.random.default_rng(0)
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj
def ci(d):
    b = np.median(rng.choice(d, (10000, len(d))), 1); return np.percentile(b, 2.5), np.percentile(b, 97.5)
def runs(pv, pt): return sorted(r for r in registered(glob.glob(f"runs/ecofix_P{pv}_T{pt}_s*")) if os.path.exists(f"{r}/done.flag"))
def P(r): return [json.loads(l) for l in open(f"{r}/probe.jsonl")]
def endp(r):
    p = P(r); g = np.array([x["gen"] for x in p], float); out = {}
    for c in CONDS:
        y = np.array([x[f"{c}.normal"]["fit"] for x in p]); last = p[-1]
        out[c] = dict(gen0=y[0], final=y[-1], AUC=float(np.trapezoid(y, g) / (g[-1] - g[0])),
                      smell=last[f"{c}.normal"]["fit"] - last[f"{c}.anosmic"]["fit"],
                      vision=last[f"{c}.normal"]["fit"] - last[f"{c}.blind"]["fit"])
    return out
E = {}
for (pv, pt), lab in CELLS.items():
    for r in runs(pv, pt): E[(lab, int(r.split("_s")[-1]))] = endp(r)
print("runs:", len(E), {lab: sum(1 for k in E if k[0] == lab) for lab in CELLS.values()})
for m in ["gen0", "AUC", "final", "smell", "vision"]:
    print(f"\n=== {m} (connectome minus control; median over 10 seeds [95% CI], p_holm over 4 controls) ===")
    for lab in list(CELLS.values()) + ["ALL"]:
        keys = [k for k in E if lab == "ALL" or k[0] == lab]
        if lab == "ALL":
            seeds = sorted({k[1] for k in E}); vals = {c: np.array([np.mean([E[(l, s)][c][m] for l in CELLS.values() if (l, s) in E]) for s in seeds]) for c in CONDS}
        else:
            vals = {c: np.array([E[k][c][m] for k in sorted(keys)]) for c in CONDS}
        ps, txt = [], []
        for n in CONDS[1:]:
            d = vals["A"] - vals[n]; ps.append(wilcoxon(d).pvalue if np.any(d != 0) else 1.0)
            lo, hi = ci(d); txt.append(f"{n} {np.median(d):+.3f}[{lo:+.2f},{hi:+.2f}] {int((d > 0).sum())}/{len(d)}")
        adj = holm(ps)
        lv = " ".join(f"{c} {np.median(vals[c]):.2f}" for c in CONDS)
        print(f" {lab:5s} | {lv} || " + " | ".join(f"{t} p{a:.3f}" for t, a in zip(txt, adj)))
print("\n=== crossing: first probe after which the control stays >= connectome at every later probe ===")
for lab in CELLS.values():
    out = []
    for n in CONDS[1:]:
        gs = []
        for k in sorted(k for k in E if k[0] == lab):
            pv, pt = [x for x, l in CELLS.items() if l == lab][0]
            r = [x for x in runs(pv, pt) if x.endswith(f"_s{k[1]}")][0]
            p = P(r); d = np.array([x[f"{n}.normal"]["fit"] - x["A.normal"]["fit"] for x in p]); g = [x["gen"] for x in p]
            i = next((i for i in range(1, len(d)) if np.all(d[i:] >= 0)), None)
            gs.append(g[i] if i is not None else np.inf)
        gs = np.array(gs, float); fin = np.isfinite(gs)
        out.append(f"{n} {fin.sum()}/{len(gs)} med {int(np.median(gs[fin])) if fin.any() else '-'}")
    print(f" {lab}: " + " | ".join(out))
