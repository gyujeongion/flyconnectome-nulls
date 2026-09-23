"""PROTOCOL_SWAP analysis: causal test of the sensory->motor shortcut.
AS = connectome + transplanted shortcuts, N2R = degree-preserving control with shortcuts removed."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon
CONDS = ["A", "AS", "N2", "N2R", "N4"]
CELLS = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0"}
rng = np.random.default_rng(0)
def ci(d):
    b = np.median(rng.choice(d, (10000, len(d))), 1); return np.percentile(b, 2.5), np.percentile(b, 97.5)
def runs(pv, pt): return sorted(r for r in glob.glob(f"runs/swap_P{pv}_T{pt}_s*") if os.path.exists(f"{r}/done.flag"))
def endp(r):
    p = [json.loads(l) for l in open(f"{r}/probe.jsonl")]; g = np.array([x["gen"] for x in p], float); out = {}
    for c in CONDS:
        y = np.array([x[f"{c}.normal"]["fit"] for x in p]); last = p[-1]
        out[c] = dict(gen0=y[0], final=y[-1], AUC=float(np.trapezoid(y, g) / (g[-1] - g[0])),
                      smell=last[f"{c}.normal"]["fit"] - last[f"{c}.anosmic"]["fit"])
    return out
E = {}
for (pv, pt), lab in CELLS.items():
    for r in runs(pv, pt): E[(lab, int(r.split("_s")[-1]))] = endp(r)
print("runs:", len(E), {lab: sum(1 for k in E if k[0] == lab) for lab in CELLS.values()})
PAIRS = [("AS", "A", "T1 transplanting shortcuts into the connectome"),
         ("N2R", "N2", "T2a removing shortcuts from the degree-preserving control"),
         ("N2R", "A", "T2b control without shortcuts vs connectome"),
         ("N2", "A", "reference: control with shortcuts vs connectome"),
         ("N4", "A", "reference: interface-preserving control vs connectome"),
         ("AS", "N2", "reference: connectome+shortcuts vs control with shortcuts")]
for m in ["gen0", "AUC", "final", "smell"]:
    print(f"\n=== {m} ===")
    for lab in list(CELLS.values()) + ["BOTH"]:
        if lab == "BOTH":
            seeds = sorted({k[1] for k in E}); V = {c: np.array([np.mean([E[(l, s)][c][m] for l in CELLS.values() if (l, s) in E]) for s in seeds]) for c in CONDS}
        else:
            ks = sorted(k for k in E if k[0] == lab); V = {c: np.array([E[k][c][m] for k in ks]) for c in CONDS}
        print(f" {lab:5s} levels: " + " ".join(f"{c} {np.median(V[c]):.2f}" for c in CONDS))
        for x, y, why in PAIRS:
            d = V[x] - V[y]; p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0; lo, hi = ci(d)
            print(f"    {x}-{y:4s} {np.median(d):+.3f} [{lo:+.2f},{hi:+.2f}] {int((d > 0).sum())}/{len(d)} p={p:.4f}   {why if lab == 'BOTH' else ''}")
