"""Summary test across the four ecologies (40 runs pooled, each run = its own seed x ecology).
Reported next to the pre-registered per-ecology tests; also shows Bonferroni x4 over ecologies."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon
COND = ["A", "N1", "N2", "N3"]
def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj
R = []
for f in sorted(glob.glob("runs/eco_P*_s*/probe.jsonl")):
    if not os.path.exists(os.path.join(os.path.dirname(f), "done.flag")): continue
    P = [json.loads(l) for l in open(f)]; g = np.array([p["gen"] for p in P], float)
    d = {}
    for c in COND:
        y = np.array([p[f"{c}.normal"]["fit"] for p in P]); last = P[-1]
        d[c] = dict(gen0=y[0], final=y[-1], AUC=np.trapezoid(y, g) / (g[-1] - g[0]),
                    abl_anosmic=last[f"{c}.normal"]["fit"] - last[f"{c}.anosmic"]["fit"])
    R.append(d)
print(f"pooled runs: {len(R)}")
for m in ["gen0", "AUC", "final", "abl_anosmic"]:
    ps, txt = [], []
    for n in COND[1:]:
        x = np.array([r["A"][m] - r[n][m] for r in R]); p = wilcoxon(x).pvalue; ps.append(p)
        txt.append(f"A-{n} median {np.median(x):+.3f} A>null {int((x > 0).sum())}/{len(x)}")
    print(f" {m:12s} " + " | ".join(f"{t} p_holm {a:.2e}" for t, a in zip(txt, holm(ps))))

# seed-level version: runs of the same seed share null realisations across ecologies -> average them first (n = 10)
print("seed-level (mean over the 4 ecologies per seed, n = 10):")
S = {}
for f in sorted(glob.glob("runs/eco_P*_s*/probe.jsonl")):
    s = int(os.path.dirname(f).split("_s")[-1]); P = [json.loads(l) for l in open(f)]
    for c in COND:
        S.setdefault(s, {}).setdefault(c, {"gen0": [], "final": []})
        S[s][c]["gen0"].append(P[0][f"{c}.normal"]["fit"]); S[s][c]["final"].append(P[-1][f"{c}.normal"]["fit"])
for m in ["gen0", "final"]:
    ps, txt = [], []
    for n in COND[1:]:
        x = np.array([np.mean(S[s]["A"][m]) - np.mean(S[s][n][m]) for s in sorted(S)]); p = wilcoxon(x).pvalue; ps.append(p)
        txt.append(f"A-{n} median {np.median(x):+.3f} A>null {int((x > 0).sum())}/{len(x)}")
    print(f" {m:12s} " + " | ".join(f"{t} p_holm {a:.4f}" for t, a in zip(txt, holm(ps))))
