"""Figure 8: causal test of the sensory->motor shortcut (transplant in / remove out)."""
import json, glob, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a86", "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6})
CONDS = ["A", "AS", "N2", "N2R", "N4"]
LBL = {"A": "connectome", "AS": "connectome + shortcuts", "N2": "degree-preserving", "N2R": "degree-preserving − shortcuts", "N4": "interface-preserving"}
COL = {"A": "#1f6feb", "AS": "#8250df", "N2": "#bf8700", "N2R": "#1a7f37", "N4": "#5c5c58"}
LS = {"A": "-", "AS": (0, (4, 1)), "N2": "-.", "N2R": (0, (3, 1, 1, 1)), "N4": ":"}
CELLS = [("0.0", "0.0", "P0T0: safe food"), ("0.2", "0.0", "P1T0: predator")]
def probes(pv, pt):
    return [[json.loads(l) for l in open(f"{r}/probe.jsonl")] for r in sorted(glob.glob(f"runs/swap_P{pv}_T{pt}_s*")) if os.path.exists(f"{r}/done.flag")]
fig = plt.figure(figsize=(10.2, 4.3))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.95], wspace=0.55, left=0.07, right=0.99, top=0.9, bottom=0.30)
for k, (pv, pt, title) in enumerate(CELLS):
    ax = fig.add_subplot(gs[0, k]); D = probes(pv, pt); gens = np.array([p["gen"] for p in D[0]]); ends = {}
    for c in CONDS:
        Y = np.array([[p[f"{c}.normal"]["fit"] for p in P] for P in D]); m = np.median(Y, 0)
        ax.fill_between(gens, np.percentile(Y, 25, 0), np.percentile(Y, 75, 0), color=COL[c], alpha=0.10, lw=0)
        ax.plot(gens, m, ls=LS[c], color=COL[c], lw=1.9, label=LBL[c]); ends[c] = m[-1]
    lo, hi = ax.get_ylim(); gap = 0.055 * (hi - lo); placed = []
    for c in sorted(ends, key=ends.get):
        y = ends[c] if not placed else max(ends[c], placed[-1] + gap); placed.append(y)
        ax.annotate(c, (gens[-1], ends[c]), xytext=(gens[-1] + 10, y), textcoords="data", va="center", fontsize=8)
    ax.set_title(f"{title}  (n={len(D)})", fontsize=9, loc="left"); ax.set_xlabel("generation"); ax.set_xlim(0, gens[-1] + 70)
    if k == 0: ax.set_ylabel("common-garden fitness")
# panel C: paired differences at generation 600, seed means over both ecologies
ax = fig.add_subplot(gs[0, 2])
S = {}
for (pv, pt, _) in CELLS:
    for r in sorted(glob.glob(f"runs/swap_P{pv}_T{pt}_s*")):
        if not os.path.exists(f"{r}/done.flag"): continue
        p = [json.loads(l) for l in open(f"{r}/probe.jsonl")][-1]; s = int(r.split("_s")[-1])
        for c in CONDS: S.setdefault(c, {}).setdefault(s, []).append(p[f"{c}.normal"]["fit"])
M = {c: np.array([np.mean(v) for _, v in sorted(S[c].items())]) for c in CONDS}
PAIRS = [("AS", "A", "shortcuts added\nto connectome"), ("N2", "A", "control with\nshortcuts"),
         ("N2R", "A", "control without\nshortcuts"), ("N4", "A", "interface-preserving\ncontrol")]
rng = np.random.default_rng(0)
for i, (x, y, lab) in enumerate(PAIRS):
    d = M[x] - M[y]; b = np.median(rng.choice(d, (10000, len(d))), 1); lo, hi = np.percentile(b, [5, 95])
    ax.errorbar(np.median(d), i, xerr=[[np.median(d) - lo], [hi - np.median(d)]], fmt="o", color=COL[x], capsize=3, lw=1.8, ms=6)
    ax.scatter(d, np.full(len(d), i) + rng.normal(0, 0.06, len(d)), s=10, color=COL[x], alpha=0.35)
ax.axvline(0, color="#8a8a86", lw=1)
ax.set_yticks(range(len(PAIRS))); ax.set_yticklabels([p[2] for p in PAIRS], fontsize=8); ax.set_ylim(-0.6, len(PAIRS) - 0.4)
ax.set_xlabel("fitness at generation 600 minus connectome\n(seed means, bars: 90% CI of the median)", fontsize=8)
ax.tick_params(axis="y", pad=1)
ax.invert_yaxis()
h, l = fig.axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
fig.savefig("figs/fig8_swap.png", dpi=200); plt.close(fig)
print("wrote figs/fig8_swap.png")
