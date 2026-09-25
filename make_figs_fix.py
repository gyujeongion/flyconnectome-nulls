"""Figure 7: corrected grid (defect 28 fixed) with shortcut-carrying and interface-preserving controls."""
import json, glob, os, numpy as np
from seeds import registered  # registered grids are seeds 0-9; see seeds.py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#8a8a86",
                     "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6})
CONDS = ["A", "N1", "N2", "N4", "N5"]
LBL = {"A": "connectome", "N1": "column shuffle", "N2": "degree-preserving", "N4": "interface-preserving swap", "N5": "interface-preserving shuffle"}
COL = {"A": "#1f6feb", "N1": "#d1242f", "N2": "#bf8700", "N4": "#1a7f37", "N5": "#8250df"}
LS = {"A": "-", "N1": "--", "N2": "-.", "N4": (0, (3, 1, 1, 1)), "N5": ":"}
CELLS = [("0.0", "0.0", "P0T0: stationary predator, safe food"), ("0.2", "0.0", "P1T0: predator"),
         ("0.0", "1.0", "P0T1: toxin"), ("0.2", "1.0", "P1T1: predator + toxin")]
def probes(pv, pt):
    return [[json.loads(l) for l in open(f"{r}/probe.jsonl")] for r in sorted(registered(glob.glob(f"runs/ecofix_P{pv}_T{pt}_s*"))) if os.path.exists(f"{r}/done.flag")]
fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), sharex=True)
for ax, (pv, pt, title) in zip(axes.flat, CELLS):
    D = probes(pv, pt); gens = np.array([p["gen"] for p in D[0]])
    ends = {}
    for c in CONDS:
        Y = np.array([[p[f"{c}.normal"]["fit"] for p in P] for P in D]); m = np.median(Y, 0)
        ax.fill_between(gens, np.percentile(Y, 25, 0), np.percentile(Y, 75, 0), color=COL[c], alpha=0.10, lw=0)
        ax.plot(gens, m, ls=LS[c], color=COL[c], lw=1.9, label=LBL[c]); ends[c] = m[-1]
    lo, hi = ax.get_ylim(); gap = 0.055 * (hi - lo); placed = []
    for c in sorted(ends, key=ends.get):
        y = ends[c] if not placed else max(ends[c], placed[-1] + gap); placed.append(y)
        ax.annotate(c, (gens[-1], ends[c]), xytext=(gens[-1] + 8, y), textcoords="data", va="center", fontsize=8, color="#2b2b2b")
    ax.set_title(f"{title}  (n={len(D)})", fontsize=9, loc="left"); ax.set_xlim(0, gens[-1] + 55)
for ax in axes[1]: ax.set_xlabel("generation")
for ax in axes[:, 0]: ax.set_ylabel("common-garden fitness")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig("figs/fig7_fix.png", dpi=200); plt.close(fig)
print("wrote figs/fig7_fix.png")
