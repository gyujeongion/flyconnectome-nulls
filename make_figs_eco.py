"""Figures 5-6 for paper-v2 (ecology grid). Run on the workstation after run_eco.sh and run_eco_assay.sh."""
import json, glob, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

os.makedirs("figs", exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a86", "axes.labelcolor": "#2b2b2b", "xtick.color": "#5c5c58",
                     "ytick.color": "#5c5c58", "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6})
COND = ["A", "N1", "N2", "N3"]
LBL = {"A": "connectome", "N1": "column shuffle", "N2": "degree-preserving", "N3": "block shuffle"}
COL = {"A": "#1f6feb", "N1": "#d1242f", "N2": "#bf8700", "N3": "#8250df"}
LS = {"A": "-", "N1": "--", "N2": "-.", "N3": ":"}          # secondary encoding, not colour alone
CELLS = [("0.0", "0.0", "no predator, no toxin"), ("0.2", "0.0", "predator"),
         ("0.0", "1.0", "toxin"), ("0.2", "1.0", "predator + toxin")]
SHORT = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}

def probes(pv, pt):
    out = []
    for r in sorted(glob.glob(f"runs/eco_P{pv}_T{pt}_s*")):
        if os.path.exists(f"{r}/done.flag"):
            out.append([json.loads(l) for l in open(f"{r}/probe.jsonl")])
    return out

# ---------- Fig 5: common-garden fitness over 600 generations, one panel per ecology ----------
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True)
for ax, (pv, pt, title) in zip(axes.flat, CELLS):
    D = probes(pv, pt); gens = np.array([p["gen"] for p in D[0]])
    ends = {}
    for c in COND:
        Y = np.array([[p[f"{c}.normal"]["fit"] for p in P] for P in D])
        m = np.median(Y, 0)
        ax.fill_between(gens, np.percentile(Y, 25, 0), np.percentile(Y, 75, 0), color=COL[c], alpha=0.12, lw=0)
        ax.plot(gens, m, LS[c], color=COL[c], lw=2, label=LBL[c])
        ends[c] = m[-1]
    # end labels, pushed apart so they never overlap
    lo, hi = ax.get_ylim(); gap = 0.06 * (hi - lo); placed = []
    for c in sorted(ends, key=ends.get):
        y = ends[c] if not placed else max(ends[c], placed[-1] + gap); placed.append(y)
        ax.annotate(c, (gens[-1], ends[c]), xytext=(gens[-1] + 8, y), textcoords="data", va="center", fontsize=8, color="#2b2b2b")
    ax.set_title(f"{SHORT[(pv, pt)]}: {title}  (n={len(D)})", fontsize=9, loc="left")
    ax.set_xlim(0, gens[-1] + 45)
for ax in axes[1]: ax.set_xlabel("generation")
for ax in axes[:, 0]: ax.set_ylabel("common-garden fitness")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig("figs/fig5_eco.png", dpi=200); plt.close(fig)

# ---------- Fig 6: cross-ecology common garden (generation-500 elites) ----------
TEST = ["eco_P0.0_T0.0", "eco_P0.2_T0.0", "eco_P0.0_T1.0", "eco_P0.2_T1.0"]
TL = ["P0T0", "P1T0", "P0T1", "P1T1"]
M = {n: np.full((4, 4), np.nan) for n in COND[1:]}
for i, (pv, pt, _) in enumerate(CELLS):
    per = []
    for r in sorted(glob.glob(f"runs/eco_P{pv}_T{pt}_s*")):
        last = {}
        if os.path.exists(f"{r}/assay.jsonl"):
            for l in open(f"{r}/assay.jsonl"):
                d = json.loads(l)
                if d["gen"] == 500 and d["assay"] in TEST: last[d["assay"]] = d["res"]
        if len(last) == 4: per.append(last)
    for n in COND[1:]:
        for j, t in enumerate(TEST):
            if per: M[n][i, j] = np.median([x[t]["A"]["fit"] - x[t][n]["fit"] for x in per])
cmap = LinearSegmentedColormap.from_list("div", ["#d1242f", "#ececea", "#1f6feb"])
lim = np.nanmax([np.nanmax(np.abs(v)) for v in M.values()])
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.3), gridspec_kw=dict(wspace=0.08))
for ax, n in zip(axes, COND[1:]):
    im = ax.imshow(M[n], cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(4):
        for j in range(4):
            v = M[n][i, j]
            ax.text(j, i, "" if np.isnan(v) else f"{v:+.2f}", ha="center", va="center", fontsize=7.5,
                    color="#ffffff" if abs(v) > 0.6 * lim else "#1b1b1b", fontweight="bold" if i == j else "normal")
    ax.set_xticks(range(4)); ax.set_xticklabels(TL, fontsize=8); ax.set_yticks(range(4))
    ax.set_yticklabels(TL if n == "N1" else [], fontsize=8)
    ax.set_xlabel("tested in"); ax.grid(False)
    ax.set_title(f"A − {n} ({LBL[n]})", fontsize=9)
axes[0].set_ylabel("evolved in")
cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02); cb.set_label("median fitness difference\n(blue: connectome ahead)", fontsize=8)
fig.savefig("figs/fig6_cross.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("wrote figs/fig5_eco.png figs/fig6_cross.png")
