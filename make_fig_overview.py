"""Figure 1 of the NeurIPS-style preprint: what each null preserves at the sensory-motor boundary (a),
and the headline comparisons (b). Every number in panel b is read from a results file."""
import json, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.rcParams.update({"font.size": 8, "font.family": "serif", "axes.spines.top": False, "axes.spines.right": False})
INK, MUTED, GRID = "#2b2b2b", "#8a8a86", "#e6e6e3"
C = {"conn": "#1f6feb", "short": "#cf222e", "bound": "#0f7b6c", "interv": "#8250df", "sham": "#1a7f37"}

# ------------------------------------------------------------------ numbers
eq = open("results/equivalence.txt").read()
def pooled(block, label):
    sec = eq.split(block, 1)[1]
    sec = sec.split("pooled (seed means)", 1)[1]
    m = re.search(re.escape(label) + r"\s+median ([+-]\d\.\d+).*?90% CI \[([+-]\d\.\d+),([+-]\d\.\d+)\]", sec)
    return tuple(map(float, m.groups()))
rows = [("connectome − N1 (column shuffle)", *pooled("=== corrected grid", "A - N1"), C["short"]),
        ("connectome − N2 (degree-preserving)", *pooled("=== corrected grid", "A - N2"), C["short"]),
        ("connectome − N4 (interior swaps)", *pooled("=== corrected grid", "A - N4"), C["bound"]),
        ("connectome − N5 (interior shuffle)", *pooled("=== corrected grid", "A - N5"), C["bound"]),
        ("shortcut transplant − connectome", *pooled("=== transplant experiment", "AS - A"), C["interv"])]
pn = json.load(open("results/paper_numbers.json"))["transplant_vs_sham_by_seed"]
rows.append(("full dose − matched sham", pn["mean_gap_all"], *pn["ci90_mean_gap"], C["interv"]))
d2 = json.load(open("results/dose_key.json"))["D2"]["pooled"]
rows.append(("matched sham − connectome", d2["mean"], *d2["ci"], C["sham"]))
bs = json.load(open("results/bsham_key.json"))                  # boundary-targeted sham (PROTOCOL_BSHAM)
rows.append(("full dose − boundary sham", bs["P"]["mean"], *bs["P"]["ci90"], C["interv"]))
rows.append(("boundary sham − connectome", bs["S3"]["mean"], *bs["S3"]["ci90"], C["sham"]))

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(5.5, 3.95))
axa = fig.add_axes([0.0, 0.52, 1.0, 0.46]); axa.set_xlim(-0.05, 3.12); axa.set_ylim(-0.30, 1.13); axa.axis("off")
axb = fig.add_axes([0.40, 0.10, 0.57, 0.37])

def net(x0, title, kind, labels=False):
    S = [(x0 + 0.20, y) for y in (0.80, 0.50, 0.20)]          # sensory: smell, taste, vision
    I = [(x0 + 0.45, 0.88), (x0 + 0.45, 0.50), (x0 + 0.45, 0.12), (x0 + 0.60, 0.69), (x0 + 0.60, 0.31)]
    M = [(x0 + 0.84, y) for y in (0.80, 0.50, 0.20)]          # descending (motor)
    def arrow(a, b, col, ls="-", lw=0.8, rad=0.0):
        axa.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=5.5, color=col, lw=lw, ls=ls,
                                      shrinkA=5, shrinkB=5.5, connectionstyle=f"arc3,rad={rad}", zorder=2))
    b_in = [(S[0], I[0]), (S[0], I[1]), (S[1], I[1]), (S[2], I[2])]
    b_out = [(I[3], M[0]), (I[1], M[1]), (I[4], M[2]), (I[2], M[2])]
    interior = {"conn": [(I[0], I[3]), (I[1], I[4]), (I[2], I[4]), (I[0], I[1])],
                "short": [(I[0], I[4]), (I[2], I[3]), (I[1], I[3])],
                "bound": [(I[0], I[4]), (I[2], I[3]), (I[1], I[0]), (I[3], I[4])]}[kind]
    for a, b in b_in + b_out:
        if kind == "short" and (a, b) in ((S[0], I[0]), (I[3], M[0])):
            continue
        arrow(a, b, INK if kind != "short" else MUTED, lw=0.9)
    for a, b in interior:
        arrow(a, b, C["bound"] if kind == "bound" else MUTED, ls=(0, (2, 1.4)) if kind == "bound" else "-")
    if kind == "short":
        arrow(S[0], M[0], C["short"], lw=1.5, rad=-0.35); arrow(S[0], M[1], C["short"], lw=1.5, rad=-0.12)
    smell_col = C["short"] if kind == "short" else INK
    axa.scatter(*zip(*S), s=46, c="white", edgecolors=[smell_col, INK, INK], linewidths=1.0, zorder=3)
    axa.scatter(*zip(*I), s=22, c="#f1f1ef", edgecolors=MUTED, linewidths=0.7, zorder=3)
    axa.scatter(*zip(*M), s=46, c="white", edgecolors=INK, linewidths=1.0, zorder=3)
    if labels:
        for (x, y), lab in zip(S, ("smell", "taste", "vision")):
            axa.text(x - 0.06, y, lab, ha="right", va="center", fontsize=6.2)
        axa.text(x0 + 0.20, 1.0, "sensory", ha="center", fontsize=5.8, color=MUTED)
        axa.text(x0 + 0.525, 1.0, "interior", ha="center", fontsize=5.8, color=MUTED)
        axa.text(x0 + 0.84, 1.0, "motor", ha="center", fontsize=5.8, color=MUTED)
    axa.text(x0 + 0.52, 1.10, title, ha="center", va="center", fontsize=7.0, fontweight="bold")
    sub = {"conn": "direct smell → motor: 0.012 %\n52 of 53 motor groups ≥ 2 synapses",
           "short": "direct smell → motor: 10.6–10.7 %\nevery motor group 1 synapse",
           "bound": "boundary edges identical\n~90 % of interior replaced"}[kind]
    axa.text(x0 + 0.52, -0.17, sub, ha="center", va="center", fontsize=6.0, color=INK, linespacing=1.2)

net(0.00, "connectome", "conn", labels=True)
net(1.02, "standard nulls (N1, N2)", "short")
net(2.04, "boundary-preserving (N4, N5)", "bound")
fig.text(0.01, 0.985, "a", fontsize=10, fontweight="bold", va="top")

y = np.arange(len(rows))[::-1]
axb.axvspan(-0.10, 0.10, color=GRID, zorder=0)
axb.axvline(0, color=MUTED, lw=0.6, zorder=1)
for yi, (lab, m, lo, hi, col) in zip(y, rows):
    axb.plot([lo, hi], [yi, yi], color=col, lw=1.6, solid_capstyle="round", zorder=2)
    axb.plot([m], [yi], "o", ms=3.6, color=col, mec="white", mew=0.6, zorder=3)
axb.set_yticks(y); axb.set_yticklabels([r[0] for r in rows], fontsize=6.4)
axb.set_xlabel("fitness difference at the last probe (90 % CI)", fontsize=6.4)
axb.tick_params(axis="x", labelsize=6.0); axb.tick_params(axis="y", length=0)
axb.set_xlim(-0.45, 0.85); axb.set_ylim(-0.6, len(rows) + 0.15)
axb.text(0.0, len(rows) - 0.25, "±0.10 bound", ha="center", fontsize=5.6, color=MUTED)
fig.text(0.01, 0.49, "b", fontsize=10, fontweight="bold", va="top")
fig.savefig("figs/fig_overview.png", dpi=300)
fig.savefig("figs/fig_overview.pdf")
print("wrote figs/fig_overview.png/.pdf")
for r in rows:
    print(f"  {r[0]:38s} {r[1]:+.3f} [{r[2]:+.3f}, {r[3]:+.3f}]")
