"""Figure 9: dose-response to the sensory->motor shortcut, with a rewiring-matched sham."""
import json, glob, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a86", "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6})

CONDS = ["A", "AS1", "AS3", "AS5", "AS10"]
SHARE = {"A": 0.012, "AS1": 1.0, "AS3": 3.0, "AS5": 5.0, "AS10": 10.0}
COL = {"A": "#1f6feb", "AS1": "#7aa7f0", "AS3": "#a97fd6", "AS5": "#8250df", "AS10": "#5a1fa8",
       "SHAM": "#1a7f37"}
CELLS = [("0.0", "0.0", "P0T0: safe food"), ("0.2", "0.0", "P1T0: predator")]


def runs(pv, pt):
    return [r for r in sorted(glob.glob(f"runs/dose_P{pv}_T{pt}_s*"))
            if os.path.exists(f"{r}/done.flag")]


def last_probe(r):
    p = None
    for line in open(f"{r}/probe.jsonl"):
        p = line
    return json.loads(p)


fig = plt.figure(figsize=(10.2, 4.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.42,
                      left=0.07, right=0.985, top=0.89, bottom=0.30)

# ---- panels A,B: final fitness against shortcut share, per ecology --------
for k, (pv, pt, title) in enumerate(CELLS):
    ax = fig.add_subplot(gs[0, k])
    R = runs(pv, pt)
    F = {c: np.array([last_probe(r)[f"{c}.normal"]["fit"] for r in R]) for c in CONDS + ["SHAM"]}
    xs = [SHARE[c] for c in CONDS]
    for r_i in range(len(R)):                     # seed trajectories
        ax.plot(xs, [F[c][r_i] for c in CONDS], color="#b9b9b4", lw=0.7, alpha=0.8, zorder=1)
    m = [F[c].mean() for c in CONDS]
    se = [F[c].std(ddof=1) / np.sqrt(len(F[c])) for c in CONDS]
    ax.errorbar(xs, m, yerr=se, color="#5a1fa8", lw=2.0, marker="o", ms=5,
                capsize=3, zorder=3, label="transplanted shortcuts")
    sm, sse = F["SHAM"].mean(), F["SHAM"].std(ddof=1) / np.sqrt(len(F["SHAM"]))
    ax.errorbar([SHARE["A"]], [sm], yerr=[sse], color=COL["SHAM"], marker="s", ms=6,
                capsize=3, zorder=4, label="sham (rewiring matched)")
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xticks([0.012, 1, 3, 5, 10])
    ax.set_xticklabels(["0.012\n(connectome)", "1", "3", "5", "10"])
    ax.set_xlabel("olfactory output on descending groups (%)")
    ax.set_title(f"{title}  (n={len(R)} seeds)", fontsize=9, loc="left")
    if k == 0:
        ax.set_ylabel("fitness at generation 600")
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")

# ---- panel C: cost of removing olfaction ---------------------------------
ax = fig.add_subplot(gs[0, 2])
for k, (pv, pt, title) in enumerate(CELLS):
    R = runs(pv, pt)
    C = {c: np.array([last_probe(r)[f"{c}.normal"]["fit"] - last_probe(r)[f"{c}.anosmic"]["fit"]
                      for r in R]) for c in CONDS}
    xs = [SHARE[c] for c in CONDS]
    m = [C[c].mean() for c in CONDS]
    se = [C[c].std(ddof=1) / np.sqrt(len(C[c])) for c in CONDS]
    ax.errorbar(xs, m, yerr=se, lw=1.9, marker="o" if k == 0 else "^", ms=5, capsize=3,
                color="#1f6feb" if k == 0 else "#bf8700",
                ls="-" if k == 0 else (0, (4, 1)), label=title.split(":")[0])
ax.set_xscale("symlog", linthresh=1.0)
ax.set_xticks([0.012, 1, 3, 5, 10])
ax.set_xticklabels(["0.012\n(connectome)", "1", "3", "5", "10"])
ax.set_xlabel("olfactory output on descending groups (%)")
ax.set_ylabel("fitness lost when olfaction is removed")
ax.set_title("dependence on smell", fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=7.5, loc="upper left")

os.makedirs("figs", exist_ok=True)
fig.savefig("figs/fig9_dose.png", dpi=200)
print("wrote figs/fig9_dose.png")
