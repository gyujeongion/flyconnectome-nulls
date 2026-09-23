"""Figures for the paper. Run on the workstation after all runs finish: ./venv/bin/python make_figs.py"""
import json, glob, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

os.makedirs("figs", exist_ok=True)
COND = ["A", "N1", "N2", "N3"]
LBL = {"A": "connectome", "N1": "column shuffle", "N2": "degree-preserving", "N3": "block shuffle"}
COL = {"A": "#1f6feb", "N1": "#d1242f", "N2": "#bf8700", "N3": "#8250df"}
SCHEMES = ["global", "pergroup", "distmatch"]
rng = np.random.default_rng(0)

def load_probe(scheme):
    D = {}
    for f in sorted(glob.glob(f"runs/paper_{scheme}_s*/probe.jsonl")):
        run = os.path.dirname(f)
        if not os.path.exists(f"{run}/done.flag"): continue
        D[int(run.split("_s")[-1])] = [json.loads(l) for l in open(f)]
    return D

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."

# ---------------- Fig 1: learning curves (common-garden probe fitness) ----------------
avail = [s for s in SCHEMES if load_probe(s)]
fig, axes = plt.subplots(1, len(avail), figsize=(5.2 * len(avail), 4), squeeze=False)
for ax, scheme in zip(axes[0], avail):
    D = load_probe(scheme); gens = np.array([p["gen"] for p in next(iter(D.values()))])
    for c in COND:
        Y = np.array([[p[f"{c}.normal"]["fit"] for p in P] for P in D.values()])
        m = np.median(Y, 0); lo = np.percentile(Y, 25, 0); hi = np.percentile(Y, 75, 0)
        ax.plot(gens, m, color=COL[c], label=f"{LBL[c]} (n={len(Y)})"); ax.fill_between(gens, lo, hi, color=COL[c], alpha=0.15)
    ax.set_xlabel("generation"); ax.set_ylabel("common-garden fitness"); ax.set_title(f"calibration: {scheme}"); ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig("figs/fig1_curves.png", dpi=150); plt.close()

# ---------------- Fig 2: pre-registered endpoints, paired differences ----------------
def endpoints(scheme):
    D = load_probe(scheme); out = {}
    for s, P in D.items():
        gens = np.array([p["gen"] for p in P])
        for c in COND:
            y = np.array([p[f"{c}.normal"]["fit"] for p in P]); last = P[-1]
            out.setdefault(c, {})[s] = dict(AUC=float(np.trapezoid(y, gens) / (gens[-1] - gens[0])), gen0=float(y[0]), final=float(y[-1]),
                                            abl_anosmic=float(last[f"{c}.normal"]["fit"] - last[f"{c}.anosmic"]["fit"]),
                                            abl_blind=float(last[f"{c}.normal"]["fit"] - last[f"{c}.blind"]["fit"]),
                                            abl_noplast=float(last[f"{c}.normal"]["fit"] - last[f"{c}.no_plasticity"]["fit"]))
    return out
METS = [("gen0", "generation 0"), ("AUC", "AUC (primary)"), ("final", "generation 275"),
        ("abl_anosmic", "olfaction ablation"), ("abl_blind", "vision ablation"), ("abl_noplast", "plasticity ablation")]
fig, axes = plt.subplots(len(avail), len(METS), figsize=(3.0 * len(METS), 3.1 * len(avail)), squeeze=False)
for r, scheme in enumerate(avail):
    E = endpoints(scheme); seeds = sorted(E["A"])
    for k, (mk, mt) in enumerate(METS):
        ax = axes[r][k]
        for i, c in enumerate(COND):
            v = np.array([E[c][s][mk] for s in seeds])
            ax.scatter(np.full(len(v), i) + rng.normal(0, 0.06, len(v)), v, s=12, color=COL[c], alpha=0.65)
            ax.hlines(np.median(v), i - 0.28, i + 0.28, color=COL[c], lw=2.5)
        A = np.array([E["A"][s][mk] for s in seeds]); top = ax.get_ylim()[1]
        for i, c in enumerate(COND[1:], start=1):
            d = A - np.array([E[c][s][mk] for s in seeds])
            p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
            ax.text(i, top, star(min(1, p * 3)), ha="center", va="top", fontsize=8)
        ax.set_xticks(range(4)); ax.set_xticklabels(["A", "N1", "N2", "N3"], fontsize=8)
        if k == 0: ax.set_ylabel(f"{scheme}\nfitness", fontsize=9)
        if r == 0: ax.set_title(mt, fontsize=10)
plt.tight_layout(); plt.savefig("figs/fig2_endpoints.png", dpi=150); plt.close()

# ---------------- Fig 3: pre-evolution circuit metrics across calibration schemes ----------------
if os.path.exists("runs/sens/circuit.jsonl"):
    R = [json.loads(l) for l in open("runs/sens/circuit.jsonl")]
    spec = lambda v: float(v.split()[0])
    MET = [("cos_KC", "KC odor similarity\n(lower = better)", lambda r: r["cos_KC"]),
           ("DAN_sugar", "dopamine response\nto sugar", lambda r: r["DAN_sugar"]),
           ("spec", "MB learning\nspecificity", lambda r: spec(r["spec_eta-0.05"])),
           ("loom", "looming escape turn\n(negative = correct)", lambda r: r["turn_loomL_minus_loomR"])]
    schemes = sorted({r["homeo"] for r in R})
    fig, axes = plt.subplots(1, len(MET), figsize=(4.0 * len(MET), 3.4))
    for ax, (key, title, f) in zip(axes, MET):
        for i, sch in enumerate(schemes):
            for j, c in enumerate(COND):
                v = np.array([f(r["res"]) for r in R if r["homeo"] == sch and r["cond"] == c])
                if not len(v): continue
                x = i + (j - 1.5) * 0.18
                ax.scatter(np.full(len(v), x) + rng.normal(0, 0.02, len(v)), v, s=10, color=COL[c], alpha=0.6)
                ax.hlines(np.median(v), x - 0.07, x + 0.07, color=COL[c], lw=2)
        ax.set_xticks(range(len(schemes))); ax.set_xticklabels(schemes, fontsize=8); ax.set_title(title, fontsize=9)
    handles = [plt.Line2D([], [], color=COL[c], marker="o", ls="", label=LBL[c]) for c in COND]
    axes[-1].legend(handles=handles, fontsize=7)
    plt.tight_layout(); plt.savefig("figs/fig3_calibration_sensitivity.png", dpi=150); plt.close()

# ---------------- Fig 4: stress grid on evolved elites ----------------
rows = [json.loads(l) for f in glob.glob("runs/paper_*/assay.jsonl") for l in open(f)]
TAGS = ["base", "stress_pred0.12", "stress_pred0.15", "stress_pred0.18", "stress_toxic0.7", "stress_toxic1.0"]
sch_av = [s for s in SCHEMES if any(r["homeo"] == s for r in rows)]
if rows and sch_av:
    fig, axes = plt.subplots(1, len(sch_av), figsize=(5.4 * len(sch_av), 3.8), squeeze=False)
    for ax, scheme in zip(axes[0], sch_av):
        D = {}
        for r in rows:
            if r["homeo"] == scheme: D.setdefault(r["seed"], {})[r["assay"]] = r["res"]
        seeds = sorted(s for s, v in D.items() if all(t in v for t in TAGS))
        for c in COND:
            m = [np.median([D[s][t][c]["fit"] for s in seeds]) for t in TAGS]
            ax.plot(range(len(TAGS)), m, "o-", color=COL[c], label=LBL[c])
        ax.set_xticks(range(len(TAGS))); ax.set_xticklabels([t.replace("stress_", "") for t in TAGS], rotation=30, fontsize=8)
        ax.set_ylabel("fitness of gen-250 elites"); ax.set_title(f"stress transfer ({scheme}, n={len(seeds)})", fontsize=10); ax.legend(fontsize=7)
    plt.tight_layout(); plt.savefig("figs/fig4_stress.png", dpi=150); plt.close()
print("figures written to figs/:", sorted(os.listdir("figs")))
