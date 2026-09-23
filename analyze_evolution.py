"""What did evolution actually change? (post hoc, uses saved elite genomes at generation 0 and 250)

Q1 Do randomised networks drift back towards connectome-like structure?
Q2 Did adaptation happen in the wiring or in the other genes (read-outs, gains, learning rates)?
Q3 Did the circuits that changed differ between connectome and nulls (visual vs olfactory pathways)?
"""
import json, glob, os, sys, numpy as np
from scipy.stats import wilcoxon

B = np.load("data/brain.npz", allow_pickle=True)
names = [str(n) for n in B["names"]]; K = len(names)
W0A = B["W0"].astype(np.float32)
side = np.array([n.split("|")[1] if "|" in n else "C" for n in names])
orn = set(int(i) for i in B["in_orn"]); sug = set(int(i) for i in B["in_sugar"]) | set(int(i) for i in B["in_bitter"])
loom = set(int(i) for i in B["in_loomL"]) | set(int(i) for i in B["in_loomR"])
mb = set(int(i) for i in B["mbon"]) | set(int(i) for i in B["dan"]); dn = set(int(i) for i in B["dn"])
pn = set(int(i) for i in B["pn"]) | set(int(i) for i in B["ln"])
PATH = {"olfactory (ORN/PN/LN)": orn | pn, "mushroom body (MBON/DAN)": mb, "visual (LC4/LPLC2)": loom,
        "taste (GRN)": sug, "descending": dn}
PATH["other"] = set(range(K)) - set().union(*PATH.values())
GENES = ["size", "alpha", "bias", "eta", "gin", "kcg", "kcl", "R", "RM", "Rb"]
COND = ["A", "N1", "N2", "N3"]

def mirror_mass(W):
    """ipsilateral vs contralateral weight mass (laterality of the evolved network)"""
    L = side == "L"; R = side == "R"
    ipsi = np.abs(W[np.ix_(L, L)]).sum() + np.abs(W[np.ix_(R, R)]).sum()
    contra = np.abs(W[np.ix_(L, R)]).sum() + np.abs(W[np.ix_(R, L)]).sum()
    return float(ipsi / (ipsi + contra + 1e-9))

def sim_to_connectome(W):
    """overlap of the support of W with the real connectome's support (Jaccard) and sign agreement"""
    a = W != 0; b = W0A != 0
    jac = float((a & b).sum() / ((a | b).sum() + 1e-9))
    both = a & b
    sgn = float((np.sign(W[both]) == np.sign(W0A[both])).mean()) if both.sum() else float("nan")
    return jac, sgn

def recip(W):
    a = W != 0
    return float((a & a.T).sum() / (a.sum() + 1e-9))

def load(scheme):
    out = {}
    for run in sorted(glob.glob(f"runs/paper_{scheme}_s*")):
        if not os.path.exists(f"{run}/done.flag"): continue
        f0, f1 = f"{run}/elite_g000000.npz", f"{run}/elite_g000250.npz"
        if not (os.path.exists(f0) and os.path.exists(f1)): continue
        out[int(run.split("_s")[-1])] = (np.load(f0), np.load(f1))
    return out

def holm(ps):
    ps = np.asarray(ps); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0
    for r, i in enumerate(o): run = max(run, (m - r) * ps[i]); adj[i] = min(1, run)
    return adj

def report(tag, vals):
    """vals: {cond: array over seeds} -> print medians and A-vs-null tests"""
    A = np.asarray(vals["A"]); ps, txt = [], []
    for n in COND[1:]:
        d = A - np.asarray(vals[n]); p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
        ps.append(p); txt.append((n, np.median(vals[n]), np.median(d), int((d > 0).sum()), len(d)))
    adj = holm(ps)
    print(f"  {tag:42s} A {np.median(A):+.4f} | " + " | ".join(
        f"{n} {m:+.4f} diff {dm:+.4f} {g}/{k} p={pa:.3f}" for (n, m, dm, g, k), pa in zip(txt, adj)))

for scheme in ["global", "pergroup", "distmatch"]:
    D = load(scheme)
    if not D: continue
    seeds = sorted(D)
    print(f"\n================ evolution trajectory, {scheme} (seeds {seeds}) ================")
    # Q2: how much each gene changed (relative to its own initial spread)
    print(" Q2. where did adaptation happen? (median |Δ| / initial SD, over elites and seeds)")
    for g in GENES + ["W"]:
        vals = {}
        for c in COND:
            v = []
            for s in seeds:
                Z0, Z1 = D[s]
                a0, a1 = Z0[f"{c}_{g}"], Z1[f"{c}_{g}"]
                sd = a0.std() + 1e-9
                v.append(float(np.abs(a1.mean(0) - a0.mean(0)).mean() / sd))
            vals[c] = np.array(v)
        report(f"Δ{g}", vals)
    # Q1: do nulls drift towards connectome structure?
    print(" Q1. structure of the evolved network (gen 250) and change from gen 0")
    for metric, f in [("support overlap with connectome (Jaccard)", lambda W: sim_to_connectome(W)[0]),
                      ("sign agreement on shared edges", lambda W: sim_to_connectome(W)[1]),
                      ("ipsilateral weight fraction", mirror_mass),
                      ("reciprocal-edge fraction", recip)]:
        for when, ix in (("gen0", 0), ("gen250", 1)):
            vals = {c: np.array([f(D[s][ix][f"{c}_W"].mean(0)) for s in seeds]) for c in COND}
            report(f"{metric} [{when}]", vals)
        vals = {c: np.array([f(D[s][1][f"{c}_W"].mean(0)) - f(D[s][0][f"{c}_W"].mean(0)) for s in seeds]) for c in COND}
        report(f"{metric} [change]", vals)
    # Q3: which pathways changed?
    print(" Q3. relative weight change by pathway. Sensory groups have no input rows (transducers),")
    print("     so for them we measure their OUTGOING weights (columns); for the rest, incoming (rows).")
    for pname, idx in PATH.items():
        ix = np.array(sorted(idx)); outgoing = pname.startswith(("olfactory", "visual", "taste"))
        vals = {}
        for c in COND:
            v = []
            for s in seeds:
                W0_, W1_ = D[s][0][f"{c}_W"].mean(0), D[s][1][f"{c}_W"].mean(0)
                A0, A1 = (W0_[:, ix], W1_[:, ix]) if outgoing else (W0_[ix], W1_[ix])
                v.append(float(np.abs(A1 - A0).sum() / (np.abs(A0).sum() + 1e-9)))
            vals[c] = np.array(v)
        report(f"{pname} [{'out' if outgoing else 'in'}] ({len(ix)})", vals)
