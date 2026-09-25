#!/usr/bin/env python3
"""Analysis for the shortcut dose-response grid, following PROTOCOL_DOSE.md.

D1  final fitness increases monotonically with shortcut share (Page's L over
    A, AS1, AS3, AS5, AS10; seed-level Spearman as a secondary read)
D2  SHAM - A lies inside the +-0.10 equivalence bound (seed-mean bootstrap)
D3  the cost of removing olfaction increases monotonically with shortcut share
D4  no difference between the AS steps at generation 0
"""
import json, glob, os, sys, itertools
from seeds import registered  # registered grids are seeds 0-9; see seeds.py
import numpy as np

CONDS = ["A", "AS1", "AS3", "AS5", "AS10"]
ALL = CONDS + ["SHAM"]
DELTA = 0.10
NBOOT = 10000
RNG = np.random.default_rng(12345)


def runs():
    out = {}
    for d in sorted(registered(glob.glob("runs/dose_*"))):
        if not os.path.exists(os.path.join(d, "done.flag")):
            continue
        name = os.path.basename(d)
        _, pred, tox, seed = name.split("_")
        out[(pred, int(seed[1:]))] = d
    return out


def final_fit(d):
    """Common-garden fitness at the last probe, per condition.

    The probe is what figures 8 and 9 plot and what the intervention analysis
    used, so the dose-response is read off the same measure rather than the
    in-run fitness in gen.jsonl.
    """
    last = None
    with open(os.path.join(d, "probe.jsonl")) as f:
        for line in f:
            last = line
    r = json.loads(last)
    return {c: r[f"{c}.normal"]["fit"] for c in ALL if f"{c}.normal" in r}, r["gen"]


def gen0_fit(d):
    with open(os.path.join(d, "gen.jsonl")) as f:
        r = json.loads(f.readline())
    return {c: r[c]["fit"] for c in ALL if c in r}


def probe_costs(d):
    """Final-generation probe: fitness drop when olfaction is ablated."""
    last = None
    with open(os.path.join(d, "probe.jsonl")) as f:
        for line in f:
            last = line
    if last is None:
        return None
    r = json.loads(last)
    keys = set(r.keys())
    ablate = None
    for cand in ("anosmic", "noodor", "nosmell", "noolf"):
        if any(k.endswith("." + cand) for k in keys):
            ablate = cand
            break
    if ablate is None:
        return None
    out = {}
    for c in ALL:
        n, a = f"{c}.normal", f"{c}.{ablate}"
        if n in r and a in r:
            out[c] = r[n]["fit"] - r[a]["fit"]
    return out


def page_L(mat):
    """Page's L for a predicted increasing trend. mat: seeds x conditions."""
    n, k = mat.shape
    ranks = np.apply_along_axis(rankdata, 1, mat)
    R = ranks.sum(axis=0)
    L = float(sum((j + 1) * R[j] for j in range(k)))
    mu = n * k * (k + 1) ** 2 / 4.0
    sd = np.sqrt(n * k ** 2 * (k + 1) * (k ** 2 - 1) / 144.0)
    z = (L - mu) / sd
    from math import erfc, sqrt
    p = 0.5 * erfc(z / sqrt(2.0))  # one-sided, increasing
    return L, z, p


def rankdata(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), float)
    r[order] = np.arange(1, len(x) + 1)
    # average ties
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    for i, c in enumerate(cnt):
        if c > 1:
            r[inv == i] = r[inv == i].mean()
    return r


def boot_ci(x, lo=5, hi=95):
    x = np.asarray(x, float)
    bs = RNG.choice(x, size=(NBOOT, len(x)), replace=True).mean(axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


def spearman(a, b):
    ra, rb = rankdata(np.asarray(a, float)), rankdata(np.asarray(b, float))
    ra, rb = ra - ra.mean(), rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def main():
    R = runs()
    ecos = sorted({p for p, _ in R})
    print(f"runs found: {len(R)}   ecologies: {ecos}")
    if not R:
        sys.exit("no completed runs")

    result = {"n_runs": len(R), "ecologies": ecos}
    SHARE = [0.012, 1.0, 3.0, 5.0, 10.0]  # nominal shortcut share, percent

    # ---- per-ecology final fitness -------------------------------------
    per_eco = {}
    for eco in ecos:
        seeds = sorted(s for p, s in R if p == eco)
        mat = np.array([[final_fit(R[(eco, s)])[0][c] for c in CONDS] for s in seeds])
        per_eco[eco] = {"seeds": seeds, "mat": mat}
        L, z, p = page_L(mat)
        rhos = [spearman(SHARE, row) for row in mat]
        print(f"\n[{eco}] n={len(seeds)} seeds, final fitness by condition")
        for j, c in enumerate(CONDS):
            m = mat[:, j].mean()
            lo, hi = boot_ci(mat[:, j])
            print(f"   {c:>5}  mean {m:6.3f}   90% CI [{lo:6.3f}, {hi:6.3f}]")
        print(f"   Page's L = {L:.0f}  z = {z:.2f}  p(one-sided, increasing) = {p:.4g}")
        print(f"   seed-level Spearman rho: mean {np.mean(rhos):+.3f}  "
              f"({sum(r > 0 for r in rhos)}/{len(rhos)} positive)")
        result[f"D1_{eco}"] = {"page_L": L, "z": z, "p": p,
                               "rho_mean": float(np.mean(rhos)),
                               "rho_pos": int(sum(r > 0 for r in rhos)),
                               "means": {c: float(mat[:, j].mean())
                                         for j, c in enumerate(CONDS)}}

    # ---- D1 pooled across ecologies (seed means) -----------------------
    seeds_common = sorted(set.intersection(*[set(per_eco[e]["seeds"]) for e in ecos]))
    pooled = np.array([[np.mean([per_eco[e]["mat"][per_eco[e]["seeds"].index(s), j]
                                 for e in ecos]) for j in range(len(CONDS))]
                       for s in seeds_common])
    L, z, p = page_L(pooled)
    rhos = [spearman(SHARE, row) for row in pooled]
    print(f"\n=== D1 (pooled seed means, n={len(seeds_common)}) ===")
    for j, c in enumerate(CONDS):
        lo, hi = boot_ci(pooled[:, j])
        print(f"   {c:>5}  mean {pooled[:, j].mean():6.3f}   90% CI [{lo:6.3f}, {hi:6.3f}]")
    print(f"   Page's L = {L:.0f}  z = {z:.2f}  p = {p:.4g}")
    print(f"   Spearman rho mean {np.mean(rhos):+.3f} ({sum(r > 0 for r in rhos)}/{len(rhos)} positive)")
    print(f"   D1 {'SUPPORTED' if p < 0.05 else 'NOT SUPPORTED'} (pre-registered p < 0.05)")
    result["D1_pooled"] = {"page_L": L, "z": z, "p": p,
                           "rho_mean": float(np.mean(rhos)),
                           "supported": bool(p < 0.05)}

    # ---- D2 SHAM vs A equivalence --------------------------------------
    print(f"\n=== D2 (SHAM - A, equivalence bound +-{DELTA}) ===")
    d2 = {}
    for eco in ecos + ["pooled"]:
        if eco == "pooled":
            diffs = np.array([np.mean([final_fit(R[(e, s)])[0]["SHAM"] -
                                       final_fit(R[(e, s)])[0]["A"] for e in ecos])
                              for s in seeds_common])
        else:
            diffs = np.array([final_fit(R[(eco, s)])[0]["SHAM"] -
                              final_fit(R[(eco, s)])[0]["A"]
                              for s in per_eco[eco]["seeds"]])
        lo, hi = boot_ci(diffs)
        inside = (lo > -DELTA) and (hi < DELTA)
        print(f"   [{eco:>5}] mean {diffs.mean():+.3f}  90% CI [{lo:+.3f}, {hi:+.3f}]  "
              f"{'INSIDE bound (equivalent)' if inside else 'NOT inside bound'}")
        d2[eco] = {"mean": float(diffs.mean()), "ci": [lo, hi], "inside": bool(inside)}
    result["D2"] = d2

    # ---- D3 olfaction dependence ---------------------------------------
    print(f"\n=== D3 (cost of ablating olfaction, by shortcut share) ===")
    costs = {}
    have = True
    for eco in ecos:
        rows = []
        for s in per_eco[eco]["seeds"]:
            c = probe_costs(R[(eco, s)])
            if c is None:
                have = False
                break
            rows.append([c[k] for k in CONDS])
        if not have:
            break
        m = np.array(rows)
        costs[eco] = m
        L, z, p = page_L(m)
        print(f"   [{eco}] " + "  ".join(f"{c}={m[:, j].mean():.3f}"
                                         for j, c in enumerate(CONDS)))
        print(f"        Page's L = {L:.0f}  z = {z:.2f}  p = {p:.4g}")
        result[f"D3_{eco}"] = {"page_L": L, "z": z, "p": p,
                               "means": {c: float(m[:, j].mean())
                                         for j, c in enumerate(CONDS)}}
    if not have:
        print("   probe file has no olfaction-ablation condition; D3 not evaluated here")
        result["D3"] = "not evaluated (no ablation probe in probe.jsonl)"

    # ---- D4 generation 0 ------------------------------------------------
    print(f"\n=== D4 (generation 0, should show no AS-step difference) ===")
    g0 = np.array([[np.mean([gen0_fit(R[(e, s)])[c] for e in ecos]) for c in CONDS]
                   for s in seeds_common])
    L, z, p = page_L(g0)
    print("   " + "  ".join(f"{c}={g0[:, j].mean():.3f}" for j, c in enumerate(CONDS)))
    print(f"   Page's L = {L:.0f}  z = {z:.2f}  p = {p:.4g}  "
          f"{'no trend (as predicted)' if p >= 0.05 else 'TREND PRESENT - D4 violated'}")
    result["D4"] = {"page_L": L, "z": z, "p": p, "supported": bool(p >= 0.05),
                    "means": {c: float(g0[:, j].mean()) for j, c in enumerate(CONDS)}}

    os.makedirs("results", exist_ok=True)
    with open("results/dose_key.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\nwrote results/dose_key.json")


if __name__ == "__main__":
    main()
