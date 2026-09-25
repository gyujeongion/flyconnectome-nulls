#!/usr/bin/env python3
"""Out-of-distribution transfer analysis, following PROTOCOL_TRANSFER.md.

Evolved elites are re-evaluated, without further evolution, in worlds they never saw. The
question is whether the populations that lean hardest on olfaction fall apart furthest when
the world changes.

    loss(c, b) = base(c) - b(c)          fitness units, larger = degraded more

The primary read is the mean loss over three blocks that move different axes of the world and
sit outside the training distribution: stress_pred0.15, transfer_toxic0.9, transfer_plume0.12.

T1  shortcut-carrying conditions lose more than the connectome (AS10-A, N1-A, N2-A)
T2  loss increases monotonically with shortcut dose (A, AS1, AS3, AS5, AS10)
T3  boundary-preserving controls (N4, N5) are inside +-0.10 of the connectome
T4  the base block reproduces the differences already reported, so T1 is not a restatement
    of a baseline gap
"""
import json, glob, os, sys
from seeds import registered  # registered grids are seeds 0-9; see seeds.py
import numpy as np

PRIMARY = ("stress_pred0.15", "transfer_toxic0.9", "transfer_plume0.12")
SECONDARY = ("stress_pred0.12", "stress_pred0.18", "stress_toxic0.7", "stress_toxic1.0",
             "reversal", "eco_P0.0_T0.0", "eco_P0.2_T0.0", "eco_P0.0_T1.0", "eco_P0.2_T1.0")
DELTA = 0.10
NBOOT = 10000
RNG = np.random.default_rng(20260923)


def read(path):
    out = {}
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        tag, res = r.get("assay"), r.get("res")
        if tag is None or not isinstance(res, dict):
            continue
        for c, v in res.items():
            if isinstance(v, dict) and "fit" in v:
                out[(tag, c)] = float(v["fit"])
    return out


def losses(prefix, conds, blocks):
    """seed -> cond -> mean loss over `blocks`, averaged across the two ecologies."""
    per = {}
    for f in sorted(registered(glob.glob(f"runs/{prefix}_P*_T0.0_s*/assay.jsonl"))):
        run = os.path.basename(os.path.dirname(f))
        seed = int(run.split("_s")[-1])
        b = read(f)
        for c in conds:
            if ("base", c) not in b:
                continue
            vals = [b[("base", c)] - b[(bl, c)] for bl in blocks if (bl, c) in b]
            if len(vals) != len(blocks):
                continue
            per.setdefault(seed, {}).setdefault(c, []).append(float(np.mean(vals)))
    return {s: {c: float(np.mean(v)) for c, v in d.items()} for s, d in per.items()}


def boot_ci(x, lo=5, hi=95):
    x = np.asarray(x, float)
    bs = RNG.choice(x, size=(NBOOT, len(x)), replace=True).mean(axis=1)
    return float(np.percentile(bs, lo)), float(np.percentile(bs, hi))


def wilcoxon_p(d):
    """Two-sided exact-ish signed-rank p by permuting signs; n is 10 so enumeration is cheap."""
    d = np.asarray([v for v in d if v != 0], float)
    n = len(d)
    if n == 0:
        return 1.0
    r = rankdata(np.abs(d))
    obs = float(np.sum(r[d > 0]))
    signs = np.array(np.meshgrid(*[[0, 1]] * n)).T.reshape(-1, n) if n <= 20 else None
    if signs is None:
        return float("nan")
    dist = signs @ r
    total = r.sum()
    lo = min(obs, total - obs)
    return float(2.0 * np.mean(dist <= lo))


def rankdata(x):
    order = np.argsort(x, kind="mergesort")
    out = np.empty(len(x), float)
    out[order] = np.arange(1, len(x) + 1)
    vals, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    for i, c in enumerate(cnt):
        if c > 1:
            out[inv == i] = out[inv == i].mean()
    return out


def holm(pairs):
    """pairs: list of (label, p). Returns dict label -> adjusted p."""
    s = sorted(pairs, key=lambda kv: kv[1])
    m, out, prev = len(s), {}, 0.0
    for i, (lab, p) in enumerate(s):
        adj = min(1.0, max(prev, (m - i) * p))
        out[lab] = adj
        prev = adj
    return out


def page_L(mat):
    n, k = mat.shape
    R = np.apply_along_axis(rankdata, 1, mat).sum(axis=0)
    L = float(sum((j + 1) * R[j] for j in range(k)))
    mu = n * k * (k + 1) ** 2 / 4.0
    sd = np.sqrt(n * k ** 2 * (k + 1) * (k ** 2 - 1) / 144.0)
    z = (L - mu) / sd
    from math import erfc, sqrt
    return L, z, 0.5 * erfc(z / sqrt(2.0))


def main():
    out = {"primary_blocks": list(PRIMARY)}
    DOSE = ["A", "AS1", "AS3", "AS5", "AS10", "SHAM"]
    FIX = ["A", "N1", "N2", "N4", "N5"]

    Ld = losses("dose", DOSE, PRIMARY)
    Lf = losses("ecofix", FIX, PRIMARY)
    seeds_d, seeds_f = sorted(Ld), sorted(Lf)
    print(f"dose seeds {len(seeds_d)}, ecofix seeds {len(seeds_f)}, primary blocks {PRIMARY}\n")

    print("mean transfer loss (seed means over both ecologies)")
    for nm, L, conds in (("dose", Ld, DOSE), ("ecofix", Lf, FIX)):
        row = {c: np.mean([L[s][c] for s in sorted(L) if c in L[s]]) for c in conds}
        print("  " + nm + ": " + "  ".join(f"{c}={row[c]:.3f}" for c in conds))
        out[f"loss_{nm}"] = {c: float(v) for c, v in row.items()}

    # ---- T1 -------------------------------------------------------------
    print("\n=== T1  shortcut-carrying conditions lose more than the connectome ===")
    tests, raw = [], {}
    for nm, L, pairs in (("dose", Ld, [("AS10", "A")]), ("ecofix", Lf, [("N1", "A"), ("N2", "A")])):
        for hi, lo in pairs:
            ss = [s for s in sorted(L) if hi in L[s] and lo in L[s]]
            d = np.array([L[s][hi] - L[s][lo] for s in ss])
            p = wilcoxon_p(d)
            lo_ci, hi_ci = boot_ci(d)
            tests.append((f"{hi}-{lo}", p))
            raw[f"{hi}-{lo}"] = (d, p, lo_ci, hi_ci, len(ss))
    adj = holm(tests)
    n_sup = 0
    for lab, (d, p, lo_ci, hi_ci, n) in raw.items():
        sup = (adj[lab] < 0.05) and (d.mean() > 0)
        n_sup += sup
        print(f"  {lab:>8}  n={n}  mean {d.mean():+.3f}  90% CI [{lo_ci:+.3f}, {hi_ci:+.3f}]  "
              f"p={p:.4f}  p_holm={adj[lab]:.4f}  {'supports' if sup else 'no'}")
        out[f"T1_{lab}"] = {"mean": float(d.mean()), "ci": [lo_ci, hi_ci],
                            "p": p, "p_holm": adj[lab], "n": n, "supports": bool(sup)}
    out["T1_supported"] = bool(n_sup >= 2)
    print(f"  T1 {'SUPPORTED' if n_sup >= 2 else 'NOT SUPPORTED'} (needs 2 of 3, got {n_sup})")

    # ---- T2 -------------------------------------------------------------
    print("\n=== T2  loss increases with shortcut dose ===")
    dc = ["A", "AS1", "AS3", "AS5", "AS10"]
    ss = [s for s in seeds_d if all(c in Ld[s] for c in dc)]
    mat = np.array([[Ld[s][c] for c in dc] for s in ss])
    L_, z, p = page_L(mat)
    print("  " + "  ".join(f"{c}={mat[:, j].mean():.3f}" for j, c in enumerate(dc)))
    print(f"  Page's L={L_:.0f}  z={z:.2f}  p={p:.4g}  {'SUPPORTED' if p < 0.05 else 'NOT SUPPORTED'}")
    out["T2"] = {"page_L": L_, "z": z, "p": p, "supported": bool(p < 0.05), "n": len(ss)}

    # ---- T3 -------------------------------------------------------------
    print(f"\n=== T3  boundary-preserving controls within +-{DELTA} of the connectome ===")
    out["T3"] = {}
    for c in ("N4", "N5"):
        ss = [s for s in seeds_f if c in Lf[s] and "A" in Lf[s]]
        d = np.array([Lf[s][c] - Lf[s]["A"] for s in ss])
        lo_ci, hi_ci = boot_ci(d)
        inside = (lo_ci > -DELTA) and (hi_ci < DELTA)
        print(f"  {c}-A  mean {d.mean():+.3f}  90% CI [{lo_ci:+.3f}, {hi_ci:+.3f}]  "
              f"{'inside' if inside else 'NOT inside'}")
        out["T3"][c] = {"mean": float(d.mean()), "ci": [lo_ci, hi_ci], "inside": bool(inside)}

    # ---- T4 -------------------------------------------------------------
    print("\n=== T4  baseline differences (transfer loss must not just restate these) ===")
    out["T4"] = {}
    for nm, prefix, conds in (("dose", "dose", DOSE), ("ecofix", "ecofix", FIX)):
        base = {}
        for f in sorted(registered(glob.glob(f"runs/{prefix}_P*_T0.0_s*/assay.jsonl"))):
            seed = int(os.path.basename(os.path.dirname(f)).split("_s")[-1])
            b = read(f)
            for c in conds:
                if ("base", c) in b:
                    base.setdefault(seed, {}).setdefault(c, []).append(b[("base", c)])
        base = {s: {c: float(np.mean(v)) for c, v in d.items()} for s, d in base.items()}
        row = {c: np.mean([base[s][c] for s in sorted(base) if c in base[s]]) for c in conds}
        print(f"  {nm} base fitness: " + "  ".join(f"{c}={row[c]:.3f}" for c in conds))
        out["T4"][nm] = {c: float(v) for c, v in row.items()}

    # ---- post hoc: the registered baseline was already out of distribution ---------------
    # The `base` block evaluates at toxicity 0.5, but every one of these populations evolved at
    # toxicity 0.0, so `base` sits 32-33 % below the fitness they reach in their own ecology and
    # the registered loss measures "already shifted -> shifted further", not "in -> out of
    # distribution". The cross-ecology blocks give clean single-axis shifts from the true
    # training baseline instead:
    #   toxicity:  eco_P{p}_T0.0  ->  eco_P{p}_T1.0      (same predator, toxicity switched on)
    #   predator:  eco_P{p}_T0.0  ->  eco_P{q}_T0.0      (same toxicity, the other predator)
    # This is not the registered analysis and is labelled post hoc wherever it is used.
    print("\n=== POST HOC: loss from the training ecology on single-axis shifts ===")
    other = {"0.0": "0.2", "0.2": "0.0"}
    offset = {}
    ph = {}
    for nm, prefix, conds, pairs in (
            ("dose", "dose", DOSE, [("AS10", "A")]),
            ("ecofix", "ecofix", FIX, [("N1", "A"), ("N2", "A"), ("N4", "A"), ("N5", "A")])):
        base_gap = []
        for axis in ("toxicity", "predator"):
            per = {}
            for f in sorted(registered(glob.glob(f"runs/{prefix}_P*_T0.0_s*/assay.jsonl"))):
                run = os.path.basename(os.path.dirname(f))
                pv = run.split("_P")[1].split("_")[0]
                s = int(run.split("_s")[-1])
                b = read(f)
                train = f"eco_P{pv}_T0.0"
                shift = f"eco_P{pv}_T1.0" if axis == "toxicity" else f"eco_P{other[pv]}_T0.0"
                for c in conds:
                    if (train, c) in b and (shift, c) in b and b[(train, c)] > 0:
                        per.setdefault(s, {}).setdefault(c, []).append(
                            (b[(train, c)] - b[(shift, c)]) / b[(train, c)])
                    if axis == "toxicity" and ("base", c) in b and (train, c) in b and b[(train, c)] > 0:
                        base_gap.append(1 - b[("base", c)] / b[(train, c)])
            per = {s: {c: float(np.mean(v)) for c, v in d.items()} for s, d in per.items()}
            means = {c: float(np.mean([per[s][c] for s in per if c in per[s]])) for c in conds}
            print(f"  [{nm} / {axis}] relative loss: " + "  ".join(f"{c}={means[c]:.3f}" for c in conds))
            res = {"relative_loss": means, "contrasts": {}}
            for hi, lo in pairs:
                ss = [s for s in sorted(per) if hi in per[s] and lo in per[s]]
                d = np.array([per[s][hi] - per[s][lo] for s in ss])
                lo_ci, hi_ci = boot_ci(d)
                npos = int((d > 0).sum())
                print(f"     {hi}-{lo}: {d.mean():+.3f}  90% CI [{lo_ci:+.3f}, {hi_ci:+.3f}]  {npos}/{len(d)} positive")
                res["contrasts"][f"{hi}-{lo}"] = {"mean": float(d.mean()), "ci": [lo_ci, hi_ci],
                                                   "n_positive": npos, "n": len(d)}
            ph[f"{nm}_{axis}"] = res
        offset[nm] = float(np.mean(base_gap))
        print(f"  [{nm}] registered `base` block sits {100 * offset[nm]:.0f} % below the training ecology")
    out["post_hoc_single_axis"] = ph
    out["post_hoc_base_block_offset"] = offset
    out["post_hoc_note"] = ("registered analysis used the `base` block (toxicity 0.5) as baseline, which is "
                            "already out of distribution for populations evolved at toxicity 0.0; "
                            "post_hoc_single_axis measures loss from the training ecology instead")

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/transfer_key.json", "w"), indent=1)
    print("\nwrote results/transfer_key.json")


if __name__ == "__main__":
    main()
