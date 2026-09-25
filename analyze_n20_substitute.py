#!/usr/bin/env python3
# HISTORICAL: this version used a substitute measure on the mis-probed extension runs (now
# runs/ecofix_pp05_*), failed its own validation, and produced results/n20_substitute_pp05.json.
# The analysis of the correctly probed extension is analyze_n20.py.
"""Corrected grid at n = 20 seeds, per ecology, on a measure that is comparable across all seeds.

Why not the registered endpoint: seeds 10-19 (the fallback extension, PROTOCOL_LEARN §3) were
probed at toxicity 0.5 because the campaign driver left --probe_ptox at the new engine's
default, while seeds 0-9 were probed at 0.0 (research log 45-2). Their last-probe fitness is
therefore not poolable with seeds 0-9.

What is used instead: the generation-500 elites of all 20 seeds re-evaluated in their own
training ecology at toxicity 0.0 (the assay's eco_P{p}_T0.0 block). This differs from the
registered endpoint in sample (8 elites vs the 16 best by training fitness) and generation
(500 vs 550), so it is reported as a robustness analysis. Its agreement with the registered
endpoint is checked on seeds 0-9, where both exist.

Estimand and interval follow the paper: connectome minus control, median over seeds, percentile
bootstrap of the median, 90 % interval for equivalence at +-0.10.
"""
import json, glob, os, re
import numpy as np

CONDS = ["A", "N1", "N2", "N4", "N5"]
NULLS = CONDS[1:]
DELTA = 0.10
RNG = np.random.default_rng(20)
NBOOT = 10000
ECO = {"0.0": "P0T0 (safe food)", "0.2": "P1T0 (predator)"}


def med_ci(d, lo=5, hi=95):
    d = np.asarray(d, float)
    b = np.median(RNG.choice(d, (NBOOT, len(d))), 1)
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))


def assay_train(run, pv):
    """fitness per condition in the run's own ecology, from the generation-500 assay"""
    f = f"runs/{run}/assay.jsonl"
    if not os.path.exists(f):
        return None
    got = None
    for line in open(f):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("assay") == f"eco_P{pv}_T0.0" and r.get("gen") == 500:
            got = {c: v["fit"] for c, v in r["res"].items()}
    return got


def probe_last(run):
    p = None
    for line in open(f"runs/{run}/probe.jsonl"):
        p = line
    p = json.loads(p)
    return {c: p[f"{c}.normal"]["fit"] for c in CONDS}


def collect():
    data = {}
    for pv in ECO:
        for s in range(20):
            run = f"ecofix_P{pv}_T0.0_s{s}"
            if not os.path.exists(f"runs/{run}/done.flag"):
                continue
            a = assay_train(run, pv)
            if a is None:
                continue
            data[(pv, s)] = {"assay": a, "probe": probe_last(run) if s < 10 else None}
    return data


def summarise(diffs):
    d = np.asarray(diffs, float)
    lo, hi = med_ci(d)
    return {"median": float(np.median(d)), "ci90": [lo, hi], "n": len(d),
            "n_positive": int((d > 0).sum()), "inside": bool(lo > -DELTA and hi < DELTA)}


def main():
    data = collect()
    out = {"measure": "gen-500 elites re-evaluated in their own ecology at toxicity 0.0 (assay eco_P*_T0.0)",
           "counts": {ECO[pv]: sum(1 for (p, _) in data if p == pv) for pv in ECO}}
    print("seeds available:", out["counts"])

    # ---- 1. is the assay measure a fair stand-in for the registered endpoint? (seeds 0-9)
    print("\n=== agreement with the registered endpoint, seeds 0-9 (connectome minus control) ===")
    agree = {}
    for n in NULLS:
        xa, xp = [], []
        for (pv, s), v in data.items():
            if s < 10 and v["probe"]:
                xa.append(v["assay"]["A"] - v["assay"][n]); xp.append(v["probe"]["A"] - v["probe"][n])
        r = float(np.corrcoef(xa, xp)[0, 1])
        same = float(np.mean(np.sign(xa) == np.sign(xp)))
        agree[n] = {"r": r, "sign_agreement": same, "n": len(xa)}
        print(f"  A-{n}: r = {r:+.2f}, same sign in {100 * same:.0f} % of {len(xa)} runs")
    out["agreement_with_registered_endpoint"] = agree

    # ---- 2. per ecology, n = 10 (seeds 0-9), n = 10 (seeds 10-19), n = 20
    print(f"\n=== per ecology: connectome minus control, median [90 % CI], equivalence at +-{DELTA} ===")
    per = {}
    for pv, lab in ECO.items():
        print(f"  {lab}")
        per[lab] = {}
        for n in NULLS:
            row = {}
            for tag, seeds in (("s0-9", range(10)), ("s10-19", range(10, 20)), ("n20", range(20))):
                d = [data[(pv, s)]["assay"]["A"] - data[(pv, s)]["assay"][n] for s in seeds if (pv, s) in data]
                if len(d) >= 5:
                    row[tag] = summarise(d)
            per[lab][n] = row
            txt = "  ".join(f"{tag} {r['median']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]"
                            f"{' eq' if r['inside'] else ''}" for tag, r in row.items())
            print(f"    A-{n}: {txt}")
    out["per_ecology"] = per

    # ---- 3. pooled over the two T0 ecologies, seed means, n = 20
    print("\n=== pooled over the two ecologies (seed means) ===")
    pooled = {}
    for n in NULLS:
        seeds = sorted({s for (_, s) in data})
        d = [np.mean([data[(pv, s)]["assay"]["A"] - data[(pv, s)]["assay"][n] for pv in ECO if (pv, s) in data])
             for s in seeds if all((pv, s) in data for pv in ECO)]
        pooled[n] = summarise(d)
        r = pooled[n]
        print(f"  A-{n}: {r['median']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]  n={r['n']}  "
              f"{r['n_positive']}/{r['n']} positive  {'inside' if r['inside'] else 'NOT inside'}")
    out["pooled_two_T0_ecologies"] = pooled

    # ---- 4. is the ecology split an artefact of averaging opposite signs?
    print("\n=== do the ecologies disagree in sign? (n = 20 per ecology) ===")
    split = {}
    for n in NULLS:
        a = per[ECO["0.0"]][n].get("n20"); b = per[ECO["0.2"]][n].get("n20")
        if a and b:
            opp = np.sign(a["median"]) != np.sign(b["median"])
            excl = (a["ci90"][0] > 0 or a["ci90"][1] < 0) and (b["ci90"][0] > 0 or b["ci90"][1] < 0)
            split[n] = {"opposite_sign": bool(opp), "both_exclude_zero": bool(excl)}
            print(f"  A-{n}: P0T0 {a['median']:+.3f}, P1T0 {b['median']:+.3f}  "
                  f"{'OPPOSITE SIGNS' if opp else 'same sign'}{' (both intervals exclude 0)' if excl else ''}")
    out["ecology_sign_split"] = split

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/n20_key.json", "w"), indent=1)
    print("\nwrote results/n20_key.json")


if __name__ == "__main__":
    main()
