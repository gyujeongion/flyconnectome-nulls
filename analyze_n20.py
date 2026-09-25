#!/usr/bin/env python3
"""Corrected grid at n = 20 seeds, per ecology, on the registered endpoint.

Seeds 10-19 of the two T0 ecologies were re-run with the corrected grid's own frozen engine
(evo_fix.py) and run_fix.sh's flags (run_followup.sh), after a first attempt was probed at the
wrong toxicity (research log 45-2; kept as runs/ecofix_pp05_*). The registered endpoint is the
common-garden fitness at the last probe (generation 550), as in analyze_fix.py.

Decision rules, written before the extension data existed (2026-09-23, research log 46):

  R1  Per-ecology equivalence at n = 20: the 90 % bootstrap interval of the median of
      (connectome - control) lies inside +-0.10.
  R2  The abstract's claim that the interior column shuffle N5 is ahead of the connectome in
      the safe-food ecology (P0T0: -0.264 at n = 10, p_holm = 0.029) counts as REPLICATED if in
      seeds 10-19 alone the median of A - N5 in P0T0 is negative and its 90 % interval excludes
      zero; NOT REPLICATED if the interval includes zero; CONTRADICTED if it is positive with the
      interval excluding zero.
  R3  The two-ecology pool (seed means over P0T0 and P1T0) is a new estimand, not the paper's
      four-ecology Table 2, and is reported as such.

Before any number is computed, every extension run's configuration must equal seed 0's apart
from run name and seed. If one differs the script stops: that check is the one that would have
caught the probe mismatch the first time.
"""
import json, glob, os, sys
import numpy as np
from scipy.stats import wilcoxon
from seeds import REGISTERED, EXTENSION

CONDS = ["A", "N1", "N2", "N4", "N5"]
NULLS = CONDS[1:]
DELTA = 0.10
RNG = np.random.default_rng(20260923)
NBOOT = 10000
ECO = {"0.0": "P0T0", "0.2": "P1T0"}


def holm(ps):
    ps = np.asarray(ps, float); o = np.argsort(ps); m = len(ps); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * ps[i]); adj[i] = min(1.0, run)
    return adj


def med_ci(d, lo=5, hi=95):
    d = np.asarray(d, float)
    b = np.median(RNG.choice(d, (NBOOT, len(d))), 1)
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))


def final(run):
    p = None
    for line in open(f"runs/{run}/probe.jsonl"):
        p = line
    p = json.loads(p)
    assert p["gen"] == 550, f"{run}: last probe at generation {p['gen']}, expected 550"
    return {c: p[f"{c}.normal"]["fit"] for c in CONDS}


def check_configs():
    ref = json.load(open("runs/ecofix_P0.0_T0.0_s0/config.json"))
    ignore = {"run", "seed"}
    bad = []
    for pv in ECO:
        for s in list(REGISTERED) + list(EXTENSION):
            run = f"ecofix_P{pv}_T0.0_s{s}"
            f = f"runs/{run}/config.json"
            if not os.path.exists(f):
                bad.append((run, "missing"))
                continue
            c = json.load(open(f))
            diff = {k: (ref["args"].get(k), c["args"].get(k)) for k in set(ref["args"]) | set(c["args"])
                    if k not in ignore and k != "pred_fixed" and ref["args"].get(k) != c["args"].get(k)}
            if c["args"].get("pred_fixed") != float(pv):
                diff["pred_fixed"] = (float(pv), c["args"].get("pred_fixed"))
            if c["E"] != ref["E"]:
                diff["E"] = "differs"
            if c.get("islands") != ref.get("islands"):
                diff["islands"] = (ref.get("islands"), c.get("islands"))
            if diff:
                bad.append((run, diff))
    return bad


def summ(d):
    d = np.asarray(d, float); lo, hi = med_ci(d)
    return {"median": float(np.median(d)), "ci90": [lo, hi], "n": len(d), "n_positive": int((d > 0).sum()),
            "inside": bool(lo > -DELTA and hi < DELTA)}


def main():
    bad = check_configs()
    if bad:
        print("CONFIG CHECK FAILED: extension runs differ from seed 0 or are missing")
        for b in bad[:10]:
            print("  ", b)
        sys.exit(1)
    print("config check: all 40 runs (2 ecologies x seeds 0-19) match seed 0 apart from name and seed")

    V = {(pv, s): final(f"ecofix_P{pv}_T0.0_s{s}") for pv in ECO for s in list(REGISTERED) + list(EXTENSION)}
    out = {"endpoint": "common-garden fitness at the last probe (generation 550)", "decision_rules": __doc__}

    print("\n=== per ecology: connectome minus control, median [90 % CI] ===")
    per = {}
    for pv, lab in ECO.items():
        print(f"  {lab}")
        per[lab] = {}
        d20_all = {n: np.array([V[(pv, s)]["A"] - V[(pv, s)][n] for s in range(20)]) for n in NULLS}
        p20 = holm([wilcoxon(d20_all[n]).pvalue if np.any(d20_all[n] != 0) else 1.0 for n in NULLS])
        for k, n in enumerate(NULLS):
            row = {}
            for tag, seeds in (("s0-9", REGISTERED), ("s10-19", EXTENSION), ("n20", range(20))):
                row[tag] = summ([V[(pv, s)]["A"] - V[(pv, s)][n] for s in seeds])
            row["n20"]["p_holm"] = float(p20[k])
            per[lab][n] = row
            print(f"    A-{n}: " + "  ".join(
                f"{t} {r['median']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]{' eq' if r['inside'] else ''}"
                for t, r in row.items()) + f"  p_holm(n20)={row['n20']['p_holm']:.3f}")
    out["per_ecology"] = per

    r2 = per["P0T0"]["N5"]["s10-19"]
    lo, hi = r2["ci90"]
    if r2["median"] < 0 and hi < 0:
        verdict = "REPLICATED"
    elif r2["median"] > 0 and lo > 0:
        verdict = "CONTRADICTED"
    else:
        verdict = "NOT REPLICATED"
    out["R2_N5_ahead_in_P0T0"] = {"seeds_10_19": r2, "verdict": verdict}
    print(f"\n=== R2: N5 ahead of the connectome in P0T0, replication in seeds 10-19 ===")
    print(f"  A-N5 {r2['median']:+.3f} [{lo:+.3f},{hi:+.3f}]  -> {verdict}")

    print("\n=== R3: two-ecology pool (seed means over P0T0 and P1T0), n = 20; not the paper's Table 2 estimand ===")
    pooled = {}
    for n in NULLS:
        d = [np.mean([V[(pv, s)]["A"] - V[(pv, s)][n] for pv in ECO]) for s in range(20)]
        pooled[n] = summ(d)
        r = pooled[n]
        print(f"  A-{n}: {r['median']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]  "
              f"{r['n_positive']}/20 positive  {'inside' if r['inside'] else 'NOT inside'} +-{DELTA}")
    out["pooled_two_T0_ecologies"] = pooled

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/n20_key.json", "w"), indent=1)
    print("\nwrote results/n20_key.json")


if __name__ == "__main__":
    main()
