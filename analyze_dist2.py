#!/usr/bin/env python3
"""Corrected grid under distribution-matched calibration (PROTOCOL_DIST2.md).

The registered corrected grid (ecofix_*, seeds 0-9) used global calibration, which matches mean
group activity but not its spread: the connectome runs with a wide activity distribution (SD
0.163, ~48 silent groups) against narrow ones for every null (SD 0.080-0.086). This grid re-runs
all four ecologies with --homeo distmatch and nothing else changed, to see whether the paper's two
comparisons survive when the distributions are matched.

Decision rules, frozen before the data existed (2026-09-23, research log 49):

  Manipulation check  per run, for N in {N4, N5}: |SD_A - SD_N| <= 0.01, |silent_A - silent_N| <= 10,
                      |sat_A - sat_N| <= 0.02. Passes if at least 36 of 40 runs pass; otherwise every
                      verdict below is reported as exploratory.
  Seed mean           per seed, the mean over the four ecologies of (A - control); n = 10. This is the
                      estimand of the paper's abstract and Table 2; the same estimand on the registered
                      global runs (seeds 0-9) is reported beside every distmatch value.
  P1 (co-primary)     distmatch seed mean A-N4 and A-N5: EQUIVALENT if the 90 % bootstrap interval of the
                      median lies inside +-0.10, otherwise NOT EQUIVALENT.
  P2 (co-primary)     distmatch seed mean A-N1 and A-N2: BEHIND if median < 0 and upper bound < 0;
                      NOT REPLICATED if the interval includes 0; REVERSED if median > 0 and lower bound > 0.
  S1 (secondary)      for N in {N1, N2, N4, N5}, per seed D = (seed mean A-N)distmatch - (seed mean A-N)global:
                      SMALL if the interval lies inside +-0.10; CHANGED if it excludes 0 and |median| > 0.10;
                      UNDETERMINED otherwise.
  S2 (secondary)      P0T0 A-N5 under distmatch: REPLICATED if median < 0 and upper < 0; NOT REPLICATED if
                      the interval includes 0; CONTRADICTED if median > 0 and lower > 0.

Before any number is computed, all 40 runs must exist and each configuration must equal the ecofix
run of the same ecology and seed apart from run name and calibration; otherwise exit 1.
"""
import ast, json, os, re, sys
import numpy as np
from scipy.stats import wilcoxon
from seeds import REGISTERED

CONDS = ["A", "N1", "N2", "N4", "N5"]
DELTA = 0.10
NBOOT = 10000
RNG = np.random.default_rng(20260924)
ECO = {("0.0", "0.0"): "P0T0", ("0.2", "0.0"): "P1T0", ("0.0", "1.0"): "P0T1", ("0.2", "1.0"): "P1T1"}
HOMEO = re.compile(r"homeostasis (\w+): (\{.*?\})")


def med_ci(d, lo=5, hi=95):
    d = np.asarray(d, float)
    b = np.median(RNG.choice(d, (NBOOT, len(d))), 1)
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))


def summ(d):
    d = np.asarray(d, float); lo, hi = med_ci(d)
    p = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
    return {"median": float(np.median(d)), "ci90": [lo, hi], "n": len(d), "n_positive": int((d > 0).sum()),
            "wilcoxon_p": p}


def fmt(r):
    return f"{r['median']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]"


def name(kind, eco, s):
    return f"eco{kind}_P{eco[0]}_T{eco[1]}_s{s}"


def final(run):
    p = None
    for line in open(f"runs/{run}/probe.jsonl"):
        p = line
    p = json.loads(p)
    assert p["gen"] == 550, f"{run}: last probe at generation {p['gen']}, expected 550"
    return {c: p[f"{c}.normal"]["fit"] for c in CONDS}


def check_configs():
    bad = []
    for eco in ECO:
        for s in REGISTERED:
            d, g = name("dist", eco, s), name("fix", eco, s)
            miss = [r for r in (d, g) if not (os.path.exists(f"runs/{r}/config.json") and os.path.exists(f"runs/{r}/done.flag"))]
            if miss:
                bad += [(r, "missing or unfinished") for r in miss]
                continue
            cd, cg = (json.load(open(f"runs/{r}/config.json")) for r in (d, g))
            diff = {k: (cg["args"].get(k), cd["args"].get(k)) for k in set(cd["args"]) | set(cg["args"])
                    if k not in {"run", "homeo"} and cd["args"].get(k) != cg["args"].get(k)}
            if cd["args"].get("homeo") != "distmatch" or cg["args"].get("homeo") != "global":
                diff["homeo"] = (cg["args"].get("homeo"), cd["args"].get("homeo"))
            if cd["E"] != cg["E"]:
                diff["E"] = "differs"
            if diff:
                bad.append((d, diff))
    return bad


def calib(run):
    out = {}
    for m in HOMEO.finditer(open(f"runs_{run}.log").read()):
        c, d = m.group(1), ast.literal_eval(m.group(2))
        out[c] = {"sd": d.get("sd", d.get("group_act_sd")), "silent": d["silent_groups"],
                  "sat": d["sat_frac"], "kc_mean": d["kc_mean"]}
    return out


def main():
    bad = check_configs()
    if bad:
        print("CONFIG CHECK FAILED: distmatch runs missing, unfinished or not matching ecofix")
        for b in bad[:10]:
            print("  ", b)
        sys.exit(1)
    n_runs = len(ECO) * len(REGISTERED)
    print(f"config check: {n_runs} distmatch runs match their ecofix counterparts apart from run name and calibration")
    out = {"endpoint": "common-garden fitness at the last probe (generation 550)", "decision_rules": __doc__}

    # ---- manipulation check ------------------------------------------------
    mc, cal = [], {}
    for eco in ECO:
        for s in REGISTERED:
            r = name("dist", eco, s); c = calib(r); cal[r] = c
            mc.append(all(abs(c["A"]["sd"] - c[n]["sd"]) <= 0.01 and abs(c["A"]["silent"] - c[n]["silent"]) <= 10
                          and abs(c["A"]["sat"] - c[n]["sat"]) <= 0.02 for n in ("N4", "N5")))
    mc_pass = sum(mc) >= 36
    out["manipulation_check"] = {"runs_passing": int(sum(mc)), "of": len(mc), "passed": mc_pass, "calibration": cal}
    print(f"\nmanipulation check: {sum(mc)}/{len(mc)} runs matched -> {'PASSED' if mc_pass else 'FAILED (verdicts exploratory)'}")
    kc_n1 = [cal[name("dist", e, s)]["N1"]["kc_mean"] for e in ECO for s in REGISTERED]
    print(f"  N1 kc_mean under distmatch: {min(kc_n1):.3f}-{max(kc_n1):.3f} (target 0.05; known limitation)")
    tag = "" if mc_pass else " (exploratory)"

    V = {(k, e, s): final(name(k, e, s)) for k in ("dist", "fix") for e in ECO for s in REGISTERED}
    diff = lambda k, e, s, n: V[(k, e, s)]["A"] - V[(k, e, s)][n]
    seedmean = lambda k, n: [np.mean([diff(k, e, s, n) for e in ECO]) for s in REGISTERED]

    print("\n=== per ecology, A minus control, median [90 % CI] (descriptive) ===")
    per = {}
    for e, lab in ECO.items():
        per[lab] = {n: {k: summ([diff(k, e, s, n) for s in REGISTERED]) for k in ("fix", "dist")} for n in CONDS[1:]}
        for n in CONDS[1:]:
            print(f"  {lab} A-{n}: global {fmt(per[lab][n]['fix'])}   distmatch {fmt(per[lab][n]['dist'])}")
    out["per_ecology"] = per

    SM = {n: {"global": summ(seedmean("fix", n)), "distmatch": summ(seedmean("dist", n))} for n in CONDS[1:]}

    print(f"\n=== P1: seed mean over four ecologies, interior-preserving controls, equivalence +-{DELTA}{tag} ===")
    for n in ("N4", "N5"):
        d = SM[n]["distmatch"]
        SM[n]["P1"] = "EQUIVALENT" if (d["ci90"][0] > -DELTA and d["ci90"][1] < DELTA) else "NOT EQUIVALENT"
        print(f"  A-{n}: global {fmt(SM[n]['global'])}   distmatch {fmt(d)}  -> {SM[n]['P1']}")

    print(f"\n=== P2: seed mean over four ecologies, shortcut-bearing controls{tag} ===")
    for n in ("N1", "N2"):
        d = SM[n]["distmatch"]; lo, hi = d["ci90"]
        SM[n]["P2"] = "BEHIND" if (d["median"] < 0 and hi < 0) else "REVERSED" if (d["median"] > 0 and lo > 0) else "NOT REPLICATED"
        print(f"  A-{n}: global {fmt(SM[n]['global'])}   distmatch {fmt(d)}  -> {SM[n]['P2']}")

    print(f"\n=== S1: change caused by calibration, seed mean (A-N)distmatch - (A-N)global{tag} ===")
    for n in CONDS[1:]:
        d = summ(np.array(seedmean("dist", n)) - np.array(seedmean("fix", n))); lo, hi = d["ci90"]
        v = ("SMALL" if (lo > -DELTA and hi < DELTA) else
             "CHANGED" if ((lo > 0 or hi < 0) and abs(d["median"]) > DELTA) else "UNDETERMINED")
        SM[n]["S1"] = {**d, "verdict": v}
        print(f"  {n}: {fmt(d)}  -> {v}")
    out["seed_mean"] = SM

    d = per["P0T0"]["N5"]["dist"]; lo, hi = d["ci90"]
    v = "REPLICATED" if (d["median"] < 0 and hi < 0) else "CONTRADICTED" if (d["median"] > 0 and lo > 0) else "NOT REPLICATED"
    out["S2_N5_ahead_P0T0"] = {"distmatch": d, "verdict": v}
    print(f"\n=== S2: N5 ahead in P0T0 under distmatch{tag} ===\n  A-N5 {fmt(d)}  -> {v}")

    both = all(SM[n]["P1"] == "EQUIVALENT" for n in ("N4", "N5")) and all(SM[n]["P2"] == "BEHIND" for n in ("N1", "N2"))
    out["interpretation"] = ("both claims hold under distribution-matched calibration" if both else
                             "at least one claim must be qualified as holding under global calibration")
    print(f"\ninterpretation{tag}: {out['interpretation']}")

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/dist2_key.json", "w"), indent=1)
    print("wrote results/dist2_key.json")


if __name__ == "__main__":
    main()
