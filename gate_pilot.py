#!/usr/bin/env python3
"""Pilot gate of PROTOCOL_LEARN.md §3: did associative learning actually evolve?

The measure is the reversal-specific benefit of plasticity, scaled by the condition's own
fitness so it does not depend on the ecology's fitness scale:

    I(c) = [ (reversal - reversal_noplast) - (base - base_noplast) ] / base

Plasticity helping in general is not enough. In the established ecology plasticity is worth a
median of 0.15 of baseline fitness (n = 220 condition x seed x ecology observations, see
results/gate_calibration.json), so a gate on that quantity would pass the very ecology whose
lack of associative learning motivated this experiment. I(c), by contrast, sits at zero there:
median -0.019 and +0.005 for the two grids, 90th percentile 0.15 to 0.18, maximum 0.43.

    G1  max_c I(c) >= 0.30      about twice the control's 90th percentile, below its maximum
    G2  I(c) > 0 in at least 2 of the 3 pilot conditions

Both thresholds were fixed from the control ecology before the candidate ecology existed.
Exit 0 = passed, run the learning grid. Exit 1 = failed, run the fallback grid.

    python3 gate_pilot.py runs/pilot_learn_s0
"""
import json, sys, os

G1_MIN = 0.30      # max relative reversal-specific plasticity benefit
G2_MIN_POS = 2     # conditions that must be positive


def blocks(path):
    """Last record per (assay block, condition) from assay.jsonl.

    The engine writes one line per block as
        {"run":..., "gen":..., "assay": "<tag>", "res": {"<cond>": {"fit":..., ...}, ...}}
    and appends, so later lines for the same tag win.
    """
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
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
            for cond, v in res.items():
                if isinstance(v, dict) and "fit" in v:
                    out[(tag, cond)] = float(v["fit"])
    return out


def main(run):
    path = os.path.join(run, "assay.jsonl")
    if not os.path.exists(path):
        print(f"GATE FAIL: no {path} (the pilot assay did not produce output)")
        return 1
    b = blocks(path)
    need = ("base", "base_noplast", "reversal", "reversal_noplast")
    have = {t for t, _ in b}
    missing = [t for t in need if t not in have]
    if missing:
        print(f"GATE FAIL: assay is missing blocks {missing}; found {sorted(have)}")
        return 1

    conds = sorted({c for (tag, c) in b if tag == "base"})
    print(f"{'condition':>10}  {'base':>8}  {'plast':>8}  {'plast|rev':>10}  {'I':>8}")
    I = {}
    for c in conds:
        try:
            base = b[("base", c)]
            if base <= 0.01:
                print(f"{c:>10}  base {base:.3f} too small to scale, skipped")
                continue
            pb = base - b[("base_noplast", c)]
            pr = b[("reversal", c)] - b[("reversal_noplast", c)]
            I[c] = (pr - pb) / base
            print(f"{c:>10}  {base:8.3f}  {pb:8.3f}  {pr:10.3f}  {I[c]:8.3f}")
        except KeyError:
            continue

    if not I:
        print("GATE FAIL: no condition had a complete set of blocks")
        return 1

    best_c = max(I, key=I.get)
    g1 = I[best_c]
    n_pos = sum(1 for v in I.values() if v > 0)
    ok1, ok2 = g1 >= G1_MIN, n_pos >= G2_MIN_POS

    print(f"\nG1  max I = {g1:+.3f} ({best_c})   needs >= {G1_MIN}   {'pass' if ok1 else 'FAIL'}")
    print(f"G2  positive in {n_pos}/{len(I)} conditions   needs >= {G2_MIN_POS}   {'pass' if ok2 else 'FAIL'}")
    print("reference, established ecology: median ~0.00, 90th pct 0.15-0.18, max 0.43 (n=220)")
    passed = ok1 and ok2
    print("\nGATE PASS: associative learning evolved, run the learning grid" if passed else
          "\nGATE FAIL: no associative learning beyond what the control ecology shows; run the fallback")

    json.dump({"I": I, "max": g1, "max_cond": best_c, "n_positive": n_pos,
               "G1_MIN": G1_MIN, "G2_MIN_POS": G2_MIN_POS, "passed": bool(passed)},
              open(os.path.join(run, "gate.json"), "w"), indent=1)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "runs/pilot_learn_s0"))
