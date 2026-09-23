#!/usr/bin/env python3
"""Pick the learning-rewarding ecology from an ecosearch log, by the rule frozen in PROTOCOL_LEARN.md §2.

The rule is deliberately mechanical: constraints, then the single best learn_margin, ties broken by
file order. Nothing here looks at wiring conditions, so the choice cannot favour the connectome.

    python3 pick_eco.py runs/ecosearch_learn/ecosearch.jsonl results/learn_eco.json
"""
import json, sys, pathlib

MIN_MARGIN = 0.3      # C5: learning must beat the best reflex by at least 3x the equivalence bound
LEARN = ("twoMem", "localSearch", "trapline")
NONLEARN = ("noEat", "reflex", "sickStop")


def main(src, dst):
    rows = []
    for line in open(src, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    print(f"candidates read: {len(rows)}")

    kept, rejected = [], {"flip0": 0, "learner_surv": 0, "reflex_surv": 0, "margin": 0}
    for i, r in enumerate(rows):
        s, cfg = r["score"], r["cfg"]
        if float(cfg.get("flip_rate", 0.0)) <= 0.0:          # C4
            rejected["flip0"] += 1; continue
        if float(s.get("learner_surv", 0.0)) < 0.5:          # C2
            rejected["learner_surv"] += 1; continue
        if float(s.get("reflex_surv", 0.0)) < 0.3:           # C3
            rejected["reflex_surv"] += 1; continue
        if float(s.get("learn_margin", -9.9)) < MIN_MARGIN:   # C5
            rejected["margin"] += 1; continue
        kept.append((i, r))

    print(f"rejected: flip_rate=0 {rejected['flip0']}, learner_surv<0.5 {rejected['learner_surv']}, "
          f"reflex_surv<0.3 {rejected['reflex_surv']}, learn_margin<{MIN_MARGIN} {rejected['margin']}")
    print(f"passing constraints: {len(kept)}")
    if not kept:
        best_seen = max((float(r["score"].get("learn_margin", -9.9)) for r in rows), default=None)
        print(f"best learn_margin anywhere in the search: {best_seen}")
        sys.exit("SEARCH FAILED: no ecology rewards learning by the frozen margin; "
                 "per PROTOCOL_LEARN the learning grid is not run")

    # single rule: highest learn_margin, ties to the earlier row
    best_i, best = min(kept, key=lambda kv: (-float(kv[1]["score"]["learn_margin"]), kv[0]))
    cfg, s = best["cfg"], best["score"]

    print("\nselected (rule: max learn_margin among constraint-passing candidates)")
    print(f"  row {best_i}  learn_margin {s['learn_margin']}  flip_rate {cfg['flip_rate']}")
    print(f"  learner_surv {s['learner_surv']}  reflex_surv {s['reflex_surv']}  nonlearn_best {s['nonlearn_best']}")
    print(f"  cfg {json.dumps(cfg, sort_keys=True)}")

    runner_up = sorted((float(r["score"]["learn_margin"]) for _, r in kept), reverse=True)[1:4]
    print(f"  next margins: {runner_up}")

    pathlib.Path(dst).parent.mkdir(parents=True, exist_ok=True)
    json.dump({"cfg": cfg, "score": s, "row": best_i, "n_candidates": len(rows),
               "n_passing": len(kept), "rule": "max learn_margin, constraints C2/C3/C4 of PROTOCOL_LEARN"},
              open(dst, "w"), indent=1)
    print(f"\nwrote {dst}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs/ecosearch_learn/ecosearch.jsonl",
         sys.argv[2] if len(sys.argv) > 2 else "results/learn_eco.json")
