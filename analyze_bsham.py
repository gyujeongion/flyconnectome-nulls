#!/usr/bin/env python3
"""Boundary-targeted, shortcut-free sham (PROTOCOL_BSHAM.md).

Decision rules, frozen before the data existed (2026-09-24, research log 54):

  Seed value   per seed, the mean over P0T0 and P1T0; n = 10. Mean over seeds with a percentile bootstrap 90 %
               interval (10,000 resamples, seed 12345), as in the dose analysis.
  P (primary)  AS10 - BS10 at the last common-garden probe (generation 550):
               SHORTCUT-SPECIFIC if the lower bound > +0.10; BOUNDARY REWIRING SUFFICES if the interval lies inside
               +-0.10; INTERMEDIATE otherwise.
  R (primary robustness)  the same contrast on the fresh-world evaluation of generation-599 elites; ROBUST if it
               falls in the same category as P, otherwise both are reported and the claim is described as not robust.
  S1  AS<k> - BS<k> increases with k (Page's L on per-seed differences, increasing), p < 0.05.
  S2  per-ecology AS10 - BS10.  S3  BS10 - A (+-0.10 reference).  S4  olfaction-ablation cost AS10 vs BS10.
  S5  AS10 - A, lower bound > 0.

Before any number is computed: all 20 runs finished, configurations equal apart from run, seed and pred_fixed,
and all 80 boundary-sham constructions complete with the shortcut share unchanged and degrees kept; else exit 1.
`--selftest` runs the decision logic on synthetic data with known effects.
"""
import json, os, re, sys, tempfile, shutil
import numpy as np
from analyze_dose import page_L  # same Page's L as the dose analysis

DOSES = ["1", "3", "5", "10"]
CONDS = ["A"] + [f"AS{d}" for d in DOSES] + [f"BS{d}" for d in DOSES]
ECO = {"0.0": "P0T0", "0.2": "P1T0"}
SEEDS = range(10)
DELTA = 0.10
BS_LINE = re.compile(r"^(BS\d+): (\d+)/(\d+) boundary rotations, smell->DN share ([0-9.]+) -> ([0-9.]+); .*in/out-degree kept (True|False); "
                     r"\|w\(Z->DN\)\|/\|w\(ORN->X\)\| median ([0-9.na]+)", re.M)


def boot(d, seed=12345, n=10000):
    d = np.asarray(d, float); b = np.random.default_rng(seed).choice(d, (n, len(d))).mean(1)
    return float(d.mean()), [float(np.percentile(b, 5)), float(np.percentile(b, 95))]


def verdict(lo, hi):
    if lo > DELTA: return "SHORTCUT-SPECIFIC"
    if lo > -DELTA and hi < DELTA: return "BOUNDARY REWIRING SUFFICES"
    return "INTERMEDIATE"


def run_name(pv, s): return f"bsham_P{pv}_T0.0_s{s}"


def checks(root):
    bad = []; ref = None; ratios = []
    for pv in ECO:
        for s in SEEDS:
            r = run_name(pv, s); d = f"{root}/runs/{r}"
            if not (os.path.exists(f"{d}/done.flag") and os.path.exists(f"{d}/config.json") and os.path.exists(f"{d}/fresh.jsonl")):
                bad.append((r, "missing or unfinished")); continue
            a = json.load(open(f"{d}/config.json"))["args"]
            if a.get("pred_fixed") != float(pv) or a.get("seed") != s:
                bad.append((r, "pred_fixed/seed mismatch"))
            if ref is None: ref = a
            diff = {k for k in set(a) | set(ref) if k not in ("run", "seed", "pred_fixed") and a.get(k) != ref.get(k)}
            if diff: bad.append((r, sorted(diff)))
            log = open(f"{root}/runs_{r}.log").read() if os.path.exists(f"{root}/runs_{r}.log") else ""
            got = BS_LINE.findall(log)
            if len(got) != 4: bad.append((r, f"{len(got)} boundary-sham construction lines")); continue
            for cond, done, want, sh0, sh1, kept, med in got:
                if done != want or kept != "True" or float(sh1) > float(sh0) + 1e-9:
                    bad.append((r, f"{cond} construction failed"))
                if med not in ("nan",): ratios.append(float(med))
    return bad, ratios


def load(root):
    last, fresh = {}, {}
    for pv in ECO:
        for s in SEEDS:
            d = f"{root}/runs/{run_name(pv, s)}"
            p = None
            for line in open(f"{d}/probe.jsonl"):
                p = line
            p = json.loads(p); assert p["gen"] == 550, f"{d}: last probe {p['gen']}"
            last[(pv, s)] = {c: (p[f"{c}.normal"]["fit"], p[f"{c}.anosmic"]["fit"]) for c in CONDS}
            fl = [json.loads(l) for l in open(f"{d}/fresh.jsonl")]
            assert len(fl) == 4, f"{d}: {len(fl)} fresh batches"
            fresh[(pv, s)] = {c: float(np.mean([f[f"{c}.normal"]["fit"] for f in fl])) for c in CONDS}
    return last, fresh


def analyse(root, quiet=False):
    bad, ratios = checks(root)
    if bad:
        print("CHECK FAILED:"); [print("  ", b) for b in bad[:10]]; sys.exit(1)
    say = (lambda *a: None) if quiet else print
    say(f"checks passed: 20 runs, 80 boundary-sham constructions complete; |w| ratio medians {min(ratios):.2f}-{max(ratios):.2f}")
    last, fresh = load(root)
    sm = lambda f: [np.mean([f(pv, s) for pv in ECO]) for s in SEEDS]
    out = {}
    m, ci = boot(sm(lambda pv, s: last[(pv, s)]["AS10"][0] - last[(pv, s)]["BS10"][0]))
    out["P"] = {"mean": m, "ci90": ci, "verdict": verdict(*ci)}
    m2, ci2 = boot(sm(lambda pv, s: fresh[(pv, s)]["AS10"] - fresh[(pv, s)]["BS10"]))
    out["R"] = {"mean": m2, "ci90": ci2, "verdict": verdict(*ci2), "robust": verdict(*ci2) == verdict(*ci)}
    mat = np.array([[np.mean([last[(pv, s)][f"AS{d}"][0] - last[(pv, s)][f"BS{d}"][0] for pv in ECO]) for d in DOSES] for s in SEEDS])
    L, z, p = page_L(mat)
    out["S1"] = {"means": dict(zip(DOSES, mat.mean(0).round(4).tolist())), "page_L": L, "z": z, "p": p, "supported": bool(p < 0.05)}
    out["S2"] = {lab: dict(zip(("mean", "ci90"), boot([last[(pv, s)]["AS10"][0] - last[(pv, s)]["BS10"][0] for s in SEEDS]))) for pv, lab in ECO.items()}
    out["S3"] = dict(zip(("mean", "ci90"), boot(sm(lambda pv, s: last[(pv, s)]["BS10"][0] - last[(pv, s)]["A"][0]))))
    cost = lambda c: (lambda pv, s: last[(pv, s)][c][0] - last[(pv, s)][c][1])
    out["S4"] = {"AS10_cost": boot(sm(cost("AS10")))[0], "BS10_cost": boot(sm(cost("BS10")))[0],
                 "diff": dict(zip(("mean", "ci90"), boot(sm(lambda pv, s: cost("AS10")(pv, s) - cost("BS10")(pv, s)))))}
    m5, ci5 = boot(sm(lambda pv, s: last[(pv, s)]["AS10"][0] - last[(pv, s)]["A"][0]))
    out["S5"] = {"mean": m5, "ci90": ci5, "replicated": ci5[0] > 0}
    say(f"P  AS10-BS10 (last probe): {m:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}] -> {out['P']['verdict']}")
    say(f"R  AS10-BS10 (fresh worlds): {m2:+.3f} [{ci2[0]:+.3f},{ci2[1]:+.3f}] -> {out['R']['verdict']} ({'ROBUST' if out['R']['robust'] else 'NOT ROBUST'})")
    say(f"S1 AS-BS by dose {out['S1']['means']}  Page z = {z:.2f}, p = {p:.3g}")
    for lab, r in out["S2"].items(): say(f"S2 {lab} AS10-BS10 {r['mean']:+.3f} [{r['ci90'][0]:+.3f},{r['ci90'][1]:+.3f}]")
    say(f"S3 BS10-A {out['S3']['mean']:+.3f} [{out['S3']['ci90'][0]:+.3f},{out['S3']['ci90'][1]:+.3f}]")
    say(f"S4 olfaction cost AS10 {out['S4']['AS10_cost']:.3f} vs BS10 {out['S4']['BS10_cost']:.3f}")
    say(f"S5 AS10-A {m5:+.3f} [{ci5[0]:+.3f},{ci5[1]:+.3f}] -> {'replicated' if out['S5']['replicated'] else 'not replicated'}")
    return out


def selftest():
    """Synthetic runs with known AS10 - BS10 effects; the verdicts must match."""
    ok = True
    # (effect, noise SD, expected): the third case has an interval that straddles +0.10
    for effect, sd, want in ((0.50, 0.02, "SHORTCUT-SPECIFIC"), (0.00, 0.02, "BOUNDARY REWIRING SUFFICES"), (0.10, 0.10, "INTERMEDIATE")):
        root = tempfile.mkdtemp(); rng = np.random.default_rng(1)
        for pv in ECO:
            for s in SEEDS:
                r = run_name(pv, s); d = f"{root}/runs/{r}"; os.makedirs(d)
                json.dump({"args": {"seed": s, "pred_fixed": float(pv), "run": r, "pop": 4608}}, open(f"{d}/config.json", "w"))
                open(f"{d}/done.flag", "w").close()
                base = 2.5 + rng.normal(0, 0.05)
                vals = {c: base for c in CONDS}
                for k, dd in enumerate(DOSES):
                    vals[f"AS{dd}"] = base + effect * (k + 1) / 4 + rng.normal(0, sd)
                    vals[f"BS{dd}"] = base + rng.normal(0, sd)
                rec = {"gen": 550, **{f"{c}.normal": {"fit": v} for c, v in vals.items()}, **{f"{c}.anosmic": {"fit": v - 0.3} for c, v in vals.items()}}
                open(f"{d}/probe.jsonl", "w").write(json.dumps({"gen": 500}) + "\n" + json.dumps(rec) + "\n")
                with open(f"{d}/fresh.jsonl", "w") as f:
                    for b in range(4): f.write(json.dumps({f"{c}.normal": {"fit": v + rng.normal(0, 0.01)} for c, v in vals.items()}) + "\n")
                open(f"{root}/runs_{r}.log", "w").write("".join(
                    f"BS{dd}: 10/10 boundary rotations, smell->DN share 0.0001 -> 0.0001; nnz 1->1; in/out-degree kept True; |w(Z->DN)|/|w(ORN->X)| median 1.00\n" for dd in DOSES))
        res = analyse(root, quiet=True)["P"]; got = res["verdict"]
        print(f"  effect {effect:+.2f} (noise {sd}): {res['mean']:+.3f} [{res['ci90'][0]:+.3f},{res['ci90'][1]:+.3f}] -> {got} (expected {want})"); ok &= got == want
        shutil.rmtree(root)
    root = tempfile.mkdtemp(); os.makedirs(f"{root}/runs")
    try:
        analyse(root, quiet=True); print("  missing data did NOT stop the analysis"); ok = False
    except SystemExit as e:
        print(f"  missing data -> exit {e.code} (expected 1)"); ok &= e.code == 1
    shutil.rmtree(root)
    print("SELFTEST", "PASSED" if ok else "FAILED"); sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    res = analyse(".")
    os.makedirs("results", exist_ok=True)
    json.dump({"decision_rules": __doc__, **res}, open("results/bsham_key.json", "w"), indent=1)
    print("wrote results/bsham_key.json")
