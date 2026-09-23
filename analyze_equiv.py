"""Council round 2 requests: paired effects per ecology with CIs, TOST-style equivalence, interaction tests,
and the change in the transplant effect from generation 0 to 600. All post hoc."""
import json, glob, os, numpy as np
from scipy.stats import wilcoxon
rng = np.random.default_rng(0)
DELTA = 0.10          # smallest effect of interest: half the shortcut artefact (|N2-A| = 0.20-0.35)
def boot_ci(d, lvl=95):
    b = np.median(rng.choice(d, (10000, len(d))), 1); a = (100 - lvl) / 2
    return np.percentile(b, a), np.percentile(b, 100 - a)
def load(pat, conds, cells):
    E = {}
    for (pv, pt), lab in cells.items():
        for r in sorted(glob.glob(f"runs/{pat}_P{pv}_T{pt}_s*")):
            if not os.path.exists(f"{r}/done.flag"): continue
            p = [json.loads(l) for l in open(f"{r}/probe.jsonl")]
            E[(lab, int(r.split('_s')[-1]))] = {c: dict(gen0=p[0][f"{c}.normal"]["fit"], final=p[-1][f"{c}.normal"]["fit"]) for c in conds}
    return E
def report(name, d):
    lo90, hi90 = boot_ci(d, 90); lo95, hi95 = boot_ci(d, 95)
    p = wilcoxon(d).pvalue if np.any(d != 0) else 1.0
    eq = f"EQUIVALENT within +-{DELTA:.2f}" if (lo90 > -DELTA and hi90 < DELTA) else f"not equivalent (90% CI leaves +-{DELTA:.2f})"
    print(f"   {name:26s} median {np.median(d):+.3f}  95% CI [{lo95:+.3f},{hi95:+.3f}]  90% CI [{lo90:+.3f},{hi90:+.3f}]  "
          f"{int((d > 0).sum())}/{len(d)}  p={p:.4f}  {eq}")

FIX = load("ecofix", ["A", "N1", "N2", "N4", "N5"], {("0.0","0.0"):"P0T0",("0.2","0.0"):"P1T0",("0.0","1.0"):"P0T1",("0.2","1.0"):"P1T1"})
SW = load("swap", ["A", "AS", "N2", "N2R", "N4"], {("0.0","0.0"):"P0T0",("0.2","0.0"):"P1T0"})
seeds = sorted({k[1] for k in FIX})
print(f"corrected grid: {len(FIX)} runs; transplant: {len(SW)} runs; delta = {DELTA}\n")
print("=== corrected grid, generation 600, connectome minus control ===")
for lab in ["P0T0", "P1T0", "P0T1", "P1T1", "pooled (seed means)"]:
    print(f" {lab}")
    for n in ["N1", "N2", "N4", "N5"]:
        if lab.startswith("pooled"):
            d = np.array([np.mean([FIX[(l, s)]["A"]["final"] - FIX[(l, s)][n]["final"] for l in ["P0T0","P1T0","P0T1","P1T1"] if (l, s) in FIX]) for s in seeds])
        else:
            d = np.array([FIX[(lab, s)]["A"]["final"] - FIX[(lab, s)][n]["final"] for s in seeds if (lab, s) in FIX])
        report(f"A - {n}", d)
print("\n=== transplant experiment, generation 600 ===")
for lab in ["P0T0", "P1T0", "pooled (seed means)"]:
    print(f" {lab}")
    for x, y in [("AS", "A"), ("N2", "A"), ("N2R", "N2"), ("N2R", "A"), ("N4", "A")]:
        if lab.startswith("pooled"):
            d = np.array([np.mean([SW[(l, s)][x]["final"] - SW[(l, s)][y]["final"] for l in ["P0T0","P1T0"] if (l, s) in SW]) for s in seeds])
        else:
            d = np.array([SW[(lab, s)][x]["final"] - SW[(lab, s)][y]["final"] for s in seeds if (lab, s) in SW])
        report(f"{x} - {y}", d)
print("\n=== interactions (paired, per seed) ===")
d = np.array([(SW[("P0T0", s)]["N2R"]["final"] - SW[("P0T0", s)]["N2"]["final"]) - (SW[("P1T0", s)]["N2R"]["final"] - SW[("P1T0", s)]["N2"]["final"]) for s in seeds])
report("removal effect: safe - pred", d)
d = np.array([np.mean([SW[(l, s)]["AS"]["final"] - SW[(l, s)]["A"]["final"] for l in ["P0T0","P1T0"]]) -
              np.mean([SW[(l, s)]["AS"]["gen0"] - SW[(l, s)]["A"]["gen0"] for l in ["P0T0","P1T0"]]) for s in seeds])
report("transplant: gen600 - gen0", d)
print(f"\n=== seeds needed for equivalence within +-{DELTA:.2f} (paired SD from the data) ===")
for name, d in [("A - N4 (pooled)", np.array([np.mean([FIX[(l, s)]["A"]["final"] - FIX[(l, s)]["N4"]["final"] for l in ["P0T0","P1T0","P0T1","P1T1"]]) for s in seeds])),
                ("A - N5 (pooled)", np.array([np.mean([FIX[(l, s)]["A"]["final"] - FIX[(l, s)]["N5"]["final"] for l in ["P0T0","P1T0","P0T1","P1T1"]]) for s in seeds]))]:
    sd = d.std(ddof=1); true = abs(np.median(d))
    n_req = next((n for n in range(5, 400) if (1.645 * sd / np.sqrt(n)) < (DELTA - true)), None) if true < DELTA else None
    print(f"   {name}: median {np.median(d):+.3f}, SD {sd:.3f} -> n ≈ {n_req if n_req else '>400 (or |effect| ≥ delta)'} seeds for a 90% CI inside +-{DELTA}")
