"""Activity-distribution diagnostics for every paper-v2 run (GPT council request): parse the homeostasis lines."""
import glob, re, ast, numpy as np
rows = {}
for f in sorted(glob.glob("runs_eco_P*_s*.log")):
    for l in open(f):
        m = re.match(r"homeostasis (A|N1|N2|N3): (\{.*?\})", l)
        if m: rows.setdefault(m.group(1), []).append(ast.literal_eval(m.group(2)))
print("condition | runs | mean activity | activity SD across groups | silent groups (of 512) | saturated fraction | KC active fraction")
for c in ["A", "N1", "N2", "N3"]:
    R = rows.get(c, [])
    q = lambda k: f"{np.median([r[k] for r in R]):.4f} [{np.min([r[k] for r in R]):.4f}-{np.max([r[k] for r in R]):.4f}]"
    print(f"{c} | {len(R)} | {q('mean')} | {q('group_act_sd')} | {q('silent_groups')} | {q('sat_frac')} | {q('kc_active')}")
