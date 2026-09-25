#!/usr/bin/env python3
"""Regenerate every manuscript number that had no saved source, and the per-seed quantities
that v4 quoted from seed 0 as if they were general.

The numerical audit of v4 (2026-09-23) found values in PAPER.md that no checked-in script
produces. This script produces them, with the paper's own conventions, and writes them to
results/paper_numbers.json so that each one has a file behind it.

Conventions follow the existing analyses: the estimand is the median across seeds of each
seed's mean difference over ecologies; intervals are percentile bootstraps of that median over
resampled seeds (analyze_swap.py, analyze_fix.py).
"""
import json, glob, os, re
from collections import deque
import numpy as np

RNG = np.random.default_rng(20260923)
NBOOT = 10000
out = {}


def med_ci(d, lo, hi):
    d = np.asarray(d, float)
    b = np.median(RNG.choice(d, (NBOOT, len(d))), 1)
    return float(np.percentile(b, lo)), float(np.percentile(b, hi))


def probes(run):
    return [json.loads(l) for l in open(f"{run}/probe.jsonl")]


def auc(ps, cond):
    g = np.array([p["gen"] for p in ps], float)
    y = np.array([p[f"{cond}.normal"]["fit"] for p in ps], float)
    return float(np.trapezoid(y, g) / (g[-1] - g[0]))


# ---------------------------------------------------------------- 1. first grid, primary endpoint
first = [r for r in sorted(glob.glob("runs/eco_P*_T*_s*"))
         if re.fullmatch(r"runs/eco_P[0-9.]+_T[0-9.]+_s\d+", r) and os.path.exists(f"{r}/done.flag")]
per_seed = {}
for r in first:
    s = int(r.split("_s")[-1]); ps = probes(r)
    per_seed.setdefault(s, []).append({c: auc(ps, c) for c in ("A", "N1", "N2")})
fg = {}
for null in ("N1", "N2"):
    seed_means = np.array([np.mean([e["A"] - e[null] for e in per_seed[s]]) for s in sorted(per_seed)])
    pooled = np.mean([e["A"] - e[null] for s in per_seed for e in per_seed[s]])
    fg[null] = {"seed_mean_median": float(np.median(seed_means)), "pooled_run_mean": float(pooled),
                "n_seeds": len(seed_means), "n_runs": sum(len(v) for v in per_seed.values())}
out["first_grid_primary_auc_A_minus_null"] = fg

# ---------------------------------------------------------------- 2. intervention, generation 0
sw = {}
for r in sorted(glob.glob("runs/swap_P*_T0.0_s*")):
    if not os.path.exists(f"{r}/done.flag"):
        continue
    s = int(r.split("_s")[-1]); p0 = probes(r)[0]
    sw.setdefault(s, []).append(p0["AS.normal"]["fit"] - p0["A.normal"]["fit"])
d0 = np.array([np.mean(v) for _, v in sorted(sw.items())])
out["intervention_gen0_AS_minus_A"] = {
    "median": float(np.median(d0)), "ci90": med_ci(d0, 5, 95), "ci95": med_ci(d0, 2.5, 97.5),
    "n_seeds": len(d0), "n_positive": int((d0 > 0).sum())}

# ---------------------------------------------------------------- 3. interior assay, smallest intact fitness
intact = []
for r in sorted(glob.glob("runs/ecofix_P*_s*")):
    f = f"{r}/assay.jsonl"
    if not os.path.exists(f):
        continue
    for l in open(f):
        try:
            d = json.loads(l)
        except json.JSONDecodeError:
            continue
        if d.get("assay") == "intact" and d.get("gen") == 500:
            intact += [v["fit"] for v in d["res"].values()]
out["interior_assay_intact_fitness"] = {"min": float(min(intact)), "median": float(np.median(intact)),
                                        "n": len(intact)} if intact else None

# ---------------------------------------------------------------- 4. olfactory -> descending paths
B = np.load("data/brain.npz", allow_pickle=True)
W0 = B["W0"]; K = W0.shape[0]
orn = sorted(set(int(i) for i in B["in_orn"])); dn = sorted(set(int(i) for i in B["dn"]))
A = np.abs(W0) > 0
dist = np.full(K, 99); q = deque()
for s_ in orn:
    dist[s_] = 0; q.append(s_)
while q:
    u = q.popleft()
    for v in np.nonzero(A[:, u])[0]:          # rows are targets: W[i, j] = j -> i
        if dist[v] == 99:
            dist[v] = dist[u] + 1; q.append(v)
hops = [int(dist[d]) for d in dn]
names = None
for key in ("group_names", "names", "labels", "gnames"):
    if key in B.files:
        names = [str(x) for x in B[key]]; break
direct = []
for i in dn:
    for j in orn:
        if W0[i, j] != 0:
            direct.append({"src": names[j] if names else int(j), "dst": names[i] if names else int(i),
                           "weight": float(W0[i, j])})
share = float(np.abs(W0[np.ix_(dn, orn)]).sum() / (np.abs(W0[:, orn]).sum() + 1e-12))
out["olfactory_to_descending"] = {
    "n_descending_groups": len(dn), "hop_counts": {str(h): hops.count(h) for h in sorted(set(hops))},
    "direct_edges": direct, "n_edges_total": int(A.sum()), "raw_shortcut_share": share,
    "names_available": names is not None}

# ---------------------------------------------------------------- 5. mutational neighbourhood intervals
txt = open("results/revision_analysis.txt", encoding="utf-8").read()
bounds, within, total = [], 0, 0
for m in re.finditer(r"CI\[([+-][0-9.]+),([+-][0-9.]+)\]\s+(NOT within|within)", txt):
    lo_, hi_ = float(m.group(1)), float(m.group(2))
    bounds += [lo_, hi_]; total += 1; within += (m.group(3) == "within")
out["mutational_neighbourhood_ci95"] = {
    "max_abs_bound_pp": round(100 * max(abs(b) for b in bounds), 2), "n_comparisons": total,
    "within_pm2pp": within, "margin_set_in_analysis_pp": 2.0}

# ---------------------------------------------------------------- 6. boundary-preserving nulls, interior overlap
ov = {"N4": [], "N5": []}
for f in sorted(glob.glob("runs_ecofix_P*_s*.log")):
    if f.endswith("_assay.log") or int(re.search(r"_s(\d+)\.log$", f).group(1)) >= 10:
        continue  # registered grid only (seeds 0-9); the n = 20 extension is reported separately
    for l in open(f, errors="ignore"):
        m = re.match(r"(N4|N5): interface .* interior overlap ([0-9.]+)", l)
        if m:
            ov[m.group(1)].append(float(m.group(2)))
out["interior_overlap_with_connectome"] = {c: {"min": min(v), "max": max(v), "n": len(v)} for c, v in ov.items() if v}

# ---------------------------------------------------------------- 7. dose grid, per-seed transplant sizes
dose = {}
for f in sorted(glob.glob("runs_dose_P0.0_T0.0_s*.log")):
    if f.endswith("_assay.log"):
        continue
    s = int(re.search(r"_s(\d+)\.log$", f).group(1))
    t = open(f, errors="ignore").read()
    for m in re.finditer(r"(AS\d+): target [0-9.]+ -> smell->DN share ([0-9.]+) after (\d+) swaps", t):
        dose.setdefault(m.group(1), []).append((s, float(m.group(2)), int(m.group(3))))
    m = re.search(r"SHAM: (\d+)/(\d+) rewiring swaps", t)
    if m:
        dose.setdefault("SHAM", []).append((s, None, int(m.group(1))))
seeds_run = sorted({s for v in dose.values() for s, _, _ in v if s < 10})
out["dose_transplant_by_seed"] = {
    c: {"share_pct": [round(100 * x, 2) for s, x, _ in v if s < 10 and x is not None],
        "swaps": [n for s, _, n in v if s < 10]}
    for c, v in dose.items()}
summ = {}
for c, v in dose.items():
    sh = [100 * x for s, x, _ in v if s < 10 and x is not None]; nn = [n for s, _, n in v if s < 10]
    summ[c] = {"share_pct_range": [round(min(sh), 2), round(max(sh), 2)] if sh else None,
               "swaps_range": [min(nn), max(nn)], "swaps_mean": round(float(np.mean(nn)), 1)}
out["dose_transplant_summary"] = summ


# ---------------------------------------------------------------- 8. transplant vs sham, per seed
def last(run):
    p = None
    for l in open(f"{run}/probe.jsonl"):
        p = l
    return json.loads(p)


as10 = {s: n for s, _, n in dose.get("AS10", []) if s < 10}
rows = []
for s in sorted(as10):
    g = np.mean([last(f"runs/dose_P{pv}_T0.0_s{s}")["AS10.normal"]["fit"]
                 - last(f"runs/dose_P{pv}_T0.0_s{s}")["SHAM.normal"]["fit"] for pv in ("0.0", "0.2")])
    rows.append((s, as10[s], float(g)))
sham_n = dose["SHAM"][0][2]
fewer = [g for s, n, g in rows if n < sham_n]; more = [g for s, n, g in rows if n >= sham_n]
out["transplant_vs_sham_by_seed"] = {
    "sham_swaps": sham_n, "rows": [{"seed": s, "as10_swaps": n, "gap": round(g, 4)} for s, n, g in rows],
    "fewer_swaps": {"n": len(fewer), "n_won": int(sum(g > 0 for g in fewer)), "mean_gap": float(np.mean(fewer))},
    "more_swaps": {"n": len(more), "n_won": int(sum(g > 0 for g in more)), "mean_gap": float(np.mean(more))},
    "corr_swapdiff_gap": float(np.corrcoef([n - sham_n for _, n, _ in rows], [g for *_, g in rows])[0, 1])}
# the paper's "full dose ahead of the sham by 0.534 (90 % CI ...)": mean over seeds of the seed-mean gap,
# percentile bootstrap of the mean, same convention and seed as analyze_dose.py (independent stream)
_g = np.array([g for *_, g in rows]); _b = np.random.default_rng(12345).choice(_g, (10000, len(_g))).mean(1)
out["transplant_vs_sham_by_seed"].update({"mean_gap_all": float(_g.mean()), "n_won_all": int((_g > 0).sum()),
                                          "ci90_mean_gap": [float(np.percentile(_b, 5)), float(np.percentile(_b, 95))]})


# ---------------------------------------------------------------- 9. numbers the NeurIPS draft quoted without a file
# (audit of 2026-09-24). No random draws here, so the values above are unchanged.
def ecol_gap(pv):
    g = [last(f"runs/dose_P{pv}_T0.0_s{s}")["AS10.normal"]["fit"] - last(f"runs/dose_P{pv}_T0.0_s{s}")["SHAM.normal"]["fit"]
         for s in sorted(as10)]
    return {"mean_gap": float(np.mean(g)), "n_won": int(sum(x > 0 for x in g)), "n": len(g)}


out["transplant_vs_sham_by_ecology"] = {"P0T0": ecol_gap("0.0"), "P1T0": ecol_gap("0.2")}

# calibration spread in the corrected grid (seeds 0-9, four ecologies), per condition: median and range over runs
HOMEO = re.compile(r"^homeostasis (\w+): (\{.*?\})", re.M)
cal = {}
for f in sorted(glob.glob("runs_ecofix_P*_T*_s*.log")):
    m = re.search(r"_s(\d+)\.log$", f)            # skips the *_assay.log files
    if not m or int(m.group(1)) >= 10:
        continue
    for c, d in HOMEO.findall(open(f).read()):
        d = eval(d, {"__builtins__": {}})
        cal.setdefault(c, []).append((d["group_act_sd"], d["silent_groups"]))
out["calibration_corrected_grid"] = {
    c: {"runs": len(v), "sd_median": float(np.median([a for a, _ in v])), "sd_range": [min(a for a, _ in v), max(a for a, _ in v)],
        "silent_median": float(np.median([b for _, b in v])), "silent_range": [min(b for _, b in v), max(b for _, b in v)]}
    for c, v in sorted(cal.items())}

# wall time of one 5-condition, 600-generation corrected-grid run: RUN -> OK in the follow-up driver log
ts = lambda l: np.datetime64(l[:19].replace(" ", "T"))
dur, start = [], {}
for l in open("followup.log"):
    m = re.search(r"(RUN|OK) (ecofix_P\S+)", l)
    if m and m.group(1) == "RUN": start[m.group(2)] = ts(l)
    elif m and m.group(2) in start: dur.append(float((ts(l) - start.pop(m.group(2))) / np.timedelta64(1, "m")))
out["run_minutes_corrected_grid"] = {"n": len(dur), "median": float(np.median(dur)), "range": [min(dur), max(dur)]} if dur else None

# size of the raw run outputs that are not redistributed
_raw = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk("runs") for f in fs)
_logs = sum(os.path.getsize(f) for f in glob.glob("runs_*.log"))
out["raw_outputs_gb"] = {"runs_dirs": round(_raw / 1e9, 2), "run_logs": round(_logs / 1e9, 2), "total": round((_raw + _logs) / 1e9, 1)}

os.makedirs("results", exist_ok=True)
json.dump(out, open("results/paper_numbers.json", "w"), indent=1)
print(json.dumps(out, indent=1, ensure_ascii=False))
