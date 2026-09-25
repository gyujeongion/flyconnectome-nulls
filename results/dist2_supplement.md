
## S15. Corrected grid under distribution-matched calibration (PROTOCOL_DIST2, D-06)

Frozen before the runs: protocol `PROTOCOL_DIST2.md` (`0124b50ad5345279`), analysis `analyze_dist2.py` (`d12b803c141faa1e`, unchanged at analysis), the corrected grid's frozen engine `evo_fix.py` and flags with `--homeo distmatch` in place of `global`. All 40 runs (4 ecologies × seeds 0–9) finished; each configuration equals its global-calibration counterpart apart from run name and calibration. Output `results/dist2_key.json`.

**Manipulation check** (registered: at least 36 of 40 runs with |SD − SD_A| ≤ 0.01, |silent − silent_A| ≤ 10 and |saturated − saturated_A| ≤ 0.02 for N4 and N5): 40/40 runs → PASSED. Calibrated activity by condition, medians over the 40 runs:

| condition | between-group SD | silent groups | Kenyon-cell mean |
|---|---|---|---|
| A | 0.118 | 7 | 0.050 |
| N1 | 0.118 | 3 | 0.050 |
| N2 | 0.118 | 6.5 | 0.050 |
| N4 | 0.118 | 4.5 | 0.050 |
| N5 | 0.118 | 4.5 | 0.050 |

Under global calibration the same quantities were SD 0.163 / 46 silent groups for the connectome against 0.080–0.087 / 8–19 for the controls (main text); distribution matching removes that difference. N1's Kenyon-cell mean ranged 0.038–0.050 (a known limitation of the scheme for the column shuffle).

**P1 and P2 — seed means over four ecologies, median with 90 % bootstrap interval:**

| contrast | global (registered grid) | distribution-matched | verdict |
|---|---|---|---|
| A − N4 | +0.002 [−0.026, +0.045] | +0.055 [+0.035, +0.097] | EQUIVALENT |
| A − N5 | −0.074 [−0.096, +0.002] | −0.009 [−0.076, +0.046] | EQUIVALENT |
| A − N1 | −0.224 [−0.337, −0.103] | −0.274 [−0.306, −0.123] | BEHIND |
| A − N2 | −0.199 [−0.321, −0.133] | −0.146 [−0.231, −0.037] | BEHIND |

**S1 — change caused by calibration, seed mean of (A − N) under distribution matching minus (A − N) under global:**

| control | change | verdict |
|---|---|---|
| N1 | +0.010 [−0.112, +0.072] | UNDETERMINED |
| N2 | +0.101 [+0.016, +0.133] | CHANGED |
| N4 | +0.062 [−0.002, +0.117] | UNDETERMINED |
| N5 | +0.076 [−0.017, +0.128] | UNDETERMINED |

**S2 — N5 ahead of the connectome in P0T0 under distribution matching:** A − N5 = −0.066 [−0.320, +0.124] → NOT REPLICATED.

**Per ecology (descriptive), A minus control, median [90 % CI]:**

| ecology | control | global | distribution-matched |
|---|---|---|---|
| P0T0 | N1 | −0.621 [−0.661, −0.469] | −0.437 [−0.563, −0.235] |
| P0T0 | N2 | −0.403 [−0.421, −0.308] | −0.478 [−0.677, −0.052] |
| P0T0 | N4 | −0.064 [−0.255, +0.002] | +0.139 [−0.128, +0.230] |
| P0T0 | N5 | −0.264 [−0.368, −0.179] | −0.066 [−0.320, +0.124] |
| P1T0 | N1 | −0.312 [−0.441, −0.197] | −0.489 [−0.664, −0.275] |
| P1T0 | N2 | −0.338 [−0.530, −0.074] | −0.164 [−0.483, −0.013] |
| P1T0 | N4 | +0.021 [−0.026, +0.038] | +0.051 [−0.120, +0.144] |
| P1T0 | N5 | +0.067 [−0.067, +0.090] | +0.041 [−0.204, +0.214] |
| P0T1 | N1 | −0.159 [−0.205, +0.010] | −0.096 [−0.182, +0.010] |
| P0T1 | N2 | −0.017 [−0.333, +0.005] | −0.020 [−0.122, +0.096] |
| P0T1 | N4 | −0.031 [−0.116, +0.085] | +0.073 [−0.066, +0.104] |
| P0T1 | N5 | −0.015 [−0.102, +0.061] | +0.036 [−0.086, +0.089] |
| P1T1 | N1 | +0.054 [−0.054, +0.127] | +0.063 [−0.040, +0.099] |
| P1T1 | N2 | +0.015 [−0.132, +0.137] | +0.104 [−0.091, +0.111] |
| P1T1 | N4 | +0.120 [−0.064, +0.283] | −0.012 [−0.099, +0.073] |
| P1T1 | N5 | +0.070 [−0.009, +0.210] | −0.027 [−0.066, +0.078] |
