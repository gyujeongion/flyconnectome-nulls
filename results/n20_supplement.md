
## S13. Extension of the corrected grid to 20 seeds (P0T0 and P1T0)

Seeds 10 to 19 of the two T0 ecologies were run with the corrected grid's frozen engine and flags after rules R1 to R3 were written (research log 46; `analyze_n20.py`). Every one of the 40 runs has the same configuration as seed 0 apart from name and seed. Endpoint: fitness at the last common-garden probe (generation 550); connectome minus control, median with 90 % bootstrap interval; p_holm from Wilcoxon signed-rank tests at n = 20, Holm-corrected across the four controls within each ecology.

| ecology | control | seeds 0–9 | seeds 10–19 | n = 20 | p_holm (n = 20) | inside ±0.10 (R1) |
|---|---|---|---|---|---|---|
| P0T0 | N1 | −0.621 [−0.661, −0.469] | −0.304 [−0.513, −0.161] | −0.491 [−0.624, −0.298] | 0.001 | no |
| P0T0 | N2 | −0.403 [−0.421, −0.308] | +0.100 [−0.068, +0.332] | −0.214 [−0.375, +0.042] | 0.248 | no |
| P0T0 | N4 | −0.064 [−0.255, +0.002] | +0.090 [−0.109, +0.217] | −0.011 [−0.114, +0.080] | 0.841 | no |
| P0T0 | N5 | −0.264 [−0.368, −0.179] | +0.031 [−0.068, +0.152] | −0.104 [−0.208, +0.016] | 0.265 | no |
| P1T0 | N1 | −0.312 [−0.441, −0.197] | −0.417 [−0.522, +0.067] | −0.333 [−0.482, −0.197] | 0.004 | no |
| P1T0 | N2 | −0.338 [−0.530, −0.074] | −0.101 [−0.335, +0.103] | −0.203 [−0.430, −0.032] | 0.046 | no |
| P1T0 | N4 | +0.021 [−0.026, +0.038] | −0.027 [−0.168, +0.050] | +0.003 [−0.042, +0.029] | 1.000 | yes |
| P1T0 | N5 | +0.067 [−0.067, +0.090] | +0.003 [−0.175, +0.144] | +0.056 [−0.097, +0.090] | 1.000 | yes |

**R2 (replication of N5's lead in P0T0 in seeds 10 to 19):** A − N5 = +0.031 [−0.068, +0.152] → **NOT REPLICATED**.

**R3 (two-ecology pool, seed means over P0T0 and P1T0, n = 20; a different estimand from Table 2's four-ecology seed means):**

| control | median [90 % CI] | seeds with connectome ahead | inside ±0.10 |
|---|---|---|---|
| N1 | −0.401 [−0.497, −0.241] | 3/20 | no |
| N2 | −0.155 [−0.273, −0.060] | 4/20 | no |
| N4 | 0.000 [−0.087, +0.083] | 10/20 | yes |
| N5 | −0.059 [−0.093, +0.018] | 7/20 | yes |
