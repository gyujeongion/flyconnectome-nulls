---
title: "Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison"
author: |
  Gyujeong Park\
  Independent Researcher\
  ORCID 0009-0006-2989-728X\
  `ionlabs2025@gmail.com`
abstract: |
  Randomised copies of a connectome are the usual baseline for asking whether measured wiring matters, and the answer depends on what the randomisation preserves. We evolved embodied foraging agents whose brains are a compressed adult *Drosophila* connectome (FlyWire v783; 512 cell-type groups and 1,000 Kenyon cells) alongside agents built on randomised wiring, in pre-registered experiments with ten seeds, four ecologies and 600 generations. Two standard randomisations, a column shuffle and degree-preserving edge swaps, route 10.6 to 10.7 % of olfactory output directly onto descending motor groups, against 0.012 % in the connectome. On the registered primary endpoint, fitness averaged over the run, no difference was detected; at the last common-garden probe the connectome was behind both controls (−0.22 and −0.20 fitness units on seed means). Against controls that keep every sensory-output and motor-input edge and rewire only the interior, the seed-mean difference lay within a ±0.10 equivalence bound (+0.002 and −0.074, unchanged under a calibration that also matches activity spread), although per ecology the interior column shuffle was ahead by 0.26 in one of four ecologies at ten seeds, a lead that ten further pre-registered seeds did not replicate. Rewiring the connectome so that it acquires the shortcut raised its fitness by 0.44 (10 of 10 seeds) and its dependence on olfaction from 0.15 to 0.99; graded doses raised both in step; at comparable swap counts the full dose was ahead of an interior-only sham by 0.53 (10 of 10 seeds); and a sham that rewired the same boundary edges without creating shortcuts matched the connectome (+0.007) while the full dose was ahead of it by 0.60. What a null preserves at the sensory-motor boundary can decide an evolutionary connectome comparison, and sensory-to-motor path statistics belong next to the degree statistics a null is said to preserve.
---

# Introduction {#sec-intro}

Whole-brain wiring diagrams of the adult fly are public [Dorkenwald et al., 2024; Schlegel et al., 2024], and models constrained by them reproduce feeding circuits [Shiu et al., 2024], visual responses [Lappalainen et al., 2024], head-direction dynamics [Duan et al., 2025] and locomotion [Vaxenburg et al., 2025]. Whether the measured wiring is better than chance is usually asked by comparing it with a randomised network that preserves some of its statistics. The verdict of such a comparison is a statement about the properties the randomisation changed [Váša & Mišić, 2022], so the choice of null decides what is being tested.

![What the null preserves at the sensory-motor boundary, and what it changes. (a) In the compressed connectome olfactory input reaches descending (motor) groups through several synapses; the standard column-shuffle and degree-preserving nulls route about a tenth of olfactory output directly onto them; the boundary-preserving nulls keep every sensory-output and motor-input edge and replace about nine in ten interior edges. (b) Differences in evolved fitness at the last common-garden probe with 90 % intervals: connectome minus each null on seed means over four ecologies (median), and the shortcut interventions over two ecologies (transplant: median; full dose minus sham and sham minus connectome: mean, as registered). The last two rows are the boundary-targeted sham (Section \ref{sec-bsham}), from separate runs with their own connectome island. The grey band is the ±0.10 equivalence bound.](figs/fig_overview.pdf){#fig:overview width=100%}

Our first pre-registered experiment found that randomised wiring out-evolves a compressed fly connectome. Following that result up showed that the two standard randomisations we had used move about a tenth of olfactory receptor output directly onto descending motor neurons, a class of connection that carries 0.012 % of calibrated olfactory output in the connectome, where 52 of 53 descending groups are two or more synapses from olfactory input. This is not an error in those nulls; it is what randomising the wiring between sensory and motor neurons does. But in an embodied task where a short reactive route is useful, it changes what the comparison measures. This paper reports what happens when the null keeps the sensory-motor boundary fixed, and what happens when the shortcut is moved into the connectome and out of the control (Figure \ref{fig:overview}).

Contributions:

- **A measurement** of how two standard randomisations change the sensory-to-motor boundary of a connectome-derived agent: the share of olfactory output reaching descending groups directly rises from 0.012 % to 10.6 to 10.7 %, and the median path length from olfactory receptors to descending groups falls from three synapses to one.
- **A pre-registered comparison** against boundary-preserving controls, which keep every edge out of a sensory group and every edge into a motor group: the connectome's end-of-run disadvantage shrinks to within a ±0.10 equivalence bound on seed means, with the per-ecology heterogeneity and the null primary endpoint reported alongside.
- **Pre-registered interventions** that move the shortcut: shortcut-creating rewiring raises the connectome's evolved fitness and transfers the controls' olfactory strategy to it; shortcut-removing rewiring costs a control its advantage where it had one; the effect rises with dose; an interior-only sham with a comparable number of swaps does not reproduce it; and neither does a sham that rewires the same boundary edges without creating shortcuts.
- **Open protocols, runs and a dated research log**, including a retraction and a correction made after pre-registration, so that each result can be traced to the protocol that preceded it.

# Related work {#sec-related}

**Connectome-constrained models.** Fly connectome models are usually evaluated on how well they reproduce measured activity or behaviour [Shiu et al., 2024; Lappalainen et al., 2024; Duan et al., 2025]. We instead ask whether measured wiring is a better starting point for evolution than randomised wiring, a question closer to Worrell et al. [2017], who found faster sensor-to-effector spreading in the fly connectome than in degree-preserving nulls.

**Null models for connectomes.** A null keeps some properties and randomises the rest, and constraining it to the property under test is standard advice [Váša & Mišić, 2022]; degree-preserving rewiring misses spatially constrained structure [Salova & Kovács, 2025]. Dhiman [2026] found that a fly visual-system model's advantage over degree-preserving nulls largely disappears under a shared random initialisation, a verdict change through initialisation and degree effects rather than path structure. Cho [2026] released a behavioural benchmark for whole-connectome simulations with shuffled and degree-preserving controls. We apply the constrain-the-null advice to the sensory-motor boundary and show a case where it changes the outcome.

**Architecture, randomisation and measurement in machine learning.** Architecture alone can encode useful behaviour [Gaier & Ha, 2019], and which structure a randomised baseline keeps determines what a comparison can attribute to it [Zhou et al., 2019]. Randomisation tests are only informative about the properties they destroy [Adebayo et al., 2018], and an apparent effect can come from the choice of measurement rather than from the system [Schaeffer et al., 2023]. With ten seeds per condition we follow the recommendation to report interval estimates over runs rather than point estimates [Agarwal et al., 2021].

# Methods {#sec-methods}

All quantities are dimensionless simulation units. Full specifications are in Appendix \ref{app-model} (brain, dynamics, world, evolution) and Appendix \ref{app-nulls} (controls).

## Agents {#sec-agents}

From FlyWire v783 we drop optic-lobe, visual-centrifugal and photoreceptor classes, group the remaining neurons by cell type × hemisphere and keep K = 512 groups, starting with a mandatory set of olfactory, gustatory, mushroom-body, dopaminergic, looming and named descending groups and filling the rest by a sensory-to-descending path score. Signed weights come from synapse counts and predicted neurotransmitters [Eckstein et al., 2024]; the matrix has 14,989 non-zero entries, and sensory rows are zeroed so that sensors act as transducers. Group rates follow leaky rectified-tanh dynamics. 1,000 Kenyon cells keep their measured projection-neuron inputs, with APL-like inhibition holding 7 to 9 % active, and a dopamine-gated plastic Kenyon-cell-to-output matrix resets at birth. Forward speed, turning and eating are read out from 24 bilateral descending-neuron types, initialised from named descending neurons (e.g. DNa01/DNa02 for turning), with a per-brain innate feeding reflex. The compression keeps named sensory and motor pathways on purpose, so the conditions compared share this boundary and differ in the wiring between its two sides.

## World, ecologies and evolution {#sec-world}

A unit square holds 24 odor-marked patches (nutritious, potentially toxic, neutral), with the odor-to-role assignment redrawn every lifetime of 400 steps and an energy budget that forces foraging. Four ecologies cross a predator that is stationary (P0) or pursuing (P1) with potentially toxic patches that are all nutritious (T0) or all toxic (T1). Each wiring condition evolves as its own island of 512 agents for 600 generations (tournament selection, 4 % elitism); mutations perturb weights, group sizes, leaks, biases, gains, learning rates and read-outs, and can add or prune edges.

## Controls and interventions {#sec-controls}

Controls are drawn per seed and shared across that seed's ecologies (Table \ref{tbl:nulls}). N1 moves each presynaptic column's non-zero weights to random non-sensory rows; N2 performs Maslov–Sneppen target swaps [Maslov & Sneppen, 2002] with 10·|E| attempts; N4 and N5 apply the same operations to *interior* edges only, whose source is not sensory and whose target is not a descending group. The engine aborts unless N4 and N5 boundary edges are bit-identical to the connectome; they replace about nine in ten interior edges (interior overlap 0.111 to 0.121 for N4 and 0.089 to 0.095 for N5 over the 40 constructions of the corrected grid). The Kenyon-cell input matrix is randomised by degree-preserving swaps in N2, N4 and N5 and by a column shuffle in N1.

\begingroup\small

Table: Randomised controls. All keep edge counts and Dale signs. {#tbl:nulls}

| code | operation | additionally preserved |
|------|------------------------|------------------------------------------------------|
| N1 | column shuffle | out-degree, per-source weight multiset |
| N2 | degree-preserving target swaps | in- and out-degree |
| N4 | N2 on interior edges only | as N2, plus every sensory-output and motor-input edge |
| N5 | N1 on interior edges only | as N1, plus every sensory-output and motor-input edge |

\endgroup

Three interventions move the shortcut directly. **AS** adds it to the connectome: pairs of edges, one from an olfactory receptor group to a non-motor target and one from a non-sensory source to a descending group, exchange targets, which keeps every degree and each weight with its source, until the olfactory shortcut share (the fraction of olfactory output weight landing directly on descending groups) reaches that of N2 (0.106). **N2R** runs the same operation in reverse on N2 until its share is 0. A **dose** series raises the connectome's share to at least 1, 3, 5 and 10 %, and a **sham** performs 101 degree-preserving swaps on edge pairs whose source is not an olfactory receptor and whose target is not a descending group, leaving the shortcut share at 0.0001. The sham does not touch the boundary edges that the dose rewires; §\ref{sec-dose} discusses what it therefore does and does not control.

## Calibration {#sec-calibration}

Every network is row-normalised and scaled by a global and a Kenyon-cell gain found by bisection so that mean group activity is 0.15 and mean Kenyon-cell activity 0.05 over a fixed stimulus set. This matches mean activity but not its spread: the connectome's median between-group activity SD is 0.163 with 46 silent groups of 512, against medians of 0.080 to 0.087 and 8 to 19 silent groups for the controls (N1 to N3 in the first grid, N4 and N5 in the corrected grid). Row normalisation also means that a boundary edge bit-identical before calibration can carry a slightly different share of weight after it; all path statistics below are computed on the calibrated matrices the runs used. A pre-registered re-run of the corrected grid under a calibration that also matches the spread of activity (between-group SD 0.118 and 3 to 7 silent groups in every condition) returned the same verdicts (Section \ref{sec-equivalence}, Appendix \ref{app-grid}).

## Evaluation and statistics {#sec-stats}

Every 50 generations the 16 agents with the highest training fitness per island are evaluated for 8 further lifetimes in a common garden: the random-number generator is reset before the batch so that all conditions and ablations face the same generated worlds. The registered **primary endpoint** is this fitness averaged over generations 0 to 550 (trapezoid rule); the registered **secondary endpoint** is the fitness at the last probe, generation 550. Olfaction, vision and plasticity are ablated in parallel copies of the batch.

Seeds are the unit of replication. Paired differences (connectome minus control) are tested with two-sided Wilcoxon signed-rank tests, Holm-corrected across controls, and summarised by medians with percentile bootstrap intervals over seeds (10,000 resamples). Runs of one seed share their controls across ecologies, so cross-ecology statements average each seed's ecologies first and test the ten seed means. For absence of a difference we require equivalence: the 90 % interval of the median paired difference must lie inside ±δ with δ = 0.10. The bound was set after the first grid, when the shortcut effect (at least 0.20) was known, and before any equivalence interval was computed; Appendix \ref{app-grid} reports how each verdict depends on δ.

## Pre-registration and a correction {#sec-prereg}

Each experiment below was run under a protocol frozen before its data (listed in Appendix \ref{app-prereg}), with engine hashes recorded. The first grid used a heading-noise term with a non-zero mean that no network could offset; it was found through a question raised in AI-assisted review of a draft and corrected, and every later experiment uses zero-mean noise. The first grid is reported because it motivated the rest; its generation-0 result is not used (Appendix \ref{app-prereg}).

# Results {#sec-results}

## Standard controls carry sensory-to-motor shortcuts and out-evolve the connectome {#sec-shortcuts}

In the connectome, exactly two of 14,989 edges run from an olfactory receptor group to a descending group, carrying 0.012 % of calibrated olfactory output; one descending group is one synapse from olfactory input, 15 are two, 36 are three and one is four. A column shuffle or degree-preserving swaps put 10.6 to 10.7 % of olfactory output directly onto descending groups and make every descending group reachable in one synapse (Table \ref{tbl:paths}). The boundary-preserving controls keep the connectome's structure; their small differences in the direct shares come from calibration. For vision, whose routes are already short, randomisation works the other way and lowers the direct share.

\begingroup\small

Table: Sensory-to-motor structure of the calibrated networks (medians over 10 seeds). "Direct": share of that sense's output weight landing on descending groups; "hops": median shortest path from that sense to a descending group. {#tbl:paths}

| condition | smell: direct | smell: hops | vision: direct | taste: direct |
|-----------------------|----------|----------|-----------|-----------|
| connectome | 0.00012 | 3 | 0.381 | 0.121 |
| N1 column shuffle | 0.107 | 1 | 0.108 | 0.113 |
| N2 degree-preserving | 0.106 | 1 | 0.128 | 0.093 |
| N4 interior swaps | 0.00010 | 3 | 0.340 | 0.125 |
| N5 interior shuffle | 0.00011 | 3 | 0.369 | 0.123 |

\endgroup

In the first pre-registered grid (4 ecologies × 10 seeds, with the uncorrected noise term), connectome populations were overtaken: at the last probe the column shuffle was ahead in all ten seeds on seed means (p_holm = 0.0059) and the degree-preserving control in eight (p_holm = 0.0195), while connectome populations depended far less on olfaction (Appendix \ref{app-first}).

## Against boundary-preserving controls the seed-mean difference lies within the equivalence bound {#sec-equivalence}

We re-ran the grid with the noise corrected and with N4 and N5 added (40 runs, 5 conditions each). On the registered primary endpoint no difference was detected against any control (seed-mean medians −0.03 for N1, −0.02 for N2, +0.03 for N4 and −0.01 for N5, none significant); the effects below emerge late in evolution and are measured on the registered secondary endpoint. At the last probe the connectome was behind the shortcut-carrying controls (p_holm = 0.029 for N1 and 0.016 for N2 on seed means), and against the boundary-preserving controls the seed-mean difference lay inside the ±0.10 bound: +0.002 for N4 and −0.074 for N5 (Table \ref{tbl:grid}). The N4 verdict holds for any bound down to 0.05; the N5 verdict holds at 0.10 but not at 0.075 and is marginal.

\begingroup\footnotesize

Table: Corrected grid, last probe (generation 550), connectome minus control: median over 10 seeds, with the 90 % interval below. Equivalence is assessed only on seed means across ecologies. {#tbl:grid}

| control | seed means | P0T0 | P1T0 | P0T1 | P1T1 |
|---------------|--------------|--------------|--------------|--------------|--------------|
| N1 (shortcut) | −0.224 | −0.621 | −0.312 | −0.159 | +0.054 |
|  | [−0.34, −0.10] | [−0.66, −0.47] | [−0.44, −0.22] | [−0.21, +0.01] | [−0.05, +0.13] |
| N2 (shortcut) | −0.199 | −0.403 | −0.338 | −0.017 | +0.015 |
|  | [−0.32, −0.13] | [−0.42, −0.31] | [−0.53, −0.07] | [−0.33, +0.01] | [−0.13, +0.14] |
| N4 (boundary) | +0.002 | −0.064 | +0.021 | −0.031 | +0.120 |
|  | [−0.026, +0.045] | [−0.26, +0.00] | [−0.03, +0.04] | [−0.12, +0.09] | [−0.06, +0.28] |
| N5 (boundary) | −0.074 | −0.264 | +0.067 | −0.015 | +0.070 |
|  | [−0.096, +0.002] | [−0.37, −0.18] | [−0.07, +0.09] | [−0.10, +0.06] | [−0.01, +0.21] |

\endgroup

The seed-mean equivalence does not hold ecology by ecology, and we do not claim it there. In the safe-food ecology (P0T0) the interior column shuffle N5 is ahead of the connectome by 0.26 (p_holm = 0.029), more than twice the bound; in the predator-plus-toxin ecology the connectome is ahead of N4 by 0.12 with an interval spanning zero. We then ran seeds 10 to 19 in the two T0 ecologies under rules fixed before those runs (Appendix \ref{app-grid}). N5's lead in P0T0 did not replicate: in the new seeds A − N5 was +0.031 (90 % CI −0.068 to +0.152), and over all 20 seeds −0.104 (−0.208 to +0.016, p_holm = 0.27). At n = 20 the per-ecology intervals for N4 and N5 lie inside ±0.10 in P1T0 (+0.003 and +0.056) but not in P0T0 (N4 −0.011, −0.114 to +0.080). N1 stays ahead in both ecologies (−0.491 and −0.333, p_holm ≤ 0.004); N2 stays ahead in P1T0 (−0.203, p_holm = 0.046), but its P0T0 lead shrank to −0.214 (−0.375 to +0.042, p_holm = 0.25), and in seeds 10 to 19 alone the connectome was ahead of N2 there (+0.100). The larger sample therefore weakens both the one per-ecology lead of a boundary-preserving null and part of the standard nulls' per-ecology advantage; the seed-mean contrast over four ecologies, which the extension did not cover, is unchanged.

The sensory strategy follows the same split. Blocking olfaction costs the connectome 0.03 against 0.68 for N1 and 0.64 for N2 (p_holm = 0.008), and 0.10 for N4 and 0.17 for N5; the connectome's differences from N4 and N5 (−0.07 and −0.11) are five- to ninefold smaller than those from N1 and N2 and not significant after correction (p_holm = 0.13 for both).

Because global calibration leaves the connectome with a wider spread of activity than the controls, we re-ran the corrected grid (40 runs) under a distribution-matched calibration that equalises mean, spread and silent fraction across conditions, with verdict rules fixed beforehand. Every run passed the manipulation check (control SD within 0.01 and silent-group count within 10 of the connectome's). The verdicts held: on seed means the boundary-preserving controls stayed inside ±0.10 (A − N4 +0.055, 90 % CI +0.035 to +0.097; A − N5 −0.009, −0.076 to +0.046) and the standard controls stayed ahead (A − N1 −0.274, −0.306 to −0.123; A − N2 −0.146, −0.231 to −0.037). The calibration shifted the N2 contrast by +0.101 (+0.016 to +0.133) and left the other three shifts undetermined, and N5's lead in P0T0 did not appear under it either (−0.066, −0.320 to +0.124; Appendix \ref{app-grid}).

## Moving the shortcut in and out {#sec-swap}

The pre-registered intervention (20 runs in the two ecologies with the clearest overtaking) rewires the shortcut in both directions while keeping degrees and Dale signs (Figure \ref{fig:swap}). Adding shortcuts to the connectome raised its fitness at the last probe in every seed, by +0.44 on seed means (10 of 10 seeds, p = 0.002; +0.445 in P0T0 and +0.475 in P1T0), and its time-averaged fitness by +0.26 (10 of 10, p = 0.002). The effect grew over evolution: +0.04 at generation 0 (p = 0.32; the interval is wider than the equivalence bound, so an initial effect is not excluded) and +0.40 from generation 0 to the last probe (10 of 10 seeds, p = 0.002). The intervention also transferred the sensory strategy, raising the cost of blocking olfaction from 0.15 to 0.99 (p = 0.002), close to the degree-preserving control's 0.84.

Removing the shortcut from the degree-preserving control cost it 0.43 fitness in the safe-food ecology (p = 0.004) and an estimated 0.00 in the predator ecology (difference between ecologies −0.36, p = 0.027); in neither ecology could the shortcut-free control be distinguished from the connectome, and its olfactory dependence fell from 0.84 to 0.58 (p = 0.027). Both interventions also rewire edges that are not shortcuts, so this experiment shows that *shortcut-creating rewiring* is sufficient for the advantage, not that the shortcut edges alone are.

![Intervention. Left and middle: fitness trajectories in the two ecologies (median and interquartile range over 10 seeds). Right: differences from the connectome at the last probe, seed means over both ecologies, with 90 % intervals. Adding shortcut-creating rewiring (AS) reproduces the control's advantage; removing it from the degree-preserving control (N2R) removes the advantage where it was clear.](figs/fig8_swap.png){#fig:swap width=100%}

## Graded doses and an interior-only sham {#sec-dose}

To ask whether the effect follows the amount of shortcut, a pre-registered dose series raised the connectome's olfactory-to-descending share to at least 1, 3, 5 and 10 % (2 ecologies × 10 seeds). Fitness rose across doses, on seed means 2.37, 2.54, 2.69, 2.75 and 2.93 from the connectome to the 10 % dose (Page's L, z = 4.68, p = 1.4 × 10⁻⁶; per-seed Spearman correlation positive in 9 of 10 seeds), in both ecologies (p = 1.5 × 10⁻⁵ and 0.004), and dependence on olfaction rose with it (p = 3.4 × 10⁻⁵ and 5.6 × 10⁻⁷; Figure \ref{fig:dose}). At generation 0 no trend was detected across doses (p = 0.98).

The full dose took 68 to 129 swaps depending on the seed (mean 98.5), comparable to the sham's 101, and was ahead of it by 0.534 on seed means (90 % CI +0.373 to +0.693, 10 of 10 seeds; +0.525 in the safe-food and +0.543 in the predator ecology). It was ahead in all seven seeds in which it used fewer swaps than the sham. The number of degree-preserving swaps therefore does not account for the effect. What the sham does not control is *where* the swaps land: the dose necessarily rewires olfactory-output and descending-input edges, while the sham avoids them, so this comparison separates shortcut-creating rewiring from interior-only rewiring, not the shortcut edges from other changes at the boundary. Section \ref{sec-bsham} closes that gap with a sham that rewires the same boundary edges. In this dose series shams were run only at the full swap count; the boundary-targeted experiment pairs a sham with every dose.

The sham itself was registered to fall inside the equivalence bound against the connectome and did not: +0.034 on seed means (90 % CI −0.062 to +0.123), ahead by 0.144 in the safe-food ecology (interval +0.072 to +0.222) and behind by 0.077 in the predator ecology (interval spanning zero). Interior-only rewiring is therefore not demonstrably inert, but its effect is small next to the dose's and does not keep its sign across ecologies (Appendix \ref{app-dose}).

![Graded shortcut transplant. Left and middle: fitness at the last probe against the share of olfactory output landing directly on descending groups, per ecology; grey lines are seeds, purple the mean with standard error, and the green square the interior-only sham (101 swaps) plotted at the connectome's own share. Right: fitness lost when olfaction is blocked, against the same axis.](figs/fig9_dose.png){#fig:dose width=100%}

## A boundary-targeted sham separates the shortcut from boundary rewiring {#sec-bsham}

A shortcut swap turns ORN→X and Y→DN into ORN→DN and Y→X. The boundary-targeted sham BS*k* rewires the same two boundary edges without creating a shortcut: it adds an interior edge Z→W and rotates the targets, giving ORN→W, Y→X and Z→DN. The olfactory weight lands in the interior instead of on a descending group, and the descending group receives Z's weight instead. Each seed performs as many rotations as its own AS*k* performed swaps. Z→W is the one of 64 random interior edges whose weight is closest in magnitude to ORN→X (median ratio 0.92 to 1.05 across constructions). Degrees, signs and the olfactory-to-descending share (0.0001) are unchanged, which the engine checks. Nine islands (A; AS1, AS3, AS5, AS10; and their shams BS1, BS3, BS5, BS10) evolved together in the two ecologies over ten seeds, under a protocol frozen before the runs (`PROTOCOL_BSHAM.md`). The primary contrast is AS10 − BS10 at the last probe, with a 90 % bootstrap interval of the seed mean: *shortcut-specific* if the lower bound exceeds +0.10, *boundary rewiring suffices* if the interval lies inside ±0.10.

The verdict was shortcut-specific. The full dose was ahead of its boundary-targeted sham by 0.600 (90 % CI +0.414 to +0.797), and by 0.568 and 0.631 in the two ecologies. The registered robustness check re-evaluated the generation-599 elites in four newly generated worlds and gave the same category (+0.639, +0.490 to +0.792). The sham itself was indistinguishable from the connectome (BS10 − A = +0.007, −0.070 to +0.086, inside ±0.10). The gap between each dose and its matched sham grew with dose (0.16, 0.28, 0.36 and 0.60 for 1, 3, 5 and 10 %; Page's L, z = 3.18, p = 7.4 × 10⁻⁴). Blocking olfaction cost the full dose 1.18 and its sham 0.35. The transplant's gain over the connectome replicated (AS10 − A = +0.606, +0.435 to +0.798). Rewiring the boundary edges the same number of times, without creating direct olfactory-to-descending connections, therefore reproduced none of the effect: it follows the shortcut edges themselves. Three runs failed with GPU out-of-memory (one of them twice), one completed attempt was deleted by a driver error and one attempt was stopped to test a memory fix; every affected run was re-run from scratch with identical flags before any outcome was inspected. Sixteen of the 20 runs used an engine copy that differs only in forcing evaluation earlier to save memory and was bit-identical to the original on the CPU (Appendix \ref{app-dose}).

## Further checks {#sec-checks}

Zeroing interior weights of evolved populations cost 18 to 23 % of fitness and permuting them 46 to 55 %, in every condition, with no detected difference between the connectome and any control (all p ≥ 0.16). The evolved agents are therefore sensitive to their interior wiring rather than running on the boundary alone; this is perturbation sensitivity, not evidence that the task requires interior computation. In the first grid, the fraction of one-step mutants beating their parent's clones was similar for connectome and control genomes (9.6 % against 9.4 and 9.8 % at generation 0), and connectome genomes lost more to large mutations (Appendix \ref{app-first}). We detected no odor-specific associative learning in any condition, and over evolution weights moved by 0.5 to 3 % of their initial spread while gains, leaks and read-outs moved by about 6 to 11 initial standard deviations (`results/evolution_analysis.txt`, the three calibration schemes of paper-v1) (Appendix \ref{app-checks}).

# Discussion {#sec-discussion}

In this model the connectome's standing depended on how the randomised control treated connections between sensors and motor neurons. The two standard nulls randomised the sensory-motor boundary along with everything else, so the comparison they support is about the boundary as much as about the interior wiring we meant to test. With the boundary held fixed, the end-of-run seed-mean difference fell within ±0.10; moving the shortcut into the connectome reproduced the controls' advantage and olfactory strategy; and the advantage followed the dose. A sham that rewired the same boundary edges the same number of times without creating shortcuts left fitness where the connectome was, so the effect is attributable to the direct olfactory-to-descending connections rather than to rewiring at the boundary as such.

Which properties a null should preserve depends on the hypothesis. A randomisation meant to test whether interior wiring matters should not change which sensors reach which effectors in one step; one meant to test the whole architecture may. Keeping every sensory-output and motor-input edge, as N4 and N5 do, is one way to hold the boundary fixed, and reporting calibrated sensory-to-motor path statistics next to degree statistics makes the choice visible. The same concern applies wherever a task rewards short reactive routes; whether it matters in static benchmarks depends on their tasks, and the statistics are cheap to report either way. Together with Dhiman [2026], who found a verdict change through initialisation and degree effects, this suggests that which null was used, and what it preserved, matters more for connectome comparisons than has been usual.

What we can say about the connectome itself is narrow. On seed means over four ecologies we found no end-of-run advantage for measured interior wiring over boundary-preserving controls, within ±0.10; that is a bounded null for one estimand, not proof of no effect, and the one per-ecology lead of the interior column shuffle did not replicate in ten further seeds. Per-ecology estimates are also less stable than the seed means: the degree-preserving control's lead in the safe-food ecology, which motivated the intervention there, shrank to an interval spanning zero at twenty seeds. Our agents are a compressed rate model in a two-dimensional world, no condition evolved detectable associative learning, and evolution adapted gains and read-outs far more than weights. A task that requires what fly circuits are built for, such as learning from experience within a lifetime, could give a different answer, and is the natural next test.

**Limitations.** One connectome, one compression and one family of two-dimensional worlds; one mutation operator, which tunes far more than it rewires; a toxin factor confounded with food abundance and a stationary rather than absent predator; global calibration in the main experiments (verdicts unchanged under distribution matching); a common-garden evaluation that reuses one generated world stream per probe; population samples rather than elites in the post hoc assays; equivalence bounded at ±0.10 on seed means with ten seeds and only where stated (twenty in two ecologies), with the bound set after the first grid; interventions that rewire more than the shortcut edges (the boundary-targeted sham controls this at the boundary, not everywhere), an interior-only sham that failed its own registered prediction; and boundary edges that mutation can add during evolution.

# Data and code availability {#sec-data .unnumbered}

Engines, analysis and figure scripts, frozen protocols, the dated research log and run outputs are at https://github.com/gyujeongion/flyconnectome-nulls; an archived copy with the manuscripts is at https://doi.org/10.5281/zenodo.22871090. Raw per-generation logs and stored genomes (about 2 GB) are not included; the analysis outputs computed from them are. Re-running a seed does not repeat a run bit for bit: after every common-garden probe, which runs under a fixed random seed, the engines reseed the evolutionary random stream from the clock, and the GPU scatter used in mutation does not fix which write survives when two mutations hit the same weight. A re-run reproduces the null constructions, calibration and generation 0 exactly and later generations only in distribution; seeds are therefore independent replicates, not reproducible trajectories. FlyWire v783 connectivity and annotations are public. Every run used one Apple M1 Ultra with MLX (about 36 minutes per 5-condition, 600-generation run).

# Use of AI tools {#sec-ai .unnumbered}

The engines, analysis and figure scripts were written, and this manuscript drafted, with an AI coding assistant (Claude Code, Anthropic); two further language models critiqued drafts. Every number is regenerated from stored run outputs by scripts in the repository, and the author directed the study and is responsible for the text (Appendix \ref{app-ai}).

# Funding and conflicts of interest {#sec-funding .unnumbered}

This research received no external funding. The author declares no conflicts of interest.
