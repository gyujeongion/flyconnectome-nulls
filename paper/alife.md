# Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison

**Gyujeong Park**, Independent Researcher. <ionlabs2025@gmail.com>

## Abstract

Whether a measured connectome outperforms a randomised copy of itself depends on which properties the randomisation preserves. We report a case in which the boundary between sensory and motor neurons decides the outcome of an embodied evolutionary comparison. Agents whose brain is a compressed adult *Drosophila* FlyWire v783 connectome evolved for 600 generations in four pre-registered ecologies alongside populations built on randomised wiring. Two randomisations in common use, a presynaptic-column shuffle and degree-preserving edge swaps, move 10.6 to 10.7 % of olfactory receptor output directly onto descending motor groups, against 0.012 % in the connectome, where 52 of 53 descending groups are two or more synapses away. Against these controls the connectome populations end behind at the last common-garden probe, generation 550 (−0.22 and −0.20 fitness units on seed means, p_holm = 0.029 and 0.016), while on the pre-registered primary endpoint no difference is detected. Against controls that preserve the boundary and randomise only the interior, the difference at that probe lies inside a ±0.10 equivalence bound on seed means (+0.002, 90 % CI −0.026 to +0.045; −0.074, 90 % CI −0.096 to +0.002), although one interior control remains ahead by 0.26 in one ecology. Adding olfactory-to-motor shortcuts to the connectome by degree-preserving swaps raises its fitness at that probe by +0.44 (10 of 10 seeds, p = 0.002) and its dependence on olfaction from 0.15 to 0.99; removing them from the control costs it 0.43 fitness where it had a clear advantage. In this model, what a null preserves at the sensory-motor boundary changes the comparison, and calibrated sensory-to-motor path statistics are worth reporting next to the degree statistics a null preserves.

## 1. Introduction

Whole-brain wiring diagrams of the adult fly are public (Dorkenwald et al., 2024; Schlegel et al., 2024) and models built on them reproduce feeding circuits (Shiu et al., 2024), visual responses (Lappalainen et al., 2024) and locomotion (Vaxenburg et al., 2025). Whether measured wiring is better than chance is usually tested against a randomised network, and the answer is a statement about the properties the randomisation did not preserve (Váša & Mišić, 2022). Degree-preserving rewiring is the usual default; for connectomes it misses structure that spatial constraints capture (Salova & Kovács, 2025), and Worrell et al. (2017) found the fly connectome spreads sensor-to-effector activity faster than such nulls. Dhiman (2026) found that switching from a naive to a degree-preserving null reversed a connectome-versus-null verdict in a fly visual model through initialisation and degree effects; Cho (2026) provides a 31-task behavioural benchmark with shuffled and degree-preserving controls. These compare networks that already exist. We asked what happens under selection, and found that in our setting the property that decides the comparison is how directly a sensor can reach a motor neuron, which the standard nulls randomise along with everything else.

## 2. Methods

**Brain.** From FlyWire v783 we drop optic-lobe, visual-centrifugal and photoreceptor classes and group the rest by cell type × hemisphere, keeping K = 512 groups (14,989 non-zero weights) chosen by a mandatory list of sensory, mushroom-body, looming and descending types plus a sensory-to-descending path score. Weights are signed synapse-count ratios with neurotransmitter-predicted signs (Eckstein et al., 2024). A separate layer of 1,000 Kenyon cells keeps its measured projection-neuron inputs, with APL-like global inhibition and dopamine-gated plasticity onto mushroom-body output neurons. Rates update as r ← (1 − α) r + α tanh(max(x, 0)) with x the sum of recurrent, Kenyon-cell, sensory and bias input. Read-outs for forward speed, turning and eating start from named descending neurons; the turn read-out subtracts a running mean of the left–right difference.

**World and evolution.** A unit square holds 24 odor-marked patches (nutritious, toxic, neutral), with odor-to-role assignment redrawn every lifetime, an energy budget that forces foraging, and a predator emitting CO₂ and looming cues. Ecologies cross a stationary harmless predator or a pursuing one with safe or toxic food. One island per wiring condition, 512 agents, tournament selection, 4 % elitism, 600 generations; mutation perturbs weights, adds and prunes edges, and tunes gains, leaks, learning rates and read-outs. Every 50 generations the top 16 agents per island are re-evaluated for 8 lifetimes each in a common garden with a fixed generator seed, intact and under four ablations.

**Controls.** N1 permutes each presynaptic column's targets; N2 applies Maslov–Sneppen degree-preserving swaps (Maslov & Sneppen, 2002); N3 permutes within super-class × hemisphere blocks. N4 and N5 apply the N2 and N1 operations only to interior edges, keeping every edge leaving a sensory group and every edge entering a descending group bit-identical to the connectome, checked at construction; this applies the standard practice of constraining a null to the property under test (Váša & Mišić, 2022) to the sensory-motor boundary. AS transplants shortcuts into the connectome by degree-preserving swaps until its olfactory shortcut share reaches 0.106; N2R removes them from N2 the same way.

**Statistics.** Seeds are the unit of replication; paired Wilcoxon tests with Holm correction across controls, percentile bootstrap intervals over seeds. Statements across ecologies average each seed's ecologies first (n = 10). For absence we require equivalence: the 90 % interval inside ±0.10 fitness units, half the smallest shortcut effect.

**A correction made after pre-registration.** The heading-noise term in the first grid was one-sided, adding 0.025 rad per step that no network could cancel. The error was found through a question raised during AI-assisted review of the first draft. Removing the drift at evaluation costs connectome populations 0.174 fitness while controls lose 0.00 to 0.04 (10 of 10 seeds, p_holm = 0.0059), so the connectome's generation-0 lead in the first grid is not used; the later overtaking was unaffected. All experiments after the first grid use zero-mean noise.

## 3. Results

**The first grid.** With N1, N2 and N3, before the noise correction, connectome populations started ahead in every ecology (12 of 12 comparisons, p_holm ≤ 0.027) and were overtaken by the last probe (N1 ahead in 10 of 10 seeds, p_holm = 0.0059; N2 in 8 of 10, p_holm = 0.0195). One explanation can be set aside: the fraction of one-step mutants beating their parent's clone distribution was the same for connectome and control genomes (9.6 % against 9.4 and 9.8 %, intervals inside ±2.5 percentage points), including later in evolution when the connectome populations were the less fit ones.

![Figure 1](figs/fig5_eco.png)

**Figure 1.** First grid, with the standard controls and before the noise correction: common-garden fitness over 600 generations in the four ecologies (median and interquartile range over 10 seeds). A, connectome; N1, column shuffle; N2, degree-preserving swaps; N3, block shuffle.

**The controls carry shortcuts.** In the compressed connectome, two of 14,989 edges run from an olfactory receptor group to a descending group, 0.012 % of calibrated olfactory output; one of 53 descending groups is one synapse away, 51 are two or three. After a column shuffle or degree-preserving swaps, 10.6 to 10.7 % of olfactory output lands directly on descending groups and every descending group is one synapse away. In the fly, olfactory receptor neurons synapse in the antennal lobe rather than onto descending neurons, and the compressed connectome keeps that structure apart from the two edges noted. The same randomisations reduce the direct visual-output share, from 0.381 in the connectome to 0.108 to 0.204.

**Table 1.** Sensory-to-motor structure of the calibrated networks (medians over 10 seeds). "Direct" is the share of that sense's output weight landing on descending groups; "hops" is the median shortest path from that sense over descending groups.

| condition | smell: direct | hops | vision: direct | taste: direct |
|---|---|---|---|---|
| connectome | 0.00012 | 3 | 0.381 | 0.121 |
| N1 column shuffle | 0.107 | 1 | 0.108 | 0.113 |
| N2 degree-preserving | 0.106 | 1 | 0.128 | 0.093 |
| N3 block shuffle | 0.065 | 1 | 0.204 | 0.067 |
| N4 interior swap | 0.00010 | 3 | 0.340 | 0.125 |
| N5 interior shuffle | 0.00011 | 3 | 0.369 | 0.123 |

The boundary edges of N4 and N5 are bit-identical to the connectome before calibration; the residual differences in the table come from row normalisation, which rescales a preserved edge when interior degrees change.

![Figure 2](figs/fig7_fix.png)

**Figure 2.** Corrected grid: connectome (A), shortcut-carrying controls (N1, N2) and boundary-preserving controls (N4, N5), median and interquartile range over 10 seeds.

**Boundary-preserving controls.** Re-running the grid with the noise term corrected and N4, N5 in place of N3 (40 runs, 5 conditions, 600 generations) gives seed-mean differences at the last probe of −0.224 [−0.34, −0.10] against N1 and −0.199 [−0.32, −0.13] against N2, both significant, but +0.002 [−0.026, +0.045] against N4 and −0.074 [−0.096, +0.002] against N5, both inside the equivalence bound. Per ecology the estimates vary and we do not claim equivalence there; in the safe-food ecology N5 remains ahead by 0.26 (p_holm = 0.029). The registered primary endpoint, time-averaged fitness, shows no detected difference against any control. The olfactory-strategy gap shrinks from 0.6 fitness units against N1 and N2 (p_holm = 0.008) to −0.07 and −0.11 against N4 and N5, neither surviving correction. With the noise term corrected and the boundary matched, the connectome no longer starts ahead either: the generation-0 seed-mean difference against N4 is −0.18 (p_holm = 0.15).

**Table 2.** Corrected grid, last probe (generation 550), connectome minus control: median over 10 seeds [90 % CI]. Equivalence (±0.10) is assessed only on seed means across ecologies.

| comparison | seed means | P0T0 | P1T0 |
|---|---|---|---|
| vs N1 | −0.224 [−0.34, −0.10] | −0.621 | −0.312 |
| vs N2 | −0.199 [−0.32, −0.13] | −0.403 | −0.338 |
| vs N4 | +0.002 [−0.03, +0.05] equivalent | −0.064 | +0.021 |
| vs N5 | −0.074 [−0.10, +0.00] equivalent | −0.264 | +0.067 |

**The intervention.** Adding shortcuts to the connectome raised fitness at the last probe by +0.44 on seed means (10 of 10 seeds, p = 0.002) and time-averaged fitness by +0.26; the effect grows over evolution (+0.40 from generation 0 to the last probe, 10 of 10 seeds), and at generation 0 no difference is detected (+0.04, p = 0.32, interval wider than the equivalence bound); its olfactory dependence rose from 0.15 to 0.99, matching the shortcut-carrying control. Removing shortcuts from N2 cost 0.43 fitness in the safe-food ecology (p = 0.004) and an estimated 0.00 in the predator ecology, after which no fitness difference from the connectome was detected. Because both interventions rewire more than the shortcut edges, path shortening is the proposed mechanism rather than an isolated variable.

![Figure 3](figs/fig8_swap.png)

**Matched rewiring.** A graded transplant raising the shortcut share to at least 1, 3, 5 and 10 % (reached 1.0 to 1.7, 3.0 to 4.6, 5.0 to 5.8 and 10.0 to 10.6 % across seeds, with 3 to 27, 11 to 48, 25 to 58 and 68 to 129 degree-preserving swaps) raises fitness in step (Page's L, z = 4.68, p = 1.4 × 10⁻⁶; per-seed Spearman positive in 9 of 10) and raises the cost of blocking olfaction with it (0.28 to 1.26 and 0.08 to 0.95 in the two ecologies), while the doses are indistinguishable at generation 0 (p = 0.98). A sham that performs 101 swaps on edges that are neither olfactory in nor descending out leaves the shortcut share at 0.0001. Rewiring comparable amounts, the full dose is ahead of the sham by 0.534 (90 % CI +0.377 to +0.700, 10 of 10 seeds), including all seven seeds in which it used fewer swaps, so the amount of rewiring does not account for the effect. The sham was registered to fall inside the ±0.10 bound against the connectome and does not (+0.034, 90 % CI −0.062 to +0.123), changing sign between ecologies, so rewiring in general also contributes something smaller and less stable.

**Figure 3.** Intervention. Left and middle: trajectories with shortcuts added to the connectome (AS) and removed from the degree-preserving control (N2R). Right: differences from the connectome at the last probe with 90 % intervals.

**The interior is used.** Generation-500 populations re-evaluated with interior weights zeroed lose 18 to 23 % of fitness and with interior weights permuted 46 to 55 %, equally in every condition, so the evolved agents are not running on the boundary alone.

**What survives from the first experiment.** The connectome's generation-0 advantage does not: it rested on the biased noise together with the boundary difference. The overtaking by N1 and N2 does, and the intervention accounts for it. Randomised networks stayed structurally distinct throughout: after 250 generations their support overlap with the connectome had moved by at most 0.003, the ipsilateral share of weight was 0.77 against 0.50, and reciprocal edges 0.43 against 0.05 to 0.12, while weights moved by 1 to 3 % of their initial magnitude and gains, leaks and read-outs by about 10 initial standard deviations.

## 4. Discussion

In this model, what the null preserved at the sensory-motor boundary decided the comparison. Which properties a null should preserve depends on the hypothesis: a randomisation meant to test interior wiring should not change which sensors reach which effectors in one step, and reporting calibrated sensory-to-motor path statistics next to degree statistics makes that choice visible. Dhiman (2026) reached a related conclusion in a fly visual model through a different mechanism, initialisation and degree effects rather than path structure. On our estimand we found no fitness advantage for measured interior wiring against boundary-preserving controls, bounded at ±0.10 fitness units, which is a bounded null rather than proof of no effect. The agents are a compressed rate model in a two-dimensional world, we detected no associative learning in any condition, and the mutation operator tunes far more than it rewires, so a task closer to fly ecology or an operator that rewires aggressively could give a different answer. Full methods, five frozen protocols, the research log and further results are in the preprint and the project repository at https://github.com/gyujeongion/flyconnectome-nulls.

## Acknowledgements

Engineering and analysis were done with AI assistance (Claude Code); drafts and results were also reviewed by two further language models, as AI-assisted review rather than peer review, and several of the controls reported here were added in response.

## References

Cho, B. (2026) flybench, v0.1.0 (software). https://github.com/brandoncho369/flybench
Dhiman, N. (2026) Topological sensitivity in connectome-constrained neural networks. arXiv:2604.04033.
Dorkenwald, S. et al. (2024) *Nature* 634, 124–138.
Eckstein, N. et al. (2024) *Cell* 187, 2574–2594.
Lappalainen, J. K. et al. (2024) *Nature* 634, 1132–1140.
Maslov, S. & Sneppen, K. (2002) *Science* 296, 910–913.
Salova, A. & Kovács, I. A. (2025) *Network Neuroscience* 9, 181–206.
Schlegel, P. et al. (2024) *Nature* 634, 139–152.
Shiu, P. K. et al. (2024) *Nature* 634, 210–219.
Váša, F. & Mišić, B. (2022) *Nature Reviews Neuroscience* 23, 493–504.
Vaxenburg, R. et al. (2025) *Nature* 643, 1312–1320.
Worrell, J. C. et al. (2017) *Network Neuroscience* 1, 415–430.
