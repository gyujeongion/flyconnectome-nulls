# Supplementary material

**A compressed fly connectome starts ahead but is overtaken by randomised wiring under evolution**, Gyujeong Park

All results in this supplement come from paper-v1 (pre-registered in `PROTOCOL.md`: 4 conditions × 10 seeds × 300 generations, one ecology with toxic probability 0.5, toxic-bite loss 0.4 and a predator pursuing at 0.1·v_max; common garden every 25 generations). Sections S2 and S3 contain the pre-registered paper-v1 endpoints; everything else is post hoc or exploratory and labelled as such. Stored genomes in paper-v1 are population samples, not elites, for the reason given in main-text §2.5; the sample text below uses the older word "elites" for them.

## S1. Calibration schemes

* **global**: W row-normalised to unit absolute input; one global gain and one Kenyon-cell gain found by bisection so mean group activity = 0.15 and mean Kenyon-cell activity = 0.05 over a fixed stimulus ensemble (16 odors × 2 concentrations, single-channel stimuli, one mixed stimulus). Exact convergence for all conditions (relative error 0).
* **pergroup**: per-group multiplicative gain feedback, 150 iterations (the exploratory-phase scheme). Residuals differ by condition (connectome 0.23, nulls 0.02–0.07).
* **rawglobal**: global gain without structural normalisation; circuit analysis only.
* **distmatch**: structural normalisation plus per-group quantile matching to a fixed log-normal reference profile (mean 0.15, SD/mean 0.8), equalising mean, SD and silent fraction. Achieved: mean 0.148–0.149, SD 0.118, sorted-profile error 0.0013 in all four conditions.

### S2 Apparent connectome advantages are calibration-dependent (exploratory)

Circuit properties of initial brains, 4 conditions × 10 seeds per scheme:

* Kenyon-cell odor decorrelation favours the connectome under **pergroup** (10/10 seeds, p_holm = 0.006) and favours the nulls under **global** and **rawglobal** (0–3/10 seeds).
* Dopaminergic responses to sugar and malaise favour the connectome under **global** (10/10, p_holm = 0.006), show no consistent difference under **pergroup**, and favour the nulls under **rawglobal**.
* Mushroom-body learning specificity shows no significant difference in any scheme.
* Looming-escape turn sign is correct in the connectome in 10/10 seeds and beats all three nulls under **global** and **pergroup** (p_holm = 0.006). Under **rawglobal** it beats only N3 (p_holm = 0.006; vs N1 p = 0.275, vs N2 p = 0.129). Laterality is thus the most calibration-robust property, but not significant against every null in every scheme.

An earlier internal result of a 43-fold learning-specificity advantage was traced to a numerically sensitive calibration and is **retracted**; under the robust schemes the ratio is 1.6–3.0 with no significant difference between conditions.

### S3 An innate head start and a different evolved sensory strategy

Pre-registered under **global** and **pergroup**; repeated post hoc under **distmatch**.

| endpoint | global | pergroup | distmatch (post hoc) |
|---|---|---|---|
| generation 0 | A > N1, N2, N3 (10/10, 10/10, 10/10), p_holm 0.006 | A > all (9/10, 10/10, 10/10), p_holm 0.006–0.008 | A > all (10/10, 9/10, 10/10), p_holm 0.006 |
| olfaction ablation | A −0.03 vs +0.31 / +0.34 / +0.30; p_holm 0.006–0.014 | A −0.08 vs +0.16 / +0.18 / +0.22; p_holm 0.006–0.008 | A −0.10 vs +0.25 / +0.32 / +0.28; p_holm 0.006–0.008 |
| vision ablation | A +0.30 vs +0.13 / +0.14 / +0.07; p_holm 0.06–0.11 | A +0.33 vs +0.09 / +0.13 / +0.17; p_holm 0.041 | A +0.47 vs +0.17 / +0.13 / +0.05; p_holm 0.012 |
| generation 275 | no difference vs N1, N2 (p_holm 1.0); vs N3 p_holm 0.059 | no difference vs N1, N2; vs N3 p_holm 0.111 | no difference vs N1, N2; vs N3 p_holm 0.252 |
| AUC (primary) | A > N1, N3 (p_holm 0.006); N2 n.s. | A > N3 only (p_holm 0.012) | A > N3 only (p_holm 0.012) |

Median generation-0 differences are +0.16 to +0.37 fitness units. The head start and the sensory-strategy divergence replicate in all three calibrations, including the one that matches the full activity distribution, so they are not explained by marginal activity statistics. The cumulative endpoint is not robust: the degree-preserving null N2 is never significantly worse on AUC or final fitness.

### S4 No evidence that the nulls' olfactory reliance is associative learning (post hoc)

On generation-250 elites the reversal × plasticity interaction is not significant in any condition or calibration (p_holm 0.39–1.0), and extra toxic bites after an odor reversal do not differ. A plasticity benefit of ≈ 0.09–0.30 exists in all conditions but is not odor-specific. We therefore find **no evidence** that the nulls' olfaction dependence is mediated by odor-specific associative learning; with 10 seeds and the observed confidence intervals a smaller interaction cannot be excluded.

Transfer: with a wider odor plume, connectome elites outperform N1 and N3 under **global** (p_holm = 0.029); under **distmatch** no transfer difference reaches significance, and higher toxicity shows none in any scheme.

### S5 The head start is amplified by, but not reducible to, the engineered reflex (post hoc)

Removing the innate feeding read-out reduces the generation-0 advantage from +0.15…+0.31 to +0.04…+0.08, still significant in all schemes (7–10/10 seeds, p 0.006–0.029). Blinding alone leaves it near its original size (+0.13…+0.30). Removing both vision and the read-out leaves +0.03…+0.07 (8–10/10 seeds, p ≤ 0.029). The random-walk control, with all read-outs zeroed, gives a difference of exactly 0.000 (0/10 seeds) in every scheme, confirming the differences arise through the behavioural pathway and not through the evaluation itself.

### S6 Convergence with the degree-preserving null is not a ceiling artifact (post hoc)

Generation-250 elites re-evaluated in harder worlds without further evolution: at predator speed 0.18 the connectome still beats N1 and N3 (global +0.228 vs N1, 9/10, p = 0.010; distmatch +0.095 vs N1, 8/10, p = 0.027 and +0.342 vs N3, 8/10, p = 0.027) but is never significantly different from N2 at any difficulty or in any scheme. At the most extreme toxicity level the connectome advantage disappears in all schemes.

### S7 Negative result: neuron-level central complex

A module of 1,041 real central-complex neurons (EPG, PEG, PEN, Delta7, PFN, hΔB, FC2, PFL, plus GLNO/LNO/ER interface neurons with measured connectivity) reproduces heading tracking of a rotating visual cue (phase-tracking ratio 1.003; 93 % of EPG variance on a two-dimensional ring). Across ≈ 1,500 pathway-gain settings from random search and an evolution strategy we found none that sustains heading in darkness, integrates angular velocity, or supports decodable path integration, for the connectome and the shuffled control alike. We report this as a limitation of our rate-model formulation, not as a statement about the connectome.

## S8. Figures from paper-v1

![Figure S1](figs/fig1_curves.png)

**Figure S1.** Common-garden fitness over 300 generations in paper-v1, per calibration scheme (median and interquartile range over 10 seeds).

![Figure S2](figs/fig2_endpoints.png)

**Figure S2.** Paper-v1 endpoints per seed: generation 0, AUC (primary), generation 275, and ablation costs.

![Figure S3](figs/fig3_calibration_sensitivity.png)

**Figure S3.** Circuit metrics before evolution under three calibration schemes.

![Figure S4](figs/fig4_stress.png)

**Figure S4.** Generation-250 samples of paper-v1 tested in harder worlds without further evolution.

## S9. Additional analyses for the ecology grid (post hoc)

Full outputs: `results/eco_analysis.txt` (registered tests), `results/eco_cross.txt` (crossing, cross-ecology, H3), `results/eco_pooled.txt` (across-ecology summaries), `results/mut_analysis.txt` (mutational neighbourhood with confidence intervals), `results/revision_analysis.txt` (sustained crossing, fitness components, absolute ablation costs, population samples against top 16), `results/eco_calib.txt` (activity diagnostics per run).


## S10. Design changes and corrections, in order

The main text reports one correction made after pre-registration (§2.7). For completeness, this section lists every design change that reached the frozen protocols or the reported results, in the order they were made; the full dated log is `RESEARCH_LOG.md` in the repository.

| when | what | effect on reported results |
|---|---|---|
| exploratory phase, before `PROTOCOL.md` | 27 changes to brain construction, calibration and world code (for example: mandatory inclusion of the gustatory-to-dopaminergic path groups; replacement of clustered Kenyon cells by individually wired ones; removal of an energy loophole that let agents survive without eating; a retracted learning-specificity result traced to a numerically unstable calibration) | none of the frozen experiments were run before these changes |
| after `PROTOCOL_ECO.md` was frozen, before any run | the setting meant to remove the predator instead enabled an adaptive predator; fixed and re-frozen before the first grid | none |
| after the first grid | stored genomes were population samples rather than elites, by an indexing error | affects only the post hoc assays of §2.8, as stated there |
| after the first grid (AI-assisted review of draft 1) | one-sided heading noise, §2.7 | the first grid's generation-0 result is not used; all later experiments corrected |
| after the first grid | the standard nulls carry olfactory-to-motor shortcuts, §3.2 | motivated the corrected grid and the intervention |

## S11. Sensitivity analysis: the sham equivalence test on an unregistered outcome

`PROTOCOL_DOSE.md` registered the common-garden probe of §2.6 as the outcome, as `PROTOCOL_ECO.md` does for every experiment in this paper. On that outcome the registered equivalence prediction for the sham fails: sham minus connectome is +0.034 on seed means, 90 % CI −0.062 to +0.123, crossing the ±0.10 bound (§3.5).

The engines also log an in-run fitness at each generation, evaluated in the training worlds rather than in the common-garden batch, and the last of those is generation 599. That outcome was not registered. On it the same difference is −0.022, 90 % CI −0.090 to +0.047, which lies inside the bound and would have supported the prediction. We report this because it is recoverable from the released logs and because the direction of the disagreement is the one that flatters us; it does not govern the registered test, which is the probe.

The disagreement is confined to the sham. Under the unregistered outcome the graded series still rises monotonically (Page's L, z = 5.12, p = 1.5 × 10⁻⁷, per-seed Spearman positive in 10 of 10 seeds) and the generation-0 null still holds (p = 0.98). The swap-matched comparison of §3.5 is likewise unaffected in sign or approximate size.

We did not investigate why the two outcomes differ for the sham specifically. The probe resets its random stream and evaluates 16 agents per island for 8 lifetimes in a fixed sequence of worlds, while the in-run number is the population's fitness in the worlds it was selected in; a condition whose advantage is partly specific to the worlds it evolved in would show the two diverging, but we have not tested that here.
