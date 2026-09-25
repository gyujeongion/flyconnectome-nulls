# Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison

Code, protocols and data for the preprint of the same name. Agents whose brain is a compressed adult *Drosophila* FlyWire v783 connectome forage and evolve in a 2D world, next to identically evolved populations built on randomised wiring. What the null model preserves at the sensory-motor boundary changes the outcome: two standard randomisations create direct olfactory-to-motor connections that the compressed connectome lacks, and boundary-preserving nulls remove the difference. Rewiring comparable amounts and varying only where the swaps land leaves most of the effect with those connections.

**Paper:** [`paper/neurips/paper_neurips.pdf`](paper/neurips/paper_neurips.pdf) (v4.2, NeurIPS-style preprint with appendix) · [`paper/biorxiv.pdf`](paper/biorxiv.pdf) (v4.2 long form) · [`paper/supplementary.pdf`](paper/supplementary.pdf). The 7-page conference version was not updated past v4.1 and is kept only in the v4.1 Zenodo record.

## What the experiments show

| | result |
|---|---|
| Standard controls (column shuffle, degree-preserving swaps) | overtake the connectome by the last probe (10 of 10 seeds, p_holm = 0.0059) |
| Olfactory output landing directly on descending motor groups | connectome 0.012 %, standard controls 10.6–10.7 % |
| Boundary-preserving controls (randomise interior only) | seed-mean difference inside a ±0.10 equivalence bound: +0.002 (90 % CI −0.026 to +0.045) and −0.074 (90 % CI −0.096 to +0.002); one interior control still ahead by 0.26 in one ecology |
| Shortcuts transplanted into the connectome | +0.44 fitness at the last probe (10 of 10 seeds, p = 0.002); olfactory dependence 0.15 → 0.99 |
| Shortcuts removed from the control | −0.43 in the ecology where it led (p = 0.004); no detected difference from the connectome afterwards |
| Graded shortcut doses (1, 3, 5, 10 %) | fitness rises in step (p = 1.4 × 10⁻⁶, 9 of 10 seeds) and so does dependence on smell (0.28 → 1.26) |
| Sham rewiring about as much as the full dose (101 swaps against 68 to 129) | full dose ahead of sham by +0.534 (90 % CI +0.373 to +0.693, 10 of 10 seeds), including all 7 seeds where it rewired less; the sham's registered equivalence prediction **failed** (+0.034, CI −0.062 to +0.123) |
| Interior wiring scrambled in evolved agents | 46–55 % of fitness lost in every condition |
| Boundary-targeted sham (same boundary edges rewired the same number of times, no shortcut created) | full dose ahead of it by +0.600 (90 % CI +0.414 to +0.797; +0.639 in fresh worlds); the sham itself matched the connectome (+0.007) |
| Corrected grid re-run under distribution-matched calibration (D-06) | verdicts unchanged: boundary-preserving controls inside ±0.10 (+0.055, −0.009), standard controls ahead (−0.274, −0.146); manipulation check 40/40 |
| Two ecologies extended to 20 seeds | N5's lead in the safe-food ecology did not replicate (seeds 10–19: +0.031, CI −0.068 to +0.152); N2's lead there shrank to an interval spanning zero |

Six pre-registration protocols were frozen before their runs, with code and data hashes: [`protocols/`](protocols/). A seventh, `PROTOCOL_DIST.md`, was frozen for an activity-distribution control that we did not run; it is kept there rather than deleted. Four more were frozen for follow-up experiments not reported in the paper: `PROTOCOL_LEARN.md` and `PROTOCOL_ETA.md` (an ecology that rewards learning; its pilots were run by mistake in a world with no toxic patches, so they did not test the hypothesis and their results are withdrawn, see the amendments), `PROTOCOL_TRANSFER.md` (out-of-distribution transfer; its registered baseline turned out to be already out of distribution, see its amendment) and `PROTOCOL_SHAM2.md` (shams paired with each intermediate dose, frozen before its runs). Two more were frozen after v4.1 and are reported in v4.2: `PROTOCOL_DIST2.md` (the corrected grid re-run with activity-distribution-matched calibration; `analyze_dist2.py`) and `PROTOCOL_BSHAM.md` (a boundary-targeted sham that rewires the same boundary edges as the shortcut transplant without creating shortcuts; `evo_bsham.py`, `analyze_bsham.py`). The n = 20 extension of two ecologies follows rules fixed in advance (`analyze_n20.py`). One registered prediction failed and is reported as failed: a control rewiring about as much as the full shortcut transplant was predicted to show no effect against the connectome, and instead crossed the equivalence bound (see §3.5 of the paper and `SUPPLEMENTARY.md` §S11). One correction was made after pre-registration (a one-sided turning-noise term in the motor model); it is described in the paper (§2.7) and in [`RESEARCH_LOG.md`](RESEARCH_LOG.md), with the re-runs that replaced the affected results. The boundary-preserving nulls apply the standard practice of constraining a null to the property under test (Váša & Mišić, 2022) to the sensory-motor boundary.

## Layout

```
evo_paper.py      evolution engine used for the first pre-registered grid
evo_fix.py        same engine with the turning-noise term corrected, plus the
                  boundary-preserving controls N4 and N5
evo_swap.py       adds the shortcut-transplant conditions AS and N2R
evo_dose.py       adds graded doses AS1/AS3/AS5/AS10 and the SHAM rewiring control
evo_bsham.py      adds the boundary-targeted shams BS1/BS3/BS5/BS10, true elites and a
                  fresh-world evaluation; evo_bsham_mem.py is the same engine with
                  evaluation forced earlier to fit 32 GB (bit-identical on the CPU)
evo_learn.py      follow-up engine: shams paired per dose, ecology search over odour
                  reversal (pick_eco.py, gate_pilot.py apply the frozen rules)
evo_assay.py      evo_mut.py  evo_noise.py  evo_int.py   read-only assay runners
build_brain.py    FlyWire v783 -> data/brain.npz (512 groups + 1,000 Kenyon cells)
analyze_*.py      the analyses reported in the paper, one file per experiment
make_figs*.py     figures
SUPPLEMENTARY.md  supplementary results (S11 sham sensitivity, S13 n = 20, S14 boundary
                  sham, S15 distribution-matched calibration)
paper/neurips/    NeurIPS-style source: main_body.md + assemble.py -> paper.md -> build.sh
make_paper_numbers.py  regenerates values quoted in the paper into results/paper_numbers.json
run_*.sh          batch drivers for a single machine
protocols/        frozen protocols: eight reported in the paper, one frozen and not run,
                  four for follow-up experiments
results/          analysis outputs quoted in the paper
figs/             figures
```

## Reproducing

Requires Python 3.11+ on Apple Silicon (the engines use [MLX](https://github.com/ml-explore/mlx)) and roughly 55 GPU-hours for the full set on one M1 Ultra.

```bash
python3 -m venv venv && ./venv/bin/pip install mlx numpy pandas scipy matplotlib
./venv/bin/python build_brain.py          # needs the FlyWire v783 files, see below
./venv/bin/python evo_fix.py --homeo global --pop 2560 --conds A,N1,N2,N4,N5 \
    --seed 0 --gens 600 --probe_every 50 --probe_reps 8 --run demo
./venv/bin/python analyze_fix.py
```

`data/brain.npz` is included, so the engines run without rebuilding the brain. To rebuild it, fetch the public FlyWire v783 connectivity table and the neuron annotations (see the header of `build_brain.py`); they are not redistributed here.

A single 600-generation run with five conditions takes about 34 minutes on an M1 Ultra. `run_fix.sh`, `run_swap.sh` and `run_dose.sh` drive the full grids from the directory they sit in.

Re-running a seed does not repeat a run bit for bit. After every common-garden probe, which uses a fixed random seed, the engines reseed the evolutionary random stream from the clock, and the GPU scatter used in mutation does not fix which write survives when two mutations hit the same weight. A re-run reproduces the null constructions, calibration and generation 0 exactly and later generations only in distribution; seeds are independent replicates, not reproducible trajectories.

## Data availability

FlyWire v783 connectivity and annotations are public and are used under their published terms (Dorkenwald et al., *Nature* 634, 124–138, 2024; Schlegel et al., *Nature* 634, 139–152, 2024). Raw run outputs (about 40 GB of per-generation logs and stored genomes) are not in this repository; the analysis outputs they produce are in `results/`.

An archived snapshot of this repository, with the manuscripts, is on Zenodo: [10.5281/zenodo.22871090](https://doi.org/10.5281/zenodo.22871090) (resolves to the latest version).

## Citation

Park, G. (2026) Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison. Preprint. Code and data: https://doi.org/10.5281/zenodo.22871090

## License

MIT, see [`LICENSE`](LICENSE). The research log and protocols are released under the same terms.
