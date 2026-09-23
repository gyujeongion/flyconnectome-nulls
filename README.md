# Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison

Code, protocols and data for the preprint of the same name. Agents whose brain is a compressed adult *Drosophila* FlyWire v783 connectome forage and evolve in a 2D world, next to identically evolved populations built on randomised wiring. What the null model preserves at the sensory-motor boundary changes the outcome: two standard randomisations create direct olfactory-to-motor connections that the compressed connectome lacks, and boundary-preserving nulls remove the difference. Rewiring comparable amounts and varying only where the swaps land leaves most of the effect with those connections.

**Paper:** [`paper/biorxiv.pdf`](paper/biorxiv.pdf) (v4.1, 18 pages) · [`paper/alife.pdf`](paper/alife.pdf) (7-page conference version) · [`paper/supplementary.pdf`](paper/supplementary.pdf)

## What the experiments show

| | result |
|---|---|
| Standard controls (column shuffle, degree-preserving swaps) | overtake the connectome by the last probe (10 of 10 seeds, p_holm = 0.0059) |
| Olfactory output landing directly on descending motor groups | connectome 0.012 %, standard controls 10.6–10.7 % |
| Boundary-preserving controls (randomise interior only) | seed-mean difference inside a ±0.10 equivalence bound: +0.002 (90 % CI −0.026 to +0.045) and −0.074 (90 % CI −0.096 to +0.002); one interior control still ahead by 0.26 in one ecology |
| Shortcuts transplanted into the connectome | +0.44 fitness at the last probe (10 of 10 seeds, p = 0.002); olfactory dependence 0.15 → 0.99 |
| Shortcuts removed from the control | −0.43 in the ecology where it led (p = 0.004); no detected difference from the connectome afterwards |
| Graded shortcut doses (1, 3, 5, 10 %) | fitness rises in step (p = 1.4 × 10⁻⁶, 9 of 10 seeds) and so does dependence on smell (0.28 → 1.26) |
| Sham rewiring about as much as the full dose (101 swaps against 68 to 129) | full dose ahead of sham by +0.534 (90 % CI +0.377 to +0.700, 10 of 10 seeds), including all 7 seeds where it rewired less; the sham's registered equivalence prediction **failed** (+0.034, CI −0.062 to +0.123) |
| Interior wiring scrambled in evolved agents | 46–55 % of fitness lost in every condition |

Six pre-registration protocols were frozen before their runs, with code and data hashes: [`protocols/`](protocols/). A seventh, `PROTOCOL_DIST.md`, was frozen for an activity-distribution control that we did not run; it is kept there rather than deleted. Four more were frozen for follow-up experiments not reported in the paper: `PROTOCOL_LEARN.md` and `PROTOCOL_ETA.md` (an ecology that rewards learning; its pilots were run by mistake in a world with no toxic patches, so they did not test the hypothesis and their results are withdrawn, see the amendments), `PROTOCOL_TRANSFER.md` (out-of-distribution transfer; its registered baseline turned out to be already out of distribution, see its amendment) and `PROTOCOL_SHAM2.md` (shams paired with each intermediate dose, frozen before its runs). One registered prediction failed and is reported as failed: a control rewiring about as much as the full shortcut transplant was predicted to show no effect against the connectome, and instead crossed the equivalence bound (see §3.5 of the paper and `SUPPLEMENTARY.md` §S11). One correction was made after pre-registration (a one-sided turning-noise term in the motor model); it is described in the paper (§2.7) and in [`RESEARCH_LOG.md`](RESEARCH_LOG.md), with the re-runs that replaced the affected results. The boundary-preserving nulls apply the standard practice of constraining a null to the property under test (Váša & Mišić, 2022) to the sensory-motor boundary.

## Layout

```
evo_paper.py      evolution engine used for the first pre-registered grid
evo_fix.py        same engine with the turning-noise term corrected, plus the
                  boundary-preserving controls N4 and N5
evo_swap.py       adds the shortcut-transplant conditions AS and N2R
evo_dose.py       adds graded doses AS1/AS3/AS5/AS10 and the SHAM rewiring control
evo_learn.py      follow-up engine: shams paired per dose, ecology search over odour
                  reversal (pick_eco.py, gate_pilot.py apply the frozen rules)
evo_assay.py      evo_mut.py  evo_noise.py  evo_int.py   read-only assay runners
build_brain.py    FlyWire v783 -> data/brain.npz (512 groups + 1,000 Kenyon cells)
analyze_*.py      the analyses reported in the paper, one file per experiment
make_figs*.py     figures
SUPPLEMENTARY.md  supplementary results, including the sensitivity analysis S11
make_paper_numbers.py  regenerates values quoted in the paper into results/paper_numbers.json
run_*.sh          batch drivers for a single machine
protocols/        frozen protocols: six reported in the paper, one frozen and not run,
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

## Data availability

FlyWire v783 connectivity and annotations are public and are used under their published terms (Dorkenwald et al., *Nature* 634, 124–138, 2024; Schlegel et al., *Nature* 634, 139–152, 2024). Raw run outputs (about 40 GB of per-generation logs and stored genomes) are not in this repository; the analysis outputs they produce are in `results/`.

An archived snapshot of this repository, with the manuscripts, is on Zenodo: [10.5281/zenodo.22871090](https://doi.org/10.5281/zenodo.22871090) (resolves to the latest version).

## Citation

Park, G. (2026) Null-model treatment of the sensory-motor boundary changes an evolutionary connectome comparison. Preprint. Code and data: https://doi.org/10.5281/zenodo.22871090

## License

MIT, see [`LICENSE`](LICENSE). The research log and protocols are released under the same terms.
