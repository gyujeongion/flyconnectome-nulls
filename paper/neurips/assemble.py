"""One-time assembly of the NeurIPS-style source paper/neurips/paper.md.

Main text: paper/neurips/main_body.md (written for this format).
References: those of PAPER.md plus the machine-learning references cited in the new sections.
Appendix: moved verbatim from PAPER.md (v4.2 draft) and SUPPLEMENTARY.md, re-headed, with section
cross-references remapped to the new labels. After this runs once, paper.md is the source to edit.
"""
import re, sys

ROOT = "../.."
paper = open(f"{ROOT}/PAPER.md").read()
supp = open(f"{ROOT}/SUPPLEMENTARY.md").read()
main = open("main_body.md").read()


def section(text, start, stop=None):
    """Text from the heading line matching `start` (exclusive) up to the next heading matching `stop`."""
    i = text.index(start)
    i = text.index("\n", i) + 1
    j = text.index(stop, i) if stop else len(text)
    return text[i:j].strip("\n -")


def fig(img, caption, label):
    return f"![{caption}]({img}){{#{label} width=100%}}"


P = lambda a, b: section(paper, a, b)
S = lambda a, b=None: section(supp, a, b)

# ---- pieces moved from PAPER.md
brain = P("### 2.1 Brain model", "### 2.2 Randomised controls")
controls = P("### 2.2 Randomised controls", "### 2.3 Calibration")
calib = P("### 2.3 Calibration", "### 2.4 World")
world = P("### 2.4 World", "### 2.5 Evolution")
evolution = P("### 2.5 Evolution", "### 2.6 Evaluation and statistics")
evalstats = P("### 2.6 Evaluation and statistics", "### 2.7 A correction")
correction = P("### 2.7 A correction made after pre-registration", "### 2.8 Assays")
assays = P("### 2.8 Assays", "### 2.9 Use of AI tools")
ai = P("### 2.9 Use of AI tools", "## 3. Results")
first = P("### 3.1 With standard controls", "### 3.2 The standard controls")
grid = P("### 3.3 Against boundary-preserving controls", "### 3.4 Adding the shortcut")
swap = P("### 3.4 Adding the shortcut", "### 3.5 At comparable rewiring")
dose = P("### 3.5 At comparable rewiring", "### 3.6 The evolved agents")
interior = P("### 3.6 The evolved agents use their interior wiring", "### 3.7 What survives")
survives = P("### 3.7 What survives from the first experiment", "## 4. Discussion")
refs = P("## References", None)


def figblock(text, old_img, caption_prefix, label):
    """Replace '![Figure N](path)\\n\\n**Figure N.** caption' with a labelled pandoc figure."""
    m = re.search(r"!\[[^\]]*\]\(" + re.escape(old_img) + r"\)\s*\n\s*\n\*\*" + re.escape(caption_prefix) + r"\*\*\s*(.+?)(?:\n\n|$)", text, re.S)
    assert m, old_img
    return text[:m.start()] + fig(old_img, m.group(1).strip(), label) + "\n\n" + text[m.end():]


first = figblock(first, "figs/fig5_eco.png", "Figure 5.", "fig:first")
grid = figblock(grid, "figs/fig7_fix.png", "Figure 7.", "fig:grid")
swap = re.sub(r"!\[Figure 8\]\(figs/fig8_swap.png\)\s*\n\s*\n\*\*Figure 8\.\*\*.*?\n\n", "", swap, flags=re.S)  # in main text
dose = re.sub(r"!\[Figure 9\]\(figs/fig9_dose.png\)\s*\n\s*\n\*\*Figure 9\.\*\*.*?(\n\n|$)", "", dose, flags=re.S)
# the corrected-grid table and the main-text paragraphs already in the main text are dropped from the appendix copy
grid = re.sub(r"\*\*Table 2\.\*\*.*?\n\n\|.*?\n\n", "", grid, flags=re.S)

# ---- pieces moved from SUPPLEMENTARY.md
s_calib = S("## S1. Calibration schemes", "### S2")
s_v1 = S("### S2 Apparent", "## S8.")
s_v1figs = S("## S8. Figures from paper-v1", "## S9.")
s_eco = S("## S9. Additional analyses", "## S10.")
s_changes = S("## S10. Design changes", "## S11.")
s_sham = S("## S11. Sensitivity analysis", "## S12.")
s_delta = S("## S12. Sensitivity of the equivalence", "## S13.")
s_n20 = S("## S13. Extension of the corrected grid", "## S14.")
s_bsham = S("## S14. Boundary-targeted sham", "## S15.")
s_dist2 = S("## S15. Corrected grid under distribution-matched", None)
for n in (1, 2, 3, 4):
    s_v1figs = figblock(s_v1figs, re.search(rf"figs/[^)]+", s_v1figs.split(f"![Figure S{n}]")[1]).group(0), f"Figure S{n}.", f"fig:v1-{n}")
s_v1 = re.sub(r"^### (S\d) ", "## ", s_v1, flags=re.M)
s_v1 = (s_v1.replace("generation-250 elites", "generation-250 population samples")
            .replace("Generation-250 elites", "Generation-250 population samples")
            .replace("connectome elites outperform", "connectome population samples outperform"))

nulls_detail = """The randomisations are implemented in `evo_fix.py` (`shuffle_control`, `null_degswap`, `null_interface_degswap`, `null_interface_colshuffle`); the dose and sham operations in `evo_dose.py`. Precisely:

- **Sensory rows.** Sensory groups receive no input in the connectome (their rows are zero). N1 places each column's weights only in non-sensory rows; degree-preserving swaps exchange the targets of two existing edges, whose targets are never sensory, so no control gives a sensory group an input.
- **Swaps.** N2 and N4 make 10·|E| swap attempts on their edge set; an attempt is rejected if the two edges share a source or a target or if either new edge already exists (for N4, existence is checked against the whole matrix, so interior swaps cannot duplicate a boundary edge). Accepted swaps are not counted by the engine; interior overlap with the connectome (0.111 to 0.121 for N4, 0.089 to 0.095 for N5) is the mixing diagnostic we report.
- **N5.** For each non-sensory source, its interior outputs are moved to randomly chosen free rows among the non-sensory, non-descending groups; its boundary outputs stay in place.
- **Kenyon-cell inputs.** The projection-neuron-to-Kenyon-cell matrix is randomised by degree-preserving swaps in N2, N4 and N5, and by shuffling which Kenyon cells receive each group's input in N1. The Kenyon-cell-to-group matrix and initial Kenyon-cell-to-output weights are shared.
- **Boundary check.** For N4 and N5 the engine compares every boundary edge with the connectome before calibration and aborts if any differs (maximum absolute difference 0 in every construction)."""

remap = [
    (r"§2\.1", r"Appendix \\ref{app-model}"), (r"§2\.2", r"Appendix \\ref{app-nulls}"),
    (r"§2\.3", r"Appendix \\ref{app-nulls}"), (r"§2\.5", r"Appendix \\ref{app-model}"),
    (r"§2\.6", r"Appendix \\ref{app-eval}"), (r"§2\.7", r"Appendix \\ref{app-prereg}"),
    (r"§2\.8", r"Appendix \\ref{app-checks}"), (r"§2\.9", r"Appendix \\ref{app-ai}"),
    (r"§3\.1", r"Appendix \\ref{app-first}"), (r"§3\.2", r"Section \\ref{sec-shortcuts}"),
    (r"§3\.3", r"Section \\ref{sec-equivalence}"), (r"§3\.4", r"Section \\ref{sec-swap}"),
    (r"§3\.5", r"Section \\ref{sec-dose}"), (r"§3\.6", r"Appendix \\ref{app-checks}"),
    (r"§3\.7", r"Appendix \\ref{app-checks}"), (r"§5\b", r"the data availability statement"),
    (r"Supplement §S11", r"Appendix \\ref{app-dose}"), (r"Supplementary S12", r"Appendix \\ref{app-grid}"),
    (r"main-text Appendix", r"Appendix"),
    (r"with the same conclusions \(supplement\)", r"with the same conclusions (Appendix \\ref{app-v1})"),
    (r"reported with their metrics in the supplement", r"reported with their metrics in Appendix \\ref{app-v1}"),
    (r"\(Fig\. 5\)", r"(Figure \\ref{fig:first})"), (r"Figure 7\b", r"Figure \\ref{fig:grid}"),
    (r"Table 1\b", r"Table \\ref{tbl:paths}"), (r"Table 2\b", r"Table \\ref{tbl:grid}"),
]

appendix = f"""
\\appendix

# Model details {{#app-model}}

## Brain model

{brain}

## World

{world}

## Evolution

{evolution}

# Controls, interventions and calibration {{#app-nulls}}

## Definitions

{controls}

## Implementation

{nulls_detail}

## Calibration

{calib}

{s_calib}

# Evaluation and statistics {{#app-eval}}

{evalstats}

# Pre-registration record and the noise correction {{#app-prereg}}

{correction}

## Design changes and corrections, in order

{s_changes}

# First grid {{#app-first}}

{first}

# Corrected grid: further results {{#app-grid}}

{grid}

## Sensitivity of the equivalence verdicts to the bound

{s_delta}

## Distribution-matched calibration (D-06)

{s_dist2}

## Extension to 20 seeds in two ecologies

{s_n20}

## Additional analyses (post hoc)

{s_eco}

# Intervention and dose-response: further details {{#app-dose}}

## Intervention

{swap}

## Dose-response and sham

{dose}

## The sham equivalence test on an unregistered outcome

{s_sham}

## Boundary-targeted sham: full results and execution record

{s_bsham}

# Further checks {{#app-checks}}

## Assays

{assays}

## Interior perturbation

{interior}

## What survives from the first experiment

{survives}

# Earlier model version (paper-v1) {{#app-v1}}

The results in this appendix come from paper-v1, the first pre-registered experiment (`PROTOCOL.md`: 4 conditions × 10 seeds × 300 generations in one ecology with toxic probability 0.5, toxic-bite loss 0.4 and a predator pursuing at 0.1·v_max; common garden every 25 generations), run before the ecology grids and before the noise correction. Stored genomes there are also population samples rather than elites. It is included for completeness; the main text does not rest on it.

{s_v1}

{s_v1figs}

# Use of AI tools {{#app-ai}}

{ai}
"""
for a, b in remap:
    appendix = re.sub(a, b, appendix)

ml_refs = """Adebayo, J., Gilmer, J., Muelly, M., Goodfellow, I., Hardt, M. & Kim, B. (2018) Sanity checks for saliency maps. *Advances in Neural Information Processing Systems* 31.
Agarwal, R., Schwarzer, M., Castro, P. S., Courville, A. & Bellemare, M. G. (2021) Deep reinforcement learning at the edge of the statistical precipice. *Advances in Neural Information Processing Systems* 34.
Duan, S., Dong, L. L. & Fiete, I. (2025) From synapses to dynamics: obtaining function from structure in a connectome constrained model of the head direction circuit. *Advances in Neural Information Processing Systems* 38.
Gaier, A. & Ha, D. (2019) Weight agnostic neural networks. *Advances in Neural Information Processing Systems* 32.
Schaeffer, R., Miranda, B. & Koyejo, S. (2023) Are emergent abilities of large language models a mirage? *Advances in Neural Information Processing Systems* 36.
Zhou, H., Lan, J., Liu, R. & Yosinski, J. (2019) Deconstructing lottery tickets: zeros, signs, and the supermask. *Advances in Neural Information Processing Systems* 32."""
_cited_text = main + appendix
all_refs = sorted([l for l in (refs + "\n" + ml_refs).splitlines()
                   if l.strip() and re.search(r"\b" + re.escape(l.split(",")[0]) + r"\b", _cited_text)], key=lambda s: s.lower())
_dropped = [l.split("(")[0].strip() for l in (refs + "\n" + ml_refs).splitlines() if l.strip() and l not in all_refs]
print("uncited references dropped:", _dropped)
references = "\n# References {.unnumbered}\n\n" + "\n\n".join(all_refs) + "\n"

out = main.rstrip() + "\n" + references + appendix
open("paper.md", "w").write(out)
left = sorted(set(re.findall(r"§[0-9S][0-9.]*", out)))
print("paper.md written;", len(out.split()), "words total; unresolved § refs:", left)
