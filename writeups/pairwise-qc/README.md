# Pairwise-QC recommendation source

This directory contains the recovered source for `pairwise-qc-recommendation.html`.
The report was developed in a Claude Desktop scratchpad on 2026-07-05 and copied
into the repository as a self-contained HTML file. The scratchpad was temporary;
these files were reconstructed from the retained session trace.

## Report files

- `qc_page.html`: editable HTML/CSS source with three figure placeholders.
- `pairwise-qc-recommendation.html`: compiled standalone report with figures embedded.
- `build_html.py`: replaces the placeholders with base64-encoded PNGs from `figs/`.
- `NOTES_theory.md`: working theory and interpretation notes.

## Analysis scripts

- `sim_qc.R`, `testability.R`: controlled simulation and directional-testability analysis.
- `sim_figs.R`, `clean_figs.R`: exploratory and final report figures.
- `a549_qc.R`, `a549_figs.R`: A549 real-data analysis and diagnostic figures.
- `extract_rd7_cell_counts.R`, `extract_rd7.R`, `analyze_rd7.R`: Replogle RD7 extraction and analysis.
- `dropcount_rd7.R`, `dropcount_a549.R`, `dropcount_survey.R`: filter drop fractions across real screens.

The scripts retain the original machine-local input data locations under
`/Users/ekatsevi/data/`. Intermediate `.rds` files and generated `figs/` are not
included; running the analysis scripts recreates them when those inputs and R
dependencies are available.

## Build

From this directory, after generating `figs/fig1_mechanism.png`,
`figs/fig2_retention.png`, and `figs/figR1_a549_entanglement.png`:

```sh
python3 build_html.py
```

The compiled report is intentionally self-contained and has no external runtime assets.
