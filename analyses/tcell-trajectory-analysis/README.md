# T-cell Trajectory Analysis of Murine Brain Tumor CAR-T scRNA-seq

## Usage

Use the current project Docker image and download the v5 data release before running the module. From this directory, run:

```bash
bash run_module.sh
```

This renders analyses 01–03, then separate CD4-like and CD8-like tradeSeq (04) and Condiments (05) reports. Stages 04 and 05 reuse cached results when available; otherwise, they can take a long time.

## Reviews

- [CD4-like shared trajectory review](docs/tcell-trajectory-cd4-review.md)
- [CD8-like shared trajectory review](docs/tcell-trajectory-cd8-review.md)
- [Condiments review: outputs and interpretation](docs/tcell-trajectory-condiments-review.md)

## Folder contents

1. `01-tcell-trajectory.Rmd`: performs the legacy combined-T-cell Slingshot analysis and generates pseudotime, lineage, marker-trend, and pseudotime-associated-gene outputs.
2. `02-tcell-cd4-cd8-annotation.Rmd`: identifies clean CD4-like and CD8-like T cells using marker expression and ScGate support, adds ProjecTILs and cluster-level annotations, calculates AUCell program scores, and saves the prepared compartment objects.
3. `03-tcell-cd4-cd8-slingshot.Rmd`: runs shared all-condition Slingshot analyses for the CD4-like and CD8-like compartments by default, caches successful trajectory inference, and generates pseudotime, lineage, gene-trend, and AUCell program outputs.
4. `04-tcell-cd4-cd8-tradeseq.Rmd`: validates the shared CD4-like and CD8-like trajectories, fits condition-aware tradeSeq models to a balanced gene set, tests lineage-associated and condition-dependent expression, validates the manually reviewed lineage labels, and generates fitted gene-trend and heatmap outputs.
5. `05-tcell-cd4-cd8-condiments.Rmd`: compares CAR constructs using topology, progression, and fate-selection tests. A global topology test refits trajectories on all tested constructs, then fits each construct's curves and compares them with permuted labels; pairwise progression and fate selection use the fixed shared trajectories for each prespecified condition pair. Separate PDFs distinguish cell-level test evidence from weighted cell distributions and mouse-level descriptive summaries.

## Analysis module directory structure

```text
.
├── 01-tcell-trajectory.Rmd
├── 01-tcell-trajectory.html
├── 02-tcell-cd4-cd8-annotation.Rmd
├── 02-tcell-cd4-cd8-annotation.html
├── 03-tcell-cd4-cd8-slingshot.Rmd
├── 03-tcell-cd4-cd8-slingshot.html
├── 04-tcell-cd4-cd8-tradeseq.Rmd
├── 04-tcell-cd4-tradeseq.html
├── 04-tcell-cd8-tradeseq.html
├── 05-tcell-cd4-cd8-condiments.Rmd
├── 05-tcell-cd4-condiments.html
├── 05-tcell-cd8-condiments.html
├── README.md
├── config
│   ├── condition-comparison-manifest.tsv
│   └── tcell-lineage-annotation-manifest.tsv
├── docs
│   ├── tcell-trajectory-cd4-review.md
│   ├── tcell-trajectory-cd8-review.md
│   └── tcell-trajectory-condiments-review.md
├── plots
│   ├── combined-tcell-trajectory
│   ├── cd4-cd8-annotation
│   ├── cd4-cd8-slingshot
│   │   ├── CD4-trajectory
│   │   └── CD8-trajectory
│   ├── cd4-cd8-tradeseq
│   └── cd4-cd8-condiments
├── results
│   ├── combined-tcell-trajectory
│   ├── cd4-cd8-annotation
│   ├── cd4-cd8-slingshot
│   │   ├── CD4-trajectory
│   │   └── CD8-trajectory
│   ├── cd4-cd8-tradeseq
│   └── cd4-cd8-condiments
├── run_module.sh
└── util
    ├── tcell_cd4_cd8_trajectory_helpers.R
    ├── tcell_condiments_helpers.R
    ├── tcell_tradeseq_helpers.R
    └── tcell_trajectory_helpers.R
```
