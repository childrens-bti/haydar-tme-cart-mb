# T-cell Trajectory Analysis of Murine Brain Tumor CAR-T scRNA-seq

## Usage

Render the existing reports for analyses 01–03:

`bash run_module.sh`

Render the CD4-like and CD8-like tradeSeq reports independently:

```bash
TRADESEQ_WORKERS=6 bash run_module.sh 04-tcell-cd4-cd8-tradeseq.Rmd --compartment CD4
TRADESEQ_WORKERS=6 bash run_module.sh 04-tcell-cd4-cd8-tradeseq.Rmd --compartment CD8
```

These commands write separate HTML reports and compartment-specific summary
tables. Valid cached `evaluateK` results are reused from the module results
directory. Fitted-model caches are checked first in the module results directory
and then in `data/v5`; a model is fitted only when neither location contains a
compatible cache.

## Reviews

- [CD4-like shared trajectory review](docs/tcell-trajectory-cd4-review.md)
- [CD8-like shared trajectory review](docs/tcell-trajectory-cd8-review.md)

## Folder contents

1. `01-tcell-trajectory.Rmd`: performs the legacy combined-T-cell Slingshot analysis and generates pseudotime, lineage, marker-trend, and pseudotime-associated-gene outputs.
2. `02-tcell-cd4-cd8-annotation.Rmd`: identifies clean CD4-like and CD8-like T cells using marker expression and ScGate support, adds ProjecTILs and cluster-level annotations, calculates AUCell program scores, and saves the prepared compartment objects.
3. `03-tcell-cd4-cd8-slingshot.Rmd`: runs condition-specific and shared Slingshot analyses for the CD4-like and CD8-like compartments, caches successful trajectory inference, and generates pseudotime, lineage, gene-trend, and AUCell program outputs.
4. `04-tcell-cd4-cd8-tradeseq.Rmd`: validates the shared CD4-like and CD8-like trajectories, fits condition-aware tradeSeq models to a balanced gene set, tests lineage-associated and condition-dependent expression, validates the manually reviewed lineage labels, and generates fitted gene-trend and heatmap outputs. The compartment-specific reports are saved as `04-tcell-cd4-tradeseq.html` and `04-tcell-cd8-tradeseq.html`.

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
├── README.md
├── config
│   ├── condition-comparison-manifest.tsv
│   └── tcell-lineage-annotation-manifest.tsv
├── docs
│   ├── tcell-trajectory-cd4-review.md
│   └── tcell-trajectory-cd8-review.md
├── plots
│   ├── combined-tcell-trajectory
│   ├── cd4-cd8-annotation
│   ├── cd4-cd8-slingshot
│   │   ├── CD4-trajectory
│   │   └── CD8-trajectory
│   └── cd4-cd8-tradeseq
├── results
│   ├── combined-tcell-trajectory
│   ├── cd4-cd8-annotation
│   ├── cd4-cd8-slingshot
│   │   ├── CD4-trajectory
│   │   └── CD8-trajectory
│   └── cd4-cd8-tradeseq
├── run_module.sh
└── util
    ├── tcell_cd4_cd8_trajectory_helpers.R
    ├── tcell_tradeseq_helpers.R
    └── tcell_trajectory_helpers.R
```
