# TLR Pathway Module Score Analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-tlr-analysis.Rmd`: performs TLR pathway module score analysis, visualization, and condition-level heatmap generation.

## Analysis module directory structure
```
tlr-analysis/
├── 01-tlr-analysis.Rmd
├── 01-tlr-analysis.html
├── README.md
├── plots
│   ├── tlr_module_scores_heatmap_by_condition.pdf
│   ├── tlr_module_scores_heatmap_by_subtype.pdf
│   ├── tlr_module_scores_heatmap_by_subtype_condition.pdf
│   ├── tlr_module_scores_heatmap_by_condition_subtype.pdf
│   ├── tlr_module_scores_umap.pdf
│   ├── tlr_module_scores_umap_by_condition.pdf
│   ├── tlr_module_scores_violin_by_condition.pdf
│   ├── tlr_module_scores_violin_by_subtype.pdf
│   └── tlr_module_scores_violin_by_subtype_condition.pdf
├── results
│   └── tlr_module_scores.rds
└── run_module.sh
```

