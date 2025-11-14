# Murine Brain Tumor CAR-T scRNA-seq Analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-seurat-processing.Rmd`: performs the standard Seurat processing and integration for 12 murine brain tumor samples with various CAR-T treatments.
2. `02-lineage-annotation.Rmd`: performs cell type annotation analysis for major cell lineages

## Analysis module directory structure

tr37_analysis/
├── results/
│   ├── cart_cluster_markers.csv
│   └── cart_lineage_markers.csv
│
├── util/
│   ├── run_doubletfinder.R
│   └── run_soupx.R
│
├── 01-seurat-processing.Rmd
├── 02-lineage-annotation.Rmd
├── 02-lineage-annotation.html
├── README.md
└── run_module.sh