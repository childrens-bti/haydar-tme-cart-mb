# Murine Brain Tumor CAR-T scRNA-seq Analysis - main lineage processing module

## Usage

`bash run_module.sh`

## Folder contents

1. `01-seurat-processing.Rmd`: performs the standard Seurat processing and integration for 12 murine brain tumor samples with various CAR-T treatments.
2. `02-lineage-annotation.Rmd`: performs cell type annotation analysis for major cell lineages.
3. `03-subclustering.Rmd`: performs subsetting and subclustering of myeloid, lymphoid, and DC lineages.
4. `plots/`: plotting outputs from the main lineage processing steps.
5. `results/`: intermediate and final artifacts (marker tables and derived outputs).
6. `util/`: helper scripts used by the module.

Approximate runtime: 5 hours on a machine with 32 cores and 128GB of RAM.

## Analysis module directory structure
```
main-lineage-processing/
├── 01-seurat-processing.Rmd
├── 01-seurat-processing.html
├── 02-lineage-annotation.Rmd
├── 02-lineage-annotation.html
├── 03-subclustering.Rmd
├── 03-subclustering.html
├── README.md
├── run_module.sh
├── plots
│   ├── cart_dotplot_lineage_proportions.png
│   ├── cart_heatmap_top_lineage_markers.pdf
│   ├── cart_umap_clusters.png
│   ├── cart_umap_lineage_annotations.png
├── results
│   ├── cart_cluster_markers.csv
│   └── cart_lineage_markers.csv
└── util
    ├── run_doubletfinder.R
    └── run_soupx.R
```