# Murine Brain Tumor CAR-T scRNA-seq Analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-seurat-processing.Rmd`: performs the standard Seurat processing and integration for 12 murine brain tumor samples with various CAR-T treatments.
2. `02-lineage-annotation.Rmd`: performs cell type annotation analysis for major cell lineages
3. `03-subclustering.Rmd`: performs subsetting and subclustering of myeloid, lymphoid, and DC lineages
4. `04-myeloid-characterization.Rmd`: performs in depth characterization of macophage and microglia subtypes present in the dataset. Includes identification of markers, annotation, composition analysis, and GSEA enrichements.
5. `05-lymphoid-characterization.Rmd`: Same as myeloid, but also performed a further round of cell filtering for T cells as there are still strong contamination from other lineages.

## Analysis module directory structure
```
murine-cart-tumor/
├── results/
│   ├── cart_cluster_markers.csv
│   ├── cart_lineage_markers.csv
│   ├── myeloid_subcluster_markers.csv
│   ├── myeloid_subtype_c5_gsea_results.csv
│   ├── myeloid_subtype_c7_gsea_results.csv
│   ├── myeloid_subtype_markers.csv
│   ├── tcell_subtype_c5_gsea_results.csv
│   ├── tcell_subtype_c7_gsea_results.csv
│   ├── tcell_subtype_hallmark_gsea_results.csv
│   └── tcell_subtype_markers.csv
│
├── util/
│   ├── run_doubletfinder.R
│   └── run_soupx.R
│
├── 01-seurat-processing.Rmd
├── 01-seurat-processing.html
├── 02-lineage-annotation.Rmd
├── 02-lineage-annotation.html
├── 03-subclustering.Rmd
├── 03-subclustering.html
├── 04-myeloid-characterization.Rmd
├── 04-myeloid-characterization.html
├── 03-subclustering.Rmd
├── 04-myeloid-characterization.Rmd
├── 04-myeloid-characterization.html
├── 05-lymphoid-characterization.Rmd
├── 05-lymphoid-characterization.html
├── README.md
└── run_module.sh
```