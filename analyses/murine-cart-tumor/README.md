# Murine Brain Tumor CAR-T scRNA-seq Analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-seurat-processing.Rmd`: performs the standard Seurat processing and integration for 12 murine brain tumor samples with various CAR-T treatments.
2. `02-lineage-annotation.Rmd`: performs cell type annotation analysis for major cell lineages
3. `03-subclustering.Rmd`: performs subsetting and subclustering of myeloid, lymphoid, and DC lineages
4. `04-myeloid-characterization.Rmd`: performs in depth characterization of macophage and microglia subtypes present in the dataset. Includes identification of markers, annotation, composition analysis, and GSEA enrichements.
5. `05-lymphoid-characterization.Rmd`: Same as myeloid, but also performed a further round of cell filtering for T cells as there are still strong contamination from other lineages.
6. `06-dendritic-characterization.Rmd`: Same as myeloid, but for dendritic cells.
7. `07-endpoint-samples-analysis.Rmd`: Focused analysis on endpoint samples to serve as validation for trends observed in samples from earlier timepoints.
8. `08-grant-umaps.Rmd`: Generates UMAP plots for inclusion in Haydar grant figures.

## Analysis module directory structure
```
murine-cart-tumor/
├── 01-seurat-processing.Rmd
├── 01-seurat-processing.html
├── 02-lineage-annotation.Rmd
├── 02-lineage-annotation.html
├── 03-subclustering.Rmd
├── 03-subclustering.html
├── 04-myeloid-characterization.Rmd
├── 04-myeloid-characterization.html
├── 05-lymphoid-characterization.Rmd
├── 05-lymphoid-characterization.html
├── 06-dendritic-characterization.Rmd
├── 06-dendritic-characterization.html
├── 07-endpoint-samples-analysis.Rmd
├── 07-endpoint-samples-analysis.html
├── 08-grant-umaps.Rmd
├── 08-grant-umaps.html
├── README.md
├── plots
│   ├── cart_dotplot_lineage_proportions.png
│   ├── cart_heatmap_top_lineage_markers.pdf
│   ├── cart_umap_clusters.png
│   ├── cart_umap_lineage_annotations.png
│   ├── dc_subtype_composition_sccomp.png
│   ├── dc_subtype_hallmark_auc_heatmap.pdf
│   ├── dc_subtype_heatmap_top_markers.pdf
│   ├── dc_subtype_proportions_barplot.png
│   ├── dc_subtype_proportions_dotplot.png
│   ├── dc_subtypes_by_condition.png
│   ├── dc_subtypes_umap.png
│   ├── myeloid_subtype_composition_sccomp.png
│   ├── myeloid_subtype_hallmark_auc_heatmap.pdf
│   ├── myeloid_subtype_heatmap_top_markers.pdf
│   ├── myeloid_subtype_proportions_barplot.png
│   ├── myeloid_subtype_proportions_dotplot.png
│   ├── myeloid_subtypes_by_condition_umap.png
│   ├── myeloid_subtypes_umap.png
│   ├── tcell_subtype_composition_sccomp_STOP.png
│   ├── tcell_subtype_composition_sccomp_tumor.png
│   ├── tcell_subtype_hallmark_auc_heatmap.pdf
│   ├── tcell_subtype_heatmap_top_markers.pdf
│   ├── tcell_subtype_proportions_barplot.png
│   ├── tcell_subtype_proportions_dotplot.png
│   ├── tcell_subtypes_by_condition.png
│   ├── tcell_subtypes_umap.png
│   ├── umap_41bb_myeloid_subtypes.png
│   ├── umap_b7h3_myeloid_subtypes.png
│   ├── umap_cd28_myeloid_subtypes.png
│   ├── umap_cd8_myeloid_subtypes.png
│   └── umap_stop_myeloid_subtypes.png
├── results
│   ├── analysis-report.Rmd
│   ├── analysis-report.html
│   ├── analysis-report.md
│   ├── analysis-report_v2.Rmd
│   ├── analysis-report_v2.html
│   ├── analysis-report_v2.md
│   ├── cart_annotated.rds
│   ├── cart_cluster_markers.csv
│   ├── cart_combined.rds
│   ├── cart_dc_subclusters.rds
│   ├── cart_dc_subtypes.rds
│   ├── cart_integrated.rds
│   ├── cart_lineage_markers.csv
│   ├── cart_lymphoid_subclusters.rds
│   ├── cart_merged.rds
│   ├── cart_myeloid_subclusters.rds
│   ├── cart_myeloid_subtypes.rds
│   ├── cart_tcell_subclusters.rds
│   ├── dc_subcluster_markers.csv
│   ├── dc_subtype_c5_gsea_results.csv
│   ├── dc_subtype_c7_gsea_results.csv
│   ├── dc_subtype_hallmark_gsea_results.csv
│   ├── dc_subtype_markers.csv
│   ├── endpoint_combined.rds
│   ├── endpoint_integrated.rds
│   ├── endpoint_merged.rds
│   ├── myeloid_subcluster_markers.csv
│   ├── myeloid_subtype_c5_gsea_results.csv
│   ├── myeloid_subtype_c7_gsea_results.csv
│   ├── myeloid_subtype_hallmark_gsea_results.csv
│   ├── myeloid_subtype_markers.csv
│   ├── nkt_subcluster_markers.csv
│   ├── tcell_subtype_c5_gsea_results.csv
│   ├── tcell_subtype_c7_gsea_results.csv
│   ├── tcell_subtype_hallmark_gsea_results.csv
│   └── tcell_subtype_markers.csv
├── run_module.sh
└── util
    ├── run_doubletfinder.R
    └── run_soupx.R
```