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
8. `08-grant-umaps.Rmd`: generates UMAP plots for inclusion in Haydar grant figures.
9. `09-lineage-plot.Rmd`: generates dot plot of top lineage markers and stacked barplot of lineage proportions across conditions.
10. `10-myeloid-functional-marker-dotplots.Rmd`: performs differential expression analysis in myeloid cells by subcluster, subtype, and treatment condition, annotates marker genes with Hallmark pathways, and generates dot plots for visualization.
11. `11-subtype-pseudobulk-gsea.Rmd`: performs Hallmark GSEA for a user-specified cell subtype across treatment conditions using a pseudobulk DESeq2 workflow.
12. `12-subtype-composition-analysis.Rmd`: performs subtype composition analysis using sccomp with 41BB-L as the baseline reference, and generates composition plots and stacked bar plots of subtype proportions across conditions for myeloid, T cell, or DC data.

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
├── 09-lineage-plot.Rmd
├── 09-lineage-plot.html
├── 10-myeloid-functional-marker-dotplots.Rmd
├── 10-myeloid-functional-marker-dotplots.html
├── 11-subtype-pseudobulk-gsea.Rmd
├── 11-subtype-pseudobulk-gsea.html
├── 12-subtype-composition-analysis-dc.html
├── 12-subtype-composition-analysis-myeloid.html
├── 12-subtype-composition-analysis-tcell.html
├── 12-subtype-composition-analysis.Rmd
├── README.md
├── plots
│   ├── cart_dotplot_lineage_proportions.png
│   ├── cart_dotplot_top_lineage_markers.pdf
│   ├── cart_heatmap_top_lineage_markers.pdf
│   ├── cart_lineage_proportions_stacked_barplot.png
│   ├── cart_umap_clusters.png
│   ├── cart_umap_lineage_annotations.png
│   ├── dc_subtype_composition_sccomp.png
│   ├── dc_subtype_composition_sccomp_41BBL_baseline.pdf
│   ├── dc_subtype_hallmark_auc_heatmap.pdf
│   ├── dc_subtype_heatmap_top_markers.pdf
│   ├── dc_subtype_proportions_barplot.png
│   ├── dc_subtype_proportions_dotplot.png
│   ├── dc_subtype_proportions_stacked_barplot.png
│   ├── dc_subtypes_by_condition.png
│   ├── dc_subtypes_umap.png
│   ├── myeloid_dotplot_subcluster_markers.pdf
│   ├── myeloid_dotplot_subcluster_markers_panelE.pdf
│   ├── myeloid_dotplot_subtype_markers_exclude_inflammatory_monocytes.pdf
│   ├── myeloid_dotplot_treatment_markers.pdf
│   ├── myeloid_dotplot_treatment_markers_panelF.pdf
│   ├── myeloid_gsea_dot_pseudobulk_41BB-L_vs_otherCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_41BB-L_vs_rest.pdf
│   ├── myeloid_gsea_dot_pseudobulk_B7H3_vs_otherCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_B7H3_vs_rest.pdf
│   ├── myeloid_gsea_dot_pseudobulk_CD28-41BB_vs_otherCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_CD28-41BB_vs_rest.pdf
│   ├── myeloid_gsea_dot_pseudobulk_CD8-41BB_vs_otherCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_CD8-41BB_vs_rest.pdf
│   ├── myeloid_gsea_dot_pseudobulk_STOP_vs_allCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_STOP_vs_rest.pdf
│   ├── myeloid_gsea_dot_pseudobulk_tumor_vs_allCAR.pdf
│   ├── myeloid_gsea_dot_pseudobulk_tumor_vs_rest.pdf
│   ├── myeloid_subtype_composition_sccomp.png
│   ├── myeloid_subtype_composition_sccomp_41BBL_baseline.pdf
│   ├── myeloid_subtype_hallmark_auc_heatmap.pdf
│   ├── myeloid_subtype_heatmap_top_markers.pdf
│   ├── myeloid_subtype_proportions_barplot.png
│   ├── myeloid_subtype_proportions_dotplot.png
│   ├── myeloid_subtype_proportions_stacked_barplot.png
│   ├── myeloid_subtypes_by_condition_umap.png
│   ├── myeloid_subtypes_umap.png
│   ├── myeloid_volcano_pseudobulk_41BB-L_vs_otherCAR.pdf
│   ├── myeloid_volcano_pseudobulk_41BB-L_vs_rest.pdf
│   ├── myeloid_volcano_pseudobulk_B7H3_vs_otherCAR.pdf
│   ├── myeloid_volcano_pseudobulk_B7H3_vs_rest.pdf
│   ├── myeloid_volcano_pseudobulk_CD28-41BB_vs_otherCAR.pdf
│   ├── myeloid_volcano_pseudobulk_CD28-41BB_vs_rest.pdf
│   ├── myeloid_volcano_pseudobulk_CD8-41BB_vs_otherCAR.pdf
│   ├── myeloid_volcano_pseudobulk_CD8-41BB_vs_rest.pdf
│   ├── myeloid_volcano_pseudobulk_STOP_vs_allCAR.pdf
│   ├── myeloid_volcano_pseudobulk_STOP_vs_rest.pdf
│   ├── myeloid_volcano_pseudobulk_tumor_vs_allCAR.pdf
│   ├── myeloid_volcano_pseudobulk_tumor_vs_rest.pdf
│   ├── tcell_subtype_composition_sccomp_41BBL_baseline.pdf
│   ├── tcell_subtype_composition_sccomp_STOP.png
│   ├── tcell_subtype_composition_sccomp_tumor.png
│   ├── tcell_subtype_hallmark_auc_heatmap.pdf
│   ├── tcell_subtype_heatmap_top_markers.pdf
│   ├── tcell_subtype_proportions_barplot.png
│   ├── tcell_subtype_proportions_dotplot.png
│   ├── tcell_subtype_proportions_stacked_barplot.png
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
│   ├── cart_cluster_markers.csv
│   ├── cart_lineage_markers.csv
│   ├── cart_lineage_proportions.tsv
│   ├── dc_subcluster_markers.csv
│   ├── dc_subtype_c5_gsea_results.csv
│   ├── dc_subtype_c7_gsea_results.csv
│   ├── dc_subtype_composition_sccomp_41BBL_baseline_results.tsv
│   ├── dc_subtype_hallmark_gsea_results.csv
│   ├── dc_subtype_markers.csv
│   ├── dc_subtype_proportions.tsv
│   ├── myeloid_DESeq2_pseudobulk_41BB-L_vs_otherCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_41BB-L_vs_rest.csv
│   ├── myeloid_DESeq2_pseudobulk_B7H3_vs_otherCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_B7H3_vs_rest.csv
│   ├── myeloid_DESeq2_pseudobulk_CD28-41BB_vs_otherCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_CD28-41BB_vs_rest.csv
│   ├── myeloid_DESeq2_pseudobulk_CD8-41BB_vs_otherCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_CD8-41BB_vs_rest.csv
│   ├── myeloid_DESeq2_pseudobulk_STOP_vs_allCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_STOP_vs_rest.csv
│   ├── myeloid_DESeq2_pseudobulk_tumor_vs_allCAR.csv
│   ├── myeloid_DESeq2_pseudobulk_tumor_vs_rest.csv
│   ├── myeloid_GSEA_pseudobulk_41BB-L_vs_otherCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_41BB-L_vs_rest_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_B7H3_vs_otherCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_B7H3_vs_rest_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_CD28-41BB_vs_otherCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_CD28-41BB_vs_rest_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_CD8-41BB_vs_otherCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_CD8-41BB_vs_rest_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_STOP_vs_allCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_STOP_vs_rest_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_tumor_vs_allCAR_hallmark.csv
│   ├── myeloid_GSEA_pseudobulk_tumor_vs_rest_hallmark.csv
│   ├── myeloid_subcluster_markers.csv
│   ├── myeloid_subtype_c5_gsea_results.csv
│   ├── myeloid_subtype_c7_gsea_results.csv
│   ├── myeloid_subtype_composition_sccomp_41BBL_baseline_results.tsv
│   ├── myeloid_subtype_hallmark_gsea_results.csv
│   ├── myeloid_subtype_markers.csv
│   ├── myeloid_subtype_proportions.tsv
│   ├── myeloid_treatment_markers.csv
│   ├── nkt_subcluster_markers.csv
│   ├── tcell_subtype_c5_gsea_results.csv
│   ├── tcell_subtype_c7_gsea_results.csv
│   ├── tcell_subtype_composition_sccomp_41BBL_baseline_results.tsv
│   ├── tcell_subtype_hallmark_gsea_results.csv
│   ├── tcell_subtype_markers.csv
│   └── tcell_subtype_proportions.tsv
├── run_module.sh
└── util
    ├── run_doubletfinder.R
    └── run_soupx.R
```