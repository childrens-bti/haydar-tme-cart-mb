# Downstream Analysis of Murine Brain Tumor CAR-T scRNA-seq

## Usage

`bash run_module.sh`

## Folder contents

1. `01-myeloid-condition-umaps.Rmd`: subsets myeloid cells by treatment condition, reprocesses each subset independently, and generates UMAP plots of myeloid subtypes for each condition.
2. `02-lineage-plot.Rmd`: generates dot plot of top lineage markers and stacked barplot of lineage proportions across conditions.
3. `03-myeloid-functional-marker-dotplots.Rmd`: performs differential expression analysis in myeloid cells by subcluster, subtype, and treatment condition, annotates marker genes with Hallmark pathways, and generates dot plots for visualization.
4. `04-myeloid-pseudobulk-gsea.Rmd`: performs pseudobulk DESeq2‑based Hallmark GSEA for myeloid cells across conditions.
5. `05-myeloid-composition-analysis.Rmd`: performs myeloid subtype composition analysis using sccomp with 41BB-L and tumor as the baseline references and generates composition plots and stacked bar plots of subtype proportions across conditions.
6. `06-tcell-reclustering.Rmd`: performs re-clustering of T cells after additional filtering and evaluates consistency with previous subtype annotations.
7. `07-dc-reclustering.Rmd`: performs re-clustering of dendritic cells after additional filtering, evaluates consistency with previous subtype annotations, and runs AUCell to evaluate Hallmark pathway activity across updated DC clusters.
8. `08-tcell-trajectory.Rmd`: performs Slingshot-based T-cell pseudotime analysis and generates UMAP visualizations of lineage trajectories.
9. `09-myeloid-milor-analysis.Rmd`: runs miloR neighborhood differential abundance testing for the myeloid population using tumor as the reference condition.
10. `10-tcell-milor-analysis.Rmd`: runs miloR neighborhood differential abundance testing for the refined T cell population using tumor as the reference condition.

## Analysis module directory structure
```
downstream-analyses/
├── 01-myeloid-condition-umaps.Rmd
├── 01-myeloid-condition-umaps.html
├── 02-lineage-plot.Rmd
├── 02-lineage-plot.html
├── 03-myeloid-functional-marker-dotplots.Rmd
├── 03-myeloid-functional-marker-dotplots.html
├── 04-myeloid-pseudobulk-gsea.Rmd
├── 04-myeloid-pseudobulk-gsea.html
├── 05-myeloid-composition-analysis.Rmd
├── 05-myeloid-composition-analysis.html
├── 06-tcell-reclustering.Rmd
├── 06-tcell-reclustering.html
├── 07-dc-reclustering.Rmd
├── 07-dc-reclustering.html
├── 08-tcell-trajectory.Rmd
├── 08-tcell-trajectory.html
├── 09-myeloid-milor-analysis.Rmd
├── 09-myeloid-milor-analysis.html
├── 10-tcell-milor-analysis.Rmd
├── 10-tcell-milor-analysis.html
├── README.md
├── input
│   ├── cart_lineage_markers.csv
│   ├── myeloid_subcluster_markers.csv
│   └── myeloid_subtype_markers.csv
├── plots
│   ├── cart_dotplot_top_lineage_markers.pdf
│   ├── cart_lineage_proportions_stacked_barplot.png
│   ├── dc
│   │   ├── filtered_dc_cluster_vs_old_subtype_heatmap.pdf
│   │   ├── filtered_dc_old_subtype_labels_on_reclustered_umap.pdf
│   │   ├── filtered_dc_old_subtype_marker_violinplots.pdf
│   │   ├── filtered_dc_reclustered_proportions_stacked_barplot.pdf
│   │   ├── filtered_dc_reclustered_umap.pdf
│   │   ├── filtered_dc_reclustered_umap_by_treatment.pdf
│   │   ├── filtered_dc_subcluster_hallmark_auc_heatmap.pdf
│   │   └── filtered_dc_top_old_label_per_new_cluster_barplot.pdf
│   ├── myeloid
│   │   ├── GSEA
│   │   │   ├── myeloid_gsea_dot_pseudobulk_41BB-L_vs_otherCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_41BB-L_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Angiogenic TAMs_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Antigen-presenting TAMs_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Antigen-presenting monocytes_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_B7H3_vs_otherCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_B7H3_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_CD28-41BB_vs_otherCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_CD28-41BB_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_CD8-41BB_vs_otherCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_CD8-41BB_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Early activated microglia_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Homeostatic microglia_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Inflammatory TAMs_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Inflammatory activated microglia_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Inflammatory monocytes_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Lipid-associated TAMs_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Perivascular macrophages (chemokine-producing)_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Perivascular macrophages_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_STOP_vs_allCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_STOP_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_Stress-response microglia_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_TAM-like phagocytic microglia_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_cDC1_vs_rest.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_tumor_vs_allCAR.pdf
│   │   │   ├── myeloid_gsea_dot_pseudobulk_tumor_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_41BB-L_vs_otherCAR.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_41BB-L_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Angiogenic TAMs_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Antigen-presenting TAMs_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Antigen-presenting monocytes_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_B7H3_vs_otherCAR.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_B7H3_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_CD28-41BB_vs_otherCAR.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_CD28-41BB_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_CD8-41BB_vs_otherCAR.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_CD8-41BB_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Early activated microglia_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Homeostatic microglia_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Inflammatory TAMs_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Inflammatory activated microglia_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Inflammatory monocytes_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Lipid-associated TAMs_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Perivascular macrophages (chemokine-producing)_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Perivascular macrophages_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_STOP_vs_allCAR.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_STOP_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_Stress-response microglia_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_TAM-like phagocytic microglia_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_cDC1_vs_rest.pdf
│   │   │   ├── myeloid_volcano_pseudobulk_tumor_vs_allCAR.pdf
│   │   │   └── myeloid_volcano_pseudobulk_tumor_vs_rest.pdf
│   │   ├── myeloid_dotplot_subcluster_markers.pdf
│   │   ├── myeloid_dotplot_subcluster_markers_panelE.pdf
│   │   ├── myeloid_dotplot_subcluster_markers_panelF.pdf
│   │   ├── myeloid_dotplot_subtype_markers.pdf
│   │   ├── myeloid_dotplot_subtype_markers_exclude_inflammatory_monocytes.pdf
│   │   ├── myeloid_dotplot_subtype_markers_panelE.pdf
│   │   ├── myeloid_dotplot_subtype_markers_panelE_exclude_inflammatory_monocytes.pdf
│   │   ├── myeloid_dotplot_subtype_markers_panelF.pdf
│   │   ├── myeloid_dotplot_subtype_markers_panelF_exclude_inflammatory_monocytes.pdf
│   │   ├── myeloid_dotplot_treatment_markers.pdf
│   │   ├── myeloid_dotplot_treatment_markers_panelE.pdf
│   │   ├── myeloid_dotplot_treatment_markers_panelF.pdf
│   │   ├── myeloid_subtype_composition_sccomp_41BBL_baseline.pdf
│   │   ├── myeloid_subtype_composition_sccomp_tumor_baseline.pdf
│   │   ├── myeloid_subtype_milor_tumor_baseline_effect_dotplot.pdf
│   │   ├── myeloid_subtype_milor_tumor_baseline_logFC_boxplot.pdf
│   │   ├── myeloid_subtype_proportions_stacked_barplot.png
│   │   ├── myeloid_umap_subcluster.pdf
│   │   ├── myeloid_umap_subtype_cluster_labels.pdf
│   │   ├── umap_41bb_myeloid_subtypes.png
│   │   ├── umap_b7h3_myeloid_subtypes.png
│   │   ├── umap_cd28_myeloid_subtypes.png
│   │   ├── umap_cd8_myeloid_subtypes.png
│   │   └── umap_stop_myeloid_subtypes.png
│   └── tcell
│       ├── figure_2E_tcell_marker_dotplot_by_condition.pdf
│       ├── supplementary_figure_5B_tcell_marker_dotplot_by_subtype.pdf
│       ├── tcell_cluster_vs_old_subtype_heatmap.pdf
│       ├── tcell_gsea_dot_pseudobulk_Activated_CD4_effector_helper_like_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Activated_effector_CD8_T_cells_stress_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Antigen_presenting_myeloid_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Cycling_proliferating_CAR_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_IFN_responsive_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_NK_like_cytotoxic_lymphocytes_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Naive_central_memory_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Quiescent_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Terminally_differentiated_cytotoxic_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Th2_like_inflammatory_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_Treg_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_gamma_delta_T_cells_vs_rest.pdf
│       ├── tcell_gsea_dot_pseudobulk_gamma_delta_Th17_like_T_cells_vs_rest.pdf
│       ├── tcell_old_subtype_labels_on_reclustered_umap.pdf
│       ├── tcell_old_subtype_marker_violinplots.pdf
│       ├── tcell_reclustered_proportions_stacked_barplot.pdf
│       ├── tcell_reclustered_umap.pdf
│       ├── tcell_reclustered_umap_by_treatment.pdf
│       ├── tcell_rename_clustered_umap.pdf
│       ├── tcell_slingshot_gene_trend_curves_all_lineages.pdf
│       ├── tcell_slingshot_gene_trend_heatmaps_all_lineages.pdf
│       ├── tcell_slingshot_lineage_curves_each_lineage_umap.pdf
│       ├── tcell_slingshot_lineage_curves_overview_umap.pdf
│       ├── tcell_slingshot_pseudotime_all_lineages_summary.pdf
│       ├── tcell_slingshot_pseudotime_umap.pdf
│       ├── tcell_subtype_milor_tumor_baseline_effect_dotplot.pdf
│       ├── tcell_subtype_milor_tumor_baseline_logFC_boxplot.pdf
│       ├── tcell_top_old_label_per_new_cluster_barplot.pdf
│       ├── tcell_volcano_pseudobulk_Activated_CD4_effector_helper_like_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Activated_effector_CD8_T_cells_stress_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Antigen_presenting_myeloid_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Cycling_proliferating_CAR_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_IFN_responsive_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_NK_like_cytotoxic_lymphocytes_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Naive_central_memory_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Quiescent_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Terminally_differentiated_cytotoxic_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Th2_like_inflammatory_T_cells_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_Treg_vs_rest.pdf
│       ├── tcell_volcano_pseudobulk_gamma_delta_T_cells_vs_rest.pdf
│       └── tcell_volcano_pseudobulk_gamma_delta_Th17_like_T_cells_vs_rest.pdf
├── results
│   ├── cart_lineage_proportions.tsv
│   ├── dc
│   │   ├── filtered_dc_reclustered_markers.csv
│   │   ├── filtered_dc_reclustered_proportions.tsv
│   │   ├── filtered_dc_reclustered_top30_markers_per_cluster.csv
│   │   └── filtered_dc_top_old_label_per_new_cluster.tsv
│   ├── myeloid
│   │   ├── GSEA
│   │   │   ├── myeloid_DESeq2_pseudobulk_41BB-L_vs_otherCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_41BB-L_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Angiogenic TAMs_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Antigen-presenting TAMs_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Antigen-presenting monocytes_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_B7H3_vs_otherCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_B7H3_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_CD28-41BB_vs_otherCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_CD28-41BB_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_CD8-41BB_vs_otherCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_CD8-41BB_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Early activated microglia_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Homeostatic microglia_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Inflammatory TAMs_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Inflammatory activated microglia_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Inflammatory monocytes_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Lipid-associated TAMs_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Perivascular macrophages (chemokine-producing)_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Perivascular macrophages_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_STOP_vs_allCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_STOP_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_Stress-response microglia_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_TAM-like phagocytic microglia_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_cDC1_vs_rest.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_tumor_vs_allCAR.csv
│   │   │   ├── myeloid_DESeq2_pseudobulk_tumor_vs_rest.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_41BB-L_vs_otherCAR_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_41BB-L_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Angiogenic TAMs_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Antigen-presenting TAMs_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Antigen-presenting monocytes_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_B7H3_vs_otherCAR_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_B7H3_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_CD28-41BB_vs_otherCAR_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_CD28-41BB_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_CD8-41BB_vs_otherCAR_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_CD8-41BB_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Early activated microglia_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Homeostatic microglia_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Inflammatory TAMs_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Inflammatory activated microglia_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Inflammatory monocytes_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Lipid-associated TAMs_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Perivascular macrophages (chemokine-producing)_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Perivascular macrophages_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_STOP_vs_allCAR_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_STOP_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_Stress-response microglia_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_TAM-like phagocytic microglia_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_cDC1_vs_rest_hallmark.csv
│   │   │   ├── myeloid_GSEA_pseudobulk_tumor_vs_allCAR_hallmark.csv
│   │   │   └── myeloid_GSEA_pseudobulk_tumor_vs_rest_hallmark.csv
│   │   ├── myeloid_subtype_composition_sccomp_41BBL_baseline_results.tsv
│   │   ├── myeloid_subtype_composition_sccomp_tumor_baseline_results.tsv
│   │   ├── myeloid_subtype_milor_tumor_baseline_effect_summary.tsv
│   │   ├── myeloid_subtype_milor_tumor_baseline_results.tsv
│   │   ├── myeloid_subtype_milor_tumor_baseline_threshold_summary.tsv
│   │   ├── myeloid_subtype_proportions.tsv
│   │   └── myeloid_treatment_markers.csv
│   └── tcell
│       ├── cart_tcell_subtypes.rds
│       ├── tcell_DESeq2_pseudobulk_Activated_CD4_effector_helper_like_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Activated_effector_CD8_T_cells_stress_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Antigen_presenting_myeloid_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Cycling_proliferating_CAR_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_IFN_responsive_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_NK_like_cytotoxic_lymphocytes_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Naive_central_memory_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Quiescent_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Terminally_differentiated_cytotoxic_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Th2_like_inflammatory_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_Treg_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_gamma_delta_T_cells_vs_rest.csv
│       ├── tcell_DESeq2_pseudobulk_gamma_delta_Th17_like_T_cells_vs_rest.csv
│       ├── tcell_GSEA_pseudobulk_Activated_CD4_effector_helper_like_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Activated_effector_CD8_T_cells_stress_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Antigen_presenting_myeloid_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Cycling_proliferating_CAR_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_IFN_responsive_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_NK_like_cytotoxic_lymphocytes_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Naive_central_memory_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Quiescent_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Terminally_differentiated_cytotoxic_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Th2_like_inflammatory_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_Treg_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_gamma_delta_T_cells_vs_rest_hallmark.csv
│       ├── tcell_GSEA_pseudobulk_gamma_delta_Th17_like_T_cells_vs_rest_hallmark.csv
│       ├── tcell_reclustered_markers.csv
│       ├── tcell_reclustered_proportions.tsv
│       ├── tcell_reclustered_top30_markers_per_cluster.csv
│       ├── tcell_slingshot_fast_spearman_pseudotime_gene_screen.csv
│       ├── tcell_subtype_milor_tumor_baseline_effect_summary.tsv
│       ├── tcell_subtype_milor_tumor_baseline_results.tsv
│       ├── tcell_subtype_milor_tumor_baseline_threshold_summary.tsv
│       └── tcell_top_old_label_per_new_cluster.tsv
├── run_module.sh
└── util
    ├── dotplot_helpers.R
    ├── milor_helpers.R
    ├── pseudobulk_gsea_helpers.R
    ├── sccomp_helpers.R
    └── trajectory_helpers.R
```