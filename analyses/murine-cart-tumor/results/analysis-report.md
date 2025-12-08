# Murine Brain Tumor CAR-T scRNA-seq Analysis Report

Prepared by: Sam Chen
Date: Dec 2025

# Overview
This document summarizes the analysis performed for the **Haydar Lab murine brain tumor CAR-T scRNA-seq**. The goal was to assess the heterogeneity of myeloid and T cells across CAR-T design conditions.

# Key takeaways
- Read depth per cell are low across samples, with high ambient RNA contamination. This may stem from the trade off between 10X flex and chromium kits and would adversely affect the analysis for smaller cell populations.
- **Interesting patterns emerged in myeloid cells**:
    - Control samples (tumor and STOP) are enriched for homeostatic and activated microglia
    - Macrophage populations expand in CAR-T treated samples
        - Lipid-associated TAMs are enriched in B7H3, 41BB-L, and CD28-41BB samples
        - B7H3 is specifically enriched for interferon responsive macrophages
        - Perivascular macrophages are depleted in B7H3 and 41BB-L samples
- **Hard to draw any conclusions from T cells because**:
    - Strong ambient RNAs are still present despite aggressive SoupX correction
    - Only < 3k "T cells" left after further filtering, which translate to ~50-150 cells per sample.
    - T cell subtype annotations are of lower confidence compared to myeloid cells
    - Therefore, results from compositional analysis and GSEA should be taken as a grain of salt.

# Summary of methods used
## Main lineages processing & annotation
- Ran SoupX and DoubletFinder algorithms to correct for ambient RNA and remove doublets, respectively
- Standard Seurat QC, processing and integration workflow for 12 10X v2 flex libraries
- Annotated main cell lineages by combining automatic algorithm (SingleR) and manual marker inspection

- Deliverables (main_dataset/)
    - UMAP plots: cart_umap_lineage_annotations.png & cart_umap_clusters.png
    - Marker gene heatmap: cart_heatmap_top_lineage_markers.pdf
    - Marker gene lists: cart_cluster_markers.csv & cart_lineage_markers.csv
    - Lineage dot plot: cart_dotplot_lineage_proportions.png

## Myeloid heterogeneity
- Subclustering of macrophage, microglia, and monocyte clusters
- Refined into interpretable functional myeloid subtypes
    - Homeostatic/stress-response microglia
    - Activated microglia
    - Phagocytic microglia
    - Lipid-associated macrophages
    - Antigen-presenting macrophages
    - Perivascular macrophages
    - Interferon-responsive macrophages
    - Monocytes
- Compositional analysis of myeloid subtypes across CAR-T design conditions using sccomp package
- Gene sets enrichment analysis using AUCell and fgsea packages, utilizing the following MSigDB gene sets:
    - Hallmark
    - C5: Ontology
    - C7: Immunology

- Deliverables (myeloid_subtypes)
    - UMAP plots: myeloid_subtypes_umap.png & myeloid_subtypes_by_condition_umap.png
    - Marker gene heatmap: myeloid_subtype_heatmap_top_markers.pdf
    - Marker gene lists: myeloid_subcluster_markers.csv & myeloid_subtype_markers.csv
    - sccomp results
        - myeloid_subtype_composition_sccomp.png
        - myeloid_subtype_proportions_barplot.png
        - myeloid_subtype_proportions_dotplot.png
    - AUCell results: myeloid_subtype_hallmark_auc_heatmap.pdf
    - GSEA results
        - myeloid_subtype_hallmark_gsea_results.csv
        - myeloid_subtype_c5_gsea_results.csv
        - myeloid_subtype_c7_gsea_results.csv

## T cell heterogeneity
- Subclustering of lymphoid cell clusters
- Further enrichment and filtering for cleaner T cell populations by in silico gating on lineage gene module scores
- Refined into interpretable T cell subtypes
    - Heat-shock/stressed cytotoxic T cells
    - Foxo1+ Satb1+ T cells (not confident)
    - CD4+ regulatory T cells
    - Cxcr6+ Bhlhe40+ T cells (not confident)
    - Early activated T cells
    - NK cells
    - Myeloid contamination
    - B cell contamination
- Compositional analysis of T cell subtypes across CAR-T design conditions
- Gene sets enrichment analysis using MSigDB gene sets.

- Deliverables (tcell_subtypes)
    - UMAP plots: tcell_subtypes_umap.png & tcell_subtypes_by_condition.png
    - Marker gene heatmap: tcell_subtype_heatmap_top_markers.pdf
    - Marker gene lists: tcell_subtype_markers.csv
    - sccomp results
        - tcell_subtype_composition_sccomp.png
        - tcell_subtype_proportions_barplot.png
        - tcell_subtype_proportions_dotplot.png
    - AUCell results: tcell_subtype_hallmark_auc_heatmap.pdf
    - GSEA results
        - tcell_subtype_hallmark_gsea_results.csv
        - tcell_subtype_c5_gsea_results.csv
        - tcell_subtype_c7_gsea_results.csv

# To do
- Assessment of dendritic cell heterogeneity, may suffer from the same problems of T cells
- Analysis of the two "endpoint" samples, which can serve as additional validation for the results reported here. 
