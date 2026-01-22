# Murine Brain Tumor CAR-T scRNA-seq Analysis Report - re-sequenced

Prepared by: Sam Chen
Date: Jan 2026

# Overview
This document summarizes the analysis performed for the **Haydar Lab murine brain tumor CAR-T scRNA-seq after resequencing to a sufficient read depth**. The goal was to assess the heterogeneity of myeloid and T cells across CAR-T design conditions.

# Key takeaways
- Read depth per cell looks good now! As a result, we observed much better separation between key cell lineages in UMAP space compared to the previous sequencing run.
- **Myeloid Cells**:
    - Control samples (tumor and STOP) are enriched for homeostatic, stress-response, and early-activated microglia
    - Macrophage subtypes with tumor-associated signatures expand in CAR-T treated samples
        - Lipid-associated, angiogenic, and inflammatory TAMs are enriched in B7H3, 41BB-L, and CD28-41BB samples
        - Perivascular macrophages are depleted in B7H3 and 41BB-L samples
    - Compared to the previous report, interferon-responsive macrophages are no longer resolved as a distinct cluster and instead co-cluster with inflammatory TAMs. This suggests they represent a transient inflammatory activation state rather than a stable macrophage subtype, and their apparent separation in the earlier analysis was likely amplified by lower sequencing depth.
- **NK/T cells**:
    - Lineage separation is much improved compared to the previous report, enabling clearer identification of T cell subtypes. 
    - However, tumor and STOP samples exhibit distinct T cell subtype compositions, potentially confounding comparisons with CAR-T treated samples.
    - Compared to tumor samples, STOP and CAR-T treated samples are:
      - enriched for CD8+ cytotoxic, CD4+ helper, and CD4+ regulatory T cells
      - depleted for central memory and iNKT cells
- **Dendritic Cells**:
    - Significant contamination from myeloid cells and lymphocytes was observed in DC subclusters, likely due to lower abundance of DCs in the dataset.
    - Despite this, we were able to identify cDC1, cDC2, migratory DC, and IL12+ inflammatory DC subtypes.
    - With lower confidence due to lower cell numbers, CAR-T treated samples appear to be enriched for IL12+ inflammatory DCs and migratory DCs compared to tumor and STOP samples.

- **Endpoint samples**:
    - Analysis of endpoint samples largely recapitulated findings from the main dataset, confirming the observed myeloid and T cell subtype compositions across CAR-T design conditions.

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
   - Homeostatic Microglia
   - Stress-Response Microglia
   - Early Activated Microglia
   - Inflammatory Activated Microglia
   - TAM-like Phagocytic Microglia
   - Angiogenic Tumor-Associated Macrophages (TAMs)
   - Antigen-Presenting TAMs
   - Inflammatory TAMs
   - Lipid-associated TAMs
   - Perivascular Macrophages
   - Perivascular Macrophages (Chemokine producing)
   - Antigen-Presenting Monocytes
   - Inflammatory Monocytes
   - cDC1 (contamination in subclustering)
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
- Subclustering of NK/T cell clusters
- Refined into interpretable T cell subtypes
   - Heat-shock activated T cells
   - CD4+ Regulatory T cells
   - CD8+ Cytotoxic T cells
   - CD4+ Helper T cells
   - Early-activated T cells
   - Naive T cells
   - Central Memory T cells
   - IFN-stimulated T cells
   - iNKT1 cells
   - iNKT2 cells
   - iNKT17 cells
   - Proliferating cells
   - Myeloid contamination
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

## DC heterogeneity
- Subclustering of dendritic cell clusters
- Refined into interpretable DC subtypes
   - cDC1
   - cDC2
   - Migratory DC
   - Inflammatory DC (IL12+)
   - Myeloid cells (contamination)
   - Lymphocytes (contamination)
- Compositional analysis of DC subtypes across CAR-T design conditions
- Gene sets enrichment analysis using MSigDB gene sets.

- Deliverables (dc_subtypes)
    - UMAP plots: dc_subtypes_umap.png & dc_subtypes_by_condition.png
    - Marker gene heatmap: dc_subtype_heatmap_top_markers.pdf
    - Marker gene lists: dc_subtype_markers.csv
    - sccomp results
        - dc_subtype_composition_sccomp.png
        - dc_subtype_proportions_barplot.png
        - dc_subtype_proportions_dotplot.png
    - AUCell results: tcell_subtype_hallmark_auc_heatmap.pdf
    - GSEA results
        - dc_subtype_hallmark_gsea_results.csv
        - dc_subtype_c5_gsea_results.csv
        - dc_subtype_c7_gsea_results.csv

## Endpoint samples analysis
- Separate analysis of endpoint samples to verify consistency with main dataset findings
- Data processing and lineage annotation performed similarly to main dataset