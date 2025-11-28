# Murine Brain Tumor CAR-T scRNA-seq Analysis Report

Prepared by: Sam Chen
Date: Nov 2025

# Overvew
This document summarizes the analysis performed for the **Haydar Lab murine brain tumor CAR-T scRNA-seq**. The goal was to assess the heterogenity of myeloid and T cells across CAR-T design conditions.

# Summary of methods used
## Main lineages processing & annotation
- Ran SoupX and DoubletFinder algorithms to correct for ambient RNA and remove doublets, respectively
- Standard Seurat QC, processing and integration workflow for 12 10X v2 flex libraries
- Annotated main cell lienages by combining automatic algorithm (SingleR) and manual marker inspection
- Deliverables
    - UMAP plots: cart_umap_lineage_annotations.png & cart_umap_clusters.png
    - Marker gene heatmap: cart_heatmap_top_lineage_markers.pdf
    - Marker gene lists: cart_cluster_markers.csv & cart_lineage_markers.csv
    - Lineage dot plot: cart_dotplot_lineage_proportions.png

## Myeloid heterogenity
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
- Deliverables
    - 

## T cell heterogenity

# Key takeaways


# To do
