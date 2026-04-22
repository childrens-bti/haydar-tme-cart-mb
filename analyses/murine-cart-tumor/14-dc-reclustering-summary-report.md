# Dendritic Cell (DC) Re-clustering Summary Report

**author: Bicna Song, Sam Chen**  
**date: 2026-04-20**

---

## What was done

- Loaded the previously annotated dendritic cell (DC) Seurat object.
- Excluded **Myeloid cells** and **Lymphocytes** before re-clustering.
- Observed that exclusion substantially reduced per-sample DC counts, particularly in: 
    - *Tumor_only_D1* (n = 27)
    - *Tumor_only_D2* (n = 45)
    - *STOP_pos_T_D1* (n = 44)
- Attempted CCA- and RPCA-based integration; both failed due to insufficient cell numbers in the smallest samples, leading to unstable anchor construction.
- Evaluated two options: 
    - Remove low-cell-count samples and re-run integration
    - Skip integration and proceed with PCA-based clustering
- Proceeded with **Option 2 (no integration)** at this stage.
- Reprocessed DCs using normalization, variable feature selection, scaling, PCA, neighbor graph construction, clustering, and UMAP **directly in PCA space**.
- Generated UMAPs colored by sample to assess batch effects.
- Generated marker genes for the new Seurat clusters.
- Compared **new Seurat clusters** to the **previous `dc_subtype` labels**.
- Summarized the top-matching old label for each new cluster.
- Identified clusters with clear old-label correspondence versus clusters that should be reviewed.
- Generated output tables and PDF plots for cluster structure, old-vs-new label agreement, marker validation, and subtype proportions.

---

## Rationale for skipping integration

- UMAP colored by sample showed **no strong sample-specific clustering**, suggesting limited batch effects within the DC compartment.
- DCs are relatively conserved transcriptionally, and PCA-based clustering is commonly sufficient for subtype resolution in rare immune populations.

**Conclusion:** Integration was not strictly necessary and was therefore omitted to ensure numerical stability and biological interpretability.

---

## Overall interpretation

Following DC-only re-clustering after exclusion of lymphocytes, the resulting clusters reflect distinct DC functional states, with clear separation of migratory and inflammatory programs, and partial mixing of classical DC lineages.

**Key observations:**

- Clusters 0, 2, 4, and 5 are dominated by Migratory DC signatures, indicating a shared transcriptional program associated with DC maturation and trafficking.
- Cluster 3 is strongly enriched for Inflammatory DC (IL12+), representing a distinct pro-inflammatory activation state.
- Cluster 1 contains a mixture of cDC1 and cDC2 cells, suggesting that these two classical DC lineages are not transcriptionally separable under the current filtering and clustering conditions.

---

## Old subtype label composition within new clusters

| New cluster | Dominant DC subtype | n | Proportion | Assessment |
|---|---|---|---|---|
| 0 | Migratory DC | 546 | 0.814 | Clear |
| 1 | cDC2 | 311 | 0.549 | Mixed |
|   | cDC1 | 251 | 0.443 | Mixed |
| 2 | Migratory DC | 532 | 0.959 | Clear |
| 3 | Inflammatory DC (IL12+) | 344 | 0.847 | Clear |
| 4 | Migratory DC| 349 | 0.936 | Clear |
| 5 | Migratory DC| 53 | 0.616 | Review |
|   | Inflammatory DC (IL12+) | 33 | 0.384 | Review |

---

## Output files generated

### Result files

- `filtered_dc_reclustered_markers.csv` — marker genes for each new Seurat cluster
- `filtered_dc_reclustered_top30_markers_per_cluster.csv` — top 30 marker genes per new Seurat cluster
- `filtered_dc_top_old_label_per_new_cluster.tsv` — top matched old label per new cluster
- `filtered_dc_reclustered_proportions.tsv` — subtype counts and proportions by condition
- `cart_dc_reclustered.rds` — final reclustered Seurat object

### Plot files

- `filtered_dc_reclustered_umap.pdf` — UMAP of new Seurat clusters
- `filtered_dc_reclustered_umap_by_treatment.pdf` — UMAP split by condition
- `filtered_dc_old_subtype_labels_on_reclustered_umap.pdf` — UMAP colored by previous subtype labels
- `filtered_dc_cluster_vs_old_subtype_heatmap.pdf` — heatmap of old-label composition within each new cluster
- `filtered_dc_old_subtype_marker_violinplots.pdf` — combined violin plots for subtype key markers across new clusters
- `filtered_dc_reclustered_proportions_stacked_barplot.pdf` — stacked barplot of subtype proportions by condition
- `filtered_dc_top_old_label_per_new_cluster_barplot.pdf` — barplot summarizing top old-label match per new cluster

---

# Summary of the DC subtypes
AI-assisted summerization of the full marker lists for each DC subtype.

### cDC1
- Key markers: Xcr1, Clec9a, Wdfy4, Mycl, Itgae (CD103)
- Specialized cross-presenting DCs optimized for antigen presentation to CD8⁺ T cells. Support cytotoxic T-cell priming and anti-tumor immunity.
- Literature support:
  - Böttcher et al., Nat. Rev. Immunol. 2018: The role of conventional type 1 dendritic cells (cDC1) in cancer immune control
  - Heger et al., Proc. Natl. Acad. Sci. U.S.A. 2023: XCR1 expression distinguishes human conventional dendritic cell type 1 with full effector functions from their immediate precursors

### cDC2
- Key markers: Ciita, H2-Ab1, H2-Oa, Mgl2, Il1b
- APCs with strong antigen presentation machinery and inflammatory bias. May represent cDC2-leaning cells responding to local inflammatory cues.
- Literature support:
  - Saito et al., Cancers (Basel) 2022: The role of type-2 conventional dendritic cells in the regulation of tumor immunity
  - Cook et al., Nat. Commun. 2025: Mgl2⁺ cDC2s coordinate fungal allergic airway type 2, but not type 17, inflammation in mice

### Migratory DC
- Key markers: Ccr7, Cd80, Fscn1, Ccl22, Il15ra
- Mature DCs undergoing activation and migration toward lymphoid niches. Exhibit enhanced co-stimulatory capacity and are poised to initiate adaptive immune responses.
- Literature support:
  - Liu et al., Cellular & Molecular Immunology. 2021: Dendritic cell migration in inflammation and immunity
  - Hong et al., Front. Pharmacol. 2022: New insights of CCR7 signaling in dendritic cell migration

### Inflammatory DC (IL12+)
- Key markers: Il12b, Ly75 (DEC205), Fscn1, Scin
- Highly activated DCs with pro-inflammatory cytokine output, particularly IL-12. Likely promote Th1 polarization and cytotoxic immune responses.
- Literature support:
  - Heufler et al., Eur. J. Immunol. 1996: Interleukin-12 is produced by dendritic cells and mediates T helper 1 development as well as interferon-gamma production by T helper 1 cells
