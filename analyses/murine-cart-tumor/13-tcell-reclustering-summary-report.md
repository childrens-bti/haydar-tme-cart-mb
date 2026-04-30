# T-cell Re-clustering Summary Report

**author: Bicna Song, Sam Chen**  
**date: 2026-04-09**

---

## What was done

- Loaded the previously annotated T-cell Seurat object.
- Excluded **Myeloid cells** and **NK cells** before re-clustering.
- Reprocessed filtered cells with normalization, variable feature selection, scaling, PCA, CCA-based integration, neighbor graph construction, clustering, and UMAP.
- Generated marker genes for the new Seurat clusters.
- Compared **new Seurat clusters** to the **previous `tcell_subtype` labels**.
- Summarized the top-matching old label for each new cluster.
- Identified clusters with clear old-label correspondence versus clusters that should be reviewed.
- Generated output tables and PDF plots for cluster structure, old-vs-new label agreement, marker validation, and subtype proportions.

---

## Overall interpretation

Re-clustering of the filtered T-cell object largely **recapitulates the previous subtype annotations**. Most new clusters show strong correspondence to a single previous `tcell_subtype` label, suggesting that **global re-annotation is not required**.

A smaller subset of clusters show only moderate label purity and should be **reviewed more carefully**, especially for possible refinement of:

- Heat-shock / early-activation states
- Proliferating states
- Central memory / helper-like states

### Bottom line

- **Likely keep existing labels for most clusters**
- **Perform targeted review for clusters 0, 2, 5, 6, and 7**
- Use marker validation plots and cluster markers to decide whether those clusters need refined labels

---

## Old subtype label match to new cluster

| New cluster | Top matched old label | Matched count | Matched proportion | Assessment |
|---|---|---:|---:|---|
| 0 | Heat-shocked T cells | 929 | 0.703 | Review |
| 1 | CD4+ Regulatory T cells | 1160 | 0.890 | Clear |
| 2 | Heat-shocked T cells | 735 | 0.621 | Review |
| 3 | iNKT1 cells | 788 | 0.956 | Clear |
| 4 | Proliferating | 790 | 0.997 | Clear |
| 5 | CD4+ Helper T cells | 493 | 0.716 | Review |
| 6 | Proliferating | 507 | 0.822 | Review |
| 7 | Central memory T cells | 403 | 0.689 | Review |
| 8 | CD8+ Cytotoxic T cells | 432 | 0.949 | Clear |
| 9 | Naive T cells | 411 | 0.909 | Clear |
| 10 | IFN-stimulated T cells | 271 | 0.900 | Clear |
| 11 | iNKT17 cells | 218 | 0.960 | Clear |
| 12 | iNKT2 cells | 122 | 0.897 | Clear |

### Clear clusters

These clusters show strong agreement with a previous subtype label and are likely fine to keep as-is unless marker genes strongly disagree.

- **1** → CD4+ Regulatory T cells
- **3** → iNKT1 cells
- **4** → Proliferating
- **8** → CD8+ Cytotoxic T cells
- **9** → Naive T cells
- **10** → IFN-stimulated T cells
- **11** → iNKT17 cells
- **12** → iNKT2 cells

### Clusters recommended for review

These clusters show moderate agreement and should be reviewed with marker plots and cluster markers before deciding whether label refinement is needed.

- **0** → Heat-shocked T cells (`matched_prop = 0.703`)
- **2** → Heat-shocked T cells (`matched_prop = 0.621`)
- **5** → CD4+ Helper T cells (`matched_prop = 0.716`)
- **6** → Proliferating (`matched_prop = 0.822`)
- **7** → Central memory T cells (`matched_prop = 0.689`)

### Notes on the review clusters

- **Clusters 0 and 2** both map to **Heat-shocked T cells**, suggesting that the original heat-shock / early-activation state may have split into multiple subclusters after re-clustering.
- **Cluster 6** maps to **Proliferating** but is less pure than cluster 4, which may indicate a secondary proliferative state or a lineage-specific cycling population.
- **Cluster 7** maps to **Central memory T cells** but with moderate purity, so marker validation is important before deciding whether to keep or refine that label.
- **Cluster 5** maps to **CD4+ Helper T cells** but should still be checked against helper, activation, and memory markers.

---

## Output files generated

### Result files

- `filtered_tcell_reclustered_markers.csv` — marker genes for each new Seurat cluster
- `filtered_tcell_reclustered_top30_markers_per_cluster.csv` — top 30 marker genes per new Seurat cluster
- `filtered_tcell_top_old_label_per_new_cluster.tsv` — top matched old label per new cluster
- `filtered_tcell_reclustered_proportions.tsv` — subtype counts and proportions by condition
- `cart_tcell_reclustered.rds` — final reclustered Seurat object

### Plot files

- `filtered_tcell_reclustered_umap.pdf` — UMAP of new Seurat clusters
- `filtered_tcell_reclustered_umap_by_treatment.pdf` — UMAP split by condition
- `filtered_tcell_old_subtype_labels_on_reclustered_umap.pdf` — UMAP colored by previous subtype labels
- `filtered_tcell_cluster_vs_old_subtype_heatmap.pdf` — heatmap of old-label composition within each new cluster
- `filtered_tcell_old_subtype_marker_violinplots.pdf` — combined violin plots for subtype key markers across new clusters
- `filtered_tcell_reclustered_proportions_stacked_barplot.pdf` — stacked barplot of subtype proportions by condition
- `filtered_tcell_top_old_label_per_new_cluster_barplot.pdf` — barplot summarizing top old-label match per new cluster

---

## Summary of the T cell subtypes
AI-assisted summerization of the full marker lists for each T cell subtype.

### Heat-shock Activated T cells
- Key markers: Hspa1a, Hspa1b, Dnajb1, Jun, Cd69, Ms4a4b, Rasal3, Btg2
- Recently activated or acutely stressed T cells, often induced by tissue dissociation, tumor microenvironment stress, or strong stimulation. Represents a state, not a stable lineage.
- Literature support:
  - Hauet-Broere et al., Ann. Rheum. Dis. 2006: Heat shock proteins induce T cell regulation of chronic inflammation
  - Calderwood &amp; Gong, Trends Biochem. Sci. 2016: Heat Shock Proteins Promote Cancer: It’s a Protection Racket

### Myeloid cells (contamination)
- Key markers: Csf1r, Tmem119, Hexb, Cst3, Ccl3
- Microglia/macrophage-like cells. Clear non-lymphoid contamination and should be excluded from T/NKT-focused analyses.Potentially central memory T cells, but not clean

### CD4+ Regulatory T cells
- Key markers: Foxp3, Ctla4, Il2ra (CD25), Ikzf2, Tnfrsf4
- Immunosuppressive CD4+ T cells involved in immune tolerance and checkpoint regulation. Likely dampen effector responses within the tissue
- Literature support:
  - Rudensky, Immunol. Rev. 2011: Regulatory T Cells and Foxp3
  - Walker, J. Autoimmun. 2013: Treg and CTLA-4: Two intertwining pathways to immune tolerance
  - Chen et al., Int. Immunopharmacol. 2016: CD4⁺CD25⁺ regulatory T cells in tumor immunity

### CD8+ Cytotoxic T cells
- Key markers: Cd8a, Cd8b1, Gzmb, Nkg7, Prf1, Pdcd1
- Activated cytotoxic effector T cells with killing capacity. PD-1 expression suggests recent activation and/or early exhaustion.
- Literature support:
  - Raskov et al., Br. J. Cancer. 2021: Cytotoxic CD8⁺ T cells in cancer and cancer immunotherapy

### CD4+ Helper T cells
- Key markers: Cd4, Icos, Tnfrsf4 (OX40), Bhlhe40
- Activated helper T cells supporting immune coordination and cytokine signaling. Less cytotoxic, more regulatory/organizing in function.
- Literature support:
  - Luckheeram et al., Clin. Dev. Immunol. 2012: CD4⁺ T Cells: Differentiation and Functions

### Early-activation T cells
- Key markers: Egr1, Jun, Fos, Atf3, Zfp36, Dusp1, Hspa1a/b, Klf2
- “Shock response” T cells—very typical after stimulation or tissue dissociation
- Literature support:
  - Malissen &amp; Bongrand, Annu. Rev. Immunol. 2015: Early T cell activation: integrating biochemical, structural, and biophysical cues
  - Chen et al., Blood. 2009: Identification of early growth response protein 1 (EGR-1) as a novel target for JUN-induced apoptosis in multiple myeloma

### Proliferating
- Key markers: Mki67, Cdk1, Ccnb1, Ccnb2, Top2a
- Actively cycling lymphocytes. Identity is secondary; represents a cell-cycle state shared across multiple T/NKT subsets.

### iNKT1 cells
- Key markers: Eomes, Tbx21, Prf1, Gzma, Nkg7
- Cytotoxic, Th1-like invariant NKT cells. Specialized for IFN-γ production and effector responses against tumors or infections.
- Literature support:
  - Daussy et al., J. Exp. Med. 2014: T-bet and Eomes instruct the development of two distinct natural killer cell lineages in the liver and in the bone marrow
  - Shimizu et al., Commun. Biol. 2019: Eomes transcription factor is required for the development and differentiation of invariant NKT cells

### NK cells
- Key markers: Ncr1, Klrb1c, Klrd1, Tyrobp, Fcgr3
- Innate cytotoxic lymphocytes lacking TCR signaling. Functionally distinct from T and iNKT cells, with strong NK receptor programs.
- Literature support:
  - Chen et al., Signal Transduction and Targeted Therapy. 2024: Comprehensive snapshots of natural killer cells functions, signaling, molecular mechanisms and clinical utilization

### Naive T cells
- Key markers: Satb1, Txk, Cd7, Tcf7, Lef1
- Quiescent or early-differentiation T cells with minimal activation or effector signatures. Likely precursors to activated or memory states.
- Literature support:
  - Willinger et al., J. Immunol. 2006: Human naive CD8 T cells down-regulate expression of the WNT pathway transcription factors lymphoid enhancer binding factor 1 and transcription factor 7 (T cell factor-1) following antigen encounter in vitro and in vivo
  - Aandahl et al., J. Immunol. 2003: CD7 is a differentiation marker that identifies multiple CD8 T cell effector subsets
  - Chen et al., Blood. 2009: Identification of early growth response protein 1 (EGR-1) as a novel target for JUN-induced apoptosis in multiple myeloma

### Central memory T cells
- Key markers: Sell (CD62L), S1pr1, Lef1, Foxp1
- Recirculating memory T cells with lymphoid-homing capacity. Less cytotoxic, optimized for long-term immune surveillance.
- Literature support:
  - Yang et al., PLoS One. 2011: The shedding of CD62L (L-selectin) regulates the acquisition of lytic activity in human tumor-reactive T lymphocytes
  - Sallusto et al., Nat. Rev. Immunol. 2004: Central memory and effector memory T cell subsets: function, generation, and maintenance
  
### IFN-stimulated T cells
- Key markers: Ifit3, Isg15, Oas3, Zbp1
- T cells responding to type I interferon signaling. Reflects cytokine exposure or viral/inflammatory sensing rather than lineage differences.
- Literature support:
  - Crouse et al., Nat. Rev. Immunol. 2015: Regulation of antiviral T cell responses by type I interferons
  - Lim et al., iScience. 2024: Type I interferon signaling regulates myeloid and T cell crosstalk in the glioblastoma tumor microenvironment

### iNKT17 cells
- Key markers: Rorc, Il23r, Il1r1, Cxcr6, Zbtb16
- Pro-inflammatory invariant NKT subset associated with IL-17 biology and tissue residency, particularly in barrier or inflamed tissues.
- Literature support:
  - Giannou et al., Immunity. 2023: Tissue resident iNKT17 cells facilitate cancer cell extravasation in liver metastasis via interleukin-22
  - Thapa et al., Sci. Rep. 2017: The differentiation of ROR-γt expressing iNKT17 cells is orchestrated by Runx1

### iNKT2 cells
- Key markers: Gata3, Il17rb, Il1rl1 (ST2), Areg
- Th2-like invariant NKT cells involved in cytokine support and tissue repair. Amphiregulin suggests regenerative or homeostatic roles.
- Literature support:
  - Krovi &amp; Gapin, Front. Immunol. 2018: Invariant Natural Killer T Cell Subsets—More Than Just Developmental Intermediates
