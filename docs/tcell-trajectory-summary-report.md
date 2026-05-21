# T-cell Trajectory Summary Report

**author: Bicna Song**  
**date: 2026-05-20**

---

## What was done

- Loaded the refined T-cell object.
- Removed non-target populations before trajectory analysis:
  - `Antigen-presenting myeloid cells`
  - `NK-like cytotoxic lymphocytes`
- Ran Slingshot pseudotime analysis using the integrated CCA embedding.
- Used `Naive / central-memory T cells` as the trajectory starting population.
- Generated UMAP visualizations for pseudotime and inferred lineage curves.
- Summarized pseudotime ordering by T-cell subtype.
- Screened variable genes for association with pseudotime along each lineage.
- Generated gene trend plots and heatmaps for top pseudotime-associated genes.

---

## Overall interpretation

The trajectory analysis supports a model in which the `Naive / central-memory T cells` population sits near the start of multiple transcriptional continua. These lineages extend toward activated, regulatory, proliferative, cytotoxic, IFN-responsive, gamma-delta-like, and inflammatory T-cell states.

The inferred lineages should be interpreted as **transcriptional state trajectories**, not as direct time-course measurements. Still, the patterns are biologically coherent and broadly match the refined T-cell subtype annotations used as input for this analysis.

### Bottom line

- Slingshot identified **7 inferred T-cell lineages**.
- The trajectory structure is consistent with differentiation/activation paths starting from the naive/central-memory compartment.
- Several lineages show strong marker programs that match expected endpoint states:
  - Treg / activated CD4 program
  - Cycling / proliferating program
  - Cytotoxic / gamma-delta-like program
  - IFN-responsive / inflammatory program
- The updated lineage UMAPs now include clearer direction arrows and direct subtype labels.

---

## How to read the results

The analysis first places cells along lineage-specific pseudotime paths. Lower pseudotime values are closer to the starting population, while higher values are closer to the inferred endpoint for that lineage.

The subtype-level pseudotime summary is useful for comparing lineages. For each lineage, subtypes with higher median pseudotime are interpreted as later states along that path. For example, `slingPseudotime_3` places IFN-responsive T cells late along lineage 3, matching the lineage endpoint seen in the UMAP.

The lineage curve UMAPs are useful for checking the direction and endpoint of each inferred path. The overview plot shows all seven paths together, while the per-lineage plots make it easier to inspect individual trajectories. These figures are saved separately as PDF outputs and are summarized below in the output file list.

The gene trend curves and heatmaps summarize genes associated with pseudotime. Trend curves show how selected genes change from early to late pseudotime, while heatmaps show broader gene programs ordered along each lineage.

---

## Input cell populations

The refined T-cell object contained **8,887 cells** before trajectory filtering.

| T-cell subtype | Cells |
|---|---:|
| Quiescent T cells | 1,322 |
| Treg | 1,303 |
| Activated effector CD8 T cells (stress) | 1,184 |
| NK-like cytotoxic lymphocytes | 824 |
| Cycling / proliferating CAR T cells | 792 |
| Activated CD4 effector / helper-like T cells | 689 |
| Antigen-presenting myeloid cells | 617 |
| Naive / central-memory T cells | 585 |
| Terminally differentiated cytotoxic T cells | 455 |
| gamma-delta T cells | 452 |
| IFN-responsive T cells | 301 |
| gamma-delta / Th17-like T cells | 227 |
| Th2-like inflammatory T cells | 136 |

The two excluded groups were removed only for the trajectory analysis, leaving the T-cell-focused state space for Slingshot.

---

## Slingshot lineages

Seven lineages were inferred from the filtered T-cell object. The number of cells with non-missing pseudotime differs by lineage because Slingshot assigns cells to lineage-specific paths.

| Lineage | Cells with pseudotime | Significant pseudotime-associated genes (`padj < 0.05`) | Top associated genes |
|---|---:|---:|---|
| `slingPseudotime_1` | 4,298 | 942 | Tnfrsf4, Ctla4, Foxp3, Il2ra, Tnfrsf9, S100a4 |
| `slingPseudotime_2` | 3,302 | 985 | S100a4, Bhlhe40, Klf2, Cxcr6, S100a11, Rbpj |
| `slingPseudotime_3` | 3,589 | 784 | Nkg7, Cxcr6, S100a4, S100a6, Isg20, Thy1 |
| `slingPseudotime_4` | 3,828 | 1,339 | Stmn1, Mki67, Top2a, Cdca3, H1f5, Ccna2 |
| `slingPseudotime_5` | 3,529 | 993 | S100a4, Bhlhe40, Thy1, S100a6, Klf2, Tnfrsf9 |
| `slingPseudotime_6` | 3,209 | 942 | Bhlhe40, Klf2, S100a4, Rbpj, Tnfrsf4, Cxcr6 |
| `slingPseudotime_7` | 2,851 | 666 | Ccl5, Nkg7, Serpinb9, Dennd4a, Il12rb2, Klrc1 |

---

## Notes on lineage programs

### Lineage 1

This lineage has strong Treg / activated CD4-associated genes, including `Foxp3`, `Ctla4`, `Il2ra`, `Tnfrsf4`, and `Ikzf2`. This is consistent with movement from the naive/central-memory compartment toward a regulatory or activated CD4-like state.

### Lineage 4

This lineage is dominated by cell-cycle genes such as `Mki67`, `Top2a`, `Cdk1`, `Nusap1`, and `Pclaf`. This likely represents a proliferative trajectory toward cycling CAR T cells.

### Lineages 3 and 7

These lineages show cytotoxic and gamma-delta/NK-like effector-associated genes such as `Nkg7`, `Ccl5`, `Xcl1`, `Ifng`, and `Il12rb2`. These paths likely capture cytotoxic activation or gamma-delta-like effector differentiation.

### Lineages 2, 5, and 6

These lineages share activation and tissue/inflammatory-state genes including `S100a4`, `Bhlhe40`, `Klf2`, `Cxcr6`, `Pdcd1`, and `Tnfrsf9`. They likely represent related but distinct activated T-cell continua.

---

## Output files generated

All output files can be found in the shared Box folder:

**/BTI-Bioinformatics/HaydarLab/10X_Yacoub/results/SR008785/tcell_reclustering/trajectory**

### Result files

- `tcell_slingshot_fast_spearman_pseudotime_gene_screen.csv`  
  Spearman screen for genes associated with pseudotime along each lineage.

### Plot files

- `tcell_slingshot_pseudotime_umap.pdf`  
  UMAP faceted by Slingshot pseudotime values for each inferred lineage.

- `tcell_slingshot_pseudotime_all_lineages_summary.pdf`  
  Dot plot summarizing median pseudotime by T-cell subtype for each lineage.

- `tcell_slingshot_lineage_curves_overview_umap.pdf`  
  Overview UMAP showing all inferred lineage curves, arrows, lineage labels, and T-cell subtype colors.

- `tcell_slingshot_lineage_curves_each_lineage_umap.pdf`  
  One UMAP per lineage. Each page highlights cells assigned to that lineage, uses the same subtype colors as the overview, includes direct subtype labels, and shows arrows indicating trajectory direction.

- `tcell_slingshot_gene_trend_curves_all_lineages.pdf`  
  Smooth gene-expression trend curves for top pseudotime-associated genes.

- `tcell_slingshot_gene_trend_heatmaps_all_lineages.pdf`  
  Heatmaps showing scaled expression of top pseudotime-associated genes ordered along pseudotime.

### Presentation deck

- `tcell-trajectory-summary-deck.pptx`  
  Slide deck summarizing the trajectory analysis for collaborator discussion.

---
