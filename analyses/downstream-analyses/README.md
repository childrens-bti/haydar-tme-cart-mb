# Downstream Analysis of Murine Brain Tumor CAR-T scRNA-seq

## Usage

`bash run_module.sh`

## Folder contents

- `01-myeloid-condition-umaps.Rmd`: subsets myeloid cells by treatment condition, reprocesses each subset independently, and generates UMAP plots of myeloid subtypes for each condition.
- `02-lineage-plot.Rmd`: generates dot plot of top lineage markers and stacked barplot of lineage proportions across conditions.
- `03-myeloid-functional-marker-dotplots.Rmd`: performs differential expression analysis in myeloid cells by subcluster, subtype, and treatment condition, annotates marker genes with Hallmark pathways, and generates dot plots for visualization.
- `04-myeloid-pseudobulk-gsea.Rmd`: performs pseudobulk DESeq2‑based Hallmark GSEA for myeloid cells across conditions.
- `05-myeloid-composition-analysis.Rmd`: performs myeloid subtype composition analysis using sccomp with 41BB-L and tumor as the baseline references and generates composition plots and stacked bar plots of subtype proportions across conditions.
- `06-tcell-reclustering.Rmd`: performs re-clustering of T cells after additional filtering and evaluates consistency with previous subtype annotations.
- `07-dc-reclustering.Rmd`: performs re-clustering of dendritic cells after additional filtering, evaluates consistency with previous subtype annotations, and runs AUCell to evaluate Hallmark pathway activity across updated DC clusters.
- `08-myeloid-milor-analysis.Rmd`: runs miloR neighborhood differential abundance testing for the myeloid population using tumor as the reference condition.
- `09-tcell-milor-analysis.Rmd`: runs miloR neighborhood differential abundance testing for the refined T cell population using tumor as the reference condition.
- `10-tcell-composition-analysis.Rmd`: performs refined T cell subtype composition analysis using sccomp with 41BB-L and tumor as the baseline references and generates composition plots.
- `11-generate-loupe-files.R`: exports the finalized all-cell, myeloid, and T-cell Seurat objects as Loupe Browser `.cloupe` files, with descriptive annotations and paired PCA/UMAP coordinate tables.

## Analysis module directory structure

```
.
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
├── 08-myeloid-milor-analysis.Rmd
├── 08-myeloid-milor-analysis.html
├── 09-tcell-milor-analysis.Rmd
├── 09-tcell-milor-analysis.html
├── 10-tcell-composition-analysis.Rmd
├── 10-tcell-composition-analysis.html
├── 11-generate-loupe-files.R
├── README.md
├── input
│   ├── cart_lineage_markers.csv
│   ├── myeloid_subcluster_markers.csv
│   └── myeloid_subtype_markers.csv
├── plots
│   ├── dc                 # Dendritic-cell reclustering plots
│   ├── myeloid
│   │   └── GSEA           # Myeloid pseudobulk GSEA and volcano plots
│   └── tcell
├── results
│   ├── dc                 # Dendritic-cell reclustering tables
│   ├── myeloid
│   │   └── GSEA           # Myeloid DESeq2 and Hallmark GSEA tables
│   └── tcell              # Reclustering, GSEA, composition, and miloR tables
├── run_module.sh
└── util
    ├── dotplot_helpers.R
    ├── milor_helpers.R
    ├── pseudobulk_gsea_helpers.R
    └── sccomp_helpers.R
```
