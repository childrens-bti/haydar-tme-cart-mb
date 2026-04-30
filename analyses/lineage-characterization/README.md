# Murine Brain Tumor CAR-T scRNA-seq Analysis - lineage characterization module

## Usage

`bash run_module.sh`

## Folder contents

1. `01-myeloid-characterization.Rmd`: characterization of myeloid/macrophage/microglia subtypes, markers, composition, and pathway activity.
2. `02-lymphoid-characterization.Rmd`: characterization of lymphoid/T cell subsets with additional filtering.
3. `03-dendritic-characterization.Rmd`: characterization of dendritic cell subsets.
4. `plots/`: output plots from lineage characterization.
5. `results/`: marker tables and GSEA outputs (excluding gitignored artifacts).

## Analysis module directory structure
```
lineage-characterization/
├── 01-myeloid-characterization.Rmd
├── 01-myeloid-characterization.html
├── 02-lymphoid-characterization.Rmd
├── 02-lymphoid-characterization.html
├── 03-dendritic-characterization.Rmd
├── 03-dendritic-characterization.html
├── README.md
├── run_module.sh
├── plots
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
│   └── tcell_subtypes_umap.png
└── results
	├── dc_subcluster_markers.csv
	├── dc_subtype_c5_gsea_results.csv
	├── dc_subtype_c7_gsea_results.csv
	├── dc_subtype_hallmark_gsea_results.csv
	├── dc_subtype_markers.csv
	├── myeloid_subcluster_markers.csv
	├── myeloid_subtype_c5_gsea_results.csv
	├── myeloid_subtype_c7_gsea_results.csv
	├── myeloid_subtype_hallmark_gsea_results.csv
	├── myeloid_subtype_markers.csv
	├── nkt_subcluster_markers.csv
	├── tcell_subtype_c5_gsea_results.csv
	├── tcell_subtype_c7_gsea_results.csv
	├── tcell_subtype_hallmark_gsea_results.csv
	└── tcell_subtype_markers.csv
```