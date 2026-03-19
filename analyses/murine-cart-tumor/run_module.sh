#!/bin/bash

set -e
set -o pipefail

# NOTE:
# The processed Seurat objects below were already included in the v4 data release
# and can be used directly for downstream analysis by loading them from data_dir
# with readRDS():
#   - cart_annotated.rds
#   - cart_myeloid_subtypes.rds
#   - cart_tcell_subtypes.rds
#   - cart_dc_subtypes.rds

# Run Seurat QC, processing, and integration on CAR-T scRNA-seq libraries
Rscript -e "rmarkdown::render('01-seurat-processing.Rmd')"
# Cell lineage annotation analysis and plots
Rscript -e "rmarkdown::render('02-lineage-annotation.Rmd')"
# Subclustering to myeloid, lymphoid, and dendritic cells
Rscript -e "rmarkdown::render('03-subclustering.Rmd')"
# Characterization of myeloid heterogeneity
Rscript -e "rmarkdown::render('04-myeloid-characterization.Rmd')"
# Characterization of lymphoid and T cell heterogeneity
Rscript -e "rmarkdown::render('05-lymphoid-characterization.Rmd')"
# Characterization of dendritic cell heterogeneity
Rscript -e "rmarkdown::render('06-dendritic-characterization.Rmd')"
# Endpoint samples analysis
Rscript -e "rmarkdown::render('07-endpoint-samples-analysis.Rmd')"
# Generate UMAPs for Haydar grant figures
Rscript -e "rmarkdown::render('08-grant-umaps.Rmd')"
# Generate lineage marker dot plots and lineage proportion stacked barplots
Rscript -e "rmarkdown::render('09-lineage-plot.Rmd')"
# Perform differential expression analysis in myeloid cells and generates pathway-annotated dot plots
Rscript -e "rmarkdown::render('10-myeloid-functional-marker-dotplots.Rmd')"

# Run Hallmark GSEA analysis for myeloid
Rscript -e "rmarkdown::render(
  '11-subtype-pseudobulk-gsea.Rmd',
  params = list(
    input_rds = 'cart_myeloid_subtypes.rds',
    subtype_label = 'myeloid'
  ),
  output_file = '11-subtype-pseudobulk-gsea-myeloid.html'
)"

# Run Hallmark GSEA analysis for tcell
Rscript -e "rmarkdown::render(
  '11-subtype-pseudobulk-gsea.Rmd',
  params = list(
    input_rds = 'cart_tcell_subtypes.rds',
    subtype_label = 'tcell'
  ),
  output_file = '11-subtype-pseudobulk-gsea-tcell.html'
)"

# Run Hallmark GSEA analysis for dc
Rscript -e "rmarkdown::render(
  '11-subtype-pseudobulk-gsea.Rmd',
  params = list(
    input_rds = 'cart_dc_subtypes.rds',
    subtype_label = 'dc'
  ),
  output_file = '11-subtype-pseudobulk-gsea-dc.html'
)"

# Run subtype composition analysis for myeloid
Rscript -e "rmarkdown::render(
  '12-subtype-composition-analysis.Rmd',
  params = list(
    input_rds = 'cart_myeloid_subtypes.rds',
    subtype_label = 'myeloid'
  ),
  output_file = '12-subtype-composition-analysis-myeloid.html'
)"

# Run subtype composition analysis for tcell
Rscript -e "rmarkdown::render(
  '12-subtype-composition-analysis.Rmd',
  params = list(
    input_rds = 'cart_tcell_subtypes.rds',
    subtype_label = 'tcell'
  ),
  output_file = '12-subtype-composition-analysis-tcell.html'
)"

# Run subtype composition analysis for dc
Rscript -e "rmarkdown::render(
  '12-subtype-composition-analysis.Rmd',
  params = list(
    input_rds = 'cart_dc_subtypes.rds',
    subtype_label = 'dc'
  ),
  output_file = '12-subtype-composition-analysis-dc.html'
)"

