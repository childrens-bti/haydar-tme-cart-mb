#!/bin/bash

set -e
set -o pipefail

# Run Seurat QC, processing, and integration on CAR-T scRNA-seq libraries
Rscript -e "rmarkdown::render('01-seurat-processing.Rmd')"
# Cell lineage annotation analysis and plots
Rscript -e "rmarkdown::render('02-lineage-annotation.Rmd')"
# Subclustering to myeloid, lymphoid, and dendritic cells
Rscript -e "rmarkdown::render('03-subclustering.Rmd')"
