#!/bin/bash

set -e
set -o pipefail

# Run Seurat QC, processing, and integration on CAR-T scRNA-seq libraries
Rscript -e "rmarkdown::render('01-seurat-processing.Rmd')"
# Cell lienage annotation analysis and plots
Rscript -e "rmarkdown::render('02-lineage-annotation.Rmd')"
