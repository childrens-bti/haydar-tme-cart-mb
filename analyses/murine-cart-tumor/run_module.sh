#!/bin/bash

set -e
set -o pipefail

# Run Seurat QC, processing, and integration on CAR-T scRNA-seq libraries
Rscript -e "rmarkdown::render('01-seurat-processing.Rmd')"
# Cell lienage annotation analysis and plots
Rscript -e "rmarkdown::render('02-lineage-annotation.Rmd')"
# Subclustering to myeloid, lymphoid, and dendritic cells
Rscript -e "rmarkdown::render('03-subclustering.Rmd')"
# Characterization of myeloid heterogenity
Rscript -e "rmarkdown::render('04-myeloid-characterization.Rmd')"
# Characterization of lymphoid heterogenity
Rscript -e "rmarkdown::render('05-lymphoid-characterization.Rmd')"
