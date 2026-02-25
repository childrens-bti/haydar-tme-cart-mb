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
# Characterization of lymphoid and T cell heterogenity
Rscript -e "rmarkdown::render('05-lymphoid-characterization.Rmd')"
# Characterization of dendritic cell heterogenity
Rscript -e "rmarkdown::render('06-dendritic-characterization.Rmd')"
# Endpoint samples analysis
Rscript -e "rmarkdown::render('07-endpoint-samples.Rmd')"
# Generate UMAPs for Haydar grant figures
Rscript -e "rmarkdown::render('08-grant-umaps.Rmd')"
# Generate lineage marker dot plots and lineage proportion stacked barplots
Rscript -e "rmarkdown::render('09-lineage-plot.Rmd')"
# Performs differential expression analysis in myeloid cells and generates pathway-annotated dot plots.
Rscript -e "rmarkdown::render('10-myeloid-differential.Rmd')"
