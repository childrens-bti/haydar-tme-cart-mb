#!/bin/bash

set -e
set -o pipefail

# Characterization of myeloid heterogeneity
Rscript -e "rmarkdown::render('01-myeloid-characterization.Rmd')"
# Characterization of lymphoid and T cell heterogeneity
Rscript -e "rmarkdown::render('02-lymphoid-characterization.Rmd')"
# Characterization of dendritic cell heterogeneity
Rscript -e "rmarkdown::render('03-dendritic-characterization.Rmd')"
