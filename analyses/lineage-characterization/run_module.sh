#!/bin/bash

set -e
set -o pipefail

# Characterization of myeloid heterogeneity
Rscript -e "rmarkdown::render('04-myeloid-characterization.Rmd')"
# Characterization of lymphoid and T cell heterogeneity
Rscript -e "rmarkdown::render('05-lymphoid-characterization.Rmd')"
# Characterization of dendritic cell heterogeneity
Rscript -e "rmarkdown::render('06-dendritic-characterization.Rmd')"
