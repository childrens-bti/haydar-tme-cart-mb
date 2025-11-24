#!/bin/bash

set -e
set -o pipefail

# Assessment of ambient RNA contamintation with SoupX visualizations
Rscript -e "rmarkdown::render('01-ambient-assessment.Rmd')"
# Testing best parameters for running SoupX with one sample
Rscript -e "rmarkdown::render('02-one-sample-test.Rmd')"
# SoupX before vs after with all samples
Rscript -e "rmarkdown::render('03-all-sample-test.Rmd')"