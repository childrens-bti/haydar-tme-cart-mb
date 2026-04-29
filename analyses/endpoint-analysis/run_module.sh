#!/bin/bash

set -e
set -o pipefail

# Endpoint samples analysis
Rscript -e "rmarkdown::render('07-endpoint-samples-analysis.Rmd')"