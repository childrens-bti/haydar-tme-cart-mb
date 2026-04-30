#!/bin/bash

set -e
set -o pipefail

# Endpoint samples analysis
Rscript -e "rmarkdown::render('01-endpoint-samples-analysis.Rmd')"