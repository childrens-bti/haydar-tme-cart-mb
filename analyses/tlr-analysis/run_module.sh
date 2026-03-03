#!/bin/bash

set -euo pipefail

# TLR pathway module score analysis and visualization
Rscript -e "rmarkdown::render('01-tlr-analysis.Rmd')"