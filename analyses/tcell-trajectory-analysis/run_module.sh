#!/bin/bash

set -e
set -o pipefail

# NOTE:
# Downloading the v5 data release is required prior to running this module.

# Analyze the combined T-cell trajectory and plot lineage and pseudotime trends.
Rscript -e "rmarkdown::render('01-tcell-trajectory.Rmd')"

# Annotate CD4-like and CD8-like cells and save their trajectory inputs.
Rscript -e "rmarkdown::render('02-tcell-cd4-cd8-annotation.Rmd')"

# Infer shared, all-condition Slingshot trajectories for CD4-like and CD8-like cells.
Rscript -e "rmarkdown::render('03-tcell-cd4-cd8-slingshot.Rmd')"

# Fit or reuse the CD4-like tradeSeq model and generate its tests and plots.
TRADESEQ_WORKERS=6 Rscript -e "rmarkdown::render('04-tcell-cd4-cd8-tradeseq.Rmd', params = list(compartment = 'CD4-like'), output_file = '04-tcell-cd4-tradeseq.html')"

# Fit or reuse the CD8-like tradeSeq model and generate its tests and plots.
TRADESEQ_WORKERS=6 Rscript -e "rmarkdown::render('04-tcell-cd4-cd8-tradeseq.Rmd', params = list(compartment = 'CD8-like'), output_file = '04-tcell-cd8-tradeseq.html')"
