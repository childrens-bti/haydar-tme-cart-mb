#!/bin/bash

set -e
set -o pipefail

# NOTE:
# Downloading the v5 data release is required prior to running this module.

# Generate myeloid UMAPs by treatment condition
Rscript -e "rmarkdown::render('01-myeloid-condition-umaps.Rmd')"

# Generate lineage marker dot plots and lineage proportion stacked barplots
Rscript -e "rmarkdown::render('02-lineage-plot.Rmd')"

# Perform differential expression analysis in myeloid cells and generate pathway-annotated dot plots
Rscript -e "rmarkdown::render('03-myeloid-functional-marker-dotplots.Rmd')"

# Run pseudobulk DESeq2-based Hallmark GSEA for myeloid cells
Rscript -e "rmarkdown::render('04-myeloid-pseudobulk-gsea.Rmd')"

# Perform subtype composition analysis for myeloid cells using sccomp
Rscript -e "rmarkdown::render('05-myeloid-composition-analysis.Rmd')"

# Perform T cell re-clustering after additional filtering
Rscript -e "rmarkdown::render('06-tcell-reclustering.Rmd')"

# Perform DC re-clustering and run AUCell on Hallmark pathways
Rscript -e "rmarkdown::render('07-dc-reclustering.Rmd')"

# Perform T cell trajectory analysis using Slingshot
Rscript -e "rmarkdown::render('08-tcell-trajectory.Rmd')"

# Run miloR neighborhood differential abundance testing for myeloid cells
Rscript -e "rmarkdown::render('09-myeloid-milor-analysis.Rmd')"

# Run miloR neighborhood differential abundance testing for refined T cells
Rscript -e "rmarkdown::render('10-tcell-milor-analysis.Rmd')"

# Perform refined T cell subtype composition analysis using sccomp
Rscript -e "rmarkdown::render('11-tcell-composition-analysis.Rmd')"
