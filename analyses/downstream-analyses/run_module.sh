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

# Prepare annotated CD4-like and CD8-like T cell objects for trajectory analysis
Rscript -e "rmarkdown::render('12-tcell-cd4-cd8-annotation.Rmd')"

# Perform CD4-like and CD8-like T cell trajectory analysis using Slingshot
Rscript -e "rmarkdown::render('13-tcell-cd4-cd8-slingshot.Rmd')"

# Compare Kdm6b expression and define Kdm6b-high and Kdm6b-low subclusters
Rscript -e "rmarkdown::render('14-kdm6b-expression-subcluster-definition.Rmd')"

# Test differential abundance of myeloid subclusters using sccomp
Rscript -e "rmarkdown::render('15-kdm6b-subcluster-abundance.Rmd', params = list(input_rds = 'cart_myeloid_subtypes.rds', ranking_tsv = 'myeloid_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'myeloid'), output_file = '15-kdm6b-myeloid-subcluster-abundance.html')"

# Test differential abundance of T-cell subclusters using sccomp
Rscript -e "rmarkdown::render('15-kdm6b-subcluster-abundance.Rmd', params = list(input_rds = 'cart_tcell_subtypes.rds', ranking_tsv = 'tcell_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'tcell'), output_file = '15-kdm6b-tcell-subcluster-abundance.html')"

# Compare transcriptional profiles of Kdm6b-high and Kdm6b-low myeloid subclusters
Rscript -e "rmarkdown::render('16-kdm6b-transcriptional-comparison.Rmd', params = list(input_rds = 'cart_myeloid_subtypes.rds', ranking_tsv = 'myeloid_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'myeloid'), output_file = '16-kdm6b-myeloid-transcriptional-comparison.html')"

# Compare transcriptional profiles of Kdm6b-high and Kdm6b-low T-cell subclusters
Rscript -e "rmarkdown::render('16-kdm6b-transcriptional-comparison.Rmd', params = list(input_rds = 'cart_tcell_subtypes.rds', ranking_tsv = 'tcell_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'tcell'), output_file = '16-kdm6b-tcell-transcriptional-comparison.html')"

# Correlate Kdm6b expression with immune-program genes across all myeloid cells
Rscript -e "rmarkdown::render('17-kdm6b-myeloid-gene-correlations.Rmd')"

# Summarize myeloid immune-program gene expression and UCell scores
Rscript -e "rmarkdown::render('18-kdm6b-immune-programs.Rmd', params = list(input_rds = 'cart_myeloid_subtypes.rds', ranking_tsv = 'myeloid_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'myeloid'), output_file = '18-kdm6b-myeloid-immune-programs.html')"

# Summarize T-cell immune-program gene expression and UCell scores
Rscript -e "rmarkdown::render('18-kdm6b-immune-programs.Rmd', params = list(input_rds = 'cart_tcell_subtypes.rds', ranking_tsv = 'tcell_kdm6b_subcluster_ranking_and_groups.tsv', subtype_label = 'tcell'), output_file = '18-kdm6b-tcell-immune-programs.html')"
