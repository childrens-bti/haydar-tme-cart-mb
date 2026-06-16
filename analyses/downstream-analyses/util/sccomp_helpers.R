# Shared sccomp utilities for downstream composition analyses.
# Assumes the Seurat object metadata contains sample, condition, and the
# requested cell-group column. The helper writes the sccomp test results to
# results_dir, then returns the model input and result table.

options(mc.cores = 1)

# Set deterministic RNG behavior before each sccomp run.
set_seed <- function(seed = 1234) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
}

# Run one condition-based sccomp composition model for a chosen reference level.
# condition_levels determines both the factor order and the baseline condition.
run_sccomp_composition <- function(obj,
                                   cell_group_col,
                                   condition_levels,
                                   analysis_label,
                                   output_suffix,
                                   reference_label,
                                   results_dir,
                                   seed = 1234) {
  set_seed(seed)

  sccomp_input <- obj@meta.data |>
    tibble::rownames_to_column("cell") |>
    dplyr::select(
      cell,
      sample = sample,
      cell_group = all_of(cell_group_col),
      condition
    ) %>%
    dplyr::mutate(
      condition = factor(condition, levels = condition_levels)
    ) %>%
    tibble::as_tibble()

  cat("Reference condition for", reference_label, "sccomp:", levels(sccomp_input$condition)[1], "\n")
  cat("Cell group column used for analysis:", cell_group_col, "\n")

  sccomp_result <- sccomp_input |>
    sccomp::sccomp_estimate(
      formula_composition = ~ condition,
      .sample = sample,
      .cell_group = cell_group,
      bimodal_mean_variability_association = TRUE,
      cores = 1,
      verbose = FALSE
    ) |>
    sccomp::sccomp_test()

  readr::write_tsv(
    sccomp_result,
    file.path(results_dir, paste0(analysis_label, "_composition_sccomp_", output_suffix, "_results.tsv"))
  )

  list(
    input = sccomp_input,
    result = sccomp_result
  )
}
