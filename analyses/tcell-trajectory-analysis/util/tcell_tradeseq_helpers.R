# Helper functions for the CD4-like and CD8-like tradeSeq analysis.
# Author: Bicna Song
# Date: 2026-09-29

# Read raw RNA counts from either a Seurat v5 layer or an older Seurat slot.
get_raw_rna_counts <- function(obj) {
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) {
      SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
    }
  )
}

# Load one shared Slingshot cache and return strictly cell-aligned tradeSeq input.
# Tumor cells are retained in the cache but excluded from the modeled conditions.
prepare_tradeseq_input <- function(cache_file, compartment, model_conditions) {
  if (!file.exists(cache_file)) {
    stop("Shared Slingshot cache was not found: ", cache_file)
  }

  cache <- readRDS(cache_file)
  required_fields <- c(
    "obj", "curve_weights", "trade_seq_pseudotime", "lineage_assignments"
  )
  missing_fields <- setdiff(required_fields, names(cache))
  if (length(missing_fields) > 0) {
    stop(
      compartment, " cache is missing required fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  obj <- cache$obj
  required_metadata <- c("condition", "sample")
  missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
  if (length(missing_metadata) > 0) {
    stop(
      compartment, " object is missing metadata: ",
      paste(missing_metadata, collapse = ", ")
    )
  }

  counts <- get_raw_rna_counts(obj)
  variable_feature_info <- SeuratObject::HVFInfo(obj, assay = "RNA")
  if (is.null(rownames(variable_feature_info)) ||
      !"variance.standardized" %in% colnames(variable_feature_info)) {
    stop(compartment, " object is missing standardized variable-feature scores")
  }
  variable_feature_score <- stats::setNames(
    variable_feature_info$variance.standardized,
    rownames(variable_feature_info)
  )
  pseudotime <- as.matrix(cache$trade_seq_pseudotime)
  cell_weights <- as.matrix(cache$curve_weights)
  cell_order <- rownames(pseudotime)

  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop(compartment, " raw RNA counts must have unique gene names")
  }
  if (is.null(cell_order) || anyDuplicated(cell_order)) {
    stop(compartment, " pseudotime must have unique cell names")
  }
  if (!identical(rownames(pseudotime), rownames(cell_weights))) {
    stop(compartment, " pseudotime and curve weights have different cell order")
  }
  if (!identical(colnames(pseudotime), colnames(cell_weights))) {
    stop(compartment, " pseudotime and curve weights have different lineages")
  }
  if (!identical(dim(pseudotime), dim(cell_weights))) {
    stop(compartment, " pseudotime and curve weights have different dimensions")
  }
  if (ncol(pseudotime) != 3) {
    stop(compartment, " shared cache must contain the expected three lineages")
  }
  if (!setequal(cell_order, colnames(counts)) ||
      !setequal(cell_order, rownames(obj@meta.data))) {
    stop(compartment, " cell names differ among counts, metadata, and Slingshot input")
  }
  if (any(!is.finite(pseudotime)) || any(!is.finite(cell_weights))) {
    stop(compartment, " pseudotime or curve weights contain non-finite values")
  }
  if (any(cell_weights < 0) || any(rowSums(cell_weights) <= 0)) {
    stop(compartment, " every modeled cell must have positive Slingshot weight")
  }

  counts <- counts[, cell_order, drop = FALSE]
  metadata <- obj@meta.data[cell_order, , drop = FALSE]
  condition <- as.character(metadata$condition)
  sample <- as.character(metadata$sample)

  if (anyNA(condition) || anyNA(sample)) {
    stop(compartment, " condition and sample metadata cannot contain missing values")
  }

  missing_conditions <- setdiff(model_conditions, unique(condition))
  if (length(missing_conditions) > 0) {
    stop(
      compartment, " cache is missing modeled conditions: ",
      paste(missing_conditions, collapse = ", ")
    )
  }

  model_cells <- cell_order[condition %in% model_conditions]
  model_condition <- factor(
    condition[match(model_cells, cell_order)],
    levels = model_conditions
  )
  model_sample <- sample[match(model_cells, cell_order)]

  list(
    compartment = compartment,
    cache_file = normalizePath(cache_file, mustWork = TRUE),
    cache_md5 = unname(tools::md5sum(cache_file)),
    counts = counts[, model_cells, drop = FALSE],
    variable_feature_score = variable_feature_score[rownames(counts)],
    pseudotime = pseudotime[model_cells, , drop = FALSE],
    cell_weights = cell_weights[model_cells, , drop = FALSE],
    condition = model_condition,
    sample = model_sample,
    cells = model_cells,
    lineage_names = colnames(pseudotime),
    n_cells_all_conditions = length(cell_order),
    excluded_conditions = sort(setdiff(unique(condition), model_conditions))
  )
}

# Select a manageable, outcome-independent gene set for the condition-aware
# model. The highest upstream variable-feature scores provide broad discovery,
# and every priority marker that passes the expression filter is added.
select_balanced_model_genes <- function(
  gene_qc,
  variable_feature_score,
  priority_programs,
  n_variable_genes = 3000L
) {
  if (n_variable_genes < 1) {
    stop("n_variable_genes must be at least 1")
  }
  if (is.null(names(variable_feature_score))) {
    stop("variable_feature_score must be named by gene")
  }

  priority_lookup <- purrr::imap_dfr(
    priority_programs,
    ~ tibble::tibble(program = .y, gene = .x)
  ) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      priority_program = paste(sort(unique(program)), collapse = "|"),
      .groups = "drop"
    )

  selection <- gene_qc %>%
    dplyr::mutate(
      variable_feature_score = as.numeric(variable_feature_score[gene])
    ) %>%
    dplyr::left_join(priority_lookup, by = "gene")

  variable_genes <- selection %>%
    dplyr::filter(passed_filter) %>%
    dplyr::arrange(
      dplyr::desc(variable_feature_score),
      dplyr::desc(n_cells_detected),
      dplyr::desc(count_variance),
      gene
    ) %>%
    dplyr::slice_head(n = n_variable_genes) %>%
    dplyr::pull(gene)

  priority_genes <- selection %>%
    dplyr::filter(passed_filter, !is.na(priority_program)) %>%
    dplyr::pull(gene)
  model_genes <- unique(c(variable_genes, priority_genes))

  selection <- selection %>%
    dplyr::mutate(
      selected_as_variable = gene %in% variable_genes,
      selected_as_priority = gene %in% priority_genes,
      selected_for_model = gene %in% model_genes,
      selection_reason = dplyr::case_when(
        selected_as_variable & selected_as_priority ~ "variable_and_priority",
        selected_as_variable ~ "variable_feature",
        selected_as_priority ~ "priority_marker",
        !passed_filter ~ "failed_expression_filter",
        TRUE ~ "not_selected"
      )
    )

  list(
    genes = model_genes,
    selection = selection,
    summary = tibble::tibble(
      n_genes_passed_filter = sum(selection$passed_filter),
      n_variable_genes_requested = as.integer(n_variable_genes),
      n_variable_genes_selected = length(variable_genes),
      n_priority_genes_added = sum(priority_genes %in% setdiff(model_genes, variable_genes)),
      n_genes_modeled = length(model_genes),
      selection_rule = paste0(
        "top ", n_variable_genes,
        " upstream standardized variable-feature scores plus all priority ",
        "markers passing the expression filter"
      )
    )
  )
}

# Filter genes deterministically by detection count and non-zero count variance.
filter_tradeseq_genes <- function(counts, min_cells = 10L) {
  if (min_cells < 1) {
    stop("min_cells must be at least 1")
  }

  n_detected <- Matrix::rowSums(counts > 0)
  mean_count <- Matrix::rowMeans(counts)
  mean_square <- Matrix::rowMeans(counts ^ 2)
  count_variance <- pmax(mean_square - mean_count ^ 2, 0)
  keep <- n_detected >= min_cells & count_variance > 0

  gene_qc <- tibble::tibble(
    gene = rownames(counts),
    n_cells_detected = as.integer(n_detected),
    count_variance = as.numeric(count_variance),
    passed_filter = keep,
    filter_reason = dplyr::case_when(
      n_detected < min_cells ~ paste0("detected_in_fewer_than_", min_cells, "_cells"),
      count_variance <= 0 ~ "zero_count_variance",
      TRUE ~ "passed"
    )
  )

  list(
    genes = gene_qc$gene[gene_qc$passed_filter],
    gene_qc = gene_qc,
    summary = tibble::tibble(
      n_genes_input = nrow(counts),
      n_genes_passed_filter = sum(keep),
      filter_rule = paste0(
        "raw RNA count > 0 in at least ", min_cells,
        " cells and non-zero count variance"
      )
    )
  )
}

# Select a fixed evaluateK gene set from the most broadly detected filtered genes.
# Ordering by detection, then variance, then gene name avoids another random sample.
select_evaluate_k_genes <- function(gene_qc, n_genes = 500L) {
  gene_qc %>%
    dplyr::filter(passed_filter) %>%
    dplyr::arrange(
      dplyr::desc(n_cells_detected),
      dplyr::desc(count_variance),
      gene
    ) %>%
    dplyr::slice_head(n = n_genes) %>%
    dplyr::pull(gene)
}

# Summarize the gene-level AIC matrix returned by tradeSeq::evaluateK().
summarize_evaluate_k <- function(aic_matrix) {
  knot_count <- as.integer(sub("^k: ", "", colnames(aic_matrix)))
  mean_aic <- colMeans(aic_matrix, na.rm = TRUE)
  se_aic <- apply(aic_matrix, 2, function(x) {
    stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  })
  finite_rows <- apply(aic_matrix, 1, function(x) all(is.finite(x)))
  best_k <- rep(NA_integer_, nrow(aic_matrix))
  best_k[finite_rows] <- knot_count[
    max.col(-aic_matrix[finite_rows, , drop = FALSE], ties.method = "first")
  ]

  tibble::tibble(
    nknots = knot_count,
    mean_aic = mean_aic,
    se_aic = se_aic,
    delta_mean_aic = mean_aic - min(mean_aic, na.rm = TRUE),
    n_genes_with_finite_aic = colSums(is.finite(aic_matrix)),
    n_genes_best = vapply(knot_count, function(k) sum(best_k == k, na.rm = TRUE), integer(1))
  )
}

# Choose the smallest knot count whose mean AIC is within two units of the best.
choose_trade_seq_k <- function(aic_summary, aic_tolerance = 2) {
  supported <- aic_summary$nknots[
    is.finite(aic_summary$delta_mean_aic) &
      aic_summary$delta_mean_aic <= aic_tolerance
  ]
  if (length(supported) == 0) {
    stop("No candidate knot count has a finite mean AIC")
  }
  min(supported)
}

# Draw a compact evaluateK selection plot from the saved AIC summary.
plot_evaluate_k_summary <- function(aic_summary, selected_k, compartment) {
  ggplot2::ggplot(aic_summary, ggplot2::aes(x = nknots, y = mean_aic)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_aic - se_aic, ymax = mean_aic + se_aic),
      width = 0.12
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(
      xintercept = selected_k,
      linetype = "dashed",
      color = "#B2182B"
    ) +
    ggplot2::scale_x_continuous(breaks = aic_summary$nknots) +
    ggplot2::labs(
      title = paste(compartment, "tradeSeq knot selection"),
      subtitle = paste("Selected", selected_k, "knots"),
      x = "Number of knots",
      y = "Mean AIC across evaluation genes"
    ) +
    ggplot2::theme_bw()
}

# Return a BiocParallel backend controlled by the TRADESEQ_WORKERS parameter.
make_tradeseq_bpparam <- function(workers, seed) {
  if (workers > 1) {
    BiocParallel::MulticoreParam(
      workers = workers,
      RNGseed = seed,
      progressbar = TRUE
    )
  } else {
    BiocParallel::SerialParam(RNGseed = seed, progressbar = TRUE)
  }
}

# Calculate the same TMM library-size offset used by tradeSeq from the complete
# count matrix. Computing it before gene subsetting preserves normalization
# while allowing fitGAM to receive only the genes that will actually be fitted.
calculate_tradeseq_offset <- function(counts) {
  norm_factors <- try(edgeR::calcNormFactors(counts), silent = TRUE)
  if (inherits(norm_factors, "try-error")) {
    message(
      "TMM normalization failed. Using unnormalized library sizes as offset."
    )
    norm_factors <- rep(1, ncol(counts))
  }

  library_sizes <- as.numeric(Matrix::colSums(counts)) * norm_factors
  offset <- log(library_sizes)
  if (any(library_sizes == 0)) {
    message("Some library sizes are zero. Offsetting these to 1.")
    offset[library_sizes == 0] <- 0
  }
  offset
}

# Compare cache contexts exactly so changed inputs or parameters trigger refitting.
tradeseq_context_matches <- function(observed, expected) {
  !is.null(observed) && identical(observed, expected)
}

# List cache-context fields that are missing or differ from the current run.
tradeseq_context_differences <- function(observed, expected) {
  if (is.null(observed)) {
    return(names(expected))
  }

  fields <- union(names(observed), names(expected))
  fields[!vapply(
    fields,
    function(field) identical(observed[[field]], expected[[field]]),
    logical(1)
  )]
}

# Run evaluateK, or reuse an AIC matrix produced from the same frozen inputs.
get_or_run_evaluate_k <- function(
  input,
  evaluation_genes,
  candidate_knots,
  seed,
  bpparam,
  output_file,
  force = FALSE
) {
  context <- list(
    cache_md5 = input$cache_md5,
    cells = input$cells,
    genes = evaluation_genes,
    conditions = levels(input$condition),
    candidate_knots = as.integer(candidate_knots),
    seed = as.integer(seed),
    tradeSeq_version = as.character(utils::packageVersion("tradeSeq"))
  )

  if (!force) {
    if (file.exists(output_file)) {
      cached <- readRDS(output_file)
      cached_context <- attr(cached, "tcell_tradeseq_context")
      if (tradeseq_context_matches(cached_context, context)) {
        message("Reusing cached evaluateK result: ", output_file)
        return(cached)
      }

      changed_fields <- tradeseq_context_differences(cached_context, context)
      message(
        "Ignoring incompatible evaluateK cache: ", output_file,
        "\nChanged context fields: ", paste(changed_fields, collapse = ", ")
      )
    } else {
      message("No evaluateK cache found: ", output_file)
    }
  } else {
    message("force_refit is TRUE; evaluateK cache reuse is disabled")
  }

  message(
    "Running evaluateK for ", input$compartment, " with ",
    length(evaluation_genes), " genes and ",
    BiocParallel::bpworkers(bpparam), " workers"
  )
  set.seed(seed)
  aic_matrix <- tradeSeq::evaluateK(
    counts = input$counts[evaluation_genes, , drop = FALSE],
    pseudotime = input$pseudotime,
    cellWeights = input$cell_weights,
    conditions = input$condition,
    k = candidate_knots,
    nGenes = length(evaluation_genes),
    plot = FALSE,
    verbose = TRUE,
    parallel = BiocParallel::bpworkers(bpparam) > 1,
    BPPARAM = bpparam
  )
  attr(aic_matrix, "tcell_tradeseq_context") <- context
  saveRDS(aic_matrix, output_file)
  aic_matrix
}

# Fit a condition-aware tradeSeq GAM, or reuse a model with an identical context.
get_or_fit_tradeseq_model <- function(
  input,
  model_genes,
  nknots,
  seed,
  bpparam,
  output_file,
  cache_files = output_file,
  force = FALSE
) {
  context <- list(
    cache_md5 = input$cache_md5,
    cells = input$cells,
    genes = model_genes,
    conditions = levels(input$condition),
    nknots = as.integer(nknots),
    seed = as.integer(seed),
    tradeSeq_version = as.character(utils::packageVersion("tradeSeq"))
  )

  if (!force) {
    cache_files <- unique(as.character(cache_files))
    for (cache_file in cache_files) {
      if (!file.exists(cache_file)) {
        message("No fitted tradeSeq model cache found: ", cache_file)
        next
      }

      cached <- readRDS(cache_file)
      cached_context <- S4Vectors::metadata(cached)$tcell_tradeseq_context
      if (tradeseq_context_matches(cached_context, context)) {
        message("Reusing cached tradeSeq model: ", cache_file)
        return(cached)
      }

      changed_fields <- tradeseq_context_differences(cached_context, context)
      message(
        "Ignoring incompatible fitted tradeSeq model cache: ", cache_file,
        "\nChanged context fields: ", paste(changed_fields, collapse = ", ")
      )
    }
  } else {
    message("force_refit is TRUE; fitted tradeSeq model cache reuse is disabled")
  }

  missing_model_genes <- setdiff(model_genes, rownames(input$counts))
  if (length(missing_model_genes) > 0) {
    stop(
      "Model genes are absent from the count matrix: ",
      paste(missing_model_genes, collapse = ", ")
    )
  }

  # tradeSeq converts its input counts to a dense matrix during fitting. Supply
  # only selected genes to avoid densifying the complete transcriptome, while
  # retaining the full-count TMM offset used by the original fit.
  model_counts <- input$counts[model_genes, , drop = FALSE]
  offset <- calculate_tradeseq_offset(input$counts)

  message(
    "Fitting ", input$compartment, " tradeSeq model with ",
    nrow(model_counts), " genes, ", ncol(model_counts), " cells, and ",
    BiocParallel::bpworkers(bpparam), " workers"
  )

  set.seed(seed)
  model <- tradeSeq::fitGAM(
    counts = model_counts,
    pseudotime = input$pseudotime,
    cellWeights = input$cell_weights,
    conditions = input$condition,
    offset = offset,
    nknots = nknots,
    verbose = TRUE,
    parallel = BiocParallel::bpworkers(bpparam) > 1,
    BPPARAM = bpparam,
    sce = TRUE
  )
  S4Vectors::metadata(model)$tcell_tradeseq_context <- context
  saveRDS(model, output_file)
  model
}

# Extract per-gene convergence flags from a fitted tradeSeq SingleCellExperiment.
get_tradeseq_convergence <- function(model, compartment) {
  converged <- SummarizedExperiment::rowData(model)$tradeSeq$converged
  tibble::tibble(
    compartment = compartment,
    gene = rownames(model),
    model_converged = as.logical(converged)
  )
}

# Stop before testing when overall convergence or priority-marker coverage is poor.
assert_tradeseq_fit_quality <- function(
  convergence,
  priority_status,
  min_convergence_fraction = 0.8,
  max_priority_failure_fraction = 0.5
) {
  convergence_fraction <- mean(convergence$model_converged, na.rm = TRUE)
  if (!is.finite(convergence_fraction) ||
      convergence_fraction < min_convergence_fraction) {
    stop(
      "Only ", round(100 * convergence_fraction, 1),
      "% of tradeSeq genes converged; statistical testing was stopped"
    )
  }

  modeled_priority <- priority_status %>% dplyr::filter(in_model)
  if (nrow(modeled_priority) > 0) {
    failure_fraction <- mean(!modeled_priority$model_converged)
    if (failure_fraction > max_priority_failure_fraction) {
      stop(
        round(100 * failure_fraction, 1),
        "% of modeled priority markers failed; statistical testing was stopped"
      )
    }
  }
  invisible(TRUE)
}

# Convert lineage numbers used by tradeSeq to the names stored by Slingshot.
make_lineage_map <- function(lineage_names) {
  tibble::tibble(
    lineage_index = seq_along(lineage_names),
    lineage = lineage_names
  )
}

# Restore condition labels after tradeSeq converts them to syntactic R names.
canonicalize_tradeseq_conditions <- function(condition, condition_levels) {
  syntactic_levels <- make.names(condition_levels)
  trade_seq_levels <- sub("^X(?=[0-9])", "", syntactic_levels, perl = TRUE)
  condition_lookup <- c(
    stats::setNames(condition_levels, condition_levels),
    stats::setNames(condition_levels, syntactic_levels),
    stats::setNames(condition_levels, trade_seq_levels)
  )
  canonical <- unname(condition_lookup[condition])
  canonical[is.na(canonical)] <- condition[is.na(canonical)]
  canonical
}

# Normalize lineage-and-condition associationTest columns to one row per test.
normalize_association_test <- function(
  result,
  compartment,
  lineage_map,
  condition_levels,
  convergence
) {
  p_columns <- grep(
    "^pvalue_lineage[0-9]+_condition",
    colnames(result),
    value = TRUE
  )
  if (length(p_columns) == 0) {
    stop("associationTest did not return lineage-by-condition p-values")
  }

  normalized <- purrr::map_dfr(p_columns, function(p_column) {
    suffix <- sub("^pvalue_", "", p_column)
    parsed <- stringr::str_match(suffix, "^lineage([0-9]+)_condition(.+)$")
    statistic_column <- paste0("waldStat_", suffix)
    mean_log_fc_global <- if ("meanLogFC" %in% colnames(result)) {
      as.numeric(result[, "meanLogFC"])
    } else {
      rep(NA_real_, nrow(result))
    }
    condition <- canonicalize_tradeseq_conditions(
      parsed[1, 3],
      condition_levels
    )

    tibble::tibble(
      compartment = compartment,
      lineage_index = as.integer(parsed[1, 2]),
      gene = rownames(result),
      test_type = "association",
      comparison = condition,
      statistic = as.numeric(result[, statistic_column]),
      effect_direction = NA_character_,
      effect_size = NA_real_,
      mean_log_fc_global = mean_log_fc_global,
      p_value = as.numeric(result[, p_column]),
      inference_level = "cell-level exploratory"
    )
  }) %>%
    dplyr::left_join(lineage_map, by = "lineage_index") %>%
    dplyr::left_join(convergence, by = c("compartment", "gene")) %>%
    dplyr::group_by(compartment, test_type, lineage, comparison) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "BH"),
      fdr_family = paste(compartment, test_type, lineage, comparison, sep = "|"),
      test_status = dplyr::case_when(
        is.na(model_converged) ~ "model_status_unknown",
        !model_converged ~ "model_not_converged",
        !is.finite(statistic) | !is.finite(p_value) ~ "not_testable",
        TRUE ~ "tested"
      ),
      test_status_note = dplyr::case_when(
        test_status == "model_status_unknown" ~
          "No model-convergence status was available",
        test_status == "model_not_converged" ~
          "The gene-level GAM did not converge",
        test_status == "not_testable" ~
          "tradeSeq did not return a finite statistic and p-value",
        TRUE ~ "Finite statistic and p-value"
      )
    ) %>%
    dplyr::ungroup()

  normalized %>%
    dplyr::select(
      compartment, lineage, gene, test_type, comparison, statistic,
      effect_direction, effect_size, mean_log_fc_global, p_value, p_adj,
      fdr_family, inference_level, model_converged, test_status,
      test_status_note
    )
}

# Normalize conditionTest pairwise columns and map numeric condition IDs to labels.
normalize_condition_test <- function(
  result,
  compartment,
  lineage_map,
  condition_levels,
  comparison_manifest,
  convergence
) {
  p_columns <- grep(
    "^pvalue_lineage[0-9]+_conds[0-9]+vs[0-9]+$",
    colnames(result),
    value = TRUE
  )
  if (length(p_columns) == 0) {
    stop("conditionTest did not return lineage-specific pairwise p-values")
  }

  normalized <- purrr::map_dfr(p_columns, function(p_column) {
    suffix <- sub("^pvalue_", "", p_column)
    parsed <- stringr::str_match(
      suffix,
      "^lineage([0-9]+)_conds([0-9]+)vs([0-9]+)$"
    )
    condition_a <- condition_levels[as.integer(parsed[1, 3])]
    condition_b <- condition_levels[as.integer(parsed[1, 4])]

    tibble::tibble(
      compartment = compartment,
      lineage_index = as.integer(parsed[1, 2]),
      gene = rownames(result),
      test_type = "condition",
      comparison = paste(condition_a, condition_b, sep = "_vs_"),
      statistic = as.numeric(result[, paste0("waldStat_", suffix)]),
      effect_direction = NA_character_,
      effect_size = NA_real_,
      mean_log_fc_global = NA_real_,
      p_value = as.numeric(result[, p_column]),
      inference_level = "cell-level exploratory"
    )
  }) %>%
    dplyr::left_join(lineage_map, by = "lineage_index") %>%
    dplyr::inner_join(
      comparison_manifest %>%
        dplyr::filter(!is.na(condition_a), !is.na(condition_b)) %>%
        dplyr::select(comparison_id, priority, planned_use),
      by = c("comparison" = "comparison_id")
    ) %>%
    dplyr::left_join(convergence, by = c("compartment", "gene")) %>%
    dplyr::group_by(compartment, test_type, lineage, comparison) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "BH"),
      fdr_family = paste(compartment, test_type, lineage, comparison, sep = "|"),
      test_status = dplyr::case_when(
        is.na(model_converged) ~ "model_status_unknown",
        !model_converged ~ "model_not_converged",
        !is.finite(statistic) | !is.finite(p_value) ~ "not_testable",
        TRUE ~ "tested"
      ),
      test_status_note = dplyr::case_when(
        test_status == "model_status_unknown" ~
          "No model-convergence status was available",
        test_status == "model_not_converged" ~
          "The gene-level GAM did not converge",
        test_status == "not_testable" ~
          "tradeSeq did not return a finite statistic and p-value",
        TRUE ~ "Finite statistic and p-value"
      )
    ) %>%
    dplyr::ungroup()

  normalized %>%
    dplyr::select(
      compartment, lineage, gene, test_type, comparison, statistic,
      effect_direction, effect_size, mean_log_fc_global, p_value, p_adj,
      fdr_family, inference_level, model_converged, test_status,
      test_status_note, priority, planned_use
    )
}

# Build one signed balanced-gene ranking per lineage from fitted start/end values.
# The score combines fitted-change magnitude and association-test strength.
build_tradeseq_rankings <- function(model, association, compartment, lineage_map) {
  predictions <- tradeSeq::predictSmooth(
    model,
    gene = rownames(model),
    nPoints = 2,
    tidy = TRUE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(lineage_index = as.integer(lineage)) %>%
    dplyr::group_by(gene, lineage_index, condition) %>%
    dplyr::summarise(
      start = yhat[which.min(time)],
      end = yhat[which.max(time)],
      .groups = "drop"
    ) %>%
    dplyr::mutate(end_minus_start = log2((end + 0.1) / (start + 0.1))) %>%
    dplyr::group_by(gene, lineage_index) %>%
    dplyr::summarise(
      end_minus_start = mean(end_minus_start, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(lineage_map, by = "lineage_index")

  strongest_association <- association %>%
    dplyr::filter(test_status == "tested", is.finite(statistic)) %>%
    dplyr::group_by(gene, lineage) %>%
    dplyr::slice_max(statistic, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      gene,
      lineage,
      strength_condition = comparison,
      wald_statistic = statistic,
      p_value,
      p_adj
    )

  predictions %>%
    dplyr::left_join(strongest_association, by = c("gene", "lineage")) %>%
    dplyr::mutate(
      compartment = compartment,
      rank_score = end_minus_start * sqrt(pmax(wald_statistic, 0)),
      rank_magnitude = abs(rank_score),
      ranking_status = dplyr::if_else(
        is.finite(rank_score),
        "ranked",
        "no_testable_association"
      ),
      ranking_method = paste0(
        "mean condition-specific log2 fitted end/start * ",
        "sqrt(max lineage-condition association Wald statistic)"
      )
    ) %>%
    dplyr::select(
      compartment, lineage, gene, rank_score, rank_magnitude,
      end_minus_start, wald_statistic, p_value, p_adj,
      strength_condition, ranking_status, ranking_method
    ) %>%
    dplyr::arrange(lineage, dplyr::desc(rank_magnitude), gene)
}

# Plot fitted priority-marker curves with conditions overlaid and lineages faceted.
plot_priority_gene_program <- function(
  model,
  genes,
  program,
  compartment,
  lineage_map,
  condition_levels,
  output_file,
  n_points = 100L
) {
  genes <- intersect(genes, rownames(model))
  if (length(genes) == 0) {
    return(invisible(FALSE))
  }

  prediction <- tradeSeq::predictSmooth(
    model,
    gene = genes,
    nPoints = n_points,
    tidy = TRUE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      gene,
      lineage_index = as.integer(lineage),
      condition = canonicalize_tradeseq_conditions(
        as.character(condition),
        condition_levels
      ),
      time,
      yhat
    ) %>%
    dplyr::left_join(lineage_map, by = "lineage_index")

  plot <- ggplot2::ggplot(
    prediction,
    ggplot2::aes(x = time, y = log1p(yhat), color = condition)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_grid(gene ~ lineage, scales = "free") +
    ggplot2::labs(
      title = paste(compartment, program, "fitted expression"),
      subtitle = "Condition-aware tradeSeq smoothers on the shared trajectory",
      x = "Pseudotime",
      y = "log1p fitted mean count",
      color = "Condition"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    output_file,
    plot,
    width = 11,
    height = max(4, 1.25 * length(genes) + 2)
  )
  invisible(TRUE)
}

# Plot condition-averaged fitted expression on a shared relative-pseudotime grid.
plot_top_dynamic_heatmaps <- function(
  model,
  rankings,
  compartment,
  lineage_map,
  output_dir,
  n_genes = 40L,
  n_points = 100L
) {
  top_genes <- rankings %>%
    dplyr::filter(
      ranking_status == "ranked",
      is.finite(p_adj),
      p_adj <= 0.05,
      is.finite(rank_magnitude)
    ) %>%
    dplyr::group_by(lineage) %>%
    dplyr::arrange(
      dplyr::desc(rank_magnitude),
      p_adj,
      gene,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = n_genes) %>%
    dplyr::ungroup()

  if (nrow(top_genes) == 0) {
    warning(compartment, " has no testable FDR-significant genes for heatmaps")
    return(invisible(tibble::tibble()))
  }

  raw_prediction <- tradeSeq::predictSmooth(
    model,
    gene = unique(top_genes$gene),
    nPoints = n_points,
    tidy = TRUE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      gene,
      lineage_index = as.integer(lineage),
      condition,
      time,
      yhat
    ) %>%
    dplyr::left_join(lineage_map, by = "lineage_index") %>%
    dplyr::inner_join(top_genes, by = c("lineage", "gene"))

  n_conditions_expected <- dplyr::n_distinct(raw_prediction$condition)
  prediction <- raw_prediction %>%
    dplyr::group_by(gene, lineage, condition) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      grid_index = dplyr::row_number(),
      relative_pseudotime = if (dplyr::n() > 1) {
        (grid_index - 1) / (dplyr::n() - 1)
      } else {
        0
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene, lineage, grid_index, relative_pseudotime) %>%
    dplyr::summarise(
      yhat = mean(yhat, na.rm = TRUE),
      n_conditions = dplyr::n_distinct(condition),
      .groups = "drop"
    ) %>%
    dplyr::group_by(gene, lineage) %>%
    dplyr::mutate(
      log_yhat = log1p(yhat),
      z_score = if (stats::sd(log_yhat) > 0) {
        (log_yhat - mean(log_yhat)) / stats::sd(log_yhat)
      } else {
        0
      },
      peak_time = relative_pseudotime[which.max(yhat)]
    ) %>%
    dplyr::ungroup()

  if (any(prediction$n_conditions != n_conditions_expected)) {
    stop("Condition averaging did not retain every modeled condition")
  }

  purrr::walk(unique(prediction$lineage), function(lineage_name) {
    plot_data <- prediction %>%
      dplyr::filter(lineage == lineage_name) %>%
      dplyr::mutate(gene = stats::reorder(gene, peak_time, FUN = mean))
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = relative_pseudotime, y = gene, fill = z_score)
    ) +
      ggplot2::geom_raster(interpolate = FALSE) +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        title = paste(compartment, lineage_name, "top dynamic genes"),
        subtitle = paste0(
          "FDR-significant genes ranked by fitted-change magnitude and test strength; ",
          "condition-averaged on a shared grid"
        ),
        x = "Relative pseudotime",
        y = NULL,
        fill = "Row z-score"
      ) +
      ggplot2::theme_bw()

    ggplot2::ggsave(
      file.path(
        output_dir,
        paste0(
          "tcell_", safe_label(compartment, lowercase = TRUE), "_",
          safe_label(lineage_name, lowercase = TRUE), "_top_dynamic_gene_heatmap.pdf"
        )
      ),
      plot,
      width = 8,
      height = 10
    )
  })
  invisible(prediction)
}

# Plot descriptive sample-specific marker trends using observed log-normalized counts.
# This sensitivity view does not convert the cell-level tradeSeq test to sample-level inference.
plot_sample_marker_sensitivity <- function(
  input,
  genes,
  compartment,
  lineage_map,
  output_file,
  min_group_cells = 30L,
  n_bins = 10L
) {
  if (min_group_cells < 2 || n_bins < 2) {
    stop("min_group_cells and n_bins must both be at least 2")
  }
  genes <- intersect(genes, rownames(input$counts))
  if (length(genes) == 0) {
    return(invisible(FALSE))
  }

  primary_lineage <- max.col(input$cell_weights, ties.method = "first")
  library_size <- Matrix::colSums(input$counts)
  expression <- as.matrix(input$counts[genes, , drop = FALSE])
  expression <- log1p(t(t(expression) / pmax(library_size, 1) * 1e4))

  plot_data <- as.data.frame(as.table(expression), stringsAsFactors = FALSE)
  colnames(plot_data) <- c("gene", "cell", "expression")
  cell_index <- match(plot_data$cell, input$cells)
  plot_data$sample <- input$sample[cell_index]
  plot_data$lineage_index <- primary_lineage[cell_index]
  plot_data$pseudotime <- input$pseudotime[
    cbind(cell_index, plot_data$lineage_index)
  ]
  plot_data <- plot_data %>%
    dplyr::left_join(lineage_map, by = "lineage_index") %>%
    dplyr::filter(is.finite(pseudotime), is.finite(expression)) %>%
    dplyr::group_by(gene, lineage, sample) %>%
    dplyr::filter(
      dplyr::n() >= min_group_cells,
      dplyr::n_distinct(pseudotime) >= 3
    ) %>%
    dplyr::mutate(
      relative_pseudotime = (pseudotime - min(pseudotime)) /
        (max(pseudotime) - min(pseudotime)),
      pseudotime_bin = dplyr::ntile(pseudotime, n_bins)
    ) %>%
    dplyr::group_by(gene, lineage, sample, pseudotime_bin) %>%
    dplyr::summarise(
      relative_pseudotime = mean(relative_pseudotime),
      mean_expression = mean(expression),
      n_cells = dplyr::n(),
      .groups = "drop"
    )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = relative_pseudotime,
      y = mean_expression,
      color = sample,
      group = sample
    )
  ) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(gene ~ lineage, scales = "free_y") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(
      title = paste(compartment, "sample-stratified marker sensitivity"),
      subtitle = paste0(
        "Descriptive equal-cell pseudotime-bin means; groups with fewer than ",
        min_group_cells, " cells omitted"
      ),
      x = "Relative pseudotime",
      y = "log1p counts per 10,000",
      color = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    output_file,
    plot,
    width = 12,
    height = max(5, 1.5 * length(genes) + 2)
  )
  invisible(TRUE)
}
