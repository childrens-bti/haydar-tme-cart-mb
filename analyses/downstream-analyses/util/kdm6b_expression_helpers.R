# ==============================================================================
# Kdm6b expression and subcluster definition helpers
# Author: Bicna Song
# Date: 2026-08-10
#
# Purpose:
#   Support 14-kdm6b-expression-subcluster-definition.Rmd by creating
#   pseudobulk summaries and defining Kdm6b-high and Kdm6b-low subcluster groups.
#
# ==============================================================================

#' Extract the D1/D2 sample number from a sample name
#'
#' Identifies a terminal sample label formatted as `_D<number>`, such as `D1`
#' in `41BB-L_pos_T_D1`. D1 and D2 denote biological samples 1 and 2 within
#' each condition; they are not matched blocks across conditions.
#'
#' @param sample_name Character vector of sample names.
#'
#' @return A character vector containing the D1/D2 label for each sample.
#'   Returns `NA` when a sample name does not end in `_D<number>`.
extract_sample_number <- function(sample_name) {
  label <- sub("^.*_(D[0-9]+)$", "\\1", sample_name)
  ifelse(label == sample_name, NA_character_, label)
}

#' Sum a sparse count matrix by a cell-level grouping variable
#'
#' @param counts Sparse genes-by-cells count matrix.
#' @param groups Character vector assigning each cell to a group.
#' @param group_levels Character vector defining output column order.
#'
#' @return Sparse genes-by-groups count matrix.
sparse_group_sum <- function(counts, groups, group_levels = unique(groups)) {
  stopifnot(length(groups) == ncol(counts))
  stopifnot(!anyDuplicated(group_levels))
  stopifnot(all(groups %in% group_levels))

  indicator <- Matrix::sparseMatrix(
    i = seq_len(ncol(counts)),
    j = match(groups, group_levels),
    x = 1,
    dims = c(ncol(counts), length(group_levels)),
    dimnames = list(colnames(counts), group_levels)
  )

  counts %*% indicator
}

#' Create Kdm6b pseudobulk summaries for one cell compartment
#'
#' Aggregates raw RNA counts from cells to biological samples using sparse
#' matrix multiplication. DESeq2 size factors are estimated from the full
#' pseudobulk transcriptome within the supplied compartment, then used to
#' report normalized Kdm6b counts.
#'
#' @param obj Seurat object containing a joined RNA `counts` layer and `sample`
#'   and `condition` metadata columns.
#' @param compartment Character label written to the returned summary table,
#'   for example `"All cells"`, `"T cells"`, or `"Myeloid cells"`.
#' @param condition_levels Character vector defining the allowed treatment
#'   conditions and their factor order.
#' @param gene Gene symbol to summarize. Defaults to `"Kdm6b"`.
#'
#' @return Tibble with one row per biological sample containing the condition,
#'   cell count, raw Kdm6b UMIs, and DESeq2-normalized Kdm6b counts.
#'
#' @details Each sample must map to exactly one condition. D1 and D2 are
#'   biological samples 1 and 2 within each condition, and each becomes a
#'   separate pseudobulk column.
aggregate_kdm6b_by_sample <- function(obj,
                                      compartment,
                                      condition_levels,
                                      gene = "Kdm6b") {
  required_metadata <- c("sample", "condition")
  stopifnot(all(required_metadata %in% colnames(obj@meta.data)))
  stopifnot(gene %in% rownames(obj))

  counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  meta <- obj@meta.data %>%
    dplyr::transmute(
      sample = as.character(sample),
      condition = as.character(condition)
    )

  sample_metadata <- meta %>%
    dplyr::distinct(sample, condition)

  if (anyDuplicated(sample_metadata$sample)) {
    stop("Each sample must map to exactly one condition.")
  }

  sample_metadata <- sample_metadata %>%
    dplyr::mutate(
      condition = factor(condition, levels = condition_levels),
      sample_number = extract_sample_number(sample)
    ) %>%
    dplyr::arrange(condition, sample)

  if (any(is.na(sample_metadata$condition))) {
    stop("At least one sample has a condition outside condition_levels.")
  }

  pseudobulk_counts <- sparse_group_sum(
    counts = counts,
    groups = meta$sample,
    group_levels = sample_metadata$sample
  )

  cell_counts <- table(meta$sample)

  coldata <- as.data.frame(sample_metadata)
  rownames(coldata) <- coldata$sample

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(pseudobulk_counts),
    colData = coldata,
    design = ~1
  )
  dds <- DESeq2::estimateSizeFactors(dds)

  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

  sample_metadata %>%
    dplyr::mutate(
      compartment = compartment,
      n_cells = as.integer(cell_counts[sample]),
      kdm6b_raw_umis = as.numeric(pseudobulk_counts[gene, sample]),
      kdm6b_deseq2_normalized_count = as.numeric(
        normalized_counts[gene, sample]
      )
    ) %>%
    dplyr::select(
      compartment, sample, sample_number, condition, n_cells, kdm6b_raw_umis,
      kdm6b_deseq2_normalized_count
    )
}

#' Aggregate Kdm6b expression within sample-by-subcluster profiles
#'
#' Creates one expression profile for every observed combination of biological
#' sample and numeric subcluster. Kdm6b is summarized by detection
#' fraction and mean log-normalized expression.
#'
#' @param obj Seurat object containing joined RNA `counts` and `data`
#'   layers plus `sample`, `condition`, `seurat_clusters`, and `subtype`
#'   metadata columns.
#' @param condition_levels Character vector defining allowed treatment
#'   conditions and their factor order.
#' @param gene Gene symbol to summarize. Defaults to `"Kdm6b"`.
#'
#' @return Tibble with one row per observed sample/subcluster combination,
#'   including cell count, detection fraction, and mean expression.
aggregate_kdm6b_by_sample_cluster <- function(obj,
                                              condition_levels,
                                              gene = "Kdm6b",
                                              subtype_col) {
  required_metadata <- c(
    "sample", "condition", "seurat_clusters", subtype_col
  )
  stopifnot(all(required_metadata %in% colnames(obj@meta.data)))
  stopifnot(gene %in% rownames(obj))

  counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  log_data <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
  condition_factor <- factor(
    as.character(obj$condition),
    levels = condition_levels
  )

  if (anyNA(condition_factor)) {
    stop("At least one cell has a condition outside condition_levels.")
  }

  obj@meta.data %>%
    dplyr::transmute(
      sample = as.character(sample),
      condition = condition_factor,
      subcluster = as.character(seurat_clusters),
      !!subtype_col := as.character(obj@meta.data[[subtype_col]]),
      kdm6b_detected = as.numeric(counts[gene, ]) > 0,
      kdm6b_log_expression = as.numeric(log_data[gene, ])
    ) %>%
    dplyr::group_by(sample, condition, subcluster, .data[[subtype_col]]) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      kdm6b_detection_fraction = mean(kdm6b_detected),
      kdm6b_mean_log_expression = mean(kdm6b_log_expression),
      .groups = "drop"
    ) %>%
    dplyr::arrange(as.integer(subcluster), condition, sample)
}

#' Rank subclusters and assign Kdm6b groups
#'
#' Ranks numeric subclusters by the median, across samples, of mean
#' log-normalized Kdm6b expression. Each qualifying sample/subcluster profile
#' contributes equally. Detection fraction is used only as a tie-breaker.
#'
#' @param sample_cluster_table Output from
#'   `aggregate_kdm6b_by_sample_cluster()`.
#' @param min_cells_per_profile Minimum cells required for a sample/subcluster
#'   profile to contribute to the ranking. Defaults to 20.
#' @param tail_fraction Fraction of subclusters assigned to each extreme
#'   tail. Defaults to 0.25, producing top-quartile high and bottom-quartile low
#'   groups while leaving the middle half intermediate.
#'
#' @return Tibble with one row per numeric subcluster containing expression
#'   summaries, sample and condition counts, rank, and `kdm6b_group`.
rank_kdm6b_subclusters <- function(sample_cluster_table,
                                   min_cells_per_profile = 20,
                                   tail_fraction = 0.25,
                                   subtype_col) {
  stopifnot(tail_fraction > 0, tail_fraction < 0.5)
  stopifnot(min_cells_per_profile > 0)

  cluster_totals <- sample_cluster_table %>%
    dplyr::group_by(subcluster, .data[[subtype_col]]) %>%
    dplyr::summarise(
      total_cells = sum(n_cells),
      .groups = "drop"
    )

  ranking <- sample_cluster_table %>%
    dplyr::filter(n_cells >= min_cells_per_profile) %>%
    dplyr::group_by(subcluster, .data[[subtype_col]]) %>%
    dplyr::summarise(
      n_sample_profiles = dplyr::n(),
      n_conditions = dplyr::n_distinct(condition),
      median_log_expression = median(
        kdm6b_mean_log_expression, na.rm = TRUE
      ),
      q25_log_expression = stats::quantile(
        kdm6b_mean_log_expression,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      ),
      q75_log_expression = stats::quantile(
        kdm6b_mean_log_expression,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      ),
      median_detection_fraction = median(
        kdm6b_detection_fraction, na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::right_join(cluster_totals, by = c("subcluster", subtype_col)) %>%
    dplyr::mutate(
      n_sample_profiles = dplyr::coalesce(n_sample_profiles, 0L),
      n_conditions = dplyr::coalesce(n_conditions, 0L)
    )

  if (any(ranking$n_sample_profiles == 0)) {
    missing_clusters <- ranking$subcluster[ranking$n_sample_profiles == 0]
    stop(
      "No sample has at least ", min_cells_per_profile,
      " cells for cluster(s): ", paste(missing_clusters, collapse = ", ")
    )
  }

  ranking <- ranking %>%
    dplyr::arrange(
      dplyr::desc(median_log_expression),
      dplyr::desc(median_detection_fraction),
      as.integer(subcluster)
    ) %>%
    dplyr::mutate(kdm6b_rank = dplyr::row_number())

  n_clusters <- nrow(ranking)
  if (n_clusters < 4) {
    stop("At least four subclusters are required for tail assignment.")
  }
  n_extreme <- max(1L, floor(n_clusters * tail_fraction))

  ranking %>%
    dplyr::mutate(
      kdm6b_group = dplyr::case_when(
        kdm6b_rank <= n_extreme ~ "Kdm6b-high",
        kdm6b_rank > n_clusters - n_extreme ~ "Kdm6b-low",
        TRUE ~ "Intermediate"
      )
    )
}
