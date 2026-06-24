# Shared miloR utilities for downstream neighborhood differential abundance analyses.
# These helpers keep common mechanics in one place while the Rmd files retain
# population-specific inputs, labels, graph settings, and output paths.

# Set deterministic RNG behavior before sampling Milo neighborhoods.
set_seed <- function(seed = 1234) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
}

# Convert labels into filesystem-safe strings for output filenames.
safe_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

# Convert a Seurat object to SingleCellExperiment and copy metadata/reductions used by miloR.
prepare_milo_sce <- function(obj,
                             cell_group_col,
                             condition_levels,
                             milo_reduction,
                             milo_plot_reduction,
                             milo_dims) {
  milo_sce <- as.SingleCellExperiment(obj)

  stopifnot(all(colnames(milo_sce) %in% rownames(obj@meta.data)))

  colData(milo_sce)$sample <- obj@meta.data[colnames(milo_sce), "sample", drop = TRUE]
  colData(milo_sce)$condition <- factor(
    obj@meta.data[colnames(milo_sce), "condition", drop = TRUE],
    levels = condition_levels
  )
  colData(milo_sce)$cell_group <- obj@meta.data[colnames(milo_sce), cell_group_col, drop = TRUE]
  colData(milo_sce)$condition_milo <- factor(
    make.names(as.character(colData(milo_sce)$condition)),
    levels = make.names(condition_levels)
  )

  if (!milo_reduction %in% reducedDimNames(milo_sce)) {
    seurat_reduction <- tolower(milo_reduction)
    stopifnot(seurat_reduction %in% names(obj@reductions))
    reducedDim(milo_sce, milo_reduction) <- Seurat::Embeddings(obj, seurat_reduction)
  }

  if (!milo_plot_reduction %in% reducedDimNames(milo_sce)) {
    seurat_plot_reduction <- tolower(milo_plot_reduction)
    if (seurat_plot_reduction %in% names(obj@reductions)) {
      reducedDim(milo_sce, milo_plot_reduction) <- Seurat::Embeddings(obj, seurat_plot_reduction)
    }
  }

  milo_d <- min(milo_dims, ncol(reducedDim(milo_sce, milo_reduction)))
  milo_layout <- if (milo_plot_reduction %in% reducedDimNames(milo_sce)) {
    milo_plot_reduction
  } else {
    milo_reduction
  }

  list(
    sce = milo_sce,
    d = milo_d,
    layout = milo_layout
  )
}

# Build the Milo graph, sample neighborhoods, and count cells per biological sample.
build_milo_neighborhoods <- function(milo_sce,
                                     milo_reduction,
                                     milo_d,
                                     milo_k,
                                     milo_prop,
                                     plot_size_hist = TRUE) {
  milo_obj <- miloR::Milo(milo_sce)
  milo_obj <- miloR::buildGraph(
    milo_obj,
    k = milo_k,
    d = milo_d,
    reduced.dim = milo_reduction
  )
  milo_obj <- miloR::makeNhoods(
    milo_obj,
    prop = milo_prop,
    k = milo_k,
    d = milo_d,
    refined = TRUE,
    reduced_dims = milo_reduction
  )

  if (isTRUE(plot_size_hist)) {
    miloR::plotNhoodSizeHist(milo_obj)
  }

  miloR::countCells(
    milo_obj,
    meta.data = data.frame(colData(milo_obj)),
    samples = "sample"
  )
}

# Build and align the sample-level design table to the Milo neighborhood count matrix.
build_milo_design <- function(milo_obj, reference_condition) {
  milo_design <- data.frame(colData(milo_obj)) %>%
    dplyr::select(sample, condition, condition_milo) %>%
    dplyr::distinct()
  rownames(milo_design) <- milo_design$sample
  milo_design <- milo_design[colnames(miloR::nhoodCounts(milo_obj)), , drop = FALSE]

  stopifnot(!any(is.na(milo_design$condition)))
  stopifnot(reference_condition %in% as.character(milo_design$condition))

  milo_design
}

# Create explicit tumor-reference model contrasts using make.names-compatible condition labels.
make_milo_contrasts <- function(condition_levels, reference_condition) {
  tibble::tibble(
    condition = condition_levels,
    condition_milo = make.names(condition_levels)
  ) %>%
    dplyr::filter(condition != reference_condition) %>%
    dplyr::mutate(
      comparison = paste(condition, "vs", reference_condition),
      contrast = paste0(
        "condition_milo",
        condition_milo,
        " - condition_milo",
        make.names(reference_condition)
      )
    )
}

# Run miloR neighborhood differential abundance testing for each condition-vs-reference contrast.
run_milo_da <- function(milo_obj,
                        milo_design,
                        condition_levels,
                        reference_condition,
                        milo_reduction,
                        milo_d,
                        cell_group_col = "cell_group") {
  milo_obj <- miloR::calcNhoodDistance(
    milo_obj,
    d = milo_d,
    reduced.dim = milo_reduction
  )

  milo_contrasts <- make_milo_contrasts(condition_levels, reference_condition)
  milo_da_results <- vector("list", nrow(milo_contrasts))

  for (i in seq_len(nrow(milo_contrasts))) {
    comparison <- milo_contrasts$comparison[i]
    contrast <- milo_contrasts$contrast[i]

    cat("Running Milo contrast:", comparison, "\n")

    milo_da <- miloR::testNhoods(
      milo_obj,
      design = ~ 0 + condition_milo,
      design.df = milo_design,
      model.contrasts = contrast,
      fdr.weighting = "graph-overlap",
      norm.method = "TMM",
      reduced.dim = milo_reduction
    )

    milo_da <- miloR::annotateNhoods(
      milo_obj,
      milo_da,
      coldata_col = cell_group_col
    )

    milo_da$comparison <- comparison
    milo_da$contrast <- contrast
    milo_da_results[[i]] <- milo_da
  }

  list(
    milo_obj = milo_obj,
    da_results = dplyr::bind_rows(milo_da_results),
    contrasts = milo_contrasts
  )
}

# Add nominal/FDR flags and summarize Milo results by contrast and annotated cell group.
summarize_milo_da <- function(milo_da_results) {
  milo_fdr_col <- if ("SpatialFDR" %in% colnames(milo_da_results)) {
    "SpatialFDR"
  } else {
    "FDR"
  }

  milo_da_results <- milo_da_results %>%
    dplyr::mutate(
      spatial_fdr = .data[[milo_fdr_col]],
      is_nominal = PValue < 0.05,
      is_fdr_significant = FDR < 0.05,
      is_spatial_fdr_significant = spatial_fdr < 0.05
    )

  threshold_summary <- milo_da_results %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_neighborhoods = dplyr::n(),
      min_p_value = min(PValue, na.rm = TRUE),
      min_fdr = min(FDR, na.rm = TRUE),
      min_spatial_fdr = min(spatial_fdr, na.rm = TRUE),
      n_nominal = sum(is_nominal, na.rm = TRUE),
      n_fdr_0_05 = sum(is_fdr_significant, na.rm = TRUE),
      n_spatial_fdr_0_05 = sum(is_spatial_fdr_significant, na.rm = TRUE),
      max_abs_logFC = max(abs(logFC), na.rm = TRUE),
      .groups = "drop"
    )

  effect_summary <- milo_da_results %>%
    dplyr::group_by(comparison, cell_group) %>%
    dplyr::summarise(
      n_neighborhoods = dplyr::n(),
      median_logFC = median(logFC, na.rm = TRUE),
      mean_logFC = mean(logFC, na.rm = TRUE),
      q25_logFC = quantile(logFC, 0.25, na.rm = TRUE),
      q75_logFC = quantile(logFC, 0.75, na.rm = TRUE),
      fraction_positive_logFC = mean(logFC > 0, na.rm = TRUE),
      n_nominal = sum(is_nominal, na.rm = TRUE),
      fraction_nominal = n_nominal / n_neighborhoods,
      median_logFC_nominal = ifelse(
        any(is_nominal, na.rm = TRUE),
        median(logFC[is_nominal], na.rm = TRUE),
        NA_real_
      ),
      min_p_value = min(PValue, na.rm = TRUE),
      min_spatial_fdr = min(spatial_fdr, na.rm = TRUE),
      median_cell_group_fraction = median(cell_group_fraction, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      neg_log10_min_p = -log10(pmax(min_p_value, .Machine$double.xmin))
    ) %>%
    dplyr::arrange(comparison, min_p_value)

  list(
    da_results = milo_da_results,
    threshold_summary = threshold_summary,
    effect_summary = effect_summary
  )
}

# Plot cell-group-level Milo trends with nominal P value as dot size and median logFC as color.
plot_milo_effect_dot <- function(effect_summary, title) {
  ggplot2::ggplot(
    effect_summary,
    ggplot2::aes(
      x = comparison,
      y = cell_group,
      size = neg_log10_min_p,
      color = median_logFC
    )
  ) +
    ggplot2::geom_point(alpha = 0.85, na.rm = TRUE) +
    ggplot2::scale_color_gradient2(
      low = "#2C7BB6",
      mid = "grey90",
      high = "#D7191C",
      midpoint = 0,
      na.value = "grey85"
    ) +
    ggplot2::scale_size_continuous(
      range = c(1, 7),
      breaks = -log10(c(0.05, 0.01, 0.001)),
      labels = c("0.05", "0.01", "0.001")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      size = "Nominal\nP value",
      color = "Median\nlogFC",
      title = title
    )
}

# Plot neighborhood logFC distributions by cell group, highlighting groups with nominal P < 0.05.
plot_milo_logfc_box <- function(milo_da_results, effect_summary, title) {
  milo_box_order <- effect_summary %>%
    dplyr::group_by(comparison) %>%
    dplyr::arrange(dplyr::desc(median_logFC), .by_group = TRUE) %>%
    dplyr::mutate(
      cell_group_facet = paste(comparison, cell_group, sep = "___"),
      nominal_status = ifelse(n_nominal > 0, "Nominal P < 0.05", "Not significant")
    ) %>%
    dplyr::ungroup()

  milo_box_levels <- rev(milo_box_order$cell_group_facet)

  plot_df <- milo_da_results %>%
    dplyr::left_join(
      dplyr::select(milo_box_order, comparison, cell_group, cell_group_facet, nominal_status),
      by = c("comparison", "cell_group")
    ) %>%
    dplyr::mutate(
      cell_group_facet = factor(cell_group_facet, levels = milo_box_levels),
      nominal_status = factor(nominal_status, levels = c("Not significant", "Nominal P < 0.05"))
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = logFC, y = cell_group_facet, fill = nominal_status)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.3, color = "grey35") +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c("Not significant" = "grey92", "Nominal P < 0.05" = "#D7191C")
    ) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::facet_wrap(~ comparison, ncol = 2, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Neighborhood logFC",
      y = NULL,
      title = title
    )
}

# Optionally save native miloR neighborhood graph plots for exploratory QC.
save_optional_milo_graphs <- function(milo_obj,
                                      milo_da_results,
                                      milo_layout,
                                      analysis_label,
                                      reference_condition,
                                      population_title,
                                      plot_dir) {
  milo_obj <- miloR::buildNhoodGraph(milo_obj)

  for (comparison_label in unique(milo_da_results$comparison)) {
    milo_da_plot <- milo_da_results %>%
      dplyr::filter(.data$comparison == comparison_label)

    p <- miloR::plotNhoodGraphDA(
      milo_obj,
      milo_da_plot,
      alpha = 0.05,
      layout = milo_layout
    ) +
      ggplot2::ggtitle(paste0(population_title, " miloR DA: ", comparison_label))

    print(p)

    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(analysis_label, "_milor_", safe_label(comparison_label), "_", safe_label(reference_condition), "_baseline.pdf")
      ),
      plot = p,
      width = 7,
      height = 6,
      device = grDevices::cairo_pdf
    )
  }
}
