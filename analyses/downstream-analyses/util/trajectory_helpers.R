# Helper functions for T cell Slingshot trajectory analysis.
#
# These functions keep the Rmd focused on the analysis flow while storing
# reusable plotting and pseudotime-gene-screen logic in one place.

# Return a log-normalized expression matrix from a Seurat object.
# Supports both Seurat v5 `layer` and older `slot` access patterns.
get_log_data <- function(seu, assay = SeuratObject::DefaultAssay(seu)) {
  tryCatch(
    SeuratObject::GetAssayData(seu, assay = assay, layer = "data"),
    error = function(e) SeuratObject::GetAssayData(seu, assay = assay, slot = "data")
  )
}

# Approximate a smooth UMAP trajectory for one Slingshot lineage.
# Cells are ordered by pseudotime, binned to reduce noise, and smoothed with
# base R splines. This is for visualization only; Slingshot itself was run in
# the chosen reduced-dimensional space, not in UMAP.
make_umap_curve_from_pseudotime <- function(umap_df, lineage_col, lineage_name, approx_points = 150) {
  
  lineage_cells <- umap_df %>%
    dplyr::filter(!is.na(.data[[lineage_col]])) %>%
    dplyr::mutate(pseudotime = .data[[lineage_col]]) %>%
    dplyr::arrange(pseudotime) %>%
    dplyr::select(cell, UMAP_1, UMAP_2, pseudotime)
  
  if (nrow(lineage_cells) < 4) {
    return(
      lineage_cells %>%
        dplyr::mutate(lineage = lineage_name)
    )
  }
  
  n_bins <- min(approx_points, nrow(lineage_cells))
  
  binned_curve <- lineage_cells %>%
    dplyr::mutate(
      pt_bin = ceiling(dplyr::row_number() / nrow(lineage_cells) * n_bins)
    ) %>%
    dplyr::group_by(pt_bin) %>%
    dplyr::summarise(
      pseudotime = median(pseudotime, na.rm = TRUE),
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(pseudotime)
  
  if (nrow(binned_curve) >= 5 && length(unique(binned_curve$pseudotime)) >= 5) {
    
    pt_grid <- seq(
      min(binned_curve$pseudotime, na.rm = TRUE),
      max(binned_curve$pseudotime, na.rm = TRUE),
      length.out = approx_points
    )
    
    smoothed_curve <- tryCatch({
      fit_x <- stats::smooth.spline(
        x = binned_curve$pseudotime,
        y = binned_curve$UMAP_1,
        spar = 0.6
      )
      
      fit_y <- stats::smooth.spline(
        x = binned_curve$pseudotime,
        y = binned_curve$UMAP_2,
        spar = 0.6
      )
      
      tibble::tibble(
        pseudotime = pt_grid,
        UMAP_1 = stats::predict(fit_x, pt_grid)$y,
        UMAP_2 = stats::predict(fit_y, pt_grid)$y,
        lineage = lineage_name
      )
    }, error = function(e) {
      binned_curve %>%
        dplyr::mutate(lineage = lineage_name)
    })
    
    return(smoothed_curve)
  }
  
  binned_curve %>%
    dplyr::mutate(lineage = lineage_name)
}

# Summarize subtype centroids for each lineage. These points are used as
# approximate lineage nodes in the UMAP visualization.
make_lineage_node_df <- function(umap_df, pseudotime_cols, min_cells = 10) {
  purrr::map2_dfr(
    pseudotime_cols,
    paste0("Lineage ", seq_along(pseudotime_cols)),
    function(lineage_col, lineage_name) {
      umap_df %>%
        dplyr::filter(!is.na(.data[[lineage_col]])) %>%
        dplyr::group_by(tcell_subtype) %>%
        dplyr::summarise(
          n = dplyr::n(),
          UMAP_1 = median(UMAP_1, na.rm = TRUE),
          UMAP_2 = median(UMAP_2, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::filter(n >= min_cells) %>%
        dplyr::mutate(lineage = lineage_name)
    }
  )
}

# Make one arrow segment per lineage, pointing from lower to higher pseudotime.
make_lineage_arrow_df <- function(curve_df, arrow_start_frac = 0.68, arrow_end_frac = 0.82) {
  curve_df %>%
    dplyr::group_by(lineage) %>%
    dplyr::arrange(pseudotime, .by_group = TRUE) %>%
    dplyr::group_modify(function(.x, .y) {
      n_points <- nrow(.x)
      
      if (n_points < 2) {
        return(tibble::tibble())
      }
      
      start_idx <- max(1, min(n_points - 1, round(n_points * arrow_start_frac)))
      end_idx <- max(start_idx + 1, min(n_points, round(n_points * arrow_end_frac)))
      
      tibble::tibble(
        UMAP_1 = .x$UMAP_1[start_idx],
        UMAP_2 = .x$UMAP_2[start_idx],
        UMAP_1_end = .x$UMAP_1[end_idx],
        UMAP_2_end = .x$UMAP_2[end_idx]
      )
    }) %>%
    dplyr::ungroup()
}

# Plot one Slingshot lineage on UMAP using the precomputed curve and node data.
plot_one_sling_lineage <- function(
  lineage_id,
  pseudotime_cols,
  umap_df,
  curve_df,
  lineage_node_df,
  subtype_colors = NULL,
  show_subtype_labels = TRUE
) {
  
  lineage_col <- pseudotime_cols[lineage_id]
  lineage_name <- paste0("Lineage ", lineage_id)
  
  lineage_umap_df <- umap_df %>%
    dplyr::filter(!is.na(.data[[lineage_col]]))
  
  crv_df <- curve_df %>%
    dplyr::filter(lineage == lineage_name)
  
  node_df <- lineage_node_df %>%
    dplyr::filter(lineage == lineage_name)
  
  arrow_df <- make_lineage_arrow_df(crv_df)
  
  label_df <- lineage_umap_df %>%
    dplyr::group_by(tcell_subtype) %>%
    dplyr::summarise(
      n = dplyr::n(),
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= 10)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = umap_df,
      ggplot2::aes(UMAP_1, UMAP_2),
      color = "grey88",
      size = 0.35,
      alpha = 0.6
    ) +
    ggplot2::geom_point(
      data = lineage_umap_df,
      ggplot2::aes(UMAP_1, UMAP_2, color = tcell_subtype),
      size = 0.55,
      alpha = 0.9
    ) +
    ggplot2::geom_path(
      data = crv_df,
      ggplot2::aes(UMAP_1, UMAP_2),
      color = "black",
      linewidth = 1.2,
      lineend = "round"
    ) +
    ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_end, yend = UMAP_2_end),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.2,
      lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(0.18, "inches"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(UMAP_1, UMAP_2),
      color = "black",
      size = 3.5
    ) +
    theme_Publication() +
    ggplot2::theme(
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background = ggplot2::element_rect(fill = "white"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(color = "T cell subtype") +
    ggplot2::ggtitle(lineage_name)
  
  if (!is.null(subtype_colors)) {
    p <- p +
      ggplot2::scale_color_manual(
        values = subtype_colors,
        breaks = names(subtype_colors),
        drop = FALSE
      )
  }
  
  if (show_subtype_labels && nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(UMAP_1, UMAP_2, label = tcell_subtype),
        inherit.aes = FALSE,
        size = 2,
        color = "black",
        fill = ggplot2::alpha("white", 0.85),
        label.size = 0.15,
        label.padding = grid::unit(0.12, "lines"),
        point.padding = grid::unit(0.25, "lines"),
        box.padding = grid::unit(0.35, "lines"),
        min.segment.length = 0,
        segment.color = "grey45",
        segment.size = 0.25,
        force = 2,
        force_pull = 0.25,
        max.iter = 10000,
        max.time = 2,
        max.overlaps = Inf,
        seed = 1234
      )
  }
  
  p
}

# Fast per-lineage pseudotime-gene screen.
# For one pseudotime column, this computes Spearman correlation between
# log-normalized expression and pseudotime for each requested gene.
fast_spearman_pseudotime_test <- function(
  seu,
  pseudotime_col,
  features = NULL,
  assay = SeuratObject::DefaultAssay(seu),
  min_cells = 50,
  min_detected_frac = 0.02
) {
  
  expr <- get_log_data(seu, assay = assay)
  
  if (is.null(features) || length(features) == 0) {
    features <- rownames(expr)
  }
  
  features <- intersect(features, rownames(expr))
  pt <- seu@meta.data[colnames(expr), pseudotime_col]
  keep_cells <- !is.na(pt)
  
  expr <- expr[features, keep_cells, drop = FALSE]
  pt <- pt[keep_cells]
  
  if (length(pt) < min_cells) {
    stop("Too few cells with non-NA pseudotime for ", pseudotime_col)
  }
  
  detected_frac <- Matrix::rowMeans(expr > 0)
  keep_genes <- detected_frac >= min_detected_frac
  expr <- expr[keep_genes, , drop = FALSE]
  detected_frac <- detected_frac[keep_genes]
  
  pt_rank <- rank(pt, ties.method = "average")
  expr_rank <- t(apply(
    as.matrix(expr),
    1,
    rank,
    ties.method = "average"
  ))
  
  # Spearman rho is Pearson correlation computed on ranked values.
  rho <- as.numeric(stats::cor(t(expr_rank), pt_rank, method = "pearson"))
  n_cells <- length(pt)
  t_stat <- rho * sqrt((n_cells - 2) / pmax(1 - rho^2, .Machine$double.eps))
  p_val <- 2 * stats::pt(abs(t_stat), df = n_cells - 2, lower.tail = FALSE)
  
  tibble::tibble(
    gene = rownames(expr),
    lineage = pseudotime_col,
    rho = rho,
    p_val = p_val,
    p_adj = p.adjust(p_val, method = "BH"),
    abs_rho = abs(rho),
    detected_frac = as.numeric(detected_frac),
    n_cells = n_cells
  ) %>%
    dplyr::arrange(p_adj, dplyr::desc(abs_rho))
}

# Select the top-ranked pseudotime genes for one lineage from the Spearman
# screen output.
get_top_pseudotime_genes <- function(pseudotime_gene_screen, lineage_name, n = 40) {
  pseudotime_gene_screen %>%
    dplyr::filter(lineage == lineage_name, !is.na(p_adj)) %>%
    dplyr::slice_min(order_by = p_adj, n = n, with_ties = FALSE) %>%
    dplyr::pull(gene) %>%
    unique()
}
