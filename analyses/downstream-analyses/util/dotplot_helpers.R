## Define dot plot helper functions
make_dotplot_with_categories <- function(obj,
                                         group_by,
                                         gene_order,
                                         category_bounds,
                                         file_name,
                                         category_labels = NULL,
                                         exclude_subtypes = NULL,
                                         ylab_text = NULL,
                                         width = 15,
                                         height = 5,
                                         y_levels = NULL,
                                         y_axis_labels = NULL,
                                         cat_label_size = 4,
                                         cat_height = 0.8,
                                         cat_y_limits = c(-0.30, 0.05),
                                         show_legend_titles = TRUE,
                                         axis_text_x_size = 10,
                                         axis_text_y_size = 10,
                                         axis_title_y_size = 11) {
  
  obj_plot <- obj
  
  if (!is.null(exclude_subtypes)) {
    obj_plot <- subset(
      obj_plot,
      subset = !myeloid_subtype %in% exclude_subtypes
    )
  }
  
  if (is.null(category_labels)) {
    category_labels <- category_bounds %>%
      mutate(category_clean = category)
  }
  
  p_main <- DotPlot(
    obj_plot,
    features = gene_order,
    cols = "RdYlBu",
    group.by = group_by,
    assay = "RNA"
  )
  
  if (!is.null(y_levels) && !is.null(y_axis_labels)) {
    p_main <- p_main + scale_y_discrete(
      limits = y_levels,
      labels = y_axis_labels
    )
  } else if (!is.null(y_levels)) {
    p_main <- p_main + scale_y_discrete(limits = y_levels)
  } else if (!is.null(y_axis_labels)) {
    p_main <- p_main + scale_y_discrete(labels = y_axis_labels)
  }
  
  p_main <- p_main +
    theme_Publication() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "plain", size = axis_title_y_size),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = axis_text_x_size),
      axis.text.y  = element_text(size = axis_text_y_size),
      legend.title = element_text(face = "plain"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylab(ylab_text) +
    geom_vline(
      data = category_bounds,
      aes(xintercept = end + 0.5),
      linetype = "dotted",
      color = "black",
      size = 0.4,
      inherit.aes = FALSE
    )
  
  if (show_legend_titles) {
    p_main <- p_main +
      labs(
        size = "Percent Expressed",
        color = "Average Expression"
      ) +
      guides(
        size = guide_legend(order = 1),
        color = guide_colorbar(order = 2)
      )
  }
  
  p_cat <- ggplot(category_labels, aes(x = mid, y = 0, label = category_clean)) +
    geom_text(
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      size = cat_label_size,
      lineheight = 1
    ) +
    scale_x_continuous(
      limits = c(1, length(gene_order)),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = cat_y_limits,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(t = 0, r = 10, b = 25, l = 10)
    )
  
  p_final <- p_main / p_cat +
    plot_layout(heights = c(12, cat_height))
  
  grDevices::cairo_pdf(
    filename = file.path(plot_dir, file_name),
    width = width,
    height = height
  )

  print(p_final)
  dev.off()
}

make_manual_dotplot <- function(obj,
                                group_by,
                                categories_list,
                                file_name,
                                exclude_subtypes = NULL,
                                ylab_text = NULL,
                                width = 15,
                                height = 5,
                                y_levels = NULL,
                                y_axis_labels = NULL,
                                cat_label_size = 5,
                                cat_height = 0.6,
                                cat_y_limits = c(-0.25, 0.05),
                                show_legend_titles = TRUE,
                                axis_text_x_size = 10,
                                axis_text_y_size = 10,
                                axis_title_y_size = 11) {
  
  plot_df <- tibble(
    category = rep(names(categories_list), lengths(categories_list)),
    gene = unname(unlist(categories_list))
  )
  
  gene_order_plot <- unname(plot_df$gene)
  
  category_bounds_plot <- plot_df %>%
    mutate(pos = match(gene, gene_order_plot)) %>%
    group_by(category) %>%
    summarise(
      start = min(pos),
      end   = max(pos),
      mid   = (start + end) / 2,
      .groups = "drop"
    ) %>%
    arrange(start)
  
  category_labels_plot <- category_bounds_plot %>%
    mutate(category_clean = category)
  
  make_dotplot_with_categories(
    obj = obj,
    group_by = group_by,
    gene_order = gene_order_plot,
    category_bounds = category_bounds_plot,
    category_labels = category_labels_plot,
    file_name = file_name,
    exclude_subtypes = exclude_subtypes,
    ylab_text = ylab_text,
    width = width,
    height = height,
    y_levels = y_levels,
    y_axis_labels = y_axis_labels,
    cat_label_size = cat_label_size,
    cat_height = cat_height,
    cat_y_limits = cat_y_limits,
    show_legend_titles = show_legend_titles,
    axis_text_x_size = axis_text_x_size,
    axis_text_y_size = axis_text_y_size,
    axis_title_y_size = axis_title_y_size
  )
}

prepare_hallmark_categories <- function(marker_df,
                                        hallmark,
                                        categories_to_remove) {
  
  category_df <- marker_df %>%
    inner_join(
      hallmark %>% select(gs_name, gene_symbol),
      by = c("gene" = "gene_symbol")
    ) %>%
    dplyr::rename(category = gs_name)
  
  category_df_clean <- category_df %>%
    distinct(gene, category, .keep_all = TRUE) %>%
    group_by(gene) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  category_df_clean <- category_df_clean %>%
    filter(!category %in% categories_to_remove) %>%
    arrange(category, gene)
  
  category_counts <- category_df_clean %>%
    group_by(category) %>%
    tally(name = "n_genes")
  
  categories_to_keep <- category_counts %>%
    filter(n_genes > 1) %>%
    pull(category)
  
  category_df_clean <- category_df_clean %>%
    filter(category %in% categories_to_keep)
  
  gene_order <- unique(category_df_clean$gene)
  
  gene_pos_df <- category_df_clean %>%
    mutate(pos = match(gene, gene_order)) %>%
    filter(!is.na(pos))
  
  category_bounds <- gene_pos_df %>%
    group_by(category) %>%
    summarize(
      start = min(pos),
      end   = max(pos),
      mid   = (start + end) / 2,
      .groups = "drop"
    )
  
  category_labels <- category_bounds %>%
    mutate(
      category_clean = category %>%
        str_remove("HALLMARK_") %>%
        str_replace_all("_", "\n")
    )
  
  list(
    gene_order = gene_order,
    category_bounds = category_bounds,
    category_labels = category_labels
  )
}

