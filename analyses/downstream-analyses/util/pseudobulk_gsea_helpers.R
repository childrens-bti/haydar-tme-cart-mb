## ## Define pseudobulk DESeq2‑based GSEA helper functions

# Write fgsea results to CSV
write_fgsea_csv <- function(df, out_path) {
  df2 <- df %>%
    mutate(
      leadingEdge = if ("leadingEdge" %in% colnames(df)) {
        vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1))
      } else {
        NA_character_
      }
    )
  write.csv(df2, out_path, row.names = FALSE)
}

# Format contrast name for plot titles
format_contrast_title <- function(out_prefix) {
  out_prefix %>%
    str_remove("^pseudobulk_") %>%
    str_replace_all("_vs_", " vs ") %>%
    str_replace_all("otherCAR", "other CARs") %>%
    str_replace_all("allCAR", "all CARs")
}

# Volcano plot
plot_volcano <- function(df,
                         gene_col = "gene",
                         lfc_col,
                         p_col,
                         title,
                         subtitle = "",
                         caption = NULL,
                         p_cutoff = 0.05,
                         lfc_cutoff = 1,
                         label_top_n = 10) {
  stopifnot(all(c(gene_col, lfc_col, p_col) %in% colnames(df)))
  
  df2 <- df %>%
    mutate(.p = .data[[p_col]], .lfc = .data[[lfc_col]]) %>%
    filter(is.finite(.p), is.finite(.lfc), !is.na(.p), !is.na(.lfc))
  
  label_df <- df2 %>%
    filter(.p <= p_cutoff, abs(.lfc) >= lfc_cutoff) %>%
    arrange(.p) %>%
    slice_head(n = label_top_n)
  
  lab_vec <- ifelse(df2[[gene_col]] %in% label_df[[gene_col]], df2[[gene_col]], "")
  
  if (is.null(caption)) {
    caption <- paste0("total = ", nrow(df2), " genes")
  }
  
  EnhancedVolcano::EnhancedVolcano(
    df2,
    lab = lab_vec,
    x = lfc_col,
    y = p_col,
    title = title,
    subtitle = subtitle,
    caption = caption,
    pCutoff = p_cutoff,
    FCcutoff = lfc_cutoff,
    legendPosition = "bottom",
    pointSize = 1.5,
    colAlpha = 0.8,
    labSize = 3,
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    colConnectors = "black",
    max.overlaps = Inf,
    titleLabSize = 10
  )
}

# GSEA dotplot (activated vs suppressed)
plot_gsea_dot <- function(gsea_df,
                          title = NULL,
                          p_cutoff = 0.05,
                          top_n = 10,
                          use_p = c("pval", "padj")) {
  use_p <- match.arg(use_p)
  stopifnot(all(c("pathway", "NES", use_p, "size") %in% colnames(gsea_df)))
  
  df <- gsea_df %>%
    mutate(
      p = .data[[use_p]],
      direction = ifelse(NES >= 0, "activated", "suppressed"),
      Count = if ("leadingEdge" %in% colnames(gsea_df)) lengths(leadingEdge) else NA_integer_,
      GeneRatio = if ("leadingEdge" %in% colnames(gsea_df)) (lengths(leadingEdge) / size) else NA_real_,
      pathway_clean = str_remove(pathway, "^HALLMARK_")
    ) %>%
    filter(is.finite(NES), is.finite(p))
  
  df2 <- df %>%
    filter(p <= p_cutoff) %>%
    group_by(direction) %>%
    arrange(p) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  if (nrow(df2) == 0) {
    message("No pathways with ", use_p, " <= ", p_cutoff, " for ", title, ".")
    return(NULL)
  }
  
  df2 <- df2 %>%
    group_by(direction) %>%
    mutate(pathway_clean = factor(pathway_clean, levels = rev(unique(pathway_clean)))) %>%
    ungroup()
  
  ggplot(df2, aes(x = GeneRatio, y = pathway_clean)) +
    geom_point(aes(size = Count, color = p), alpha = 0.9) +
    facet_grid(. ~ direction, scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
    labs(
      x = "GeneRatio",
      y = NULL,
      size = "Count",
      color = use_p,
      title = title
    ) +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# Run pseudobulk DESeq2 + fgsea for one contrast
run_pseudobulk_fgsea <- function(pb_counts,
                                 coldata,
                                 cond2_levels,
                                 hallmark_pathways,
                                 out_prefix,
                                 subtype_label,
                                 subtype_title,
                                 results_dir,
                                 plot_dir,
                                 minSize = 15,
                                 maxSize = 500,
                                 volcano_lfc_cutoff = 1,
                                 volcano_p_cutoff = 0.05,
                                 gsea_p_cutoff = 0.05) {
  coldata$cond2 <- factor(coldata$cond2, levels = cond2_levels)
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = coldata,
    design = ~ cond2
  )
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("cond2", cond2_levels[2], cond2_levels[1]))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(stat), is.finite(stat)) %>%
    arrange(padj)
  
  write.csv(
    res_df,
    file.path(results_dir, paste0(subtype_label, "_DESeq2_", out_prefix, ".csv")),
    row.names = FALSE
  )
  
  ranks <- setNames(res_df$stat, res_df$gene)
  ranks <- sort(ranks, decreasing = TRUE)
  
  set.seed(1234)
  gsea <- fgsea(
    pathways = hallmark_pathways,
    stats = ranks,
    minSize = minSize,
    maxSize = maxSize
  ) %>%
    arrange(pval)
  
  write_fgsea_csv(
    gsea,
    file.path(results_dir, paste0(subtype_label, "_GSEA_", out_prefix, "_hallmark.csv"))
  )
  
  plot_title <- paste(subtype_title, "-", format_contrast_title(out_prefix))
  
  pvol <- plot_volcano(
    df = res_df,
    gene_col = "gene",
    lfc_col = "log2FoldChange",
    p_col = "pvalue",
    title = plot_title,
    p_cutoff = volcano_p_cutoff,
    lfc_cutoff = volcano_lfc_cutoff,
    label_top_n = 10
  )
  ggsave(
    file.path(plot_dir, paste0(subtype_label, "_volcano_", out_prefix, ".pdf")),
    pvol, width = 9, height = 7
  )
  
  pgsea <- plot_gsea_dot(
    gsea_df = gsea,
    title = paste0("GSEA: ", plot_title),
    p_cutoff = gsea_p_cutoff,
    top_n = 10,
    use_p = "pval"
  )
  
  if (!is.null(pgsea)) {
    ggsave(
      file.path(plot_dir, paste0(subtype_label, "_gsea_dot_", out_prefix, ".pdf")),
      pgsea, width = 10, height = 6
    )
  }
  
  list(de = res_df, gsea = gsea, volcano = pvol, gsea_dot = pgsea)
}

# Generic pseudobulk contrast wrapper
run_pseudobulk_contrast <- function(pb_counts,
                                    coldata,
                                    case_condition,
                                    ref_conditions,
                                    ref_label,
                                    hallmark_pathways,
                                    out_prefix,
                                    subtype_label,
                                    subtype_title,
                                    results_dir,
                                    plot_dir,
                                    minSize,
                                    maxSize,
                                    volcano_lfc_cutoff,
                                    volcano_p_cutoff,
                                    gsea_p_cutoff) {
  keep_conditions <- c(case_condition, ref_conditions)
  
  keep_samples <- rownames(coldata)[coldata$condition %in% keep_conditions]
  coldata_sub <- coldata[keep_samples, , drop = FALSE]
  pb_sub <- pb_counts[, keep_samples, drop = FALSE]
  
  coldata_sub$cond2 <- ifelse(coldata_sub$condition == case_condition, case_condition, ref_label)
  
  run_pseudobulk_fgsea(
    pb_counts = pb_sub,
    coldata = coldata_sub,
    cond2_levels = c(ref_label, case_condition),
    hallmark_pathways = hallmark_pathways,
    out_prefix = out_prefix,
    subtype_label = subtype_label,
    subtype_title = subtype_title,
    results_dir = results_dir,
    plot_dir = plot_dir,
    minSize = minSize,
    maxSize = maxSize,
    volcano_lfc_cutoff = volcano_lfc_cutoff,
    volcano_p_cutoff = volcano_p_cutoff,
    gsea_p_cutoff = gsea_p_cutoff
  )
}
