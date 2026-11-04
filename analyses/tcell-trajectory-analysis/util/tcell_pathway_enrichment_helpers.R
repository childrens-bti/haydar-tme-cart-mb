# Helpers for pathway enrichment of fitted T-cell lineage changes.
# Purpose: describe which biological programs change from each shared lineage's
# start to its endpoint, as one line of evidence for the reviewed lineage labels.
# Author: Bicna Song
# Date: 2026-10-25

# Validate the complete modeled-gene table and retain converged, finite endpoint
# changes for genes with a testable association in at least one condition of the
# lineage. Genes that are untestable in every condition have unsupported fitted
# endpoints (e.g. +173 log2) that would dominate weighted enrichment scores.
# Association p-values do not filter this input.
read_pathway_rankings <- function(file, qc_file, compartment, annotations) {
  if (!file.exists(file) || !file.exists(qc_file)) {
    stop("Missing lineage-ranking or gene-fit QC table for ", compartment)
  }
  ranks <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  qc <- read.delim(qc_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("compartment", "lineage", "gene", "end_minus_start", "wald_statistic") %in% names(ranks)) ||
      !all(c("compartment", "gene", "model_converged") %in% names(qc))) {
    stop("Incomplete lineage-ranking or gene-fit QC columns for ", compartment)
  }
  expected <- annotations$lineage[annotations$compartment == compartment]
  if (length(expected) == 0L || anyDuplicated(expected) ||
      anyNA(ranks[c("compartment", "lineage", "gene")]) ||
      anyNA(qc[c("compartment", "gene")]) ||
      any(ranks$compartment != compartment) || any(qc$compartment != compartment) ||
      any(!nzchar(ranks$gene)) || any(!nzchar(qc$gene)) ||
      anyDuplicated(ranks[c("lineage", "gene")]) || anyDuplicated(qc$gene) ||
      !setequal(unique(ranks$lineage), expected) ||
      !all(vapply(split(ranks$gene, ranks$lineage), setequal, logical(1), qc$gene))) {
    stop("Compartment, lineage or modeled-gene keys are inconsistent: ", file)
  }
  if (!is.numeric(ranks$end_minus_start) ||
      !(is.numeric(ranks$wald_statistic) || all(is.na(ranks$wald_statistic))) ||
      !is.logical(qc$model_converged)) {
    stop("Endpoint changes must be numeric and convergence must be logical")
  }
  ranks$model_converged <- qc$model_converged[match(ranks$gene, qc$gene)]
  ranks$analysis_status <- ifelse(is.na(ranks$model_converged) | !ranks$model_converged,
    "not_converged", ifelse(!is.finite(ranks$end_minus_start), "nonfinite_endpoint",
      ifelse(is.finite(ranks$wald_statistic), "used", "no_testable_association")))
  # This is the saved mean of available condition-specific log2 fitted
  # (end + 0.1)/(start + 0.1) ratios; upstream averaging omits missing values.
  # It is not a CAR-versus-CAR contrast.
  ranks$score <- ranks$end_minus_start
  summary <- ranks |>
    dplyr::group_by(.data$compartment, .data$lineage) |>
    dplyr::summarise(
      n_modeled = dplyr::n(), n_used = sum(.data$analysis_status == "used"),
      n_not_converged = sum(.data$analysis_status == "not_converged"),
      n_nonfinite_endpoint = sum(.data$analysis_status == "nonfinite_endpoint"),
      n_no_testable_association = sum(.data$analysis_status == "no_testable_association"),
      n_positive = sum(.data$score[.data$analysis_status == "used"] > 0),
      n_negative = sum(.data$score[.data$analysis_status == "used"] < 0),
      n_zero = sum(.data$score[.data$analysis_status == "used"] == 0),
      n_ties = sum(duplicated(.data$score[.data$analysis_status == "used"])),
      .groups = "drop"
    )
  if (any(summary$n_used < 15L) || any(summary$n_positive == 0L) || any(summary$n_negative == 0L)) {
    stop("Every lineage needs at least 15 finite genes and scores of both signs")
  }
  list(ranks = ranks[ranks$analysis_status == "used", ], summary = summary,
    excluded = ranks[ranks$analysis_status != "used",
      c("compartment", "lineage", "gene", "analysis_status")])
}

# Freeze exact mouse gene-set memberships as a compressed, versioned TSV. Subsequent
# renders read this file offline instead of querying a changing database.
load_mouse_pathways <- function(file, db_version = "2025.1.Mm") {
  if (file.exists(file)) {
    message("Using saved mouse gene sets: ", file)
    members <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    message("Downloading mouse MSigDB on first use; internet access is required")
    hallmark <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
    go_bp <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus",
                              collection = "M5", subcollection = "GO:BP")
    members <- dplyr::bind_rows(
      dplyr::transmute(hallmark, collection = "Hallmark", pathway = .data$gs_name,
                        gene = .data$gene_symbol, db_version = .data$db_version),
      dplyr::transmute(go_bp, collection = "GO_BP", pathway = .data$gs_name,
                        gene = .data$gene_symbol, db_version = .data$db_version)
    ) |>
      dplyr::filter(!is.na(.data$gene), nzchar(.data$gene)) |>
      dplyr::distinct() |>
      dplyr::arrange(.data$collection, .data$pathway, .data$gene)
  }
  needed <- c("collection", "pathway", "gene", "db_version")
  if (!all(needed %in% names(members)) || nrow(members) == 0L ||
      anyNA(members[needed]) || any(!nzchar(members$gene)) ||
      any(!nzchar(members$pathway)) ||
      !setequal(members$collection, c("Hallmark", "GO_BP")) ||
      !identical(unique(members$db_version), db_version) ||
      anyDuplicated(members[c("collection", "pathway", "gene")])) {
    stop("Invalid gene-set snapshot or wrong MSigDB release: ", file)
  }
  if (!file.exists(file)) {
    connection <- gzfile(file, open = "wt")
    tryCatch(write.table(members, connection, sep = "\t", quote = FALSE, row.names = FALSE),
      finally = close(connection))
    message("Saved mouse gene-set snapshot: ", file)
  }
  members
}

# Read the four CD8 T-cell signatures of Miller et al. 2019 (PMC6673650,
# supplementary table 1): genes up-regulated in each tumor-infiltrating CD8
# cluster, one sheet per cluster. Mouse symbols are used as given.
load_cd8_signatures <- function(file, collection = "Miller2019_CD8") {
  sheets <- c(`Progenitor Exh.` = "MILLER2019_PROGENITOR_EXHAUSTED",
              `Effector-like` = "MILLER2019_EFFECTOR_LIKE",
              `Terminally Exh.` = "MILLER2019_TERMINALLY_EXHAUSTED",
              Proliferating = "MILLER2019_PROLIFERATING")
  if (!file.exists(file) || !setequal(readxl::excel_sheets(file), names(sheets))) {
    stop("Missing CD8 signature table or unexpected sheets: ", file)
  }
  members <- purrr::imap_dfr(sheets, function(pathway, sheet) {
    genes <- readxl::read_excel(file, sheet = sheet)$gene
    tibble::tibble(collection = collection, pathway = pathway,
      gene = unique(genes[!is.na(genes) & nzchar(genes)]), db_version = "PMC6673650_SuppTable1")
  })
  if (any(table(members$pathway) < 15L)) stop("A CD8 signature has fewer than 15 genes")
  members
}

# Rank once, record all pathway overlaps, and test only sets represented by
# min_size..max_size genes. Size limits refer to the modeled ranking, not the
# full database set. Gene symbols remain case-sensitive mouse symbols. Fix
# ordering and the RNG algorithm so earlier R commands cannot change a job.
enrich_lineage_pathways <- function(ranks, members, compartment, lineage,
                                    collection, min_size, max_size, seed, workers, eps,
                                    n_perm_simple = 100000L) {
  ranks <- ranks[ranks$lineage == lineage, ]
  ranks <- ranks[order(-ranks$score, ranks$gene, method = "radix"), ]
  scores <- stats::setNames(ranks$score, ranks$gene)
  sets <- members[members$collection == collection, ]
  pathways <- split(sets$gene, sets$pathway)
  pathways <- pathways[order(names(pathways), method = "radix")]
  coverage <- tibble::tibble(
    pathway = names(pathways), n_genes_in_set = lengths(pathways),
    n_genes_in_ranking = vapply(pathways, function(x) sum(x %in% names(scores)), integer(1))
  ) |>
    dplyr::mutate(coverage_fraction = .data$n_genes_in_ranking / .data$n_genes_in_set,
      eligible = .data$n_genes_in_ranking >= min_size & .data$n_genes_in_ranking <= max_size)
  fit <- tibble::tibble(pathway = character(), ES = double(), NES = double(),
    p_value = double(), log2err = double(), leading_edge = character())
  if (any(coverage$eligible)) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
      sample.kind = "Rejection")
    bpparam <- if (workers == 1L) BiocParallel::SerialParam(RNGseed = seed) else
      BiocParallel::SnowParam(workers = workers, RNGseed = seed, type = "SOCK")
    raw <- fgsea::fgseaMultilevel(
      pathways = pathways[coverage$eligible], stats = scores,
      minSize = min_size, maxSize = max_size, eps = eps, sampleSize = 101L,
      scoreType = "std", gseaParam = 1, nPermSimple = n_perm_simple, BPPARAM = bpparam
    )
    fit <- tibble::tibble(pathway = raw$pathway, ES = raw$ES, NES = raw$NES,
      p_value = raw$pval, log2err = raw$log2err,
      leading_edge = vapply(raw$leadingEdge, paste, character(1), collapse = ";"))
  }
  coverage |>
    dplyr::left_join(fit, by = "pathway") |>
    dplyr::mutate(
      compartment = compartment, lineage = lineage, collection = collection,
      db_version = unique(sets$db_version), n_ranked_genes = length(scores),
      ranking_method = "mean_condition_log2_fitted_end_start_testable_genes",
      seed = seed,
      test_status = dplyr::case_when(!.data$eligible ~ "excluded_size",
        is.finite(.data$p_value) & is.finite(.data$NES) ~ "estimated",
        TRUE ~ "not_estimated"),
      direction = dplyr::case_when(.data$NES > 0 ~ "endpoint_higher",
        .data$NES < 0 ~ "endpoint_lower", TRUE ~ NA_character_)
    )
}

# One BH family spans all attempted eligible pathways AND all three lineages
# within a compartment/collection. Excluded-size sets are not extra tests;
# eligible but unestimated p-values stay NA and count conservatively in n.
adjust_pathway_fdr <- function(results) {
  results |>
    dplyr::group_by(.data$compartment, .data$collection) |>
    dplyr::group_modify(function(data, key) {
      data$p_adj <- NA_real_
      keep <- which(data$eligible)
      data$p_adj[keep] <- stats::p.adjust(data$p_value[keep], method = "BH", n = length(keep))
      data$n_tests_fdr_family <- length(keep)
      data
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(fdr_family = paste(.data$compartment, .data$collection, sep = "|"))
}

# Select supported review candidates only; never fill a sparse result with
# nominally significant or nonsignificant pathways and never relabel a lineage.
# Within each lineage/collection/direction, a pathway is skipped when its modeled
# genes overlap an already selected pathway at Jaccard >= max_jaccard, so related
# GO BP terms (e.g. cell killing and its regulation) do not fill the panel.
select_pathway_review <- function(results, members, ranked_genes, fdr = 0.05, n = 5L,
                                  max_jaccard = 0.5) {
  results |>
    dplyr::filter(.data$test_status == "estimated", is.finite(.data$p_adj),
      .data$p_adj < fdr, !is.na(.data$direction)) |>
    dplyr::arrange(.data$p_adj, dplyr::desc(abs(.data$NES)), .data$pathway) |>
    dplyr::group_by(.data$compartment, .data$lineage, .data$collection, .data$direction) |>
    dplyr::group_modify(function(data, key) {
      keep <- drop_redundant_pathways(data$pathway, members, ranked_genes, max_jaccard, n = n)
      data[data$pathway %in% keep, ]
    }) |>
    dplyr::ungroup()
}

# Greedily keep pathways in the given priority order, skipping any whose modeled
# genes overlap an already kept pathway at Jaccard >= max_jaccard.
drop_redundant_pathways <- function(pathways, members, ranked_genes, max_jaccard = 0.5,
                                    n = Inf) {
  members <- members[members$gene %in% ranked_genes & members$pathway %in% pathways, ]
  set_genes <- split(members$gene, members$pathway)
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  keep <- character()
  for (pathway in pathways) {
    redundant <- any(vapply(keep, function(k)
      jaccard(set_genes[[pathway]], set_genes[[k]]) >= max_jaccard, logical(1)))
    if (!redundant) keep <- c(keep, pathway)
    if (length(keep) == n) break
  }
  keep
}

# Summarise each lineage by its own selected pathways, without reference to any
# expected program, so the reviewer can judge whether they fit the assigned label.
# Counts cover every FDR-supported pathway; listed pathways are the non-redundant
# figure selection, ordered by FDR. Nothing here assigns or changes a lineage label.
summarise_lineage_pathways <- function(results, review, fdr = 0.05) {
  describe <- function(data) {
    if (!nrow(data)) return("")
    data <- dplyr::arrange(data, .data$p_adj)
    paste0(short_pathway_name(data$pathway), " (NES ", sprintf("%+.2f", data$NES),
      ", FDR ", signif(data$p_adj, 2), ")", collapse = "; ")
  }
  supported <- results[which(results$test_status == "estimated" & results$p_adj < fdr), ]
  dplyr::distinct(results, .data$compartment, .data$lineage, .data$assigned_program) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_fdr_higher = sum(supported$lineage == .data$lineage & supported$NES > 0),
      n_fdr_lower = sum(supported$lineage == .data$lineage & supported$NES < 0),
      top = list(review[review$lineage == .data$lineage, ])) |>
    dplyr::mutate(
      hallmark_higher = describe(top[top$collection == "Hallmark" & top$NES > 0, ]),
      hallmark_lower = describe(top[top$collection == "Hallmark" & top$NES < 0, ]),
      go_bp_higher = describe(top[top$collection == "GO_BP" & top$NES > 0, ]),
      go_bp_lower = describe(top[top$collection == "GO_BP" & top$NES < 0, ])) |>
    dplyr::ungroup() |>
    dplyr::select(-"top") |>
    dplyr::arrange(.data$compartment, .data$lineage)
}

short_pathway_name <- function(pathway) {
  gsub("_", " ", sub("^(HALLMARK|GOBP|MILLER2019)_", "", pathway))
}

lineage_label <- function(lineage, program) {
  paste0(sub("slingPseudotime_", "L", lineage, fixed = TRUE), "\n",
    stringr::str_wrap(gsub("/", "/ ", program, fixed = TRUE), width = 16))
}

# Build matrix rows for one compartment and collection: the selected review
# pathways, grouped by the lineage that selected them. Every row shows all three
# lineages, so a pathway selected by several lineages is drawn once, under the
# first of them (L1, L2, L3), and a pathway redundant with an earlier row
# (Jaccard >= max_jaccard) is not repeated.
pathway_matrix_rows <- function(results, review, compartment, collection,
                                members, ranked_genes, max_jaccard = 0.5) {
  results <- results[results$compartment == compartment & results$collection == collection, ]
  selected <- review[review$compartment == compartment & review$collection == collection, ] |>
    dplyr::arrange(.data$lineage, .data$p_adj, dplyr::desc(abs(.data$NES))) |>
    dplyr::distinct(.data$pathway, .keep_all = TRUE)
  keep <- drop_redundant_pathways(selected$pathway, members[members$collection == collection, ],
    ranked_genes, max_jaccard)
  selected <- selected[selected$pathway %in% keep, ]
  programs <- dplyr::distinct(results, .data$lineage, .data$assigned_program)
  selected$program_category <- paste0("Top in ",
    sub("slingPseudotime_", "L", selected$lineage, fixed = TRUE), ": ",
    programs$assigned_program[match(selected$lineage, programs$lineage)])
  rows <- dplyr::select(selected, "program_category", "pathway") |>
    dplyr::mutate(row_order = dplyr::row_number())
  rows |>
    dplyr::inner_join(results, by = "pathway", relationship = "one-to-many") |>
    dplyr::mutate(
      lineage_label = lineage_label(.data$lineage, .data$assigned_program),
      program_category = factor(.data$program_category, levels = unique(rows$program_category)),
      pathway_label = factor(stringr::str_wrap(short_pathway_name(.data$pathway), 45),
        levels = rev(unique(stringr::str_wrap(short_pathway_name(rows$pathway), 45)))))
}

# Matrix rows for a prespecified collection: every set in the given order,
# whether or not it is FDR-supported, in one row group.
signature_matrix_rows <- function(results, compartment, collection, pathways, group) {
  results <- results[results$compartment == compartment & results$collection == collection, ]
  rows <- tibble::tibble(program_category = group, pathway = pathways,
    row_order = seq_along(pathways))
  rows |>
    dplyr::inner_join(results, by = "pathway", relationship = "one-to-many") |>
    dplyr::mutate(
      lineage_label = lineage_label(.data$lineage, .data$assigned_program),
      program_category = factor(.data$program_category, levels = group),
      pathway_label = factor(short_pathway_name(.data$pathway),
        levels = rev(short_pathway_name(pathways))))
}

# Pathway-by-lineage dot matrix. Filled dots are FDR-supported (color = NES,
# size = FDR); open grey dots were tested but not supported; blank cells were
# not tested at the size limits. NES is clamped to a fixed range so colors are
# comparable across lineages, compartments and collections. compact = TRUE fixes
# the panel size per lineage and row, so pages with few rows are not stretched
# to fill the page.
plot_pathway_matrix <- function(data, title, fdr = 0.05, nes_limit = 2.5, compact = FALSE) {
  data <- data |>
    dplyr::filter(.data$test_status == "estimated") |>
    dplyr::mutate(
      significant = is.finite(.data$p_adj) & .data$p_adj < fdr,
      fdr_score = pmin(-log10(pmax(.data$p_adj, 1e-8)), 8),
      NES_plot = pmax(pmin(.data$NES, nes_limit), -nes_limit))
  panel_size <- if (compact) {
    rows_per_group <- tapply(as.character(data$pathway_label), data$program_category,
      dplyr::n_distinct)
    ggplot2::theme(panel.widths = grid::unit(dplyr::n_distinct(data$lineage_label), "in"),
      # A minimum height leaves room for a wrapped row-group label.
      panel.heights = grid::unit(pmax(0.3 * rows_per_group[!is.na(rows_per_group)], 0.55), "in"))
  }
  ggplot2::ggplot(data, ggplot2::aes(.data$lineage_label, .data$pathway_label)) +
    # Train the axes on all rows first so row order follows the factor levels,
    # not the order in which the two dot layers first meet each pathway.
    ggplot2::geom_blank() +
    ggplot2::geom_point(data = data[!data$significant, ],
      shape = 1, size = 2, color = "grey65", stroke = 0.5) +
    ggplot2::geom_point(data = data[data$significant, ],
      ggplot2::aes(size = .data$fdr_score, fill = .data$NES_plot),
      shape = 21, color = "grey20", stroke = 0.3) +
    ggplot2::facet_grid(program_category ~ ., scales = "free_y", space = "free_y",
      labeller = ggplot2::label_wrap_gen(width = 22)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948",
      midpoint = 0, limits = c(-nes_limit, nes_limit), name = "NES") +
    ggplot2::scale_size_continuous(range = c(2.5, 7.5), limits = c(-log10(fdr), 8),
      breaks = c(2, 4, 8), labels = c("0.01", "1e-4", "\u2264 1e-8"), name = "FDR") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 8.5),
      axis.text.x.top = ggplot2::element_text(size = 9, lineheight = 0.9),
      axis.text.y = ggplot2::element_text(size = 8.5, lineheight = 0.85),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", color = NA)) +
    panel_size +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}
