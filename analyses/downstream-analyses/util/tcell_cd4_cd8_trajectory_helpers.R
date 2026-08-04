# Helper functions for CD4/CD8 T-cell trajectory analysis.
#
# These functions support the CD4/CD8 annotation and Slingshot reports by
# keeping marker scoring, annotation summaries, clustering, and trajectory
# execution outside the main reports.

# Return features that are present in the Seurat object.
available_features <- function(obj, features) {
  intersect(features, rownames(obj))
}

# Convert a free-text label into a file-safe label.
safe_label <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# Assign stable colors from the complete set of cluster-level T-cell labels.
make_tcell_subtype_palette <- function(labels) {
  label_levels <- sort(unique(stats::na.omit(as.character(labels))))

  if (length(label_levels) == 0) {
    return(stats::setNames(character(), character()))
  }

  stats::setNames(scales::hue_pal()(length(label_levels)), label_levels)
}

# Calculate an average log-normalized expression score for one marker set.
score_marker_set <- function(obj, features, assay = SeuratObject::DefaultAssay(obj)) {
  features <- available_features(obj, features)

  if (length(features) == 0) {
    return(rep(NA_real_, ncol(obj)))
  }

  mat <- get_log_data(obj, assay = assay)[features, , drop = FALSE]
  Matrix::colMeans(mat)
}

# Return the metadata column used for one marker-program score.
marker_score_column <- function(program_name) {
  if (program_name %in% c("CD4", "CD8")) {
    return(paste0(tolower(program_name), "_score"))
  }

  paste0("state_", program_name, "_score")
}

aucell_score_column <- function(program_name) {
  paste0("aucell_", program_name, "_score")
}

# Add marker-program scores used for CD4/CD8 calling and functional-state summaries.
add_marker_program_scores <- function(obj) {
  for (program_name in names(marker_programs)) {
    obj[[marker_score_column(program_name)]] <- score_marker_set(
      obj,
      marker_programs[[program_name]]
    )
  }

  obj
}

add_aucell_program_scores <- function(obj, aucell_programs) {
  aucell_programs <- purrr::map(aucell_programs, ~ available_features(obj, .x))
  aucell_programs <- aucell_programs[lengths(aucell_programs) > 0]

  if (length(aucell_programs) == 0) {
    return(obj)
  }

  expr <- get_log_data(obj)
  rankings <- AUCell::AUCell_buildRankings(
    expr,
    plotStats = FALSE,
    verbose = FALSE
  )
  auc <- AUCell::AUCell_calcAUC(
    aucell_programs,
    rankings,
    aucMaxRank = ceiling(0.05 * nrow(rankings)),
    verbose = FALSE
  )
  auc_mat <- SummarizedExperiment::assay(auc)

  for (program_name in rownames(auc_mat)) {
    obj[[aucell_score_column(program_name)]] <- as.numeric(auc_mat[program_name, colnames(obj)])
  }

  obj
}

# Return marker-program names used for functional-state annotation.
functional_state_program_names <- function() {
  setdiff(names(marker_programs), c("CD4", "CD8"))
}

# Convert the top scoring state program into a readable CD4/CD8 state label.
format_functional_state_label <- function(program_name, compartment) {
  if (is.na(program_name)) {
    return(paste(compartment, "T cells"))
  }

  cd4_labels <- c(
    naive_memory = "Naive / memory CD4 T cells",
    activated = "Activated CD4 T cells",
    treg = "Treg",
    cytotoxic = "Cytotoxic-like CD4 T cells",
    exhausted = "Exhausted/dysfunctional CD4 T cells",
    proliferating = "Proliferating CD4 T cells",
    ifn_response = "IFN-responsive CD4 T cells",
    stress_response = "Stress-response CD4 T cells"
  )

  cd8_labels <- c(
    naive_memory = "Naive / memory CD8 T cells",
    activated = "Activated CD8 T cells",
    treg = "Treg-marker-high CD8-like T cells",
    cytotoxic = "Cytotoxic CD8 T cells",
    exhausted = "Exhausted/dysfunctional CD8 T cells",
    proliferating = "Proliferating CD8 T cells",
    ifn_response = "IFN-responsive CD8 T cells",
    stress_response = "Stress-response CD8 T cells"
  )

  label_map <- if (compartment == "CD4-like") cd4_labels else cd8_labels
  label <- unname(label_map[[program_name]])

  if (is.null(label) || is.na(label)) {
    return(paste(compartment, "T cells"))
  }

  label
}

# Classify cells as CD4-like, CD8-like, or ambiguous using marker scores.
classify_cd4_cd8 <- function(obj) {
  obj <- add_marker_program_scores(obj)

  obj$cd4_cd8_marker_group <- dplyr::case_when(
    !is.na(obj$cd4_score) &
      obj$cd4_score >= min_marker_score &
      obj$cd4_score - obj$cd8_score >= min_score_delta ~ "CD4-like",
    !is.na(obj$cd8_score) &
      obj$cd8_score >= min_marker_score &
      obj$cd8_score - obj$cd4_score >= min_score_delta ~ "CD8-like",
    TRUE ~ "Ambiguous"
  )

  obj$cd4_cd8_marker_group <- factor(
    obj$cd4_cd8_marker_group,
    levels = c("CD4-like", "CD8-like", "Ambiguous")
  )

  obj
}

# Get the generic mouse CD4T or CD8T ScGate model.
get_generic_mouse_scgate_model <- function(model_type = c("CD4", "CD8")) {
  model_type <- match.arg(model_type)
  model_name <- paste0(model_type, "T")

  model_db <- scGate::get_scGateDB()

  if (is.null(model_db) || !"mouse" %in% names(model_db)) {
    stop("The ScGate model database does not contain mouse models.")
  }

  if (!"generic" %in% names(model_db$mouse)) {
    stop("The mouse ScGate database does not contain a generic model group.")
  }

  if (!model_name %in% names(model_db$mouse$generic)) {
    stop("The mouse generic ScGate database does not contain model: ", model_name)
  }

  model_db$mouse$generic[[model_name]]
}

# Run ScGate CD4 and CD8 models as supporting CD4/CD8 annotation.
# This intentionally requires the generic mouse ScGate CD4T/CD8T models.
run_scgate_cd4_cd8 <- function(obj) {
  if (!requireNamespace("scGate", quietly = TRUE)) {
    stop("scGate is not installed. Install scGate before running CD4/CD8 ScGate support.")
  }

  cd4_model <- get_generic_mouse_scgate_model("CD4")
  cd8_model <- get_generic_mouse_scgate_model("CD8")
  obj$scgate_cd4_model <- "mouse_generic_CD4T"
  obj$scgate_cd8_model <- "mouse_generic_CD8T"

  obj$scgate_cd4 <- NA_character_
  obj$scgate_cd8 <- NA_character_

  obj_cd4 <- scGate::scGate(data = obj, model = cd4_model)

  if ("is.pure" %in% colnames(obj_cd4@meta.data)) {
    obj$scgate_cd4 <- as.character(obj_cd4$is.pure)
  }

  obj_cd8 <- scGate::scGate(data = obj, model = cd8_model)

  if ("is.pure" %in% colnames(obj_cd8@meta.data)) {
    obj$scgate_cd8 <- as.character(obj_cd8$is.pure)
  }

  obj$scgate_cd4_cd8_group <- dplyr::case_when(
    obj$scgate_cd4 == "Pure" & obj$scgate_cd8 != "Pure" ~ "CD4-like",
    obj$scgate_cd8 == "Pure" & obj$scgate_cd4 != "Pure" ~ "CD8-like",
    TRUE ~ "Ambiguous"
  )

  obj
}

# Run SingleR for one mouse reference at the cluster level.
run_singler_reference <- function(obj, ref, ref_name, label_col = "label.fine") {
  sce <- Seurat::as.SingleCellExperiment(obj)
  clusters <- as.character(obj$seurat_clusters)

  pred <- tryCatch(
    SingleR::SingleR(
      test = sce,
      ref = ref,
      labels = ref[[label_col]],
      clusters = clusters
    ),
    error = function(e) {
      warning("SingleR failed for ", ref_name, ": ", conditionMessage(e))
      NULL
    }
  )

  label_column <- paste0("singler_", ref_name, "_label")
  pruned_column <- paste0("singler_", ref_name, "_pruned_label")

  obj@meta.data[[label_column]] <- NA_character_
  obj@meta.data[[pruned_column]] <- NA_character_

  if (is.null(pred)) {
    return(obj)
  }

  pred_df <- as.data.frame(pred) %>%
    tibble::rownames_to_column("seurat_clusters") %>%
    dplyr::select(seurat_clusters, labels, pruned.labels)

  label_map <- stats::setNames(pred_df$labels, pred_df$seurat_clusters)
  pruned_map <- stats::setNames(pred_df$pruned.labels, pred_df$seurat_clusters)

  obj@meta.data[[label_column]] <- unname(label_map[clusters])
  obj@meta.data[[pruned_column]] <- unname(pruned_map[clusters])

  obj
}

# Run SingleR with both mouse references and store cluster-level labels.
run_singler_annotations <- function(obj) {
  immgen_ref <- celldex::ImmGenData()
  mouse_rnaseq_ref <- celldex::MouseRNAseqData()

  obj <- run_singler_reference(
    obj,
    ref = immgen_ref,
    ref_name = "immgen",
    label_col = "label.fine"
  )

  obj <- run_singler_reference(
    obj,
    ref = mouse_rnaseq_ref,
    ref_name = "mousernaseq",
    label_col = "label.fine"
  )

  obj
}

# Run ProjecTILs classifier with the default mouse TIL reference.
run_projectils_annotation <- function(obj) {
  obj$projectils_label <- NA_character_
  obj$projectils_confidence <- NA_real_
  obj$projectils_label_column <- NA_character_
  obj$projectils_confidence_column <- NA_character_

  if (!requireNamespace("ProjecTILs", quietly = TRUE)) {
    warning("ProjecTILs is not installed. Skipping ProjecTILs annotation.")
    return(obj)
  }

  if ("RNA" %in% names(obj@assays)) {
    SeuratObject::DefaultAssay(obj) <- "RNA"
  }

  projectils_data <- new.env(parent = emptyenv())
  suppressWarnings(
    utils::data(
      "Hs2Mm.convert.table",
      package = "ProjecTILs",
      envir = projectils_data
    )
  )

  if (!exists("Hs2Mm.convert.table", envir = projectils_data)) {
    projectils_data_path <- file.path(
      system.file(package = "ProjecTILs"),
      "data",
      "Hs2Mm.convert.table.RData"
    )

    if (file.exists(projectils_data_path)) {
      load(projectils_data_path, envir = projectils_data)
    }
  }

  if (!exists("Hs2Mm.convert.table", envir = projectils_data)) {
    warning("ProjecTILs ortholog table Hs2Mm.convert.table is unavailable. Skipping ProjecTILs annotation.")
    return(obj)
  }

  obj_projectils <- tryCatch(
    {
      options(timeout = max(getOption("timeout"), 3000))
      ref <- ProjecTILs::load.reference.map()
      ProjecTILs::ProjecTILs.classifier(
        query = obj,
        ref = ref,
        filter.cells = FALSE,
        labels.col = "functional.cluster",
        ncores = 1,
        ortholog_table = projectils_data$Hs2Mm.convert.table
      )
    },
    error = function(e) {
      warning("ProjecTILs annotation failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(obj_projectils)) {
    return(obj)
  }

  common_cells <- intersect(colnames(obj), colnames(obj_projectils))
  if (length(common_cells) == 0) {
    warning("ProjecTILs returned an object with no cell names in common with the input object.")
    return(obj)
  }

  projectils_meta_cols <- colnames(obj_projectils@meta.data)
  label_candidates <- c(
    "functional.cluster",
    "functional.cluster.label",
    "functional_cluster",
    "ProjecTILs.functional.cluster",
    "projectils_label",
    "predicted.celltype",
    "predicted.id",
    "celltype",
    "cell_type"
  )
  label_col <- label_candidates[label_candidates %in% projectils_meta_cols][1]

  if (is.na(label_col)) {
    label_col <- projectils_meta_cols[
      stringr::str_detect(
        tolower(projectils_meta_cols),
        "functional.*cluster|cluster.*functional|projectils.*label|predicted|cell.*type"
      )
    ][1]
  }

  if (!is.na(label_col)) {
    projectils_label_map <- stats::setNames(
      as.character(obj_projectils@meta.data[[label_col]]),
      colnames(obj_projectils)
    )
    obj$projectils_label[common_cells] <- unname(projectils_label_map[common_cells])
    obj$projectils_label_column[common_cells] <- label_col
  } else {
    warning(
      "ProjecTILs annotation completed, but no label column was found. Available metadata columns: ",
      paste(projectils_meta_cols, collapse = ", ")
    )
  }

  confidence_candidates <- c(
    "functional.cluster.conf",
    "functional.cluster.confidence",
    "functional_cluster_conf",
    "ProjecTILs.functional.cluster.conf",
    "projectils_confidence",
    "prediction.score.max",
    "predicted.score",
    "confidence"
  )
  confidence_col <- confidence_candidates[confidence_candidates %in% projectils_meta_cols][1]

  if (is.na(confidence_col)) {
    confidence_col <- projectils_meta_cols[
      stringr::str_detect(tolower(projectils_meta_cols), "conf|score|prob") &
        vapply(obj_projectils@meta.data, is.numeric, logical(1))
    ][1]
  }

  if (!is.na(confidence_col)) {
    projectils_confidence_map <- stats::setNames(
      as.numeric(obj_projectils@meta.data[[confidence_col]]),
      colnames(obj_projectils)
    )
    obj$projectils_confidence[common_cells] <- unname(projectils_confidence_map[common_cells])
    obj$projectils_confidence_column[common_cells] <- confidence_col
  }

  if (all(is.na(obj$projectils_label))) {
    warning("ProjecTILs annotation did not assign any non-NA labels.")
  }

  obj
}

# Summarize marker, ProjecTILs, SingleR, and ScGate evidence for each cluster annotation.
summarize_cluster_annotation <- function(obj, compartment) {
  marker_score_cols <- purrr::map_chr(names(marker_programs), marker_score_column)
  marker_score_cols <- intersect(marker_score_cols, colnames(obj@meta.data))
  state_score_cols <- purrr::map_chr(functional_state_program_names(), marker_score_column)
  state_score_cols <- intersect(state_score_cols, colnames(obj@meta.data))

  marker_summary <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(
        dplyr::all_of(marker_score_cols),
        ~ median(.x, na.rm = TRUE),
        .names = "median_{.col}"
      ),
      .groups = "drop"
    )

  marker_long <- marker_summary %>%
    tidyr::pivot_longer(
      cols = starts_with("median_"),
      names_to = "marker_program",
      values_to = "median_score"
    ) %>%
    dplyr::mutate(
      marker_program = stringr::str_remove(marker_program, "^median_"),
      marker_program = stringr::str_remove(marker_program, "_score$")
    ) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = median_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(seurat_clusters, top_marker_program = marker_program, top_marker_score = median_score)

  state_summary <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(state_score_cols),
        ~ median(.x, na.rm = TRUE),
        .names = "median_{.col}"
      ),
      .groups = "drop"
    )

  state_long <- state_summary %>%
    tidyr::pivot_longer(
      cols = starts_with("median_state_"),
      names_to = "functional_state_program",
      values_to = "top_functional_state_score"
    ) %>%
    dplyr::mutate(
      functional_state_program = stringr::str_remove(functional_state_program, "^median_state_"),
      functional_state_program = stringr::str_remove(functional_state_program, "_score$")
    ) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = top_functional_state_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      functional_state_label = purrr::map_chr(
        functional_state_program,
        format_functional_state_label,
        compartment = compartment
      )
    ) %>%
    dplyr::rename(top_functional_state_program = functional_state_program)

  previous_label_summary <- tibble::tibble(
    seurat_clusters = character(),
    top_previous_tcell_subtype = character(),
    top_previous_tcell_subtype_n = integer(),
    top_previous_tcell_subtype_prop = numeric()
  )

  if ("previous_tcell_subtype" %in% colnames(obj@meta.data)) {
    previous_label_summary <- obj@meta.data %>%
      dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
      dplyr::filter(!is.na(previous_tcell_subtype)) %>%
      dplyr::count(seurat_clusters, previous_tcell_subtype, name = "previous_label_n") %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::mutate(previous_label_prop = previous_label_n / sum(previous_label_n)) %>%
      dplyr::slice_max(order_by = previous_label_n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        seurat_clusters,
        top_previous_tcell_subtype = previous_tcell_subtype,
        top_previous_tcell_subtype_n = previous_label_n,
        top_previous_tcell_subtype_prop = previous_label_prop
      )
  }

  singler_cols <- intersect(
    c("singler_immgen_label", "singler_mousernaseq_label"),
    colnames(obj@meta.data)
  )

  singler_summary <- purrr::map_dfr(
    singler_cols,
    function(label_col) {
      obj@meta.data %>%
        dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
        dplyr::filter(!is.na(.data[[label_col]])) %>%
        dplyr::distinct(seurat_clusters, .data[[label_col]]) %>%
        dplyr::transmute(
          seurat_clusters,
          reference = label_col,
          cluster_label = .data[[label_col]]
        )
    }
  )

  if (nrow(singler_summary) > 0) {
    singler_summary <- singler_summary %>%
      tidyr::pivot_wider(
        names_from = reference,
        values_from = cluster_label,
        names_glue = "{reference}"
      )
  } else {
    singler_summary <- tibble::tibble(seurat_clusters = character())
  }

  projectils_summary <- tibble::tibble(
    seurat_clusters = character(),
    projectils_label = character(),
    projectils_label_n = integer(),
    projectils_label_prop = numeric(),
    median_projectils_confidence = numeric()
  )

  if ("projectils_label" %in% colnames(obj@meta.data)) {
    projectils_summary <- obj@meta.data %>%
      dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
      dplyr::filter(!is.na(projectils_label)) %>%
      dplyr::count(seurat_clusters, projectils_label, name = "projectils_label_n") %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::mutate(projectils_label_prop = projectils_label_n / sum(projectils_label_n)) %>%
      dplyr::slice_max(order_by = projectils_label_n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
        obj@meta.data %>%
          dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
          dplyr::group_by(seurat_clusters) %>%
          dplyr::summarise(
            median_projectils_confidence = median(projectils_confidence, na.rm = TRUE),
            .groups = "drop"
          ),
        by = "seurat_clusters"
      )
  }

  scgate_summary <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::count(seurat_clusters, scgate_cd4_cd8_group, name = "scgate_n") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = scgate_n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(seurat_clusters, top_scgate_group = scgate_cd4_cd8_group, top_scgate_n = scgate_n)

  marker_summary %>%
    dplyr::left_join(marker_long, by = "seurat_clusters") %>%
    dplyr::left_join(state_long, by = "seurat_clusters") %>%
    dplyr::left_join(previous_label_summary, by = "seurat_clusters") %>%
    dplyr::left_join(projectils_summary, by = "seurat_clusters") %>%
    dplyr::left_join(singler_summary, by = "seurat_clusters") %>%
    dplyr::left_join(scgate_summary, by = "seurat_clusters")
}

# Add annotation metadata to a clustered object and summarize clusters.
annotate_clustered_object <- function(obj, compartment) {
  obj <- add_marker_program_scores(obj)
  obj <- run_singler_annotations(obj)

  if (!"scgate_cd4_cd8_group" %in% colnames(obj@meta.data)) {
    obj <- run_scgate_cd4_cd8(obj)
  }

  cluster_annotation <- summarize_cluster_annotation(
    obj,
    compartment = compartment
  )

  obj$cd4_cd8_all_cluster <- as.character(obj$seurat_clusters)

  functional_state_map <- stats::setNames(
    cluster_annotation$functional_state_label,
    cluster_annotation$seurat_clusters
  )

  obj$functional_state_label <- unname(
    functional_state_map[as.character(obj$seurat_clusters)]
  )
  projectils_label_map <- stats::setNames(
    cluster_annotation$projectils_label,
    cluster_annotation$seurat_clusters
  )
  projectils_label_prop_map <- stats::setNames(
    cluster_annotation$projectils_label_prop,
    cluster_annotation$seurat_clusters
  )

  obj$cluster_projectils_label <- unname(
    projectils_label_map[as.character(obj$seurat_clusters)]
  )
  obj$cluster_projectils_label_prop <- unname(
    projectils_label_prop_map[as.character(obj$seurat_clusters)]
  )

  previous_subtype_map <- stats::setNames(
    cluster_annotation$top_previous_tcell_subtype,
    cluster_annotation$seurat_clusters
  )
  previous_subtype_prop_map <- stats::setNames(
    cluster_annotation$top_previous_tcell_subtype_prop,
    cluster_annotation$seurat_clusters
  )

  obj$top_previous_tcell_subtype <- unname(
    previous_subtype_map[as.character(obj$seurat_clusters)]
  )
  obj$top_previous_tcell_subtype_prop <- unname(
    previous_subtype_prop_map[as.character(obj$seurat_clusters)]
  )

  list(
    obj = obj,
    cluster_annotation = cluster_annotation
  )
}

# Cluster one CD4-like or CD8-like object before annotation and trajectory.
cluster_subset <- function(obj, label, resolution = 0.4, dims = 20) {
  DefaultAssay(obj) <- "RNA"
  dims <- min(dims, ncol(obj) - 1)
  pca_npcs <- min(max(30, dims), ncol(obj) - 1)

  if (dims < 2) {
    stop("At least 3 cells are required for PCA-based clustering.")
  }

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = pca_npcs, verbose = FALSE)

  reduction_name <- "pca"
  obj <- FindNeighbors(obj, reduction = reduction_name, dims = seq_len(dims), verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction_name, dims = seq_len(dims), verbose = FALSE)
  obj$cd4_cd8_condition_label <- label
  obj@misc$cd4_cd8_cluster_dims <- dims
  obj@misc$cd4_cd8_cluster_reduction <- reduction_name
  obj
}

# Choose the Slingshot root cluster using ProjecTILs naive-like labels first,
# then the naive/memory marker score if ProjecTILs support is unavailable.
get_root_cluster_info <- function(obj) {
  if (!"cluster_projectils_label" %in% colnames(obj@meta.data)) {
    obj$cluster_projectils_label <- NA_character_
  }

  cluster_scores <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      median_naive_memory_score = median(state_naive_memory_score, na.rm = TRUE),
      projectils_label = dplyr::first(na.omit(cluster_projectils_label)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      projectils_naive_like = stringr::str_detect(
        tolower(dplyr::coalesce(projectils_label, "")),
        "naive"
      ),
      has_naive_memory_score = is.finite(median_naive_memory_score) &
        median_naive_memory_score > 0
    )

  projectils_candidates <- cluster_scores %>%
    dplyr::filter(projectils_naive_like) %>%
    dplyr::arrange(
      dplyr::desc(median_naive_memory_score),
      seurat_clusters
    )

  if (nrow(projectils_candidates) > 0) {
    root <- projectils_candidates[1, , drop = FALSE]
    return(tibble::tibble(
      root_cluster = root$seurat_clusters,
      root_cluster_selection = "ProjecTILs naive/memory label",
      root_cluster_projectils_label = root$projectils_label,
      root_cluster_naive_memory_score = root$median_naive_memory_score
    ))
  }

  score_candidates <- cluster_scores %>%
    dplyr::filter(has_naive_memory_score) %>%
    dplyr::arrange(
      dplyr::desc(median_naive_memory_score),
      seurat_clusters
    )

  if (nrow(score_candidates) > 0) {
    root <- score_candidates[1, , drop = FALSE]
    return(tibble::tibble(
      root_cluster = root$seurat_clusters,
      root_cluster_selection = "naive/memory marker score",
      root_cluster_projectils_label = root$projectils_label,
      root_cluster_naive_memory_score = root$median_naive_memory_score
    ))
  }

  tibble::tibble(
    root_cluster = NA_character_,
    root_cluster_selection = NA_character_,
    root_cluster_projectils_label = NA_character_,
    root_cluster_naive_memory_score = NA_real_
  )
}

# Mark the selected Slingshot root cluster in the per-cluster QC table.
add_root_cluster_qc <- function(cluster_counts, root_cluster, root_projectils_label) {
  cluster_counts %>%
    dplyr::mutate(
      is_root_cluster = seurat_clusters == root_cluster,
      root_cluster_projectils_label = dplyr::if_else(
        is_root_cluster,
        root_projectils_label,
        NA_character_
      )
    )
}

# Extract cell-aligned pseudotime and curve weights for tradeSeq.
get_slingshot_trade_seq_inputs <- function(obj, plotting_pseudotime, reduction_name = "SLING") {
  slingshot_result <- obj@misc$slingshot[[reduction_name]]

  if (is.null(slingshot_result$SlingshotDataSet)) {
    stop("SlingshotDataSet is unavailable for reduction ", reduction_name)
  }

  trade_seq_pseudotime <- as.matrix(
    slingshot::slingPseudotime(slingshot_result$SlingshotDataSet, na = FALSE)
  )
  curve_weights <- as.matrix(
    slingshot::slingCurveWeights(slingshot_result$SlingshotDataSet)
  )
  cell_order <- rownames(plotting_pseudotime)

  align_cells <- function(x, value_name) {
    if (!is.null(rownames(x)) && all(cell_order %in% rownames(x))) {
      return(x[cell_order, , drop = FALSE])
    }
    if (nrow(x) == length(cell_order)) {
      rownames(x) <- cell_order
      return(x)
    }
    stop(value_name, " and plotting pseudotime have different cell counts")
  }

  trade_seq_pseudotime <- align_cells(trade_seq_pseudotime, "tradeSeq pseudotime")
  curve_weights <- align_cells(curve_weights, "Slingshot curve weights")

  if (ncol(trade_seq_pseudotime) != ncol(plotting_pseudotime) ||
      ncol(curve_weights) != ncol(plotting_pseudotime)) {
    stop("Slingshot pseudotime and curve weights have different lineage counts")
  }

  colnames(trade_seq_pseudotime) <- colnames(plotting_pseudotime)
  colnames(curve_weights) <- colnames(plotting_pseudotime)

  list(
    pseudotime = trade_seq_pseudotime,
    cell_weights = curve_weights
  )
}

# Assign each cell to the lineage with its largest Slingshot curve weight.
get_slingshot_lineage_assignments <- function(curve_weights) {
  weight_matrix <- curve_weights
  weight_matrix[is.na(weight_matrix)] <- 0
  has_lineage <- rowSums(weight_matrix) > 0
  primary_lineage <- rep(NA_character_, nrow(weight_matrix))
  primary_lineage[has_lineage] <- colnames(weight_matrix)[
    max.col(weight_matrix[has_lineage, , drop = FALSE], ties.method = "first")
  ]

  tibble::tibble(
    cell = rownames(weight_matrix),
    primary_lineage = primary_lineage
  )
}

# Keep Seurat clusters with enough cells and summarize Slingshot eligibility.
prepare_slingshot_subset <- function(obj, compartment, condition) {
  obj$seurat_clusters <- droplevels(factor(obj$seurat_clusters))

  cluster_counts <- as.data.frame(table(obj$seurat_clusters))
  colnames(cluster_counts) <- c("seurat_clusters", "n")

  cluster_counts <- cluster_counts %>%
    dplyr::mutate(
      seurat_clusters = as.character(seurat_clusters),
      retained_for_slingshot = n >= min_cells_per_cluster,
      compartment = compartment,
      condition = condition
    )

  retained_groups <- cluster_counts %>%
    dplyr::filter(retained_for_slingshot) %>%
    dplyr::pull(seurat_clusters)

  retained_cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% retained_groups]
  n_retained_cells <- length(retained_cells)
  n_retained_groups <- length(retained_groups)
  slingshot_run <- n_retained_cells >= min_cells_for_slingshot &&
    n_retained_groups >= min_clusters_for_slingshot
  skip_reason <- dplyr::case_when(
    n_retained_cells < min_cells_for_slingshot ~ "Too few cells after removing small clusters",
    n_retained_groups < min_clusters_for_slingshot ~ "Too few retained clusters",
    TRUE ~ NA_character_
  )

  eligibility <- tibble::tibble(
    compartment = compartment,
    condition = condition,
    total_slingshot_clusters = dplyr::n_distinct(as.character(obj$seurat_clusters)),
    min_cells_for_slingshot = min_cells_for_slingshot,
    min_clusters_for_slingshot = min_clusters_for_slingshot,
    min_cells_per_cluster = min_cells_per_cluster,
    skip_reason = skip_reason,
    root_cluster = NA_character_,
    root_cluster_selection = NA_character_,
    root_cluster_naive_memory_score = NA_real_
  )

  if (!slingshot_run) {
    return(list(
      obj = obj,
      cluster_counts = cluster_counts,
      eligibility = eligibility,
      slingshot_run = slingshot_run
    ))
  }

  obj_filtered <- subset(obj, cells = retained_cells)
  obj_filtered$seurat_clusters <- droplevels(factor(obj_filtered$seurat_clusters))

  list(
    obj = obj_filtered,
    cluster_counts = cluster_counts,
    eligibility = eligibility,
    slingshot_run = slingshot_run
  )
}

# Run Slingshot on one condition subset or all-condition compartment using Seurat clusters.
run_condition_slingshot <- function(
  obj,
  compartment,
  condition,
  gene_trend_genes,
  aucell_programs = NULL,
  plot_aucell_pseudotime = FALSE,
  subtype_colors = NULL,
  cache_file = NULL,
  cache_context = list(),
  force_rerun_slingshot = FALSE,
  output_plot_dir = plot_dir,
  output_results_dir = results_dir
) {
  label <- safe_label(paste(compartment, condition, sep = "_"))
  plot_label <- paste(compartment, condition)
  message("Processing ", label)

  expected_cache_metadata <- c(
    list(
      cache_format_version = 2L,
      compartment = compartment,
      condition = condition,
      min_cells_for_slingshot = min_cells_for_slingshot,
      min_clusters_for_slingshot = min_clusters_for_slingshot,
      min_cells_per_cluster = min_cells_per_cluster,
      slingshot_package_version = as.character(utils::packageVersion("slingshot")),
      seurat_extend_package_version = as.character(utils::packageVersion("SeuratExtend"))
    ),
    cache_context
  )

  cached_inference <- NULL
  if (!force_rerun_slingshot && !is.null(cache_file) && file.exists(cache_file)) {
    cached_inference <- tryCatch(
      readRDS(cache_file),
      error = function(e) {
        warning("Could not read Slingshot cache for ", label, ": ", conditionMessage(e))
        NULL
      }
    )

    cache_is_valid <- !is.null(cached_inference) &&
      is.list(cached_inference) &&
      identical(cached_inference$cache_metadata, expected_cache_metadata) &&
      all(c("obj", "cluster_counts", "eligibility", "sling", "curve_weights", "trade_seq_pseudotime") %in% names(cached_inference))

    if (!cache_is_valid) {
      message("Ignoring stale or incomplete Slingshot cache for ", label)
      cached_inference <- NULL
    }
  }

  used_cached_slingshot <- !is.null(cached_inference)
  if (used_cached_slingshot) {
    message("Using cached Slingshot result for ", label)
    obj <- cached_inference$obj
    cluster_counts <- cached_inference$cluster_counts
    eligibility <- cached_inference$eligibility
    sling <- cached_inference$sling
    curve_weights <- cached_inference$curve_weights
    trade_seq_pseudotime <- cached_inference$trade_seq_pseudotime
    lineage_assignments <- cached_inference$lineage_assignments

    if (is.null(lineage_assignments)) {
      lineage_assignments <- get_slingshot_lineage_assignments(curve_weights)
      cached_inference$lineage_assignments <- lineage_assignments
      saveRDS(cached_inference, cache_file)
      message("Added lineage assignments to Slingshot cache: ", cache_file)
    }

    cached_root_cluster <- cached_inference$root_cluster
    if (is.null(cached_root_cluster) && "root_cluster" %in% colnames(eligibility)) {
      cached_root_cluster <- eligibility$root_cluster
    }

    if ((is.null(cached_root_cluster) || all(is.na(cached_root_cluster))) &&
        "is_root_cluster" %in% colnames(cluster_counts)) {
      cached_root_cluster <- cluster_counts %>%
        dplyr::filter(!is.na(is_root_cluster) & is_root_cluster) %>%
        dplyr::pull(seurat_clusters)
    }

    if (length(cached_root_cluster) == 0 || all(is.na(cached_root_cluster))) {
      stop("Cached Slingshot result is missing its root cluster: ", cache_file)
    }

    root_cluster <- as.character(cached_root_cluster[[1]])
  } else {
    slingshot_subset <- prepare_slingshot_subset(
      obj,
      compartment = compartment,
      condition = condition
    )

    cluster_counts <- slingshot_subset$cluster_counts
    eligibility <- slingshot_subset$eligibility

    if (!slingshot_subset$slingshot_run) {
      message("Skipping Slingshot for ", label, ": not enough cells or clusters.")
      cluster_counts <- cluster_counts %>%
        dplyr::mutate(
          is_root_cluster = NA,
          root_cluster_projectils_label = NA_character_
        )

      return(list(
        obj = obj,
        cluster_counts = cluster_counts,
        eligibility = eligibility,
        sling_summary = NULL
      ))
    }

    obj <- slingshot_subset$obj
    root_info <- get_root_cluster_info(obj)
    root_cluster <- root_info$root_cluster[[1]]

    if (is.na(root_cluster)) {
      message("Skipping Slingshot for ", label, ": no naive-like root cluster was identified.")
      eligibility$skip_reason <- "No naive-like root cluster identified"
      cluster_counts <- cluster_counts %>%
        dplyr::mutate(
          is_root_cluster = NA,
          root_cluster_projectils_label = NA_character_
        )

      return(list(
        obj = obj,
        cluster_counts = cluster_counts,
        eligibility = eligibility,
        sling_summary = NULL
      ))
    }

    root_projectils_label <- root_info$root_cluster_projectils_label[[1]]
    eligibility$root_cluster <- root_cluster
    eligibility$root_cluster_selection <- root_info$root_cluster_selection[[1]]
    eligibility$root_cluster_naive_memory_score <- root_info$root_cluster_naive_memory_score[[1]]
    cluster_counts <- add_root_cluster_qc(cluster_counts, root_cluster, root_projectils_label)

    stored_reduction <- obj@misc$cd4_cd8_cluster_reduction

    if (is.null(stored_reduction) ||
        length(stored_reduction) == 0 ||
        !stored_reduction %in% names(obj@reductions)) {
      stored_reduction <- "pca"
    }

    stored_dims <- obj@misc$cd4_cd8_cluster_dims

    if (is.null(stored_dims) || length(stored_dims) == 0) {
      stored_dims <- ncol(Embeddings(obj, stored_reduction))
    }

    sling_dims <- min(stored_dims, ncol(Embeddings(obj, stored_reduction)))

    obj[["SLING"]] <- SeuratObject::CreateDimReducObject(
      embeddings = Embeddings(obj, stored_reduction)[, seq_len(sling_dims), drop = FALSE],
      key = "SLING_",
      assay = DefaultAssay(obj)
    )

    obj <- RunSlingshot(
      obj,
      group.by = "seurat_clusters",
      reducedDim = "SLING",
      start.clus = root_cluster
    )

    sling <- obj@misc$slingshot$SLING$SlingPseudotime
    obj@meta.data[, colnames(sling)] <- as.data.frame(sling)
    trade_seq_inputs <- get_slingshot_trade_seq_inputs(obj, sling)
    trade_seq_pseudotime <- trade_seq_inputs$pseudotime
    curve_weights <- trade_seq_inputs$cell_weights
    lineage_assignments <- get_slingshot_lineage_assignments(curve_weights)

    if (!is.null(cache_file)) {
      dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(
        list(
          cache_metadata = expected_cache_metadata,
          obj = obj,
          cluster_counts = cluster_counts,
          eligibility = eligibility,
          sling = sling,
          curve_weights = curve_weights,
          trade_seq_pseudotime = trade_seq_pseudotime,
          lineage_assignments = lineage_assignments,
          root_cluster = root_cluster
        ),
        cache_file
      )
      message("Saved Slingshot result cache: ", cache_file)
    }
  }

  if (!"functional_state_label" %in% colnames(obj@meta.data)) {
    obj$functional_state_label <- NA_character_
  }

  if (!"cluster_projectils_label" %in% colnames(obj@meta.data)) {
    obj$cluster_projectils_label <- NA_character_
  }

  if (!"cluster_projectils_label_prop" %in% colnames(obj@meta.data)) {
    obj$cluster_projectils_label_prop <- NA_real_
  }

  if (!"previous_tcell_subtype" %in% colnames(obj@meta.data)) {
    obj$previous_tcell_subtype <- NA_character_
  }

  cluster_labels <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::select(
      seurat_clusters,
      dplyr::any_of(c(
        "top_previous_tcell_subtype",
        "top_previous_tcell_subtype_prop",
        "cluster_projectils_label",
        "cluster_projectils_label_prop"
      ))
    ) %>%
    dplyr::distinct()

  sling_df <- as.data.frame(sling) %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(
      cols = starts_with("slingPseudotime"),
      names_to = "lineage",
      values_to = "pseudotime"
    ) %>%
    dplyr::filter(!is.na(pseudotime)) %>%
    dplyr::left_join(
      obj@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::mutate(naive_memory_score = state_naive_memory_score) %>%
        dplyr::select(
          cell,
          seurat_clusters,
          dplyr::any_of("functional_state_label"),
          dplyr::any_of(c("cluster_projectils_label", "cluster_projectils_label_prop")),
          dplyr::any_of(c("singler_immgen_label", "singler_mousernaseq_label")),
          naive_memory_score
        ),
      by = "cell"
  ) %>%
    dplyr::group_by(
      lineage,
      seurat_clusters,
      dplyr::across(dplyr::any_of("functional_state_label")),
      dplyr::across(dplyr::any_of(c("cluster_projectils_label", "cluster_projectils_label_prop"))),
      dplyr::across(dplyr::any_of(c("singler_immgen_label", "singler_mousernaseq_label")))
    ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_pseudotime = median(pseudotime),
      mean_pseudotime = mean(pseudotime),
      median_naive_memory_score = median(naive_memory_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      functional_state_label = if ("functional_state_label" %in% names(.)) {
        functional_state_label
      } else {
        NA_character_
      },
      projectils_label = dplyr::coalesce(cluster_projectils_label, NA_character_),
      projectils_label_prop = cluster_projectils_label_prop
    ) %>%
    dplyr::left_join(cluster_labels, by = "seurat_clusters") %>%
    dplyr::mutate(
      root_cluster = .env$root_cluster,
      seurat_cluster_top_previous_tcell_subtype = paste0(
        as.character(seurat_clusters),
        ": ",
        dplyr::coalesce(top_previous_tcell_subtype, "Unannotated")
      ),
      seurat_cluster_projectils_label = paste0(
        as.character(seurat_clusters),
        ": ",
        dplyr::coalesce(projectils_label, top_previous_tcell_subtype, "Unannotated")
      )
    )

  lineage_cols <- colnames(sling)
  lineage_display_names <- stats::setNames(
    paste0("lineage ", seq_along(lineage_cols)),
    lineage_cols
  )

  pseudotime_plots <- lapply(
    lineage_cols,
    function(lineage_col) {
      DimPlot2(
      obj,
      reduction = "umap",
      features = lineage_col,
      cols = "C",
      theme = NoAxes()
    ) +
      ggplot2::ggtitle(paste0(plot_label, ": ", lineage_display_names[[lineage_col]])) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8)
      )
    }
  )

  grDevices::cairo_pdf(
    filename = file.path(output_plot_dir, paste0("tcell_", label, "_slingshot_pseudotime_umap.pdf")),
    width = 6,
    height = 5,
    onefile = TRUE,
    family = "sans"
  )

  for (p in pseudotime_plots) {
    print(p)
  }

  grDevices::dev.off()

  summary_plot_height <- max(
    3.5,
    0.28 * max(table(sling_df$lineage))
  )

  summary_plots <- lapply(
    lineage_cols,
    function(lineage_name) {
      sling_df %>%
        dplyr::filter(lineage == lineage_name) %>%
        dplyr::mutate(
          cluster_label = seurat_cluster_projectils_label
        ) %>%
        ggplot2::ggplot(
          ggplot2::aes(
            x = median_pseudotime,
            y = reorder(cluster_label, median_pseudotime)
          )
        ) +
        ggplot2::geom_point(ggplot2::aes(size = n), color = "#2166AC") +
        theme_Publication() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_text(size = 10),
          axis.text.x = ggplot2::element_text(size = 8),
          axis.text.y = ggplot2::element_text(size = 8),
          legend.title = ggplot2::element_text(size = 9),
          legend.text = ggplot2::element_text(size = 8)
        ) +
        ggplot2::labs(
          x = "Median Slingshot pseudotime",
          y = "Cluster",
          size = "Cells",
          title = paste0(plot_label, ": ", lineage_display_names[[lineage_name]])
        )
    }
  )

  grDevices::cairo_pdf(
    filename = file.path(output_plot_dir, paste0("tcell_", label, "_slingshot_pseudotime_all_lineages_summary.pdf")),
    width = 7,
    height = summary_plot_height,
    onefile = TRUE,
    family = "sans"
  )

  for (p in summary_plots) {
    print(p)
  }

  grDevices::dev.off()

  umap_df <- Embeddings(obj, "umap") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")

  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

  umap_df <- umap_df %>%
    dplyr::left_join(
      obj@meta.data %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::transmute(
          cell,
          seurat_clusters,
          condition,
          dplyr::across(dplyr::all_of(colnames(sling)))
        ),
      by = "cell"
    ) %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::left_join(cluster_labels, by = "seurat_clusters") %>%
    dplyr::mutate(
      tcell_subtype = paste0(
        seurat_clusters,
        ": ",
        dplyr::coalesce(cluster_projectils_label, top_previous_tcell_subtype, "Unannotated")
      )
    )

  subtype_levels <- sort(unique(as.character(umap_df$tcell_subtype)))
  if (is.null(subtype_colors)) {
    subtype_colors <- make_tcell_subtype_palette(subtype_levels)
  }

  curve_df <- purrr::map2_dfr(
    lineage_cols,
    paste0("Lineage ", seq_along(lineage_cols)),
    ~ make_umap_curve_from_pseudotime(
      umap_df = umap_df,
      lineage_col = .x,
      lineage_name = .y,
      approx_points = 150
    )
  )

  lineage_node_df <- make_lineage_node_df(umap_df, lineage_cols, min_cells = 10)
  lineage_plots <- lapply(
    seq_along(lineage_cols),
    plot_one_sling_lineage,
    pseudotime_cols = lineage_cols,
    umap_df = umap_df,
    curve_df = curve_df,
    lineage_node_df = lineage_node_df,
    subtype_colors = subtype_colors,
    show_subtype_labels = TRUE
  )

  grDevices::cairo_pdf(
    filename = file.path(output_plot_dir, paste0("tcell_", label, "_slingshot_lineage_curves_each_lineage_umap.pdf")),
    width = 6,
    height = 5.5,
    onefile = TRUE,
    family = "sans"
  )

  for (p in lineage_plots) {
    print(p)
  }

  grDevices::dev.off()

  gene_trends <- plot_slingshot_gene_trends(
    obj = obj,
    sling = sling,
    genes = gene_trend_genes,
    analysis_label = label,
    plot_dir = output_plot_dir
  )

  aucell_pseudotime <- NULL
  if (plot_aucell_pseudotime && !is.null(aucell_programs)) {
    aucell_pseudotime <- plot_aucell_pseudotime_tracks(
      obj = obj,
      sling = sling,
      aucell_programs = aucell_programs,
      analysis_label = label,
      plot_dir = output_plot_dir,
      results_dir = output_results_dir
    )
  }

  list(
    obj = obj,
    cluster_counts = cluster_counts,
    eligibility = eligibility,
    sling_summary = sling_df,
    gene_trends = gene_trends,
    aucell_pseudotime = aucell_pseudotime,
    cache_file = cache_file,
    lineage_assignments = lineage_assignments,
    used_cached_slingshot = used_cached_slingshot
  )
}
