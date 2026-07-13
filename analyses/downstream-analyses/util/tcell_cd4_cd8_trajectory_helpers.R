# Helper functions for CD4/CD8 T-cell trajectory analysis.
#
# These functions support `12-tcell-cd4-cd8-trajectory.Rmd` by keeping marker
# scoring, SingleR/ScGate annotation summaries, clustering, and Slingshot
# execution outside the main report.

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
  
  paste0(program_name, "_score")
}

# Add marker-program scores used for CD4/CD8 calling and cluster annotation.
add_marker_program_scores <- function(obj) {
  for (program_name in names(marker_programs)) {
    obj[[marker_score_column(program_name)]] <- score_marker_set(
      obj,
      marker_programs[[program_name]]
    )
  }
  
  obj
}

# Return the metadata column used for one functional-state score.
state_score_column <- function(program_name) {
  paste0("state_", program_name, "_score")
}

# Add functional-state scores used to annotate CD4-like and CD8-like clusters.
add_functional_state_scores <- function(obj) {
  for (program_name in names(functional_state_programs)) {
    obj[[state_score_column(program_name)]] <- score_marker_set(
      obj,
      functional_state_programs[[program_name]]
    )
  }
  
  obj
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
    treg = "Treg-like CD8 T cells",
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
  
  model_db <- tryCatch(
    scGate::get_scGateDB(),
    error = function(e) {
      warning("Could not load ScGate model database: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(model_db) || !"mouse" %in% names(model_db)) {
    return(NULL)
  }
  
  if (!"generic" %in% names(model_db$mouse)) {
    warning("The mouse ScGate database does not contain a generic model group.")
    return(NULL)
  }
  
  if (!model_name %in% names(model_db$mouse$generic)) {
    warning("The mouse generic ScGate database does not contain model: ", model_name)
    return(NULL)
  }
  
  model_db$mouse$generic[[model_name]]
}

# Create simple marker-based CD4 and CD8 models if database models are unavailable.
make_fallback_scgate_models <- function() {
  cd4_model <- scGate::gating_model(
    name = "CD4T",
    signature = c("Cd4", "Cd8a-", "Cd8b1-")
  )
  
  cd8_model <- scGate::gating_model(
    name = "CD8T",
    signature = c("Cd8a", "Cd8b1", "Cd4-")
  )
  
  list(CD4 = cd4_model, CD8 = cd8_model)
}

# Run ScGate CD4 and CD8 models as supporting CD4/CD8 annotation.
# The function first tries the generic mouse ScGate models and falls back to
# simple marker-based models if generic models are not available.
run_scgate_cd4_cd8 <- function(obj) {
  cd4_model <- get_generic_mouse_scgate_model("CD4")
  cd8_model <- get_generic_mouse_scgate_model("CD8")
  obj$scgate_cd4_model <- ifelse(is.null(cd4_model), "fallback_marker_model", "mouse_generic_CD4T")
  obj$scgate_cd8_model <- ifelse(is.null(cd8_model), "fallback_marker_model", "mouse_generic_CD8T")
  
  if (is.null(cd4_model) || is.null(cd8_model)) {
    fallback_models <- make_fallback_scgate_models()
    
    if (is.null(cd4_model)) {
      cd4_model <- fallback_models$CD4
    }
    
    if (is.null(cd8_model)) {
      cd8_model <- fallback_models$CD8
    }
  }
  
  obj$scgate_cd4 <- NA_character_
  obj$scgate_cd8 <- NA_character_
  
  obj_cd4 <- tryCatch(
    scGate::scGate(data = obj, model = cd4_model),
    error = function(e) {
      warning("scGate CD4 model failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(obj_cd4) && "is.pure" %in% colnames(obj_cd4@meta.data)) {
    obj$scgate_cd4 <- as.character(obj_cd4$is.pure)
  }
  
  obj_cd8 <- tryCatch(
    scGate::scGate(data = obj, model = cd8_model),
    error = function(e) {
      warning("scGate CD8 model failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(obj_cd8) && "is.pure" %in% colnames(obj_cd8@meta.data)) {
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

# Summarize marker, SingleR, and ScGate evidence for each cluster annotation.
summarize_cluster_annotation <- function(obj, compartment) {
  marker_score_cols <- purrr::map_chr(names(marker_programs), marker_score_column)
  marker_score_cols <- intersect(marker_score_cols, colnames(obj@meta.data))
  state_score_cols <- purrr::map_chr(names(functional_state_programs), state_score_column)
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
    dplyr::left_join(singler_summary, by = "seurat_clusters") %>%
    dplyr::left_join(scgate_summary, by = "seurat_clusters")
}

# Add annotation metadata to a clustered object and summarize clusters.
annotate_clustered_object <- function(obj, compartment) {
  obj <- add_marker_program_scores(obj)
  obj <- add_functional_state_scores(obj)
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
  
  if (dims < 2) {
    stop("At least 3 cells are required for PCA-based clustering.")
  }
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = seq_len(dims), verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = seq_len(dims), verbose = FALSE)
  obj$cd4_cd8_condition_label <- label
  obj@misc$cd4_cd8_cluster_dims <- dims
  obj
}

# Choose the Slingshot cluster with the highest naive/memory score as the root.
get_root_cluster <- function(obj) {
  cluster_scores <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
      median_naive_memory_score = median(naive_memory_score, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(median_naive_memory_score), dplyr::desc(n))
  
  cluster_scores$seurat_clusters[[1]]
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
    total_cells = ncol(obj),
    total_slingshot_clusters = dplyr::n_distinct(as.character(obj$seurat_clusters)),
    retained_cells = n_retained_cells,
    retained_slingshot_clusters = n_retained_groups,
    min_cells_for_slingshot = min_cells_for_slingshot,
    min_clusters_for_slingshot = min_clusters_for_slingshot,
    min_cells_per_cluster = min_cells_per_cluster,
    slingshot_run = slingshot_run,
    skip_reason = skip_reason,
    root_cluster = NA_character_,
    root_cluster_median_naive_memory_score = NA_real_
  )
  
  if (!eligibility$slingshot_run) {
    return(list(
      obj = obj,
      cluster_counts = cluster_counts,
      eligibility = eligibility
    ))
  }
  
  obj_filtered <- subset(obj, cells = retained_cells)
  obj_filtered$seurat_clusters <- droplevels(factor(obj_filtered$seurat_clusters))
  
  list(
    obj = obj_filtered,
    cluster_counts = cluster_counts,
    eligibility = eligibility
  )
}

# Run Slingshot on one condition subset using Seurat clusters.
run_condition_slingshot <- function(obj, compartment, condition) {
  label <- safe_label(paste(compartment, condition, sep = "_"))
  plot_label <- paste(compartment, condition)
  message("Processing ", label)
  
  slingshot_subset <- prepare_slingshot_subset(
    obj,
    compartment = compartment,
    condition = condition
  )
  
  cluster_counts <- slingshot_subset$cluster_counts
  eligibility <- slingshot_subset$eligibility
  
  if (!eligibility$slingshot_run) {
    message("Skipping Slingshot for ", label, ": not enough cells or clusters.")
    return(list(
      obj = obj,
      cluster_counts = cluster_counts,
      eligibility = eligibility,
      sling_summary = NULL
    ))
  }
  
  obj <- slingshot_subset$obj
  root_cluster <- get_root_cluster(obj)
  root_cluster_median_naive_memory_score <- median(
    obj$naive_memory_score[as.character(obj$seurat_clusters) == root_cluster],
    na.rm = TRUE
  )
  eligibility$root_cluster <- root_cluster
  eligibility$root_cluster_median_naive_memory_score <- root_cluster_median_naive_memory_score
  
  stored_pca_dims <- obj@misc$cd4_cd8_cluster_dims
  
  if (is.null(stored_pca_dims) || length(stored_pca_dims) == 0) {
    stored_pca_dims <- ncol(Embeddings(obj, "pca"))
  }
  
  pca_dims <- min(stored_pca_dims, ncol(Embeddings(obj, "pca")))
  
  obj[["PCA20"]] <- SeuratObject::CreateDimReducObject(
    embeddings = Embeddings(obj, "pca")[, seq_len(pca_dims), drop = FALSE],
    key = "PCA20_",
    assay = DefaultAssay(obj)
  )
  
  obj <- RunSlingshot(
    obj,
    group.by = "seurat_clusters",
    reducedDim = "PCA20",
    start.clus = root_cluster
  )
  
  sling <- obj@misc$slingshot$PCA20$SlingPseudotime
  obj@meta.data[, colnames(sling)] <- as.data.frame(sling)
  
  if (!"functional_state_label" %in% colnames(obj@meta.data)) {
    obj$functional_state_label <- NA_character_
  }
  
  if (!"previous_tcell_subtype" %in% colnames(obj@meta.data)) {
    obj$previous_tcell_subtype <- NA_character_
  }
  
  cluster_previous_labels <- obj@meta.data %>%
    dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    dplyr::select(
      seurat_clusters,
      dplyr::any_of(c("top_previous_tcell_subtype", "top_previous_tcell_subtype_prop"))
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
        dplyr::select(
          cell,
          seurat_clusters,
          dplyr::any_of("functional_state_label"),
          dplyr::any_of(c("singler_immgen_label", "singler_mousernaseq_label")),
          naive_memory_score
        ),
      by = "cell"
  ) %>%
    dplyr::group_by(
      lineage,
      seurat_clusters,
      dplyr::across(dplyr::any_of("functional_state_label")),
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
      }
    ) %>%
    dplyr::left_join(cluster_previous_labels, by = "seurat_clusters") %>%
    dplyr::mutate(
      root_cluster = root_cluster,
      seurat_cluster_top_previous_tcell_subtype = paste0(
        as.character(seurat_clusters),
        ":",
        dplyr::coalesce(top_previous_tcell_subtype, "Unannotated")
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
    filename = file.path(plot_dir, paste0("tcell_", label, "_slingshot_pseudotime_umap.pdf")),
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
          cluster_label = seurat_cluster_top_previous_tcell_subtype
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
    filename = file.path(plot_dir, paste0("tcell_", label, "_slingshot_pseudotime_all_lineages_summary.pdf")),
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
    dplyr::left_join(cluster_previous_labels, by = "seurat_clusters") %>%
    dplyr::mutate(
      tcell_subtype = paste0(
        seurat_clusters,
        ":",
        dplyr::coalesce(top_previous_tcell_subtype, "Unannotated")
      )
    )
  
  subtype_levels <- sort(unique(as.character(umap_df$tcell_subtype)))
  subtype_colors <- stats::setNames(
    scales::hue_pal()(length(subtype_levels)),
    subtype_levels
  )
  
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
    filename = file.path(plot_dir, paste0("tcell_", label, "_slingshot_lineage_curves_each_lineage_umap.pdf")),
    width = 6,
    height = 5.5,
    onefile = TRUE,
    family = "sans"
  )
  
  for (p in lineage_plots) {
    print(p)
  }
  
  grDevices::dev.off()
  
  list(
    obj = obj,
    cluster_counts = cluster_counts,
    eligibility = eligibility,
    sling_summary = sling_df
  )
}
