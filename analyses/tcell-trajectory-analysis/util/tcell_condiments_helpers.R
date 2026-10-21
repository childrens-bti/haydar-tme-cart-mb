# Helpers for condition comparisons on shared CD4-like and CD8-like trajectories.
# Author: Bicna Song
# Date: 2026-10-13

# Configure the workers used for topology permutation fits.
make_condiments_bpparam <- function(workers, seed) {
  if (workers == 1L) {
    return(BiocParallel::SerialParam(RNGseed = seed))
  }
  if (.Platform$OS.type == "windows") {
    return(BiocParallel::SnowParam(workers = workers, RNGseed = seed))
  }
  BiocParallel::MulticoreParam(workers = workers, RNGseed = seed)
}

# Load the frozen all-condition trajectory and align every cell-level component.
prepare_condiments_input <- function(cache_file, compartment) {
  if (!file.exists(cache_file)) {
    stop("Shared Slingshot cache was not found: ", cache_file)
  }
  cache <- readRDS(cache_file)
  needed <- c("obj", "trade_seq_pseudotime", "curve_weights")
  if (!all(needed %in% names(cache))) {
    stop(compartment, " cache is missing: ",
         paste(setdiff(needed, names(cache)), collapse = ", "))
  }

  obj <- cache$obj
  pseudotime <- as.matrix(cache$trade_seq_pseudotime)
  weights <- as.matrix(cache$curve_weights)
  cells <- rownames(pseudotime)
  if (is.null(cells) || anyDuplicated(cells) ||
      !identical(cells, rownames(weights)) ||
      !identical(colnames(pseudotime), colnames(weights)) ||
      !identical(dim(pseudotime), dim(weights)) ||
      !identical(cells, colnames(obj)) ||
      !identical(cells, rownames(obj@meta.data))) {
    stop(compartment, " cache has inconsistent cell or lineage order")
  }
  if (ncol(pseudotime) < 2L || any(!is.finite(pseudotime)) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      any(rowSums(weights) <= 0)) {
    stop(compartment, " cache has invalid pseudotime or curve weights")
  }
  # Shared cells can have full weight on several lineages; use probabilities
  # rather than assigning each cell to just one of them.
  weights <- weights / rowSums(weights)
  if (!all(c("condition", "sample") %in%
           colnames(obj@meta.data))) {
    stop(compartment, " cache lacks condition or sample metadata")
  }
  condition <- as.character(obj@meta.data$condition)
  sample <- as.character(obj@meta.data$sample)
  if (anyNA(condition) || anyNA(sample) || any(!nzchar(condition)) ||
      any(!nzchar(sample))) {
    stop(compartment, " cache contains missing condition or sample IDs")
  }
  if (dplyr::n_distinct(condition) < 2L ||
      dplyr::n_distinct(sample) < 2L) {
    stop(compartment, " requires at least two conditions and samples")
  }
  sample_conditions <- tapply(condition, sample, function(x) length(unique(x)))
  if (any(sample_conditions != 1L)) {
    stop(compartment, " has a sample assigned to multiple conditions")
  }
  # Inputs for topologyTest: the original reduced space, clusters and root.
  if (!"SLING" %in% names(obj@reductions)) {
    stop(compartment, " cache lacks the SLING reduction used by Slingshot")
  }
  sling_embeddings <- SeuratObject::Embeddings(obj, "SLING")[cells, , drop = FALSE]
  root_cluster <- cache$root_cluster
  if (is.null(root_cluster) && !is.null(cache$eligibility$root_cluster)) {
    root_cluster <- cache$eligibility$root_cluster[[1]]
  }
  if (is.null(root_cluster) || is.na(root_cluster)) {
    stop(compartment, " cache lacks the Slingshot root cluster")
  }

  list(
    compartment = compartment, cache_file = cache_file,
    cache_md5 = unname(tools::md5sum(cache_file)),
    obj = obj, cells = cells,
    condition = condition, sample = sample,
    pseudotime = pseudotime, weights = weights,
    lineage_names = colnames(pseudotime),
    clusters = as.character(obj@meta.data$seurat_clusters),
    root_cluster = as.character(root_cluster),
    sling_embeddings = sling_embeddings
  )
}

# Short mouse labels (D1, D2) for plots; fall back to numbering within condition.
make_mouse_labels <- function(sample, condition) {
  mouse <- stringr::str_extract(sample, "D[0-9]+$")
  missing <- is.na(mouse)
  if (any(missing)) {
    ranks <- stats::ave(seq_along(sample), condition, FUN = function(i) {
      match(sample[i], sort(unique(sample[i])))
    })
    mouse[missing] <- paste0("M", ranks[missing])
  }
  mouse
}

# Summarize each mouse's position on each lineage using Slingshot weights.
summarize_sample_pseudotime <- function(input, annotations) {
  labels <- annotations[annotations$compartment == input$compartment, ]
  purrr::map_dfr(seq_along(input$lineage_names), function(j) {
    lineage <- input$lineage_names[[j]]
    purrr::map_dfr(unique(input$sample), function(s) {
      cells <- which(input$sample == s)
      w <- input$weights[cells, j]
      tibble::tibble(
        compartment = input$compartment,
        condition = unique(input$condition[cells]), sample = s,
        lineage = lineage,
        lineage_program = labels$assigned_program[match(lineage, labels$lineage)],
        n_cells = length(cells),
        weight_sum = sum(w),
        weighted_mean_pseudotime = if (sum(w) > 0)
          stats::weighted.mean(input$pseudotime[cells, j], w)
        else NA_real_
      )
    })
  })
}

# Summarize each mouse's normalized lineage weights without hard assignment.
summarize_sample_lineage_probabilities <- function(input, annotations) {
  labels <- annotations[annotations$compartment == input$compartment, ]
  purrr::map_dfr(unique(input$sample), function(s) {
    cells <- which(input$sample == s)
    purrr::map_dfr(seq_along(input$lineage_names), function(j) {
      lineage <- input$lineage_names[[j]]
      w <- input$weights[cells, j]
      tibble::tibble(
        compartment = input$compartment,
        condition = unique(input$condition[cells]), sample = s,
        lineage = lineage,
        lineage_program = labels$assigned_program[match(lineage, labels$lineage)],
        n_cells = length(cells), weight_sum = sum(w),
        mean_lineage_probability = mean(w)
      )
    })
  })
}

# Use a cache only when its complete input and parameter context matches.
load_condiments_result_cache <- function(file, context, force = FALSE) {
  if (!force && file.exists(file)) {
    cached <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.list(cached) && identical(cached$context, context) &&
        is.list(cached$results)) {
      message("Reusing condiments cache: ", file)
      return(cached)
    }
    message("Ignoring incompatible condiments cache: ", file)
  }
  list(context = context, results = list())
}

# Save each completed comparison immediately; failed comparisons retry next time.
get_or_run_condiments_comparison <- function(cache, key, run, file) {
  if (key %in% names(cache$results)) {
    message("Reusing condiments comparison: ", key)
    return(list(cache = cache, value = cache$results[[key]]))
  }
  message("Running condiments comparison: ", key)
  value <- tryCatch(run(), error = function(e) e)
  if (!inherits(value, "error")) {
    cache$results[[key]] <- value
    saveRDS(cache, file)
  } else {
    warning("Condiments comparison failed for ", key, ": ",
            conditionMessage(value))
  }
  list(cache = cache, value = value)
}

# Restrict each test to its prespecified two conditions in frozen cell order.
subset_condiments_comparison <- function(input, condition_a, condition_b) {
  keep <- which(input$condition %in% c(condition_a, condition_b))
  condition <- factor(input$condition[keep], levels = c(condition_a, condition_b))
  if (length(keep) < 2L || any(table(condition) == 0L)) {
    stop(input$compartment, " lacks cells for ", condition_a, " versus ",
         condition_b)
  }
  list(
    compartment = input$compartment,
    condition = condition,
    pseudotime = input$pseudotime[keep, , drop = FALSE],
    weights = input$weights[keep, , drop = FALSE],
    n_cells = length(keep), n_samples = dplyr::n_distinct(input$sample[keep])
  )
}

# Compare the frozen shared pseudotime distributions within each lineage.
run_condiments_progression <- function(pair_data) {
  condiments::progressionTest(
    pseudotime = pair_data$pseudotime,
    cellWeights = pair_data$weights,
    conditions = pair_data$condition,
    global = TRUE, lineages = TRUE, method = "KS"
  )
}

# Compare soft lineage allocations using the package's MMD test.
run_condiments_fate <- function(pair_data) {
  condiments::fateSelectionTest(
    cellWeights = pair_data$weights,
    conditions = pair_data$condition,
    global = TRUE, pairwise = TRUE, method = "mmd"
  )
}

# Give progression and fate tests the same output columns and explicit failures,
# and map package lineage numbers back to the shared lineage IDs and labels.
normalize_condiments_test <- function(value, test, compartment, comparison,
                                     n_samples, n_cells, lineage_names,
                                     annotations) {
  test_type <- match.arg(test, c("progression", "fate"))
  empty <- tibble::tibble(
    compartment = compartment, test = test_type, lineage = NA_character_,
    lineage_program = NA_character_, comparison = comparison,
    method = NA_character_, statistic = NA_real_, p_value = NA_real_,
    n_samples = n_samples, n_cells = n_cells,
    status = "failed", status_note = if (inherits(value, "error"))
      conditionMessage(value) else "Test returned no result"
  )
  if (inherits(value, "error") || is.null(value) || length(value) == 0L) {
    return(empty)
  }
  result <- as.data.frame(value)
  if (nrow(result) == 0L ||
      !all(c("statistic", "p.value") %in% names(result))) {
    return(empty)
  }
  row_ids <- if (test_type == "progression") result$lineage else result$pair
  if (length(row_ids) != nrow(result)) {
    return(empty)
  }
  row_ids <- as.character(row_ids)
  labels <- annotations[annotations$compartment == compartment, , drop = FALSE]
  lineage <- rep(NA_character_, length(row_ids))
  program <- rep(NA_character_, length(row_ids))
  if (test_type == "progression") {
    idx <- suppressWarnings(as.integer(row_ids))
    valid <- !is.na(idx) & idx >= 1L & idx <= length(lineage_names)
    lineage[valid] <- lineage_names[idx[valid]]
    program[valid] <- labels$assigned_program[
      match(lineage[valid], labels$lineage)
    ]
    lineage[row_ids == "All"] <- "All"
    program[row_ids == "All"] <- "All lineages"
  } else {
    lineage[row_ids == "All"] <- "All"
    program[row_ids == "All"] <- "All lineage pairs"
    for (i in which(row_ids != "All")) {
      pair <- strsplit(row_ids[[i]], "vs", fixed = TRUE)[[1]]
      idx <- suppressWarnings(as.integer(pair))
      if (length(idx) == 2L && all(!is.na(idx)) &&
          all(idx >= 1L & idx <= length(lineage_names))) {
        ids <- lineage_names[idx]
        lineage[[i]] <- paste(ids, collapse = " vs ")
        program[[i]] <- paste(labels$assigned_program[
          match(ids, labels$lineage)
        ], collapse = " vs ")
      }
    }
  }
  # Resolve the method before tibble creates a multirow column named 'test'.
  test_method <- if ("method" %in% names(result)) {
    as.character(result$method)
  } else if (test_type == "progression") {
    "KS"
  } else {
    "mmd"
  }
  tibble::tibble(
    compartment = compartment, test = test_type,
    lineage = lineage, lineage_program = program,
    comparison = comparison,
    method = test_method,
    statistic = as.numeric(result$statistic),
    p_value = as.numeric(result$p.value),
    n_samples = n_samples, n_cells = n_cells,
    status = ifelse(is.finite(as.numeric(result$p.value)), "ok", "untestable"),
    status_note = ifelse(is.finite(as.numeric(result$p.value)), NA_character_,
                         "Package returned a non-finite p-value")
  )
}


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

# Build a process-local compatibility fix for condiments 1.14 / Slingshot 2.14.
# Missing-cluster cleanup can leave Lineage1 and Lineage3 as the two remaining
# paths. Slingshot's shrinkage code must locate these by name, not assume that
# the suffix of Lineage3 means the third curve. Original IDs stay unchanged.
# Copies of the affected functions are used; no package namespace is modified.
# Returns the environment holding the patched copies.
make_condiments_topology_env <- function() {
  if (as.character(utils::packageVersion("condiments")) != "1.14.0" ||
      as.character(utils::packageVersion("slingshot")) != "2.14.0") {
    stop("The topology curve-index compatibility fix requires condiments 1.14.0 and slingshot 2.14.0")
  }

  # Replace only the known expressions; fail if the installed code has changed.
  replacements <- 0L
  replace_expression <- function(x, old, new) {
    if (identical(x, old)) {
      replacements <<- replacements + 1L
      return(new)
    }
    if (is.function(x)) {
      body(x) <- replace_expression(body(x), old, new)
      return(x)
    }
    if (is.call(x)) {
      x <- as.call(lapply(x, replace_expression, old = old, new = new))
    }
    x
  }
  method <- methods::selectMethod("getCurves", "PseudotimeOrdering",
                                  fdef = slingshot::getCurves)
  fixed_body <- replace_expression(
    body(method), quote(as.numeric(gsub("Lineage", "", n))),
    quote(match(n, names(lineages)))
  )
  if (replacements != 3L) {
    stop("Unexpected Slingshot curve-index implementation: ", replacements,
         " matching expressions (expected 3)")
  }
  get_curves <- eval(call("function", formals(method), fixed_body),
                     envir = environment(method))

  # Retain the package's permutation generation, missing-lineage extrapolation,
  # and statistical tests. Only condition-specific getCurves calls are rerouted.
  local_env <- new.env(parent = asNamespace("condiments"))
  for (name in c(".topologyTest", ".generate_permutations_curves",
                 ".condition_sling", ".sling_cond")) {
    fun <- get(name, envir = asNamespace("condiments"))
    environment(fun) <- local_env
    assign(name, fun, envir = local_env)
  }
  replacements <- 0L
  body(local_env$.sling_cond) <- replace_expression(
    body(local_env$.sling_cond), quote(slingshot::getCurves), quote(get_curves)
  )
  if (replacements != 1L) stop("Unexpected condiments condition-fit implementation")
  local_env$get_curves <- get_curves
  condition_fit <- function(sds, conditions, approx_points = 100,
                            adjust_skeleton = FALSE, verbose = TRUE, ...) {
    if (adjust_skeleton) stop("Topology fitting must retain the shared skeleton")
    if (!methods::is(sds, "PseudotimeOrdering")) {
      sds <- slingshot::as.PseudotimeOrdering(sds)
    }
    .sling_cond(sds, conditions, approx_points = approx_points,
                verbose = verbose, ...)
  }
  environment(condition_fit) <- local_env
  local_env$slingshot_conditions <- condition_fit
  local_env
}

# The patched topologyTest() used for every topology run.
make_condiments_topology_test <- function() {
  make_condiments_topology_env()$.topologyTest
}

# Refit Slingshot on all tested conditions together, using the same reduced
# space, clusters and root cluster as the shared fit. The topology test and
# the topology-curve figure both start from this fit.
refit_condiments_topology <- function(input, conditions) {
  index <- which(input$condition %in% conditions)
  if (dplyr::n_distinct(input$condition[index]) != length(conditions)) {
    stop(input$compartment, " lacks cells for some topology conditions")
  }
  clusters <- input$clusters[index]
  if (!input$root_cluster %in% clusters) {
    stop("Root cluster ", input$root_cluster, " has no cells in the tested conditions")
  }
  fit <- slingshot::slingshot(
    input$sling_embeddings[index, , drop = FALSE],
    clusterLabels = clusters,
    start.clus = input$root_cluster
  )
  list(index = index, fit = fit)
}

# Test whether fitting each condition separately on the all-condition refit
# changes the trajectories more than permuted labels do. With more than two
# conditions, condiments uses the classifier test.
run_condiments_topology <- function(input, conditions, method = "Classifier",
                                    rep = 100L, threshs = 0.01,
                                    bpparam = BiocParallel::SerialParam()) {
  refit_data <- refit_condiments_topology(input, conditions)
  refit <- refit_data$fit
  topology_test <- make_condiments_topology_test()
  result <- topology_test(
    sds = refit,
    conditions = input$condition[refit_data$index],
    rep = rep,
    threshs = threshs,
    methods = method,
    parallel = BiocParallel::bpworkers(bpparam) > 1L,
    BPPARAM = bpparam
  )
  attr(result, "n_refit_lineages") <- ncol(slingshot::slingPseudotime(refit))
  result
}

# Return the global topology result in the shared result format (one row per
# compartment).
normalize_condiments_topology <- function(value, compartment, n_samples, n_cells) {
  base <- tibble::tibble(
    compartment = compartment, test = "topology", lineage = "All",
    lineage_program = "Trajectory structure", comparison = "all_conditions",
    n_samples = n_samples, n_cells = n_cells
  )
  if (inherits(value, "error") || is.null(value) || NROW(value) == 0L) {
    return(dplyr::mutate(
      base, method = NA_character_, statistic = NA_real_, p_value = NA_real_,
      status = "failed",
      status_note = if (inherits(value, "error")) conditionMessage(value)
      else "Test returned no result"
    ))
  }
  result <- as.data.frame(value)[1, , drop = FALSE]
  p_value <- as.numeric(result$p.value)
  n_lineages <- attr(value, "n_refit_lineages")
  dplyr::mutate(
    base,
    method = as.character(result$method),
    statistic = as.numeric(result$statistic),
    p_value = p_value,
    status = ifelse(is.finite(p_value), "ok", "untestable"),
    status_note = if (is.null(n_lineages)) NA_character_ else
      paste0("All-condition refit has ", n_lineages, " lineages")
  )
}

# Fit each condition's curves on the all-condition skeleton, exactly as the
# topology test does for the real labels (patched curve fitting included).
# Refit lineages are matched to the shared lineages by their end clusters.
# Returns the fitted principal-curve coordinates (first four SLING dimensions,
# or fewer if the space has fewer) for each condition and for the
# all-condition reference fit.
fit_condiments_condition_curves <- function(input, conditions) {
  refit_data <- refit_condiments_topology(input, conditions)
  sds <- refit_data$fit
  if (!methods::is(sds, "PseudotimeOrdering")) {
    sds <- slingshot::as.PseudotimeOrdering(sds)
  }
  condition_fits <- make_condiments_topology_env()$.sling_cond(
    sds, input$condition[refit_data$index],
    approx_points = slingshot::slingParams(sds)$approx_points, verbose = FALSE
  )

  end_cluster <- function(lineages) {
    vapply(lineages, function(path) as.character(utils::tail(path, 1)), character(1))
  }
  shared <- input$obj@misc$slingshot$SLING$SlingshotDataSet
  shared_ends <- if (is.null(shared)) character() else
    end_cluster(slingshot::slingLineages(shared))
  shared_ids <- if (length(shared_ends) == length(input$lineage_names)) {
    stats::setNames(input$lineage_names, shared_ends)
  } else {
    character()
  }
  refit_ends <- end_cluster(slingshot::slingLineages(sds))
  lineage_map <- stats::setNames(
    dplyr::coalesce(unname(shared_ids[refit_ends]),
                    paste0("unmatched ", names(refit_ends))),
    names(refit_ends)
  )

  curve_coordinates <- function(fit, group) {
    purrr::imap_dfr(slingshot::slingCurves(fit), function(curve, lin) {
      points <- curve$s[curve$ord, seq_len(min(4L, ncol(curve$s))), drop = FALSE]
      colnames(points) <- paste0("dim_", seq_len(ncol(points)))
      lineage <- if (lin %in% names(lineage_map)) lineage_map[[lin]] else
        paste("unmatched", lin)
      dplyr::bind_cols(
        tibble::tibble(group = group, lineage = lineage, step = seq_len(nrow(points))),
        tibble::as_tibble(points)
      )
    })
  }
  curves <- dplyr::bind_rows(
    curve_coordinates(sds, "All conditions"),
    purrr::imap_dfr(condition_fits, curve_coordinates)
  )

  list(curves = curves, lineage_map = lineage_map,
       refit_end_clusters = refit_ends)
}

# Draw the curves compared by the topology test, in two dimensions of the
# SLING space where Slingshot fitted them (the condiments vignette view).
# These are the fitted principal curves themselves, not smoothed cell paths.
# `dims` selects the pair of SLING dimensions; lineages that overlap in
# dimensions 1-2 can separate in dimensions 3-4. Page 1 compares constructs
# within each lineage; page 2 shows each construct's own curves over its cells.
plot_condiments_topology_curves <- function(input, curve_fit, annotations,
                                            conditions, dims = c(1L, 2L)) {
  dim_cols <- paste0("dim_", dims)
  if (!all(dim_cols %in% names(curve_fit$curves))) {
    stop("Curve coordinates are unavailable for SLING dimensions ",
         paste(dims, collapse = " and "))
  }
  # Plot the selected dimensions under fixed column names.
  cells_df <- tibble::tibble(dim_1 = input$sling_embeddings[, dims[[1]]],
                             dim_2 = input$sling_embeddings[, dims[[2]]],
                             condition = input$condition) |>
    dplyr::filter(.data$condition %in% conditions)
  curves <- curve_fit$curves |>
    dplyr::transmute(.data$group, .data$lineage, .data$step,
                     dim_1 = .data[[dim_cols[[1]]]], dim_2 = .data[[dim_cols[[2]]]])
  dims_label <- paste0("SLING dimensions ", dims[[1]], "-", dims[[2]])
  axis_x <- paste("SLING dimension", dims[[1]])
  axis_y <- paste("SLING dimension", dims[[2]])

  labels <- condiments_lineage_labels(annotations, input$compartment)
  extra <- setdiff(unique(curves$lineage), names(labels))
  labels[extra] <- paste0(extra, ": no matching shared lineage")
  short_id <- function(x) sub("slingPseudotime_", "L", x, fixed = TRUE)
  construct_levels <- display_conditions(conditions)
  reference <- dplyr::filter(curves, .data$group == "All conditions") |>
    dplyr::select(-"group")
  construct_curves <- dplyr::filter(curves, .data$group != "All conditions") |>
    dplyr::mutate(group = factor(.data$group, levels = construct_levels))
  curve_ends <- construct_curves |>
    dplyr::group_by(.data$group, .data$lineage) |>
    dplyr::slice_max(.data$step, n = 1) |>
    dplyr::ungroup()
  background <- dplyr::select(cells_df, "dim_1", "dim_2")
  arrow <- grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
  subtitle <- "Curves compared by the global topology test | Descriptive"
  theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 10),
                   plot.margin = ggplot2::margin(12, 18, 12, 12))

  by_lineage <- ggplot2::ggplot() +
    ggplot2::geom_point(data = background, ggplot2::aes(.data$dim_1, .data$dim_2),
                        colour = "grey88", size = 0.4) +
    ggplot2::geom_path(data = reference, ggplot2::aes(.data$dim_1, .data$dim_2),
                       colour = "grey20", linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_path(data = construct_curves,
                       ggplot2::aes(.data$dim_1, .data$dim_2, colour = .data$group),
                       linewidth = 0.9, arrow = arrow) +
    ggplot2::facet_wrap(~ lineage, nrow = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_colour_manual(values = condition_colors, labels = condition_plot_label,
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
    theme +
    ggplot2::labs(
      title = paste0(input$compartment, " Topology: each construct's curves on the same lineage (",
                     dims_label, ")"),
      subtitle = subtitle,
      caption = paste(
        "Colored lines: each construct's curves; arrows point to the lineage end. Dashed: all-condition curve. Grey points: all tested cells.",
        "A missing curve means that construct lost the lineage.",
        sep = "\n"),
      x = axis_x, y = axis_y)

  construct_cells <- dplyr::mutate(cells_df,
    group = factor(.data$condition, levels = construct_levels))
  by_construct <- ggplot2::ggplot() +
    ggplot2::geom_point(data = background, ggplot2::aes(.data$dim_1, .data$dim_2),
                        colour = "grey90", size = 0.4) +
    ggplot2::geom_point(data = construct_cells,
                        ggplot2::aes(.data$dim_1, .data$dim_2, colour = .data$group),
                        size = 0.5, alpha = 0.5, show.legend = FALSE) +
    ggplot2::geom_path(data = reference,
                       ggplot2::aes(.data$dim_1, .data$dim_2, group = .data$lineage),
                       colour = "grey45", linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_path(data = construct_curves,
                       ggplot2::aes(.data$dim_1, .data$dim_2, group = .data$lineage),
                       colour = "grey10", linewidth = 0.8, arrow = arrow) +
    ggplot2::geom_text(data = curve_ends,
                       ggplot2::aes(.data$dim_1, .data$dim_2,
                                    label = short_id(.data$lineage)),
                       size = 3.4, fontface = "bold", vjust = -0.6) +
    ggplot2::facet_wrap(~ group, nrow = 1, labeller = ggplot2::as_labeller(
      function(x) condition_label_string(x), default = ggplot2::label_parsed)) +
    ggplot2::scale_colour_manual(values = condition_colors) +
    theme +
    ggplot2::labs(
      title = paste0(input$compartment, " Topology: each construct's own curves (",
                     dims_label, ")"),
      subtitle = subtitle,
      caption = paste(
        "Colored points: that construct's cells. Solid lines: its curves. Dashed: all-condition curves.",
        paste(strwrap(paste(gsub("\n", " ", labels, fixed = TRUE), collapse = "; "),
                      width = 150), collapse = "\n"),
        sep = "\n"),
      x = axis_x, y = axis_y)

  list(by_lineage, by_construct)
}

# Order the conditions consistently across all figures.
display_conditions <- function(condition, exclude_conditions = character()) {
  present <- setdiff(unique(as.character(condition)), exclude_conditions)
  c(intersect(condition_display_order, present),
    setdiff(sort(present), condition_display_order))
}

# Keep technical lineage IDs explicit beside their reviewed biological labels.
condiments_lineage_labels <- function(annotations, compartment) {
  labels <- annotations[annotations$compartment == compartment, ]
  stats::setNames(vapply(seq_len(nrow(labels)), function(i) {
    id <- sub("slingPseudotime_", "L", labels$lineage[[i]], fixed = TRUE)
    paste(strwrap(paste0(id, ": ", labels$assigned_program[[i]]), width = 30),
          collapse = "\n")
  }, character(1)), labels$lineage)
}

# Display one pairwise test at a time; FDR is statistical evidence, not effect
# magnitude. Rows show global and lineage-specific results, separating
# CAR-vs-CAR from CAR-vs-control comparisons. The global topology test has one
# result per compartment and is reported as a table, not plotted.
plot_condiments_test_results <- function(results, compartment, test, annotations) {
  test_type <- match.arg(test, c("progression", "fate"))
  data <- results |>
    dplyr::filter(.data$compartment == .env$compartment,
                  .data$test == .env$test_type) |>
    dplyr::mutate(
      score = pmin(-log10(pmax(.data$p_adj, .Machine$double.xmin)), 6),
      label = vapply(seq_len(dplyr::n()), function(i) {
        p <- .data$p_adj[[i]]
        if (is.na(p)) return(.data$status[[i]])
        formatC(p, format = if (p < 0.01) "e" else "f",
                digits = if (p < 0.01) 1 else 3)
      }, character(1)),
      text_color = ifelse(!is.na(.data$score) & .data$score > 3, "white", "grey15")
    )
  conditions <- display_conditions(c(data$condition_a, data$condition_b))
  data <- data |>
    dplyr::arrange(.data$comparison_set != "primary_pairwise",
                   match(.data$condition_a, conditions),
                   match(.data$condition_b, conditions)) |>
    dplyr::mutate(
      row = paste0(condition_label_string(.data$condition_a), ' ~ "vs" ~ ',
                    condition_label_string(.data$condition_b)),
      column = dplyr::if_else(.data$lineage == "All", "All",
        gsub("slingPseudotime_", "L", .data$lineage, fixed = TRUE)),
      set = factor(ifelse(.data$comparison_set == "primary_pairwise",
                           "CAR vs CAR", "CAR vs control"),
                    levels = c("CAR vs CAR", "CAR vs control"))
    )
  data$row <- factor(data$row, levels = rev(unique(data$row)))
  data$column <- factor(data$column, levels = c("All", setdiff(unique(data$column), "All")))
  p <- ggplot2::ggplot(data, ggplot2::aes(.data$column, .data$row)) +
    ggplot2::facet_grid(set ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) parse(text = x)) +
    ggplot2::scale_x_discrete(position = "top")
  caption <- if (test_type == "progression")
    "All = global test; L1-L3 = shared lineages." else
    "All = global test; other columns = lineage pairs."
  caption <- paste(caption,
    "FDR is adjusted within each comparison across its four tests.", sep = " ")
  key <- paste(gsub("\n", " ", condiments_lineage_labels(annotations, compartment),
                     fixed = TRUE), collapse = "; ")
  caption <- paste(caption, paste(strwrap(key, width = 145), collapse = "\n"),
                    sep = "\n")
  p + ggplot2::geom_tile(ggplot2::aes(fill = .data$score),
                         color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, color = .data$text_color),
                       size = 3.6, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(low = "white", high = "#2166AC",
                                 limits = c(0, 6), na.value = "grey85",
                                 name = "-log10 FDR\n(capped at 6)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title.position = "plot", plot.caption.position = "plot",
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.margin = ggplot2::margin(12, 18, 12, 12)) +
    ggplot2::labs(
      title = paste(compartment, if (test_type == "fate")
        "Fate-selection test overview" else "Progression-test overview"),
      subtitle = "Cell-level test | Two mice per construct",
      caption = caption, x = NULL, y = NULL
    )
}

# Label cached FDRs without confusing an unavailable test with a negative one.
condiments_fdr_label <- function(p, status) {
  vapply(seq_along(p), function(i) {
    if (is.na(status[[i]]) || status[[i]] != "ok" || !is.finite(p[[i]])) {
      return("Test unavailable")
    }
    value <- formatC(p[[i]], format = if (p[[i]] < 0.01) "e" else "f",
                       digits = if (p[[i]] < 0.01) 1 else 3)
    paste0("FDR = ", value, "\n", if (p[[i]] < 0.05)
      "Difference detected" else "No significant difference")
  }, character(1))
}

# Focus on one prespecified condition pair. Weighted empirical CDFs show the
# same pseudotimes/soft weights used by the KS test. Means show direction only;
# the cached two-sided test can also detect crossing or shape differences.
plot_condiments_progression_distribution <- function(input, annotations,
                                                      comparison, results) {
  conditions <- c(comparison$condition_a[[1]], comparison$condition_b[[1]])
  data <- purrr::map_dfr(input$lineage_names, function(lineage) {
    purrr::map_dfr(conditions, function(condition) {
      keep <- input$condition == condition & input$weights[, lineage] > 0
      tibble::tibble(pseudotime = input$pseudotime[keep, lineage],
                     weight = input$weights[keep, lineage]) |>
        dplyr::group_by(.data$pseudotime) |>
        dplyr::summarise(weight = sum(.data$weight), .groups = "drop") |>
        dplyr::arrange(.data$pseudotime) |>
        dplyr::mutate(cumulative = cumsum(.data$weight) / sum(.data$weight),
                      mean_pseudotime = stats::weighted.mean(.data$pseudotime, .data$weight),
                      lineage = lineage, condition = condition)
    })
  })
  data$condition <- factor(data$condition, levels = conditions)
  labels <- condiments_lineage_labels(annotations, input$compartment)
  tests <- results |>
    dplyr::filter(.data$compartment == .env$input$compartment,
      .data$comparison == .env$comparison$comparison_id[[1]], .data$test == "progression")
  tests <- tests[match(names(labels), tests$lineage), ]
  labels[] <- paste(labels, condiments_fdr_label(tests$p_adj, tests$status), sep = "\n")
  means <- dplyr::distinct(data, lineage, condition, mean_pseudotime)
  ggplot2::ggplot(data, ggplot2::aes(.data$pseudotime, .data$cumulative,
                                      color = .data$condition)) +
    ggplot2::geom_hline(yintercept = 0.5, color = "grey85", linewidth = 0.4) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::geom_vline(data = means, ggplot2::aes(xintercept = .data$mean_pseudotime,
      color = .data$condition), linetype = "dashed", linewidth = 0.55,
      show.legend = FALSE) +
    ggplot2::geom_text(data = means, ggplot2::aes(x = Inf,
      y = ifelse(.data$condition == conditions[[1]], 0.12, 0.04),
      label = sprintf("mean: %.1f", .data$mean_pseudotime), color = .data$condition),
      inherit.aes = FALSE, hjust = 1.05, size = 3.4, show.legend = FALSE) +
    ggplot2::facet_wrap(~ lineage, nrow = 1, scales = "free_x",
                        labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_color_manual(values = condition_colors,
                                labels = condition_plot_label) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank(),
                   panel.spacing.x = grid::unit(1.2, "lines"),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 10),
                   plot.margin = ggplot2::margin(12, 18, 12, 12)) +
    ggplot2::labs(title = paste(input$compartment, "Progression: are cells later along the same lineage?"),
      subtitle = paste("Focused comparison:", comparison$comparison_id[[1]],
        "| Cell-level distributions and cached KS FDRs"),
      caption = paste(
        "Rightward = later shared pseudotime. Dashed lines and numbers: pooled cell-weighted means.",
        "The KS test detects any distribution difference, not specifically 'later'.",
        sep = "\n"),
      x = "Shared pseudotime", y = "Cumulative lineage-weighted fraction", color = NULL)
}

# Compare the per-cell relative weights actually used in pairwise fate tests.
# For L1 vs L3, x = w3 / (w1 + w3); cells with zero weight in both are excluded.
# No cell is assigned to a single lineage. Draw all eligible cells rather than
# reproducing the test's random condition-size balancing/subsampling.
plot_condiments_fate_distribution <- function(input, annotations, comparison, results) {
  conditions <- c(comparison$condition_a[[1]], comparison$condition_b[[1]])
  pairs <- utils::combn(input$lineage_names, 2, simplify = FALSE)
  labels <- annotations[annotations$compartment == input$compartment, ]
  tests <- results |>
    dplyr::filter(.data$compartment == .env$input$compartment,
      .data$comparison == .env$comparison$comparison_id[[1]], .data$test == "fate")
  panel_labels <- vapply(pairs, function(pair) {
    pair_id <- paste(pair, collapse = " vs ")
    ids <- sub("slingPseudotime_", "L", pair, fixed = TRUE)
    programs <- sub(" trajectory$", "", labels$assigned_program[match(pair, labels$lineage)])
    test <- tests[match(pair_id, tests$lineage), ]
    paste(
      paste(strwrap(paste0(ids[[2]], ": ", programs[[2]]), width = 36), collapse = "\n"),
      paste(strwrap(paste0("relative to ", ids[[1]], ": ", programs[[1]]), width = 36), collapse = "\n"),
      condiments_fdr_label(test$p_adj, test$status), sep = "\n")
  }, character(1))
  names(panel_labels) <- vapply(pairs, paste, character(1), collapse = " vs ")
  data <- purrr::map_dfr(pairs, function(pair) {
    pair_id <- paste(pair, collapse = " vs ")
    pair_weight <- rowSums(input$weights[, pair, drop = FALSE])
    purrr::map_dfr(conditions, function(condition) {
      keep <- input$condition == condition & pair_weight > 0
      relative_weight <- input$weights[keep, pair[[2]]] / pair_weight[keep]
      tibble::tibble(relative_weight = relative_weight) |>
        dplyr::count(relative_weight, name = "n_cells") |>
        dplyr::arrange(.data$relative_weight) |>
        dplyr::mutate(cumulative = cumsum(.data$n_cells) / sum(.data$n_cells),
          mean_relative_weight = stats::weighted.mean(.data$relative_weight, .data$n_cells),
          lineage = pair_id, condition = condition)
    })
  })
  data$condition <- factor(data$condition, levels = conditions)
  means <- dplyr::distinct(data, lineage, condition, mean_relative_weight)
  ggplot2::ggplot(data, ggplot2::aes(.data$relative_weight, .data$cumulative,
                                      color = .data$condition)) +
    ggplot2::geom_hline(yintercept = 0.5, color = "grey85", linewidth = 0.4) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::geom_vline(data = means, ggplot2::aes(xintercept = .data$mean_relative_weight,
      color = .data$condition), linetype = "dashed", linewidth = 0.55,
      show.legend = FALSE) +
    ggplot2::geom_text(data = means, ggplot2::aes(x = Inf,
      y = ifelse(.data$condition == conditions[[1]], 0.12, 0.04),
      label = sprintf("mean: %.2f", .data$mean_relative_weight), color = .data$condition),
      inherit.aes = FALSE, hjust = 1.05, size = 3.4, show.legend = FALSE) +
    ggplot2::facet_wrap(~ lineage, nrow = 1,
      labeller = ggplot2::as_labeller(panel_labels)) +
    ggplot2::scale_color_manual(values = condition_colors, labels = condition_plot_label) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
      labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 10),
      panel.spacing.x = grid::unit(1.2, "lines"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10),
      plot.margin = ggplot2::margin(12, 18, 12, 12)) +
    ggplot2::labs(title = paste(input$compartment, "Fate selection: which lineage has stronger relative membership?"),
      subtitle = paste("Focused comparison:", comparison$comparison_id[[1]],
        "| Cell-level distributions and cached MMD FDRs"),
      caption = paste(
        "For L1 vs L3, x = L3 weight / (L1 + L3 weight): 0 favors L1, 1 favors L3. Dashed lines and numbers: cell means.",
        "The MMD test detects any distribution difference, not its direction.",
        sep = "\n"),
      x = "Relative weight of the named lineage within the pair", y = "Cumulative fraction of cells", color = NULL)
}

# Show every mouse as a separate bar; do not imply that D1/D2 are paired across
# treatment groups or that two mice support a precise population-level estimate.
plot_condiments_progression_mice <- function(summary, annotations, compartment,
                                             exclude_conditions) {
  conditions <- display_conditions(summary$condition, exclude_conditions)
  data <- summary |>
    dplyr::filter(.data$compartment == .env$compartment,
                  .data$condition %in% conditions) |>
    dplyr::mutate(condition = factor(.data$condition, levels = rev(conditions)),
                  mouse = make_mouse_labels(.data$sample, as.character(.data$condition)))
  # Horizontal dodging stacks the first level at the bottom; reverse the levels
  # so D1 is the upper bar in each pair, matching the legend's reading order.
  data$mouse <- factor(data$mouse, levels = rev(sort(unique(data$mouse))))
  labels <- condiments_lineage_labels(annotations, compartment)
  ggplot2::ggplot(data, ggplot2::aes(.data$weighted_mean_pseudotime,
                                      .data$condition, fill = .data$mouse)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                      width = 0.65, orientation = "y") +
    ggplot2::facet_wrap(~ lineage, nrow = 1, scales = "free_x",
                        labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_y_discrete(labels = condition_plot_label) +
    ggplot2::scale_fill_manual(values = c(D1 = "#32688E", D2 = "#A7BCCB")) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank(),
                   panel.spacing.x = grid::unit(1.2, "lines"),
                   plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(title = paste(compartment, "Progression by mouse"),
      subtitle = "Descriptive | Each bar is one mouse",
      caption = "Each bar is one mouse's lineage-weighted mean pseudotime. Compare constructs within one panel; pseudotime scales differ between panels.",
      x = "Lineage-weighted mean shared pseudotime", y = NULL, fill = "Mouse")
}

# Mean normalized lineage weights give a soft composition per mouse: they sum
# to one but are not counts of cells assigned to one lineage or terminal fate.
plot_condiments_fate_mice <- function(summary, annotations, compartment,
                                      exclude_conditions) {
  conditions <- display_conditions(summary$condition, exclude_conditions)
  data <- summary |>
    dplyr::filter(.data$compartment == .env$compartment,
                  .data$condition %in% conditions) |>
    dplyr::mutate(condition = factor(.data$condition, levels = conditions),
                  mouse = make_mouse_labels(.data$sample, as.character(.data$condition)),
                  lineage = factor(.data$lineage,
                    levels = annotations$lineage[annotations$compartment == compartment]))
  totals <- dplyr::distinct(data, condition, mouse, sample, n_cells)
  labels <- condiments_lineage_labels(annotations, compartment)
  colors <- stats::setNames(c("#0072B2", "#D55E00", "#009E73"), names(labels))
  ggplot2::ggplot(data, ggplot2::aes(.data$mouse, .data$mean_lineage_probability,
                                      fill = .data$lineage)) +
    ggplot2::geom_col(width = 0.7, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(data = totals, ggplot2::aes(.data$mouse, 1.04,
                       label = paste0("n=", .data$n_cells)),
                       inherit.aes = FALSE, size = 3.2) +
    ggplot2::facet_wrap(~ condition, nrow = 1,
                        labeller = ggplot2::as_labeller(function(x)
                          condition_label_string(x), default = ggplot2::label_parsed)) +
    ggplot2::scale_fill_manual(values = colors, labels = labels) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25),
                                labels = function(x) paste0(round(100 * x), "%"),
                                limits = c(0, 1.1), expand = c(0, 0)) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", strip.background = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(title = paste(compartment, "Lineage membership by mouse"),
      subtitle = "Descriptive | Each bar is one mouse; n = cells",
      caption = "Segments are mean normalized lineage weights, not counts of cells assigned to one lineage.",
      x = "Mouse", y = "Mean lineage membership", fill = NULL)
}

# Keep a topic's test evidence and descriptive views on separate PDF pages.
save_condiments_pages <- function(file, pages) {
  grDevices::pdf(file, width = 12, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page in pages) print(page)
}
