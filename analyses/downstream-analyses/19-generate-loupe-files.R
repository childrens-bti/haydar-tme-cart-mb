#!/usr/bin/env Rscript

# Export the finalized all-cell, myeloid, and T-cell Seurat objects for
# interactive review in 10x Genomics Loupe Browser.  Descriptive biological
# annotations, rather than numeric Seurat clusters, are the primary identities
# to make the exported views interpretable without the analysis notebooks.

suppressPackageStartupMessages({
  library(rprojroot)
  library(Seurat)
  library(loupeR)
})

root_dir <- find_root(has_dir(".git"))
setwd(root_dir)
set.seed(1234)

output_dir <- file.path(root_dir, "analyses", "downstream-analyses", "results", "loupe")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

export_specs <- list(
  list(
    export_name = "all_cells",
    object_path = file.path(root_dir, "data", "v5", "cart_annotated.rds"),
    metadata_columns = c(
      "sample", "condition", "cell_type", "general_label", "immgen_label",
      "seurat_clusters"
    ),
    identity_column = "cell_type",
    output_stem = "cart_all_cells"
  ),
  list(
    export_name = "myeloid",
    object_path = file.path(root_dir, "data", "v5", "cart_myeloid_subtypes.rds"),
    metadata_columns = c(
      "sample", "condition", "cell_type", "general_label", "immgen_label",
      "seurat_clusters", "myeloid_subtype", "kdm6b_group"
    ),
    identity_column = "myeloid_subtype",
    kdm6b_ranking_path = file.path(
      root_dir, "analyses", "downstream-analyses", "results", "myeloid",
      "kdm6b", "myeloid_kdm6b_subcluster_ranking_and_groups.tsv"
    ),
    output_stem = "cart_myeloid_subtypes"
  ),
  list(
    export_name = "tcell",
    object_path = file.path(root_dir, "data", "v5", "cart_tcell_subtypes.rds"),
    metadata_columns = c(
      "sample", "condition", "cell_type", "general_label", "immgen_label",
      "seurat_clusters", "tcell_subtype", "kdm6b_group"
    ),
    identity_column = "tcell_subtype",
    kdm6b_ranking_path = file.path(
      root_dir, "analyses", "downstream-analyses", "results", "tcell",
      "kdm6b", "tcell_kdm6b_subcluster_ranking_and_groups.tsv"
    ),
    output_stem = "cart_tcell_subtypes"
  )
)

select_export_specs <- function(export_specs) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) == 0) {
    return(export_specs)
  }
  if (length(arguments) != 2 || arguments[[1]] != "--exports") {
    stop("Usage: 19-generate-loupe-files.R [--exports name1,name2]")
  }

  requested_exports <- trimws(strsplit(arguments[[2]], ",", fixed = TRUE)[[1]])
  available_exports <- vapply(export_specs, `[[`, character(1), "export_name")
  unknown_exports <- setdiff(requested_exports, available_exports)
  if (length(unknown_exports) > 0) {
    stop(
      "Unknown export(s): ", paste(unknown_exports, collapse = ", "),
      ". Available exports: ", paste(available_exports, collapse = ", ")
    )
  }

  Filter(
    function(specification) specification$export_name %in% requested_exports,
    export_specs
  )
}

export_specs <- select_export_specs(export_specs)

validate_reduction <- function(object, reduction_name, minimum_dimensions) {
  if (!reduction_name %in% Reductions(object)) {
    stop("Missing ", reduction_name, " reduction")
  }

  embeddings <- Embeddings(object, reduction = reduction_name)
  if (ncol(embeddings) < minimum_dimensions) {
    stop(
      "Reduction '", reduction_name, "' has ", ncol(embeddings),
      " dimensions; at least ", minimum_dimensions, " are required"
    )
  }

  embeddings
}

add_kdm6b_group <- function(object, specification) {
  if (is.null(specification$kdm6b_ranking_path)) {
    return(object)
  }
  if (!file.exists(specification$kdm6b_ranking_path)) {
    stop("Missing Kdm6b ranking table: ", specification$kdm6b_ranking_path)
  }
  if (!"seurat_clusters" %in% colnames(object[[]])) {
    stop("Missing seurat_clusters metadata required for Kdm6b group assignment")
  }

  kdm6b_ranking <- read.delim(
    specification$kdm6b_ranking_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_columns <- c("subcluster", "kdm6b_group")
  missing_columns <- setdiff(required_columns, colnames(kdm6b_ranking))
  if (length(missing_columns) > 0) {
    stop(
      "Kdm6b ranking table is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  kdm6b_ranking$subcluster <- as.character(kdm6b_ranking$subcluster)
  if (anyDuplicated(kdm6b_ranking$subcluster)) {
    stop("Kdm6b ranking table contains duplicate subcluster assignments")
  }
  kdm6b_by_cluster <- setNames(
    as.character(kdm6b_ranking$kdm6b_group),
    kdm6b_ranking$subcluster
  )
  object_clusters <- as.character(object[["seurat_clusters"]][, 1])
  kdm6b_group <- unname(kdm6b_by_cluster[object_clusters])
  if (anyNA(kdm6b_group)) {
    missing_clusters <- sort(unique(object_clusters[is.na(kdm6b_group)]))
    stop(
      "Kdm6b ranking table has no assignment for Seurat cluster(s): ",
      paste(missing_clusters, collapse = ", ")
    )
  }

  object[["kdm6b_group"]] <- factor(
    kdm6b_group,
    levels = c("Kdm6b-high", "Intermediate", "Kdm6b-low")
  )
  object
}

export_loupe <- function(specification) {
  if (!file.exists(specification$object_path)) {
    stop("Missing input Seurat object: ", specification$object_path)
  }

  object <- readRDS(specification$object_path)
  if (!inherits(object, "Seurat")) {
    stop("Input is not a Seurat object: ", specification$object_path)
  }
  if (!"RNA" %in% Assays(object)) {
    stop("Missing RNA assay with raw UMI counts: ", specification$object_path)
  }
  if (is.null(object[["RNA"]]$counts) || ncol(object[["RNA"]]$counts) == 0) {
    stop("RNA assay has no counts layer: ", specification$object_path)
  }
  # loupeR exports the active assay's counts layer.  Force the raw RNA assay
  # here so integrated, normalized, and scaled values are never exported.
  DefaultAssay(object) <- "RNA"
  object <- add_kdm6b_group(object, specification)
  missing_metadata_columns <- setdiff(
    specification$metadata_columns,
    colnames(object[[]])
  )
  if (length(missing_metadata_columns) > 0) {
    stop(
      "Missing metadata column(s) in ", specification$object_path, ": ",
      paste(missing_metadata_columns, collapse = ", ")
    )
  }

  metadata <- object[[]][, specification$metadata_columns, drop = FALSE]
  metadata[] <- lapply(metadata, as.character)
  invalid_metadata_columns <- names(metadata)[vapply(
    metadata,
    function(x) anyNA(x) || any(!nzchar(x)),
    logical(1)
  )]
  if (length(invalid_metadata_columns) > 0) {
    stop(
      "Metadata column(s) contain missing or empty labels: ",
      paste(invalid_metadata_columns, collapse = ", ")
    )
  }
  primary_annotation_column <- specification$identity_column
  annotation <- metadata[[primary_annotation_column]]

  pca_embeddings <- validate_reduction(object, "pca", minimum_dimensions = 2)
  umap_embeddings <- validate_reduction(object, "umap", minimum_dimensions = 2)
  if (!identical(sort(rownames(pca_embeddings)), sort(rownames(umap_embeddings)))) {
    stop("PCA and UMAP reductions do not contain the same cells")
  }

  loupe_path <- file.path(output_dir, paste0(specification$output_stem, ".cloupe"))
  if (file.exists(loupe_path)) {
    stop("Loupe file already exists; remove it before rerunning: ", loupe_path)
  }

  coordinate_table <- data.frame(
    cell_id = rownames(pca_embeddings),
    pca_1 = pca_embeddings[, 1],
    pca_2 = pca_embeddings[, 2],
    umap_1 = umap_embeddings[rownames(pca_embeddings), 1],
    umap_2 = umap_embeddings[rownames(pca_embeddings), 2],
    metadata[match(rownames(pca_embeddings), rownames(metadata)), , drop = FALSE],
    check.names = FALSE
  )
  coordinate_path <- file.path(
    output_dir, paste0(specification$output_stem, "_pca_umap_coordinates.tsv")
  )
  if (file.exists(coordinate_path)) {
    stop("Coordinate table already exists; remove it before rerunning: ", coordinate_path)
  }
  write.table(
    coordinate_table, coordinate_path, sep = "\t", quote = FALSE,
    row.names = FALSE, na = ""
  )

  # PCA remains available in the coordinate table, but is omitted from the
  # Cloupe file so the subtype-specific UMAP is the only projection on open.
  object[["pca"]] <- NULL
  # The inherited merged UMAP is not valid after subset-specific reprocessing.
  # Keep only the recomputed `umap` reduction so Loupe has one unambiguous UMAP
  # projection for each myeloid or T-cell export.
  if ("merged_umap" %in% Reductions(object)) {
    object[["merged_umap"]] <- NULL
  }
  Idents(object) <- factor(annotation)
  for (metadata_column in specification$metadata_columns) {
    object[[metadata_column]] <- factor(metadata[[metadata_column]])
  }

  create_loupe_from_seurat(
    obj = object,
    output_dir = output_dir,
    output_name = specification$output_stem,
    metadata_cols = specification$metadata_columns,
    dedup_clusters = FALSE
  )

  message("Created: ", loupe_path)
  message("Created: ", coordinate_path)
}

invisible(lapply(export_specs, export_loupe))
sessionInfo()
