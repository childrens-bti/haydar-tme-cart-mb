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
    object_path = file.path(root_dir, "data", "v5", "cart_annotated.rds"),
    annotation_column = "cell_type",
    output_stem = "cart_all_cells"
  ),
  list(
    object_path = file.path(root_dir, "data", "v5", "cart_myeloid_subtypes.rds"),
    annotation_column = "myeloid_subtype",
    output_stem = "cart_myeloid_subtypes"
  ),
  list(
    object_path = file.path(root_dir, "data", "v5", "cart_tcell_subtypes.rds"),
    annotation_column = "tcell_subtype",
    output_stem = "cart_tcell_subtypes"
  )
)

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
  if (!specification$annotation_column %in% colnames(object[[]])) {
    stop(
      "Missing descriptive annotation column '",
      specification$annotation_column, "' in ", specification$object_path
    )
  }

  annotation <- as.character(object[[specification$annotation_column]][, 1])
  if (anyNA(annotation) || any(!nzchar(annotation))) {
    stop(
      "Annotation column '", specification$annotation_column,
      "' contains missing or empty labels"
    )
  }

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
    cluster = annotation[match(rownames(pca_embeddings), rownames(object[[]]))],
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

  # loupeR serializes only two-dimensional reductions.  Retain PC1/PC2 under
  # the standard 'pca' name in this in-memory export copy and leave the source
  # RDS unchanged; UMAP is already two dimensional.
  object[["pca"]] <- CreateDimReducObject(
    embeddings = pca_embeddings[, 1:2, drop = FALSE],
    key = "PC_",
    assay = DefaultAssay(object)
  )
  Idents(object) <- factor(annotation)
  object[[specification$annotation_column]] <- factor(annotation)

  create_loupe_from_seurat(
    obj = object,
    output_dir = output_dir,
    output_name = specification$output_stem,
    metadata_cols = specification$annotation_column,
    dedup_clusters = TRUE
  )

  message("Created: ", loupe_path)
  message("Created: ", coordinate_path)
}

invisible(lapply(export_specs, export_loupe))
sessionInfo()
