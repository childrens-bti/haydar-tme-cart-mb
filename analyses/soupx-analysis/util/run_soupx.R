# Function to run Soupx for background correction of scRNA-seq data
library(SoupX)
run_soupx <- function(input_dir, sample_name, guide_genes, use_clustering = TRUE,
                      min_cells = 1, min_features = 50) {
  # Load the raw data
  raw_mat <- Read10X(data.dir = file.path(input_dir, sample_name, "raw_bc"))
  filt_mat <- Read10X(data.dir = file.path(input_dir, sample_name, "filtered_bc"))
  
  # intersect genes
  common_genes <- intersect(rownames(raw_mat), rownames(filt_mat))
  # subset and reorder them identically
  raw_mat  <- raw_mat[common_genes, ]
  filt_mat <- filt_mat[common_genes, ]
  
  # Create a SoupChannel object
  sc <- SoupChannel(tod = raw_mat, toc = filt_mat)
  
  # basic seurat workflow from filtered data for clustering
  sobj <- CreateSeuratObject(counts = filt_mat, project = sample_name, min.cells = min_cells, min.features = min_features)
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj)
  sobj <- FindNeighbors(sobj, dims = 1:20)
  sobj <- FindClusters(sobj)
  
  sc <- setClusters(sc, sobj$seurat_clusters)
  
  # manually estimating contamination fraction using marker genes without clusters
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = guide_genes, clusters = use_clustering)
  sc <- calculateContaminationFraction(sc, nonExpressedGeneList = guide_genes, useToEst = useToEst,
                                       forceAccept = TRUE)
  
  if (sc$metaData$rho[1] > 0.35) {
    warning(paste0("Setting contamination fraction to 0.35 for sample ", sample_name, 
                   " since estimated value is too high: ", round(sc$metaData$rho[1], 3)))
    sc = setContaminationFraction(sc, 0.35)
  }
  
  # Adjust counts to remove contamination
  clean_counts <- adjustCounts(sc, roundToInt = TRUE)
  
  # build final corrected Seurat object
  clean_sobj <- CreateSeuratObject(counts = clean_counts, project = sample_name,
                                   min.cells = min_cells, min.features = min_features)
  
  return(clean_sobj)
}