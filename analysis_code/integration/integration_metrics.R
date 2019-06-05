# Run metrics 
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

method <- Seurat:::ExtractField(
	string = Seurat:::ExtractField(string = basename(args[1]), field = 5, delim = "_"),
	field = 1,
	delim = "\\."
)
dataset.name <- Seurat:::ExtractField(
  string = Seurat:::ExtractField(string = basename(args[1]), field = 2, delim = "_"),
  field = 1,
  delim = "\\."
)

dataset <- readRDS(file = args[1])
if (method == "none") {
  # Calculate metrics on uncorrected data
  genes.use <- dataset@misc$integration.features
  DefaultAssay(object = dataset) <- "RNA"
  dataset[["integrated"]] <- NULL
  dataset <- ScaleData(object = dataset, features = genes.use, verbose = FALSE)
  dataset <- RunPCA(object = dataset, features = genes.use, npcs = 30, verbose = FALSE)
}

reduction <- "pca"
dims <- 1:30
if (method == "seuratV2") {
  reduction <- "cca.aligned"
#  dims <- 1:10
#  if (dataset.name == "bipolar") {
#    dims <- 1:20
#  }
}

# silhouette metric
library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = dataset[[reduction]])[, dims])
clusters <- dataset$celltype
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
dataset$sil <- sil[, 3]

# mixing metric
max.k <- 300
mm <- max.k - MixingMetric(object = dataset, grouping.var = "replicate", reduction = reduction, dims = dims, max.k = max.k)

DefaultAssay(object = dataset) <- "RNA"
# Local structure preservation
ls <- LocalStruct(object = dataset, grouping.var = "replicate", reduction = reduction, reduced.dims = dims, orig.dims = 1:30)
ls <- unname(obj = unlist(x = ls))

all.metrics <- list(
  silhouette = dataset$sil, 
  mixing.metric = mm,
  local.struct = ls
)

saveRDS(object = all.metrics, file = args[2])
