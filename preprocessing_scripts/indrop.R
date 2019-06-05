suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

input.file <- paste0(getwd(), "/raw_data/pancreas/indrop.csv.gz")

raw.data <- suppressWarnings(expr = fread(
  cmd = paste0("gzip -dc ", input.file),
  data.table = FALSE,
  fill = FALSE,
  showProgress = FALSE)
)
rownames(raw.data) <- raw.data[, 1]
raw.data <- raw.data[, -1]

# second column corresponds to the assigned cluster ID
clusters <- raw.data[, 2]
names(x = clusters) <- rownames(x = raw.data)
# exclude barcode and cluster ID from expression matrix
# transpose to fit with expected Seurat input format
raw.data <- raw.data[, c(-1, -2)]
raw.data <- Matrix(data = as.matrix(x = raw.data), sparse = TRUE)
raw.data <- t(x = raw.data)

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = raw.data, project = "INDROP")
# Not filtering here, all cells are annotated with cluster info
seurat.object <- AddMetaData(object = seurat.object, metadata = clusters, col.name = "assigned_cluster")
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object[["tech"]] <- "indrop"
seurat.object[["replicate"]] <- seurat.object[["orig.ident"]]
Idents(object = seurat.object) <- "tech"
object.list <- SplitObject(object = seurat.object, split.by = "replicate")
for (i in 1:length(x = object.list)) {
  object.list[[i]] <- FindVariableFeatures(
    object = object.list[[i]], 
    selection.method = "vst", 
    nfeatures = 2000,
    verbose = FALSE
  )
}
# Save rds file
saveRDS(object = object.list, file = args[2])
