suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

input.file <- paste0(getwd(), "/raw_data/pancreas/celseq.csv.gz")

raw.data <- suppressWarnings(expr = fread(
  cmd = paste0("gzip -dc ", input.file), 
  data.table = FALSE,
  fill = FALSE,
  showProgress = FALSE)
)
rownames(raw.data) <- raw.data[, 1]
raw.data <- raw.data[, -1]

# remove chromosome annotation and enforce unique gene names
rownames(raw.data) <- make.unique(sapply(X = rownames(raw.data), FUN = function(x) {
  Seurat:::ExtractField(string = x, field = 1, delim = "__")
}))

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = raw.data, project = "CELSEQ")
seurat.object <- subset(x = seurat.object, subset = nFeature_RNA > 1750)
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(object = seurat.object, verbose = FALSE, 
                                      selection.method = "vst", nfeatures = 2000)
seurat.object[["tech"]] <- "celseq"
seurat.object[["replicate"]] <- "celseq"
Idents(object = seurat.object) <- "tech"

# Save rds file
saveRDS(object = seurat.object, file = args[2])
