suppressMessages(library(data.table))
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])

input.file <- paste0(getwd(), "/raw_data/pancreas/smartseq2.txt")

raw.data <- suppressWarnings(expr = fread(
  input = input.file,
  data.table = FALSE,
  fill = FALSE)
)
cell.names <- colnames(x = raw.data)[2:3515]

# weird file structure here. According to the provided README:
# First line contains sample IDS (starting with #samples)
# Columns 1:3514 correspond to rpkm values
# Columns 3515:7028 correspond to read counts
# Rows 1:26179 correspond to data for RefSeq genes
# Rows 26180:26271 correspond to data for the 92 external RNA spike-in controls (ERCCs)
# Row 26272 (last) contains data for ‘eGFP’

# Keep only RefSeq genes, ERCCs start at 26179 here due to indexing off by one from README
raw.data <- raw.data[1:26178, ]
rownames(raw.data) <- make.unique(names = raw.data[, 1])
raw.data <- raw.data[, 3517:7030]
colnames(raw.data) <- cell.names
raw.data <- Matrix(data = as.matrix(x = raw.data), sparse = TRUE)

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = raw.data, project = "SMARTSEQ2")
seurat.object <- subset(x = seurat.object, subset = nFeature_RNA > 2500)
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(object = seurat.object, verbose = FALSE,
                                      selection.method = "vst", nfeatures = 2000)
seurat.object[["tech"]] <- "smartseq2"
seurat.object[["replicate"]] <- "smartseq2"
Idents(object = seurat.object) <- "tech"

# Save rds file
saveRDS(object = seurat.object, file = args[2])
