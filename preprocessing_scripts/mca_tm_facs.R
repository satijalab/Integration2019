# combine matrix files single matrix
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)

input.dir <- paste0(getwd(), "/raw_data/mca/FACS/")

data.files <- paste0(input.dir, list.files(input.dir))
metadata <- read.csv(paste0(getwd(), "/raw_data/mca/annotations_FACS.csv"), row.names = 1)

data.matrices <- list()

message("Reading in TM FACS matrices")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.files))
for(i in 1:length(data.files)) {
  data.mat <- suppressWarnings(fread(input = data.files[i], showProgress = FALSE))
  rows <- as.vector(as.matrix(data.mat[, 1]))
  data.mat[, 1] <- NULL
  data.matrices[[i]] <- as.matrix(x = data.mat)
  rownames(data.matrices[[i]]) <- rows
  # only keep cells for which there is metadata
  cells.to.keep <- intersect(colnames(data.matrices[[i]]), rownames(metadata))
  data.matrices[[i]] <- data.matrices[[i]][]
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

message("Merging TM FACS matrices")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.matrices), min = 2)
final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = data.matrices[[1]], mat2 = data.matrices[[2]])
for(i in 3:length(x = data.matrices)) {
  if (ncol(x = data.matrices[[i]]) == 0) next
  final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = final.mat, mat2 = data.matrices[[i]])
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

metadata1 <- read.csv(file = paste0(getwd(), "/raw_data/mca/metadata_FACS.csv"), row.names = 1)
metadata2 <- metadata

metadata <- lapply(X = 1:nrow(metadata2), FUN = function(x) {
  id <- Seurat:::ExtractField(string = rownames(metadata2)[x], field = 2, delim = "\\.")
  metadata1[id, ]
})
metadata <- do.call(rbind, metadata)
rownames(metadata) <- rownames(metadata2)
metadata <- cbind(metadata, metadata2)

# updated github annotations
annotations <- read.csv(file = paste0(getwd(), "/raw_data/mca/annotations_facs.csv"), header = T)
rownames(x = annotations) <- annotations$cell

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = final.mat, meta.data = metadata, project = "MCA_TM_FACS")
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object <- subset(x = seurat.object, subset = nFeature_RNA > 1000)

seurat.object[["tech"]] <- "FACS"
Idents(object = seurat.object) <- "tech"

seurat.object <- AddMetaData(object = seurat.object, metadata = annotations)

# Save rds file
saveRDS(object = seurat.object, file = args[2])
