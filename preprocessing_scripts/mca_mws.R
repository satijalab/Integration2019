# combine matrix files single matrix
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

input.dir <- paste0(getwd(), "/raw_data/mca/500more_dge/")

data.files <- paste0(input.dir, list.files(input.dir))
assignments <- read.csv(paste0(getwd(), "/raw_data/mca/MCA_CellAssignments.csv"), row.names = 2)

data.matrices <- list()

message("Reading in MWS MCA files - this will take a while")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.files))
for(i in 1:length(x = data.files)) {
  data.mat <- suppressWarnings(fread(cmd = paste0("gzip -dc ", data.files[i]), showProgress = FALSE))
  rows <- as.vector(as.matrix(data.mat[, 1]))
  data.mat[, 1] <- NULL
  data.matrices[[i]] <- as.matrix(x = data.mat)
  rownames(data.matrices[[i]]) <- rows
  cells.to.keep <- intersect(x = colnames(x = data.matrices[[i]]), y = rownames(x = assignments))
  data.matrices[[i]] <- as(object = data.matrices[[i]][, cells.to.keep], Class = "dgCMatrix")
  invisible(gc())
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

message("Merging MCA files")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.matrices), min = 2)
final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = data.matrices[[1]], mat2 = data.matrices[[2]])
for(i in 3:length(x = data.matrices)) {
  if (ncol(x = data.matrices[[i]]) == 0) next
  final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = final.mat, mat2 = data.matrices[[i]])
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = final.mat, meta.data = assignments, project = "MCA_MICROWELL_SEQ")
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)

seurat.object[["tech"]] <- "uWell"
Idents(object = seurat.object) <- "tech"

# Save rds file
saveRDS(object = seurat.object, file = args[1])
