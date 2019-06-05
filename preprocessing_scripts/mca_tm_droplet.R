# combine matrix files single matrix
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

input.dir <- paste0(getwd(), "/raw_data/mca/droplet/")

data.files <- paste0(input.dir, list.files(input.dir))
metadata <- read.csv(paste0(getwd(), "/raw_data/mca/annotations_droplets.csv"), row.names = 1)

data.matrices <- list()

message("Reading in TM droplet matrices")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.files))
for(i in 1:length(data.files)) {
  data.matrices[[i]] <- Read10X(data.dir = data.files[i])
  colnames(data.matrices[[i]]) <- paste0(
    Seurat:::ExtractField(string = Seurat:::ExtractField(string = basename(data.files[i]), field = 1:3), field = 2, delim = "-"),
    "_",
    colnames(data.matrices[[i]])
  )
  cells.to.keep <- intersect(x = colnames(x = data.matrices[[i]]), y = rownames(x = metadata))
  # only keep cells for which there is metadata
  data.matrices[[i]] <- data.matrices[[i]][, cells.to.keep]
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

message("Merging TM droplet matrices")
pb <- txtProgressBar(char = '=', style = 3, max = length(data.matrices), min = 2)
final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = data.matrices[[1]], mat2 = data.matrices[[2]])
for(i in 3:length(x = data.matrices)) {
  if (ncol(x = data.matrices[[i]]) == 0) next
  final.mat <- Seurat:::RowMergeSparseMatrices(mat1 = final.mat, mat2 = data.matrices[[i]])
  setTxtProgressBar(pb = pb, value = i)
}
close(pb)

metadata1 <- read.csv(file = paste0(getwd(), "/raw_data/mca/metadata_droplet.csv"), row.names = 1)
metadata2 <- metadata

metadata <- lapply(X = 1:nrow(metadata2), FUN = function(x) {
  channel <- Seurat:::ExtractField(string = rownames(metadata2)[x], field = 1:3)
  metadata1[channel, ]
})
metadata <- do.call(rbind, metadata)
rownames(metadata) <- rownames(metadata2)
metadata <- cbind(metadata, metadata2)

# updated github annotations
annotations <- read.csv(file = paste0(getwd(), "/raw_data/mca/annotations_droplet.csv"), row.names = 1,header = T)

# Basic Seurat object setup and preprocessing
seurat.object <- CreateSeuratObject(counts = final.mat, meta.data = metadata, project = "MCA_TM_DROPLET")
seurat.object <- NormalizeData(object = seurat.object, verbose = FALSE)
seurat.object <- subset(x = seurat.object, subset = nCount_RNA > 1000)

seurat.object[["tech"]] <- "10X"
Idents(object = seurat.object) <- "tech"

seurat.object <- AddMetaData(object = seurat.object, metadata = annotations)

# Save rds file
saveRDS(object = seurat.object, file = args[2])
