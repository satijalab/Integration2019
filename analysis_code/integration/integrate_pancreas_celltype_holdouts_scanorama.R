suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

ad.loc <- paste0(getwd(), "/analysis_data/pancreas/")
corrected.matrices <- paste0(
  ad.loc,
  list.files(path = ad.loc)[grepl(pattern = "scanorama_corrected", x = list.files(path = ad.loc))]
)

data.matrices <- list()
for (i in 1:length(x = corrected.matrices)) {
  data <- fread(input = corrected.matrices[i], data.table = FALSE, header = TRUE)
  rownames(x = data) <- data$Genes
  data <- data[, -1]
  # get the cell names from the txt file
  txt.filename <- paste0(ad.loc, Seurat:::ExtractField(string = basename(corrected.matrices[i]), field = 1, delim = ".scanorama"), ".txt")
  cell.names <- names(x = read.table(file = txt.filename, nrows = 1, header = TRUE, row.names = 1))
  if (!grepl(pattern = "indrop", x = basename(corrected.matrices[i]))){
    cell.names <- gsub(pattern = "\\.", replacement = "-", x = cell.names)
  }
  cell.names <- gsub(pattern = "X", replacement = "", x = cell.names)
  colnames(data) <- cell.names
  data.matrices[[i]] <- data
}

combined.mat <- do.call(what = cbind, args = data.matrices)
pancreas.integrated <- readRDS(file = args[1])
genes.use <- pancreas.integrated@misc$integration.features
combined.mat <- as.matrix(combined.mat[, colnames(x = pancreas.integrated)])

sc.assay <- CreateAssayObject(data = combined.mat)
pancreas.integrated[["integrated"]] <- sc.assay
DefaultAssay(object = pancreas.integrated) <- "integrated"

pancreas.integrated <- ScaleData(object = pancreas.integrated, features = genes.use)
pancreas.integrated <- RunPCA(object = pancreas.integrated, features = genes.use, npcs = 30, verbose = FALSE)
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- RunUMAP(object = pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- RunTSNE(object = pancreas.integrated, reduction = "pca", dims = 1:30, tsne.method = "FIt-SNE")

saveRDS(object = pancreas.integrated, file = args[2])
