suppressMessages(library(data.table))
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

ad.loc <- paste0(getwd(), "/analysis_data/bipolar/")
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
  colnames(data) <- cell.names
  data.matrices[[i]] <- data
}

combined.mat <- do.call(what = cbind, args = data.matrices)
bipolar.integrated <- readRDS(file = args[1])
genes.use <- bipolar.integrated@misc$integration.features
combined.mat <- as.matrix(combined.mat[, colnames(x = bipolar.integrated)])

sc.assay <- CreateAssayObject(data = combined.mat)
bipolar.integrated[["integrated"]] <- sc.assay
DefaultAssay(object = bipolar.integrated) <- "integrated"
bipolar.integrated <- ScaleData(object = bipolar.integrated, features = genes.use)
bipolar.integrated <- RunPCA(object = bipolar.integrated, features = genes.use, npcs = 30, verbose = FALSE)
bipolar.integrated <- FindNeighbors(object = bipolar.integrated, reduction = "pca", dims = 1:30)
bipolar.integrated <- FindClusters(object = bipolar.integrated, resolution = 0.5)
bipolar.integrated <- RunUMAP(object = bipolar.integrated, reduction = "pca", dims = 1:30)
bipolar.integrated <- RunTSNE(object = bipolar.integrated, reduction = "pca", dims = 1:30, tsne.method = "FIt-SNE")

saveRDS(object = bipolar.integrated, file = args[2])
