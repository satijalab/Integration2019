# weird conda issue with installing R data packages from bioconductor
# just reinstalls here

if(!"GenomeInfoDbData" %in% rownames(installed.packages())){
    source("https://bioconductor.org/biocLite.R")
    biocLite("GenomeInfoDbData")
}

suppressMessages(library(scran))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

set.seed(42)

bipolar.integrated <- readRDS(file = args[1])
DefaultAssay(object = bipolar.integrated) <- "RNA"
bipolar.integrated[["integrated"]] <- NULL

genes.use <- bipolar.integrated@misc$integration.features
Idents(object = bipolar.integrated) <- "replicate"

rep1 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 1)]
rep2 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 2)]
rep3 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 3)]
rep4 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 4)]
rep5 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 5)]
rep6 <- GetAssayData(object = bipolar.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = bipolar.integrated, idents = 6)]

mnnc.results <- mnnCorrect(
  as.matrix(x = rep1),
  as.matrix(x = rep2),
  as.matrix(x = rep3),
  as.matrix(x = rep4),
  as.matrix(x = rep5),
  as.matrix(x = rep6),
  pc.approx = TRUE
)

mnnc.combined <- do.call(what = cbind, args = mnnc.results$corrected)
colnames(x = mnnc.combined) <- c(colnames(x = rep1), colnames(x = rep2), 
                                 colnames(x = rep3), colnames(x = rep4), 
                                 colnames(x = rep5), colnames(x = rep6))
mnnc.combined <- mnnc.combined[, Cells(x = bipolar.integrated)]
mnn.combined <- CreateAssayObject(data = mnnc.combined)

bipolar.integrated[["integrated"]] <- mnn.combined
DefaultAssay(object = bipolar.integrated) <- "integrated"

VariableFeatures(object = bipolar.integrated) <- genes.use
bipolar.integrated <- ScaleData(object = bipolar.integrated)
bipolar.integrated <- RunPCA(object = bipolar.integrated, features = genes.use, npcs = 30, verbose = FALSE)
bipolar.integrated <- FindNeighbors(object = bipolar.integrated, reduction = "pca", dims = 1:30)
bipolar.integrated <- FindClusters(object = bipolar.integrated, resolution = 0.5)
bipolar.integrated <- RunTSNE(object = bipolar.integrated, dims = 1:30, reduction = "pca", tsne.method = "FIt-SNE")
bipolar.integrated <- RunUMAP(object = bipolar.integrated, dims = 1:30, reduction = "pca")

saveRDS(object = bipolar.integrated, file = args[2])




