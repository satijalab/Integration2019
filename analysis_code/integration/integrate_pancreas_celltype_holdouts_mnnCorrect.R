# weird conda issue with installing R data packages from bioconductor
# just reinstalls here

if(!"GenomeInfoDbData" %in% rownames(installed.packages())){
    source("https://bioconductor.org/biocLite.R")
    biocLite("GenomeInfoDbData")
}

suppressMessages(library(Matrix))
suppressMessages(library(scran))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

set.seed(42)

pancreas.integrated <- readRDS(file = args[1])
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated[["integrated"]] <- NULL
genes.use <- pancreas.integrated@misc$integration.features
Idents(object = pancreas.integrated) <- "replicate"

celseq <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "celseq")]
celseq2 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "celseq2")]
smartseq2 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "smartseq2")]
fluidigmc1 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "fluidigmc1")]
indrop1 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "indrop1")]
indrop2 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "indrop2")]
indrop3 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "indrop3")]
indrop4 <- GetAssayData(object = pancreas.integrated, assay = "RNA", slot = "data")[genes.use, WhichCells(object = pancreas.integrated, idents = "indrop4")]

mnnc.results <- mnnCorrect(
  as.matrix(celseq),
  as.matrix(celseq2),
  as.matrix(smartseq2),
  as.matrix(fluidigmc1),
  as.matrix(indrop1),
  as.matrix(indrop2),
  as.matrix(indrop3),
  as.matrix(indrop4),
  pc.approx = TRUE
)

mnnc.combined <- do.call(what = cbind, args = mnnc.results$corrected)
colnames(x = mnnc.combined) <- c(colnames(x = celseq), colnames(x = celseq2), 
                                 colnames(x = smartseq2), colnames(x = fluidigmc1), 
                                 colnames(x = indrop1), colnames(x = indrop2), colnames(x = indrop3),
                                 colnames(x = indrop4))
mnnc.combined <- mnnc.combined[, Cells(x = pancreas.integrated)]
mnn.combined <- CreateAssayObject(data = mnnc.combined)

pancreas.integrated[["integrated"]] <- mnn.combined
DefaultAssay(object = pancreas.integrated) <- "integrated"
VariableFeatures(object = pancreas.integrated) <- genes.use
pancreas.integrated <- ScaleData(object = pancreas.integrated)
pancreas.integrated <- RunPCA(object = pancreas.integrated, features = genes.use, npcs = 30, verbose = FALSE)
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- RunTSNE(object = pancreas.integrated, dims = 1:30, reduction = "pca", tsne.method = "FIt-SNE")
pancreas.integrated <- RunUMAP(object = pancreas.integrated, dims = 1:30, reduction = "pca")

saveRDS(object = pancreas.integrated, file = args[2])




