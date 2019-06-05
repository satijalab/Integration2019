args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[1])
seurat.object <- readRDS(file = args[2])

DefaultAssay(object = seurat.object) <- "RNA"
seurat.object[["integrated"]] <- NULL

if (grepl(pattern = "integrated_mca.rds", x = args[1])) {
  dims.use <- 1:100
} else {
  dims.use <- 1:30
}

seurat.object <- FindVariableFeatures(object = seurat.object, selection.method = "vst", nfeatures = 2000)
seurat.object <- ScaleData(object = seurat.object, features = VariableFeatures(object = seurat.object))
seurat.object <- RunPCA(object = seurat.object, features = VariableFeatures(object = seurat.object), npcs = max(dims.use) + 10, verbose = FALSE)
seurat.object <- RunUMAP(object = seurat.object, reduction = "pca", dims = dims.use)
seurat.object <- RunTSNE(object = seurat.object, reduction = "pca", dims = dims.use, tsne.method = 'FIt-SNE')

saveRDS(object = seurat.object, file = args[3])
