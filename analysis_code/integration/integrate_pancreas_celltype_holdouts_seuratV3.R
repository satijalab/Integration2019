suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])
#args <- c("~/Projects/muir/seurat_objects/integrated_pancreas_celltype_holdouts_none.rds")

pancreas.integrated <- readRDS(file = args[1])

DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated[["integrated"]] <- NULL

all.pancreas <- SplitObject(object = pancreas.integrated, split.by = "replicate")
genes.use <- pancreas.integrated@misc$integration.features

pancreas.anchors <- FindIntegrationAnchors(
  object.list = all.pancreas, 
  anchor.features = genes.use, 
  scale = TRUE, 
  l2.norm = TRUE, 
  dims = 1:30, 
  k.anchor = 5, 
  k.filter = 200, 
  k.score = 30, 
  max.features = 200, 
  eps = 0, 
  verbose = TRUE
)

pancreas.integrated <- IntegrateData(
  anchorset = pancreas.anchors, 
  new.assay.name = "integrated", 
  dims = 1:30, 
  k.weight = 100, 
  sd.weight = 1, 
  eps = 0, 
  verbose = TRUE
)


# take through clustering
pancreas.integrated[["RNA"]]@key <- "rna_"
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(object = pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(object = pancreas.integrated, verbose = FALSE, npcs = 30)
pancreas.integrated <- RunUMAP(object = pancreas.integrated, dims = 1:30)
pancreas.integrated <- RunTSNE(object = pancreas.integrated, dims = 1:30, tsne.method = "FIt-SNE")
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated$metric_clusters <- pancreas.integrated[["integrated_snn_res.0.3"]]

saveRDS(object = pancreas.integrated, file = args[2])
