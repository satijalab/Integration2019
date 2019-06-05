args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])

# args <- c("~/Projects/muir/seurat_objects/integrated_bipolar_celltype_holdouts_none.rds")

bipolar.integrated <- readRDS(file = args[1])
DefaultAssay(object = bipolar.integrated) <- "RNA"
bipolar.integrated[["integrated"]] <- NULL
all.bipolar <- SplitObject(object = bipolar.integrated, split.by = "replicate")
genes.use <- bipolar.integrated@misc$integration.features

suppressMessages(library(future))
plan(strategy = 'multicore', workers = 6)
options(future.globals.maxSize = 10^12)

bipolar.anchors <- FindIntegrationAnchors(
  object.list = all.bipolar, 
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

bipolar.integrated <- IntegrateData(
  anchorset = bipolar.anchors, 
  new.assay.name = "integrated", 
  dims = 1:30, 
  k.weight = 100, 
  sd.weight = 1, 
  eps = 0, 
  verbose = TRUE
)


# take through clustering
bipolar.integrated[["RNA"]]@key <- "rna_"
DefaultAssay(object = bipolar.integrated) <- "integrated"
bipolar.integrated <- ScaleData(object = bipolar.integrated, verbose = FALSE)
bipolar.integrated <- RunPCA(object = bipolar.integrated, verbose = FALSE, npcs = 30)
bipolar.integrated <- RunUMAP(object = bipolar.integrated, dims = 1:30)
bipolar.integrated <- RunTSNE(object = bipolar.integrated, dims = 1:30, tsne.method = "FIt-SNE")
bipolar.integrated <- FindNeighbors(object = bipolar.integrated, reduction = "pca", dims = 1:30)
bipolar.integrated <- FindClusters(object = bipolar.integrated, resolution = 0.3)
bipolar.integrated$metric_clusters <- bipolar.integrated[["integrated_snn_res.0.3"]]

saveRDS(object = bipolar.integrated, file = args[2])
