args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[3])
#args <- c("~/Projects/muir/seurat_objects/bipolar.rds")

bipolar <- readRDS(file = args[1])
all.bipolar <- SplitObject(object = bipolar, split.by = "replicate")
for(i in 1:length(x = all.bipolar)) {
  all.bipolar[[i]] <- FindVariableFeatures(
    object = all.bipolar[[i]], 
    selection.method = "vst", 
    nfeatures = 2000, 
    verbose = FALSE)
}

suppressMessages(library(future))
plan(strategy = 'multicore', workers = 6)
options(future.globals.maxSize = 10^12)

bipolar.anchors <- FindIntegrationAnchors(
  object.list = all.bipolar, 
  anchor.features = 2000, 
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

# take through clustering and annotate
DefaultAssay(object = bipolar.integrated) <- "integrated"
bipolar.integrated <- ScaleData(
  object = bipolar.integrated, 
  features = VariableFeatures(object = bipolar.integrated), 
  verbose = FALSE
)
bipolar.integrated <- RunPCA(
  object = bipolar.integrated,
  features = VariableFeatures(object = bipolar.integrated),
  verbose = FALSE, npcs = 30
)
bipolar.integrated <- RunTSNE(
  object = bipolar.integrated,
  reduction = "pca",
  dims = 1:30, 
  tsne.method = "FIt-SNE",
  verbose = FALSE
)
bipolar.integrated <- RunUMAP(
  object = bipolar.integrated,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE,
  nneighbors = 5
)
bipolar.integrated <- FindNeighbors(
  object = bipolar.integrated, 
  reduction = "pca", 
  dims = 1:30, 
  verbose = FALSE
)
bipolar.integrated <- FindClusters(object = bipolar.integrated, resolution = 0.4)

saveRDS(object = bipolar.integrated, file = args[2])
