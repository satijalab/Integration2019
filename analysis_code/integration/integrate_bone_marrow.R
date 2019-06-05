library(Seurat)
suppressMessages(library(future))
set.seed(1234)

plan(strategy = 'multicore', workers = 8)
options(future.globals.maxSize = 10^12)

args <- commandArgs(trailingOnly = TRUE)

bm <- readRDS(file = args[1])
umap_coords <- readRDS("raw_data/ica/ica_umap_nn200.rds")

# remove doublets from each batch
bm.list.full <- SplitObject(object = bm)
accept.cells <- c()
for (i in 1:length(x = bm.list.full)) {
  dname <- names(x = bm.list.full)[[i]]
  scrublet <- read.table(file = paste0('raw_data/ica/ICA_scrublet_results/', dname, '_doublets.tsv'), sep = '\t')
  dubs <- scrublet[, 1] > 0.4
  accept.cells <- c(accept.cells, colnames(x = bm.list.full[[i]])[!dubs])
}
rm(bm.list.full)

bm <- bm[, accept.cells]
bm.list <- SplitObject(object = bm)
rm(bm)
gc()
bm.list <- sapply(X = bm.list, FUN = FindVariableFeatures, nfeatures = 2000)

bm.anchors <- FindIntegrationAnchors(
  object.list = bm.list,
  anchor.features = 2000,
  dims = 1:60,
  eps = 1,
  verbose = TRUE
)

saveRDS(object = bm.anchors, file = "analysis_data/bm.anchors.rds")

bm.integrated <- IntegrateData(
  anchorset = bm.anchors, 
  dims = 1:60,
  eps = 1
)

DefaultAssay(object = bm.integrated) <- "integrated"
bm.integrated <- ScaleData(object = bm.integrated)
bm.integrated <- RunPCA(object = bm.integrated, npcs = 60)

saveRDS(object = bm.integrated, file = args[2])