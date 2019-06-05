library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

activity.scores <- readRDS(file = "raw_data/atac/activity_scores.quantitative.rds")
cell.md <- read.table(file = "raw_data/atac/cell_metadata.txt", header = TRUE, row.names = 1, sep = "\t")

atac <- CreateSeuratObject(counts = activity.scores)
atac <- AddMetaData(object = atac, metadata = cell.md)
pfc <- subset(x = atac, tissue == "PreFrontalCortex")

atac.pks <- readRDS(file = "raw_data/atac/atac_matrix.binary.qc_filtered.rds")

pfc.pks <- atac.pks[, colnames(x = pfc)]
pfc[['peaks']] <- CreateAssayObject(counts = pfc.pks)
pfc <- subset(x = pfc, subset = nCount_peaks > 5000)
pfc <- subset(x = pfc, subset = cell_label %in% c('Collisions', 'Unknown'), invert = TRUE)
pfc <- NormalizeData(object = pfc, scale.factor = 1e6)
DefaultAssay(pfc) <- 'peaks'
VariableFeatures(pfc) <- names(which(Matrix::rowSums(pfc) > 100))
pfc <- RunLSI(object = pfc, features = VariableFeatures(pfc), n = 50, scale.max = NULL)
pfc <- FindNeighbors(object = pfc, dims = 1:30, reduction = "lsi")
pfc <- FindClusters(object = pfc, graph = "peaks_snn", resolution = 1)

saveRDS(object = pfc, file = args[2])
