source("analysis_code/projection/downsampling_functions.R")

indrop <- readRDS("seurat_objects/inDrop.rds")
indrop <- SplitObject(indrop, attribute.1 = 'replicate', subset.raw = TRUE)
indrop <- indrop[[1]]
pancreas.integrated <- readRDS("seurat_objects/pancreas_integrated_no_indrop1.rds")
celltypes <- readRDS("analysis_data/2018_06_22_integrated_pancreas_clusterIDs.rds")
pancreas.integrated <- AddMetaData(pancreas.integrated, celltypes, col.name = 'celltype')

sce <- as.SingleCellExperiment(pancreas.integrated)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
assays(sce)$counts <- as.matrix(assays(sce)$counts)
assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)

umi_range <- seq(100, 10000, by = 100)
results <- lapply(umi_range, runDownsampling, reference=pancreas.integrated, sce_reference=sce, query=indrop)

df <- as.data.frame(results[[1]])
for(i in 2:length(results)){
  df[i,] <- results[[i]]
}

write.table(df, "analysis_data/downsampling_pancreas.tsv", quote = FALSE, sep = "\t")