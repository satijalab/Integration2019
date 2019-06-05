library(methods)

source("analysis_code/projection/downsampling_functions.R")

bipolar <- readRDS("seurat_objects/bipolar.rds")
metadata <- bipolar@meta.data
bipolar_list <- SplitObject(bipolar, attribute.1 = 'orig.ident', subset.raw = TRUE)
bipolar6 <- bipolar_list$Bipolar6
bipolar <- readRDS("seurat_objects/bipolar_integrated_no_rep6.rds")
celltypes <- metadata$celltype
names(celltypes) <- rownames(metadata)

# subset cell types in query
celltype_df <- data_frame('celltype' = celltypes, 'cell' = names(celltypes)) %>% 
	  filter(!(celltype %in% c("Unknown", "Doublets"))) %>% 
	  group_by(celltype) %>% 
	  top_n(200)

sampled_celltypes <- celltype_df$celltype
names(sampled_celltypes) <- celltype_df$cell

bipolar6 <- SubsetData(bipolar6, cells.use = names(sampled_celltypes), subset.raw = TRUE)

sce <- as.SingleCellExperiment(bipolar)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
assays(sce)$counts <- as.matrix(assays(sce)$counts)
assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)

umi_range <- seq(100, 2000, by = 100)

results <- lapply(umi_range, runDownsampling, reference=bipolar, sce_reference=sce, query=bipolar6, celltypes=sampled_celltypes)
df <- as.data.frame(results[[1]])
for(i in 2:length(results)){
	df[i,] <- results[[i]]
}

write.table(df, "analysis_data/downsampling_bipolar_subsampled_cells.tsv", quote = FALSE, sep = "\t")
