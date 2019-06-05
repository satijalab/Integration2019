library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(alluvial)
library(dplyr)
library(methods)

bipolar <- readRDS("seurat_objects/bipolar.rds")
metadata <- bipolar@meta.data
bipolar_list <- SplitObject(bipolar, attribute.1 = 'orig.ident', subset.raw = TRUE)
bipolar6 <- bipolar_list$Bipolar6
bipolar_list$Bipolar6 <- NULL

bipolar <- MultiIntegrateData(object.list = bipolar_list, var.genes = bipolar_list[[1]]@var.genes,
                              order.by = 'similarity',
                              score.mnn = TRUE, k = 500, k.mnn = 8, score.mnn.k = 10,
                              eps = 0, reduction.use = 'cca.cosine', do.pairwise = TRUE,
                              dims.use = 1:20, sd = 1/2, do.cpp = TRUE,
                              verbose = TRUE, keep.sparse = TRUE)

bipolar <- FindVariableGenes(bipolar, do.plot = FALSE, selection.method = 'dispersion')
bipolar <- ScaleData(bipolar, genes.use = bipolar@var.genes)
bipolar <- RunPCA(bipolar, do.print = FALSE, pcs.compute = 30)

saveRDS(bipolar, file = "seurat_objects/bipolar_integrated_no_rep6.rds")

projections <- ProjectCells(reference = bipolar, query = bipolar6, label.name = 'celltype',
                           correct.query = TRUE, diet = TRUE, eps = 0,
                           do.cpp = TRUE, verbose = TRUE, sd = 1, k = 100,
                           use.cosine = TRUE, dims.use = 1:30, normalize.with.lm = TRUE)

projections$true_label <- metadata[rownames(projections),]$celltype
all_projections <- projections
projections[projections$prediction.score < 0.5 , 'predicted.id'] <- 'unassigned'

seurat_projection <- projections %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_accuracy <- sum(projections$predicted.id == projections$true_label) / nrow(projections) * 100

# run scmap
sce <- as.SingleCellExperiment(bipolar)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
assays(sce)$counts <- as.matrix(assays(sce)$counts)
assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = FALSE)

bipolar6.sce <- as.SingleCellExperiment(bipolar6)
rowData(bipolar6.sce)$feature_symbol <- rownames(bipolar6.sce)
assays(bipolar6.sce)$counts <- as.matrix(assays(bipolar6.sce)$counts)
assays(bipolar6.sce)$logcounts <- as.matrix(assays(bipolar6.sce)$logcounts)

# scmap-cluster
sce <- indexCluster(sce, cluster_col = 'celltype')
scmapCluster_results <- scmapCluster(projection = bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cluster_index))

# scmap-cell
sce <- indexCell(sce)
scmapCell_results <- scmapCell(projection = bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)))

# evaluate
celltypes <- metadata$celltype
names(celltypes) <- rownames(metadata)
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
rownames(scmap.cluster.results) <- rownames(colData(bipolar6.sce))
scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])

scmap_cluster <- scmap.cluster.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster <- as.data.frame(scmap_cluster)
scmap_cluster$yan <- as.character(scmap_cluster$yan)
scmap_cluster$true_label <- as.character(scmap_cluster$true_label)

scmap_cluster_accuracy <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
rownames(scmap.cell.results) <- rownames(colData(bipolar6.sce))
scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])

scmap_cell <- scmap.cell.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cell <- as.data.frame(scmap_cell)

scmap_cell_accuracy <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

pdf("figures/bipolar_projection.pdf", height = 8, width = 15)
par(mfrow=c(1,3))
alluvial(seurat_projection[,1:2], freq=seurat_projection$count)
mtext(paste0("Seurat:", round(seurat_accuracy, 2)), 3, line=3)
alluvial(scmap_cluster[,1:2], freq=scmap_cluster$count)
mtext(paste0("scmap-cluster: ", round(scmap_cluster_accuracy, 2)), 3, line=3)
alluvial(scmap_cell[,1:2], freq=scmap_cell$count)
mtext(paste0("scmap-cell: ", round(scmap_cell_accuracy, 2)), 3, line=3)
dev.off()

##### Dropping thresholds for scmap
# scmap-cluster
scmapCluster_results <- scmapCluster(projection = bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cluster_index), threshold = 0)

# scmap-cell
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)), w=1, threshold = 0)

# evaluate
celltypes <- metadata$celltype
names(celltypes) <- rownames(metadata)
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
rownames(scmap.cluster.results) <- rownames(colData(bipolar6.sce))
scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])

scmap_cluster <- scmap.cluster.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster <- as.data.frame(scmap_cluster)
scmap_cluster$yan <- as.character(scmap_cluster$yan)
scmap_cluster$true_label <- as.character(scmap_cluster$true_label)

scmap_cluster_accuracy <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
rownames(scmap.cell.results) <- rownames(colData(bipolar6.sce))
scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])

scmap_cell <- scmap.cell.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cell <- as.data.frame(scmap_cell)
scmap_cell_accuracy <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

seurat_projection <- all_projections %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_accuracy <- sum(all_projections$predicted.id == all_projections$true_label) / nrow(all_projections) * 100

pdf("figures/bipolar_projection_low_threshold.pdf", height = 8, width = 15)
par(mfrow=c(1,3))
alluvial(seurat_projection[,1:2], freq=seurat_projection$count)
mtext(paste0("Seurat:", round(seurat_accuracy, 2)), 3, line=3)
alluvial(scmap_cluster[,1:2], freq=scmap_cluster$count)
mtext(paste0("scmap-cluster: ", round(scmap_cluster_accuracy, 2)), 3, line=3)
alluvial(scmap_cell[,1:2], freq=scmap_cell$count)
mtext(paste0("scmap-cell: ", round(scmap_cell_accuracy, 2)), 3, line=3)
dev.off()

### Downsample query UMIs
subsample_bipolar6 <- SampleUMI(bipolar6@raw.data[,bipolar6@cell.names], max.umi = 500)
dimnames(subsample_bipolar6) <- dimnames(bipolar6@data)
subsample_bipolar6 <- CreateSeuratObject(subsample_bipolar6, min.cells = 0, min.genes = 0, project = 'bipolar6')
subsample_bipolar6 <- NormalizeData(subsample_bipolar6)

# seurat projection
projection <- ProjectCells(reference = bipolar, query = subsample_bipolar6, label.name = 'celltype',
                           correct.query = TRUE, diet = TRUE, eps = 0,
                           do.cpp = TRUE, verbose = TRUE, sd = 1, k = 100,
                           use.cosine = TRUE, dims.use = 1:30, normalize.with.lm = TRUE)

projection$true_label <- as.character(celltypes[rownames(projection)])
all_projections <- projection
projection[projection$prediction.score < 0.5, 'predicted.id'] <- 'unassigned'

seurat_projection <- projection %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_accuracy <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100

subsample_bipolar6.sce <- as.SingleCellExperiment(subsample_bipolar6)
rowData(subsample_bipolar6.sce)$feature_symbol <- rownames(subsample_bipolar6.sce)
assays(subsample_bipolar6.sce)$counts <- as.matrix(assays(subsample_bipolar6.sce)$counts)
assays(subsample_bipolar6.sce)$logcounts <- as.matrix(assays(subsample_bipolar6.sce)$logcounts)

# scmap-cluster
scmapCluster_results <- scmapCluster(projection = subsample_bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cluster_index))

# scmap-cell
scmapCell_results <- scmapCell(projection = bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)))

# evaluate
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
rownames(scmap.cluster.results) <- rownames(colData(subsample_bipolar6.sce))
scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])

scmap_cluster <- scmap.cluster.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster <- as.data.frame(scmap_cluster)
scmap_cluster$yan <- as.character(scmap_cluster$yan)
scmap_cluster$true_label <- as.character(scmap_cluster$true_label)

scmap_cluster_accuracy <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
rownames(scmap.cell.results) <- rownames(colData(subsample_bipolar6.sce))
scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])

scmap_cell <- scmap.cell.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cell <- as.data.frame(scmap_cell)
scmap_cell_accuracy <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

pdf("figures/bipolar_projection_500umi.pdf", height = 8, width = 15)
par(mfrow=c(1,3))
alluvial(seurat_projection[,1:2], freq=seurat_projection$count)
mtext(paste0("Seurat:", round(seurat_accuracy, 2)), 3, line=3)
alluvial(scmap_cluster[,1:2], freq=scmap_cluster$count)
mtext(paste0("scmap-cluster: ", round(scmap_cluster_accuracy, 2)), 3, line=3)
alluvial(scmap_cell[,1:2], freq=scmap_cell$count)
mtext(paste0("scmap-cell: ", round(scmap_cell_accuracy, 2)), 3, line=3)
dev.off()

### Drop thresholds
# scmap-cluster
scmapCluster_results <- scmapCluster(projection = subsample_bipolar6.sce, index_list = list(yan = metadata(sce)$scmap_cluster_index), threshold = 0)

# scmap-cell
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)), w = 1, threshold = 0)

# evaluate
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
rownames(scmap.cluster.results) <- rownames(colData(subsample_bipolar6.sce))
scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])

scmap_cluster <- scmap.cluster.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster <- as.data.frame(scmap_cluster)
scmap_cluster$yan <- as.character(scmap_cluster$yan)
scmap_cluster$true_label <- as.character(scmap_cluster$true_label)

scmap_cluster_accuracy <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
rownames(scmap.cell.results) <- rownames(colData(subsample_bipolar6.sce))
scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])

scmap_cell <- scmap.cell.results %>% 
  group_by(yan, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cell <- as.data.frame(scmap_cell)
scmap_cell_accuracy <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

seurat_projection <- all_projections %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_accuracy <- sum(all_projections$predicted.id == all_projections$true_label) / nrow(all_projections) * 100

pdf("figures/bipolar_projection_500umi_low_threshold.pdf", height = 8, width = 15)
par(mfrow=c(1,3))
alluvial(seurat_projection[,1:2], freq=seurat_projection$count)
mtext(paste0("Seurat:", round(seurat_accuracy, 2)), 3, line=3)
alluvial(scmap_cluster[,1:2], freq=scmap_cluster$count)
mtext(paste0("scmap-cluster: ", round(scmap_cluster_accuracy, 2)), 3, line=3)
alluvial(scmap_cell[,1:2], freq=scmap_cell$count)
mtext(paste0("scmap-cell: ", round(scmap_cell_accuracy, 2)), 3, line=3)
dev.off()
