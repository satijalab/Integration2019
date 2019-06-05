library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(dplyr)
library(methods)
library(ggalluvial)

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

bipolar6.sce <- as.SingleCellExperiment(bipolar6)
rowData(bipolar6.sce)$feature_symbol <- rownames(bipolar6.sce)
assays(bipolar6.sce)$counts <- as.matrix(assays(bipolar6.sce)$counts)
assays(bipolar6.sce)$logcounts <- as.matrix(assays(bipolar6.sce)$logcounts)

# Seurat
projection <- ProjectCells(reference = bipolar, query = bipolar6, label.name = 'celltype',
                           correct.query = TRUE, diet = TRUE, eps = 1,
                           do.cpp = TRUE, verbose = TRUE, sd = 1, k = 50,
                           use.cosine = TRUE, dims.use = 1:20, normalize.with.lm = TRUE)

projection$true_label <- as.character(sampled_celltypes[rownames(projection)])
projection_full <- projection
accuracy_full <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100
projection[projection$prediction.score < 0.5, 'predicted.id'] <- 'unassigned'
accuracy_threshold <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100

seurat_projection <- projection %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_projection_full <- projection_full %>% 
  group_by(true_label, predicted.id) %>%
  mutate(count = n()) %>%
  select(-prediction.score) %>% 
  arrange(-count) %>% 
  unique()

seurat_threshold <- ggplot(as.data.frame(seurat_projection), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("Seurat: ",round(accuracy_threshold,2),"%")) +
  theme(plot.title = element_text(size=20))

seurat_full <- ggplot(as.data.frame(seurat_projection_full), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("Seurat: ",round(accuracy_full,2),"%")) +
  theme(plot.title = element_text(size=20))

# Scmap
# cluster
sce <- indexCluster(sce, cluster_col = 'celltype')
scmapCluster_results <- scmapCluster(projection = bipolar6.sce, index_list = list(predicted.id = metadata(sce)$scmap_cluster_index))
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$predicted.id <- as.character(scmap.cluster.results$predicted.id)
rownames(scmap.cluster.results) <- rownames(colData(bipolar6.sce))
scmap.cluster.results$true_label <- as.character(sampled_celltypes[rownames(scmap.cluster.results)])
scmap_cluster_accuracy_threshold <- sum(scmap.cluster.results$predicted.id == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap_cluster_threshold <- scmap.cluster.results %>% 
  group_by(predicted.id, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster_threshold_plot <- ggplot(as.data.frame(scmap_cluster_threshold), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("scmap-cluster: ",round(scmap_cluster_accuracy_threshold,2),"%")) +
  theme(plot.title = element_text(size=20))

scmapCluster_results <- scmapCluster(projection = bipolar6.sce, index_list = list(predicted.id = metadata(sce)$scmap_cluster_index), threshold = 0)
scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
scmap.cluster.results$predicted.id <- as.character(scmap.cluster.results$predicted.id)
rownames(scmap.cluster.results) <- rownames(colData(bipolar6.sce))
scmap.cluster.results$true_label <- as.character(sampled_celltypes[rownames(scmap.cluster.results)])
scmap_cluster_accuracy_no_threshold <- sum(scmap.cluster.results$predicted.id == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100

scmap_cluster_full <- scmap.cluster.results %>% 
  group_by(predicted.id, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cluster_full_plot <- ggplot(as.data.frame(scmap_cluster_full), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("scmap-cluster: ",round(scmap_cluster_accuracy_no_threshold,2),"%")) +
  theme(plot.title = element_text(size=20))

# cell
sce <- indexCell(sce)
scmapCell_results <- scmapCell(projection = bipolar6.sce, index_list = list(predicted.id = metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)))
scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$predicted.id <- as.character(scmap.cell.results$predicted.id)
rownames(scmap.cell.results) <- rownames(colData(bipolar6.sce))
scmap.cell.results$true_label <- as.character(sampled_celltypes[rownames(scmap.cell.results)])
scmap_cell_accuracy <- sum(scmap.cell.results$predicted.id == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

scmap_cell_threshold <- scmap.cell.results %>% 
  group_by(predicted.id, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()dshudx23


scmap_cell_threshold_plot <- ggplot(as.data.frame(scmap_cell_threshold), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("scmap-cell: ",round(scmap_cell_accuracy,2),"%")) +
  theme(plot.title = element_text(size=20))

scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)), w = 1, threshold = 0)
scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
scmap.cell.results$predicted.id <- as.character(scmap.cell.results$predicted.id)
rownames(scmap.cell.results) <- rownames(colData(bipolar6.sce))
scmap.cell.results$true_label <- as.character(sampled_celltypes[rownames(scmap.cell.results)])
scmap_cell_accuracy_no_threshold <- sum(scmap.cell.results$predicted.id == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100

scmap_cell_no_threshold <- scmap.cell.results %>% 
  group_by(predicted.id, true_label) %>%
  mutate(count = n()) %>%
  arrange(-count) %>% 
  unique()

scmap_cell_full_plot <- ggplot(as.data.frame(scmap_cell_no_threshold), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
  geom_alluvium(aes(fill = true_label), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
  theme_void() +
  theme(legend.position = 'none') +
  ggtitle(paste0("scmap-cell: ",round(scmap_cell_accuracy_no_threshold,2),"%")) +
  theme(plot.title = element_text(size=20))

plot_grid(seurat_threshold, scmap_cluster_threshold_plot, scmap_cell_threshold_plot, nrow=1) + ggsave("figures/bipolar_projection_thresholds.pdf")
plot_grid(seurat_full, scmap_cluster_full_plot, scmap_cell_full_plot, nrow=1) + ggsave("figures/bipolar_projection_full.pdf")
