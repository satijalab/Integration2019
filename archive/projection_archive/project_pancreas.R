library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(ggalluvial)
library(dplyr)
library(methods)

celseq <- readRDS("seurat_objects/celseq.rds")
celseq2 <- readRDS("seurat_objects/celseq2.rds")
fluidigmc1 <- readRDS("seurat_objects/fluidigmc1.rds")
smartseq2 <- readRDS("seurat_objects/smartseq2.rds")
indrop <- readRDS("seurat_objects/inDrop.rds")

genes.use <- c(celseq@var.genes, celseq2@var.genes, smartseq2@var.genes,
               fluidigmc1@var.genes, indrop@var.genes)
genes.use <- names(which(table(genes.use) > 1))

#### Integration keeping out indrop1 ####
indrop <- SplitObject(indrop, attribute.1 = 'replicate', subset.raw = TRUE)
all.pancreas <- c(list(celseq, smartseq2, celseq2, fluidigmc1), indrop[2:4])

pancreas.integrated <- MultiIntegrateData(object.list = all.pancreas,
                                          var.genes = genes.use,
                                          order.by = "similarity",
                                          score.mnn = TRUE, k = 300,
                                          k.mnn = 10, score.mnn.k = 50,
                                          eps = 0, nn.dims.use = 1:20,
                                          do.pairwise = TRUE,
                                          reduction.use = "cca.cosine",
                                          dims.use = 1:20,
                                          nn.reduction.use = "pca",
                                          sd = 1/5, keep.sparse = TRUE,
                                          do.cpp = TRUE, verbose = TRUE)

pancreas.integrated <- FindVariableGenes(object = pancreas.integrated,
                                        do.plot = FALSE,
                                        display.progress = FALSE,
                                        selection.method = "dispersion")
pancreas.integrated <- ScaleData(object = pancreas.integrated,
                                genes.use = pancreas.integrated@var.genes,
                                display.progress = FALSE)
pancreas.integrated <- RunPCA(object = pancreas.integrated,
                             pc.genes = pancreas.integrated@var.genes,
                             do.print = FALSE, pcs.compute = 30)

saveRDS(pancreas.integrated, file = "seurat_objects/pancreas_integrated_no_indrop1.rds")
indrop1 <- indrop[[1]]

# transfer celltypes
celltypes <- readRDS("analysis_data/2018_06_22_integrated_pancreas_clusterIDs.rds")
pancreas.integrated <- AddMetaData(pancreas.integrated, celltypes, col.name = 'celltype')

runSeuratProject <- function(reference, query, label.name, celltypes, threshold=0.5){
  projection <- ProjectCells(reference = reference, query = query, label.name = label.name,
                             correct.query = TRUE, diet = TRUE, eps = 0,
                             do.cpp = TRUE, verbose = TRUE, sd = 1, k = 100,
                             use.cosine = TRUE, dims.use = 1:30, normalize.with.lm = TRUE)
  
  projection$true_label <- as.character(celltypes[rownames(projection)])
  all_projections <- projection
  projection[projection$prediction.score < threshold, 'predicted.id'] <- 'unassigned'
  
  seurat_projection <- projection %>% 
    group_by(true_label, predicted.id) %>%
    mutate(count = n()) %>%
    select(-prediction.score) %>% 
    arrange(-count) %>% 
    unique()
  
  seurat_accuracy <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100
  
  p1 <- ggplot(as.data.frame(seurat_projection), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
    geom_alluvium(aes(fill = true_label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
    theme_void() +
    theme(legend.position = 'none') +
    ggtitle(paste0("Seurat: ",round(seurat_accuracy,2),"%")) +
    theme(plot.title = element_text(size=20))
  return(p1)
}

runScmapCluster <- function(reference, query, label.name, celltypes, threshold=0.7){
  reference <- indexCluster(reference, cluster_col = label.name)
  scmapCluster_results <- scmapCluster(projection = query, index_list = list(predicted.id = metadata(reference)$scmap_cluster_index), threshold = threshold)
  scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
  scmap.cluster.results$predicted.id <- as.character(scmap.cluster.results$predicted.id)
  rownames(scmap.cluster.results) <- rownames(colData(query))
  scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])
  
  scmap_cluster <- scmap.cluster.results %>% 
    group_by(predicted.id, true_label) %>%
    mutate(count = n()) %>%
    arrange(-count) %>% 
    unique()
  
  scmap_cluster <- as.data.frame(scmap_cluster)
  scmap_cluster$predicted.id <- as.character(scmap_cluster$predicted.id)
  scmap_cluster$true_label <- as.character(scmap_cluster$true_label)
  
  scmap_cluster_accuracy <- sum(scmap.cluster.results$predicted.id == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100
  
  p1 <- ggplot(as.data.frame(scmap_cluster), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
    geom_alluvium(aes(fill = true_label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
    theme_void() +
    theme(legend.position = 'none') +
    ggtitle(paste0("scmap-cluster: ",round(scmap_cluster_accuracy,2),"%")) +
    theme(plot.title = element_text(size=20))
  return(p1)
}

runScmapCell <- function(reference, query, celltypes, threshold=0.5, w=3){
  reference <- indexCell(reference)
  scmapCell_results <- scmapCell(projection = query, index_list = list(predicted.id = metadata(reference)$scmap_cell_index))
  scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(reference)$celltype)), w = w, threshold = threshold)
  
  scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
  scmap.cell.results$predicted.id <- as.character(scmap.cell.results$predicted.id)
  rownames(scmap.cell.results) <- rownames(colData(indrop.sce))
  scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])
  
  scmap_cell <- scmap.cell.results %>% 
    group_by(predicted.id, true_label) %>%
    mutate(count = n()) %>%
    arrange(-count) %>% 
    unique()
  
  scmap_cell <- as.data.frame(scmap_cell)
  scmap_cell_accuracy <- sum(scmap.cell.results$predicted.id == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100
  
  p1 <- ggplot(as.data.frame(scmap_cell), aes(weight = count, axis1 = true_label, axis2 = predicted.id)) +
    geom_alluvium(aes(fill = true_label), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("True Label", "Predicted")) +
    theme_void() +
    theme(legend.position = 'none') +
    ggtitle(paste0("scmap-cell: ",round(scmap_cell_accuracy,2),"%")) +
    theme(plot.title = element_text(size=20))
  return(p1)
}

# run scmap
sce <- as.SingleCellExperiment(pancreas.integrated)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
assays(sce)$counts <- as.matrix(assays(sce)$counts)
assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)

indrop.sce <- as.SingleCellExperiment(indrop1)
rowData(indrop.sce)$feature_symbol <- rownames(indrop.sce)
assays(indrop.sce)$counts <- as.matrix(assays(indrop.sce)$counts)
assays(indrop.sce)$logcounts <- as.matrix(assays(indrop.sce)$logcounts)

seurat_plot <- runSeuratProject(reference = pancreas.integrated, query = indrop1, label.name = 'celltype', celltypes=celltypes, threshold = 0.5)
scmap_cluster_plot <- runScmapCluster(reference = sce, query = indrop.sce, label.name = 'celltype', threshold = 0.7, celltypes=celltypes)
scmap_cell_plot <- runScmapCell(reference = sce, query = indrop.sce, celltypes=celltypes, w = 3, threshold = 0.5)

cowplot::plot_grid(seurat_plot, scmap_cluster_plot, scmap_cell_plot, nrow=1) +
  ggsave("figures/pancreas_projection.pdf", height = 10, width = 15)

### Downsample query UMIs
subsample_indrop <- SampleUMI(indrop1@raw.data[,indrop1@cell.names], max.umi = 500)
dimnames(subsample_indrop) <- dimnames(indrop1@data)
subsample_indrop <- CreateSeuratObject(subsample_indrop, min.cells = 0, min.genes = 0, project = 'indrop')
subsample_indrop <- NormalizeData(subsample_indrop)

subsample_indrop.sce <- as.SingleCellExperiment(subsample_indrop)
rowData(subsample_indrop.sce)$feature_symbol <- rownames(subsample_indrop.sce)
assays(subsample_indrop.sce)$counts <- as.matrix(assays(subsample_indrop.sce)$counts)
assays(subsample_indrop.sce)$logcounts <- as.matrix(assays(subsample_indrop.sce)$logcounts)

seurat_plot <- runSeuratProject(reference = pancreas.integrated, query = subsample_indrop, label.name = 'celltype', celltypes=celltypes, threshold = 0.5)
scmap_cluster_plot <- runScmapCluster(reference = sce, query = subsample_indrop.sce, label.name = 'celltype', threshold = 0.7, celltypes=celltypes)
scmap_cell_plot <- runScmapCell(reference = sce, query = subsample_indrop.sce, celltypes=celltypes, w = 3, threshold = 0.5)

cowplot::plot_grid(seurat_plot, scmap_cluster_plot, scmap_cell_plot, nrow=1) +
  ggsave("figures/pancreas_projection_downsample.pdf")

##### Dropping thresholds for scmap
seurat_plot <- runSeuratProject(reference = pancreas.integrated, query = indrop1, label.name = 'celltype', celltypes=celltypes, threshold = 0)
scmap_cluster_plot <- runScmapCluster(reference = sce, query = indrop.sce, label.name = 'celltype', threshold = 0, celltypes=celltypes)
scmap_cell_plot <- runScmapCell(reference = sce, query = indrop.sce, celltypes=celltypes, w = 1, threshold = 0)

cowplot::plot_grid(seurat_plot, scmap_cluster_plot, scmap_cell_plot, nrow=1) +
  ggsave("figures/pancreas_projection_downsample_low_threshold.pdf")
