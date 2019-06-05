library(SingleCellExperiment)
library(scmap)
library(Seurat)
library(alluvial)
library(dplyr)
library(methods)

downsampleObject <- function(object, numi){
  subsample_matrix <- SampleUMI(object@raw.data[,object@cell.names], max.umi = numi)
  dimnames(subsample_matrix) <- dimnames(object@data)
  obj <- CreateSeuratObject(subsample_matrix, min.cells = 0, min.genes = 0, project = 'projection')
  obj <- NormalizeData(obj)
  return(obj)
}

findAccuracySeurat <- function(reference, query, ground_truth, label.name = 'celltype'){
  projection <- ProjectCells(reference = reference, query = query, label.name = label.name,
                             correct.query = TRUE, diet = TRUE, eps = 0,
                             do.cpp = TRUE, verbose = FALSE, sd = 1, k = 100,
                             use.cosine = TRUE, dims.use = 1:30, normalize.with.lm = TRUE)
  projection$true_label <- as.character(ground_truth[rownames(projection)])
  accuracy_full <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100
  projection[projection$prediction.score < 0.5, 'predicted.id'] <- 'unassigned'
  accuracy_threshold <- sum(projection$predicted.id == projection$true_label) / nrow(projection) * 100
  return(list(accuracy_full, accuracy_threshold))
}

findAccuracySCMAP <- function(reference, query, ground_truth){
  # cluster
  reference <- indexCluster(reference, cluster_col = 'celltype')
  scmapCluster_results <- scmapCluster(projection = query, index_list = list(yan = metadata(reference)$scmap_cluster_index))
  scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
  scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
  rownames(scmap.cluster.results) <- rownames(colData(query))
  scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])
  scmap_cluster_accuracy <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100
  scmapCluster_results <- scmapCluster(projection = query, index_list = list(yan = metadata(reference)$scmap_cluster_index), threshold = 0)
  scmap.cluster.results <- as.data.frame(scmapCluster_results$scmap_cluster_labs)
  scmap.cluster.results$yan <- as.character(scmap.cluster.results$yan)
  rownames(scmap.cluster.results) <- rownames(colData(query))
  scmap.cluster.results$true_label <- as.character(celltypes[rownames(scmap.cluster.results)])
  scmap_cluster_accuracy_no_threshold <- sum(scmap.cluster.results$yan == scmap.cluster.results$true_label) / nrow(scmap.cluster.results) * 100
  
  # cell
  reference <- indexCell(reference)
  scmapCell_results <- scmapCell(projection = query, index_list = list(yan = metadata(reference)$scmap_cell_index))
  scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(reference)$celltype)))
  scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
  scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
  rownames(scmap.cell.results) <- rownames(colData(query))
  scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])
  scmap_cell_accuracy <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100
  scmapCell_results <- scmapCell(projection = query, index_list = list(yan = metadata(reference)$scmap_cell_index))
  scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(reference)$celltype)), w = 1, threshold = 0)
  scmap.cell.results <- as.data.frame(scmapCell_clusters$scmap_cluster_labs)
  scmap.cell.results$yan <- as.character(scmap.cell.results$yan)
  rownames(scmap.cell.results) <- rownames(colData(query))
  scmap.cell.results$true_label <- as.character(celltypes[rownames(scmap.cell.results)])
  scmap_cell_accuracy_no_threshold <- sum(scmap.cell.results$yan == scmap.cell.results$true_label) / nrow(scmap.cell.results) * 100
  
  return(list(scmap_cluster_accuracy_no_threshold, scmap_cluster_accuracy, scmap_cell_accuracy_no_threshold, scmap_cell_accuracy))
}

runDownsampling <- function(numi, sce_reference, reference, query, celltypes){
  dsobj <- downsampleObject(query, numi)
  scmap_ds <- as.SingleCellExperiment(dsobj)
  rowData(scmap_ds)$feature_symbol <- rownames(scmap_ds)
  assays(scmap_ds)$counts <- as.matrix(assays(scmap_ds)$counts)
  assays(scmap_ds)$logcounts <- as.matrix(assays(scmap_ds)$logcounts)
  seurat_accuracy <- findAccuracySeurat(reference, dsobj, celltypes, label.name = 'celltype')
  scmap_accuracy <- findAccuracySCMAP(sce_reference, scmap_ds, celltypes)
  return(list('umi' = numi, 'seurat_full' = seurat_accuracy[[1]], 'seurat_threshold' = seurat_accuracy[[2]],
              'scmap_cluster_full' = scmap_accuracy[[1]], 'scmap_cluster_threshold' = scmap_accuracy[[2]],
              'scmap_cell_full' = scmap_accuracy[[3]], 'scmap_cell_threshold' = scmap_accuracy[[4]]))
}