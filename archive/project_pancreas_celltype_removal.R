suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scmap))
suppressMessages(library(ggalluvial))
suppressMessages(library(cowplot))
args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/celseq.rds", "~/Projects/muir/seurat_objects/celseq2.rds",
#           "~/Projects/muir/seurat_objects/smartseq2.rds", "~/Projects/muir/seurat_objects/fluidigmc1.rds",
#           "~/Projects/muir/seurat_objects/inDrop.rds", "~/Projects/muir/seurat_objects/integrated_pancreas.rds",
#           "celseq", "quiescent_stellate")

celseq <- readRDS(file = args[1])
celseq2 <- readRDS(file = args[2])
smartseq2 <- readRDS(file = args[3])
fluidigmc1 <- readRDS(file = args[4])
indrop <- readRDS(file = args[5])
pancreas.integrated <- readRDS(file = args[6])
cell.type.removed <- args[7]
dataset.names <- c("celseq", "celseq2", "smartseq2", "fluidigmc1", "indrop1", "indrop2", "indrop3", "indrop4")

ref.master <- c(celseq, celseq2, smartseq2, fluidigmc1, indrop)

for(i in 1:length(ref.master)){
  ref.master[[i]]["celltype"] <- Idents(object = pancreas.integrated)
  Idents(ref.master[[i]]) <- "celltype"
}

# leave one dataset out of reference, remove given cell type from all but query

accuracy_results <- data.frame(query.dataset = character(), celltype.removed = character(), method = character(), accuracy = numeric())

for(i in 1:length(ref.master)) {
  ref <- ref.master
  query <- ref.master[[i]]
  if (table(query["celltype"])[cell.type.removed] < 10  || is.na(table(query["celltype"])[cell.type.removed] )) next
  ref[[i]] <- NULL
  for(j in 1:length(ref)) {
    if(cell.type.removed %in% unique(ref[[j]]["celltype"])$celltype) {
      ref[[j]] <- SubsetData(object = ref[[j]], ident.remove = cell.type.removed)
    }
  }
  # gene selection
  for(j in 1:length(x = ref)){
    ref[[j]] <- FindVariableFeatures(
      object = ref[[j]], 
      selection.method = "vst", 
      nfeatures = 2000, 
      verbose = FALSE)
  }
  genes.use <- SelectIntegrationFeatures(obj.list = ref, nfeatures = 2000)
  for(j in 1:length(x = ref)) {
    ref[[j]] <- ScaleData(object = ref[[j]], features = genes.use, verbose = FALSE)
    ref[[j]] <- RunPCA(object = ref[[j ]], npcs = 40, features = genes.use, verbose = FALSE)
  }
  # integrate the reference datasets
  ref.integrated <- MultiIntegrateData(
    obj.list = ref, 
    assay = "RNA",
    reduction  = "cca.cosine", 
    k.nn = 300,
    k.mnn = 10, 
    sd = 1, 
    do.cpp = TRUE, 
    verbose = TRUE, 
    dims = 1:30,
    features = genes.use,
    features.integrate = genes.use
  )
  DefaultAssay(object = ref.integrated) <- "integrated"
  VariableFeatures(ref.integrated) <- genes.use
  ref.integrated <- ScaleData(object = ref.integrated, verbose = FALSE)
  ref.integrated <- RunPCA(object = ref.integrated,
                           features.use = VariableFeatures(ref.integrated),
                           verbose = FALSE, npcs = 50)
  
  # downsample query to 100 cells (max) per ident with removed pop being 20%
  Idents(object = query) <- "celltype"
  query1 <- SubsetData(object = query, max.cells.per.ident = 100, ident.remove = cell.type.removed)
  newCells <- round(ncol(query1)/4)
  query2 <- SubsetData(query, ident.use = cell.type.removed, max.cells.per.ident = newCells)
  query <- SubsetData(query, cells = c(colnames(query1), colnames(query2)))
  
  # perform projection 
  seurat.projection <- ProjectCells(
    reference = ref.integrated, query = query, label.name = 'celltype', eps = 0,
    do.cpp = TRUE, verbose = TRUE, sd = 1, k.neighbors = 100, k.mnn = 10, k.filter = 50, 
    k.weights = 50, use.cosine = TRUE, dims = 1:30, normalize.with.lm = TRUE
  )
  
  seurat.projection$true.label <- query["celltype", , drop = TRUE]
  seurat.projection <- seurat.projection[order(seurat.projection$prediction.score), ]
  seurat.projection$predicted.id.unassigned <- seurat.projection$predicted.id
  seurat.projection[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id"] <- cell.type.removed
  seurat.projection[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id.unassigned"] <- "Unassigned"
  
  seurat.accuracy <- table(seurat.projection$predicted.id == seurat.projection$true.label)[2] / 
    sum(table(seurat.projection$predicted.id == seurat.projection$true.label) )
  seurat.results <- as.data.frame(table(seurat.projection$predicted.id.unassigned, seurat.projection$true.label))
  colnames(seurat.results) <- c("predicted", "actual", "count")
  p1 <- ggplot(seurat.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
    geom_alluvium(aes(fill = actual), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
    ggtitle(paste0("Seurat: ",round(seurat.accuracy,2) * 100,"%")) +
    theme_cowplot() +
    NoAxes() +
    theme(legend.position = 'none')
  ggsave(p1, filename = paste0("figures/pancreas_seurat_projection_", dataset.names[i], "_", cell.type.removed, "_removed.pdf"), width = 10, height = 8)
  
  # predict using scMAP
  sce <- SingleCellExperiment(assays = list(logcounts = ref.integrated[["RNA"]]@data, counts = ref.integrated[["RNA"]]@counts))
  colData(sce) <- S4Vectors::DataFrame(ref.integrated[])
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  assays(sce)$counts <- as.matrix(assays(sce)$counts)
  assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
  sce <- selectFeatures(sce, n_features = 500, suppress_plot = TRUE)
  
  sce.query <- SingleCellExperiment(assays = list(logcounts = query[["RNA"]]@data, counts = query[["RNA"]]@counts))
  colData(sce.query) <- S4Vectors::DataFrame(query[])
  rowData(sce.query)$feature_symbol <- rownames(sce.query)
  assays(sce.query)$counts <- as.matrix(assays(sce.query)$counts)
  assays(sce.query)$logcounts <- as.matrix(assays(sce.query)$logcounts)
  
  # scmap-cluster
  sce <- indexCluster(sce, cluster_col = 'celltype')
  scmapcluster.results <- scmapCluster(
    projection = sce.query, 
    index_list = list(yan = metadata(sce)$scmap_cluster_index), 
    threshold = -Inf  # force assignments 
  )
  
  # scmap-cell
  sce <- indexCell(sce)
  scmapcell.results <- scmapCell(projection = sce.query, index_list = list(yan = metadata(sce)$scmap_cell_index))
  scmapcell.clusters <- scmapCell2Cluster(scmapcell.results, list(as.character(colData(sce)$celltype)), threshold = -Inf, w = 1)
  scmapcell.results <- data.frame(predicted.id = as.vector(scmapcell.clusters$combined_labs), siml = unname(scmapcell.clusters$scmap_cluster_siml), 
                                  true.label = colData(sce.query)$celltype, stringsAsFactors = FALSE, row.names = rownames(colData(sce.query)))
  scmapcell.results$predicted.id.unassigned <- scmapcell.results$predicted.id
  scmapcell.results$actual.prediction <- scmapcell.results$predicted.id
  scmapcell.results <- scmapcell.results[order(scmapcell.results$siml), ]
  scmapcell.results[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id"] <- cell.type.removed
  scmapcell.results[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id.unassigned"] <- "Unassigned"
  
  scmapcluster.results <- data.frame(predicted.id = as.vector(scmapcluster.results$combined_labs), siml = unname(scmapcluster.results$scmap_cluster_siml), 
                                     true.label = colData(sce.query)$celltype, stringsAsFactors = FALSE, row.names = rownames(colData(sce.query)))
  scmapcluster.results$predicted.id.unassigned <- scmapcluster.results$predicted.id
  scmapcluster.results$actual.prediction <- scmapcluster.results$predicted.id
  scmapcluster.results <- scmapcluster.results[order(scmapcluster.results$siml), ]
  scmapcluster.results[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id"] <- cell.type.removed
  scmapcluster.results[1:(length(which(query["celltype"] == cell.type.removed))), "predicted.id.unassigned"] <- "Unassigned"
  
  scmapcell.accuracy <- table(scmapcell.results$predicted.id == scmapcell.results$true.label)[2] / 
    sum(table(scmapcell.results$predicted.id == scmapcell.results$true.label))
  scmapcell.results.plot <- as.data.frame(table(scmapcell.results$predicted.id.unassigned, scmapcell.results$true.label))
  colnames(scmapcell.results.plot) <- c("predicted", "actual", "count")
  p2 <- ggplot(scmapcell.results.plot, aes(y = count, axis1 = actual, axis2 = predicted)) +
    geom_alluvium(aes(fill = actual), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
    ggtitle(paste0("scmap_cell: ",round(scmapcell.accuracy,2) * 100,"%")) +
    theme_cowplot() +
    NoAxes() +
    theme(legend.position = 'none')
  ggsave(p2, filename = paste0("figures/pancreas_scmapcell_projection_", dataset.names[i], "_", cell.type.removed, "_removed.pdf"), width = 10, height = 8)
  
  scmapcluster.accuracy <- table(scmapcluster.results$predicted.id == scmapcluster.results$true.label)[2] / 
    sum(table(scmapcluster.results$predicted.id == scmapcluster.results$true.label))
  scmapcluster.results.plot <- as.data.frame(table(scmapcluster.results$predicted.id.unassigned, scmapcluster.results$true.label))
  colnames(scmapcluster.results.plot) <- c("predicted", "actual", "count")
  p3 <- ggplot(scmapcluster.results.plot, aes(y = count, axis1 = actual, axis2 = predicted)) +
    geom_alluvium(aes(fill = actual), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
    scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
    ggtitle(paste0("scmap_cluster: ",round(scmapcluster.accuracy,2) * 100,"%")) +
    theme_cowplot() +
    NoAxes() +
    theme(legend.position = 'none')
  ggsave(p3, filename = paste0("figures/pancreas_scmapcluster_projection_", dataset.names[i], "_", cell.type.removed, "_removed.pdf"), width = 10, height = 8)
  
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = dataset.names[i], 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "seurat", 
                                                         accuracy = seurat.accuracy))
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = dataset.names[i], 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "scmap_cell", 
                                                         accuracy = scmapcell.accuracy))
  accuracy_results <- rbind(accuracy_results, data.frame(query.dataset = dataset.names[i], 
                                                         celltype.removed = cell.type.removed, 
                                                         method = "scmap_cluster", 
                                                         accuracy = scmapcluster.accuracy))
  save(ref.integrated, query, seurat.projection, scmapcluster.results, scmapcell.results, accuracy_results,
       file = paste0("~/data/pancreas_projections/intermediate_data_", dataset.names[i], "_", cell.type.removed, "_removed.Rda"))
}

saveRDS(accuracy_results, file = paste0("analysis_data/", cell.type.removed, "_accuracy.rds"))

