suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scmap))
suppressMessages(library(ggalluvial))
suppressMessages(library(cowplot))
args <- commandArgs(trailingOnly = TRUE)

celseq <- readRDS(file = args[1])
celseq2 <- readRDS(file = args[2])
smartseq2 <- readRDS(file = args[3])
fluidigmc1 <- readRDS(file = args[4])
indrop <- readRDS(file = args[5])
dataset.name <- args[7]
dataset.names <- c("celseq", "celseq2", "smartseq2", "fluidigmc1", "indrop1", "indrop2", "indrop3", "indrop4")

# set up reference set and query
ref <- c(celseq, celseq2, smartseq2, fluidigmc1, indrop)
query <- ref[[which(dataset.name == dataset.names)]]
ref[[which(dataset.name == dataset.names)]] <- NULL

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
  k.mnn = 5, 
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
pancreas.integrated <- readRDS(file = args[6])

# add cell type labels from the integrated pancreas object
query["celltype"] <- Idents(object = pancreas.integrated)
ref.integrated["celltype"] <- Idents(object = pancreas.integrated)

# downsample query to 100 cells (max) per ident
Idents(object = query) <- "celltype"
query <- SubsetData(object = query, max.cells.per.ident = 100)

# perform projection 
seurat.projection <- ProjectCells(
  reference = ref.integrated, query = query, label.name = 'celltype', eps = 0,
  do.cpp = TRUE, verbose = TRUE, sd = 1, k.neighbors = 100, k.mnn = 10, k.filter = 50, 
  k.weights = 50, use.cosine = TRUE, dims = 1:30, normalize.with.lm = TRUE
)
seurat.projection$true.label <- query["celltype", , drop = TRUE]

# predict using scMAP
sce <- SingleCellExperiment(assays = list(logcounts = ref.integrated[["integrated"]]@data, counts = ref.integrated[["RNA"]]@counts[rownames(ref.integrated),]))
colData(sce) <- S4Vectors::DataFrame(ref.integrated[])
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
assays(sce)$counts <- as.matrix(assays(sce)$counts)
assays(sce)$logcounts <- as.matrix(assays(sce)$logcounts)
# not sure if this is strictly correct (using corrected and uncorrected for model building to select features)
sce <- selectFeatures(sce, n_features = 500, suppress_plot = FALSE)

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
scmapCell_results <- scmapCell(projection = sce.query, index_list = list(yan = metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)), threshold = -Inf, w = 1)
scmap.results <- data.frame(scmap_cluster = scmapcluster.results$combined_labs, scmap_cell = scmapCell_clusters$combined_labs, true.label = colData(sce.query)$celltype)
rownames(scmap.results) <- colnames(x = sce.query)
results <- list(seurat.projection = seurat.projection, scmap.results = scmap.results)

# Make the alluvial plots
seurat.accuracy <- table(seurat.projection$predicted.id == seurat.projection$true.label)[2] / 
  sum(table(seurat.projection$predicted.id == seurat.projection$true.label) )
seurat.results <- as.data.frame(table(seurat.projection$predicted.id, seurat.projection$true.label))
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
ggsave(p1, filename = paste0("figures/pancreas_seurat_projection_", dataset.name , ".pdf"), width = 10, height = 8)

scmapcell.results <- as.data.frame(table(scmap.results$scmap_cell, scmap.results$true.label))
scmapcell.accuracy <- table(as.vector(scmap.results$scmap_cell) == scmap.results$true.label)[2] / 
  sum(table(as.vector(scmap.results$scmap_cell) == scmap.results$true.label))
colnames(scmapcell.results) <- c("predicted", "actual", "count")
p2 <- ggplot(scmapcell.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
  geom_alluvium(aes(fill = actual), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
  ggtitle(paste0("scmap_cell: ",round(scmapcell.accuracy,2)*100,"%")) +
  theme_cowplot() +
  NoAxes() +
  theme(legend.position = 'none')
ggsave(p2, filename = paste0("figures/pancreas_scmapcell_projection_", dataset.name , ".pdf"), width = 10, height = 8)

scmapcluster.results <- as.data.frame(table(scmap.results$scmap_cluster, scmap.results$true.label))
scmapcluster.accuracy <- table(as.vector(scmap.results$scmap_cluster) == scmap.results$true.label)[2] / 
  sum(table(as.vector(scmap.results$scmap_cluster) == scmap.results$true.label))
colnames(scmapcluster.results) <- c("predicted", "actual", "count")
p3 <- ggplot(scmapcluster.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
  geom_alluvium(aes(fill = actual), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
  ggtitle(paste0("scmap_cluster: ",round(scmapcluster.accuracy, 2) * 100,"%")) +
  theme_cowplot() +
  NoAxes() +
  theme(legend.position = 'none')
ggsave(p3, filename = paste0("figures/pancreas_scmapcluster_projection_", dataset.name , ".pdf"), width = 10, height = 8)
