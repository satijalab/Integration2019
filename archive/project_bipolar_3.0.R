suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scmap))
suppressMessages(library(ggalluvial))
suppressMessages(library(cowplot))
args <- commandArgs(trailingOnly = TRUE)

bipolar <- readRDS(file = args[1])

Idents(object = bipolar) <- "celltype"
bipolar <- SubsetData(object = bipolar, ident.remove = c("Unknown", "Doublets"))
# remove Unknown and doublets

Idents(object = bipolar) <- "orig.ident"
bipolar.reps <- SplitObject(object = bipolar, attribute.1 = 'orig.ident')

# remove one, integrate the rest, downsample removed, project
rep.to.remove <- as.numeric(args[2])

message(paste0("Processing with bipolar ", rep.to.remove, " removed."))
ref <- bipolar.reps[-rep.to.remove]
query <- unname(bipolar.reps[rep.to.remove])[[1]]

# gene selection
for(j in 1:length(x = ref)){
  ref[[j]] <- FindVariableFeatures(
    object = ref[[j]], 
    selection.method = "vst", 
    num.features = 2000, 
    verbose = FALSE)
}
genes.use <- SelectIntegrationFeatures(obj.list = ref, num.features = 2000)
for(j in 1:length(x = ref)) {
  ref[[j]] <- ScaleData(object = ref[[j]], features.use = genes.use, verbose = FALSE)
  ref[[j]] <- RunPCA(object = ref[[j ]], compute.dims = 40, features.use = genes.use, verbose = FALSE)
}
ref.integrated <- MultiIntegrateData(
  obj.list = ref, 
  assay.use = "RNA",
  reduction.use = "cca.cosine", 
  k.nn = 300,
  k.mnn = 5, 
  sd = 1, 
  do.cpp = TRUE, 
  verbose = TRUE, 
  dims.use = 1:30,
  features.use = genes.use
)
DefaultAssay(object = ref.integrated) <- "integrated"
ref.integrated <- FindVariableFeatures(object = ref.integrated, selection.method = "dispersion", num.features = 2000, verbose = FALSE)
ref.integrated <- ScaleData(object = ref.integrated, verbose = FALSE)
ref.integrated <- RunPCA(object = ref.integrated,
                             features.use = VariableFeatures(ref.integrated),
                             verbose = FALSE, compute.dims = 50)
saveRDS(object = ref.integrated, file = paste0("~/data/bipolar_projection/bipolar_norep", rep.to.remove, ".rds"))  

# downsample query
Idents(object = query) <- "celltype"
query <- SubsetData(object = query, max.cells.per.ident = 50)
saveRDS(object = query, file = paste0("~/data/bipolar_projection/bipolar_query_norep", rep.to.remove, ".rds"))

seurat.projection <- ProjectCells(
  reference = ref.integrated, query = query, label.name = 'celltype', eps = 0,
  do.cpp = TRUE, verbose = TRUE, sd = 1, k.neighbors = 300, k.mnn = 10, k.filter = 300, 
  k.weights = 300, use.cosine = TRUE, dims.use = 1:30, normalize.with.lm = TRUE
)

seurat.projection$true.label <- query["celltype", , drop = TRUE]
saveRDS(object = seurat.projection, file = paste0("~/data/bipolar_projection/seuratpred_norep", rep.to.remove, ".rds"))

# predict using scMAP
sce <- SingleCellExperiment(assays = list(logcounts = ref.integrated[["integrated"]]@data, counts = ref.integrated[["RNA"]]@counts))
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
scmapCluster.results <- scmapCluster(
  projection = sce.query, 
  index_list = list(yan = metadata(sce)$scmap_cluster_index), 
  threshold = -Inf  # force assignments 
)

# scmap-cell
sce <- indexCell(sce)
scmapCell_results <- scmapCell(projection = sce.query, index_list = list(yan = metadata(sce)$scmap_cell_index))
scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(sce)$celltype)), threshold = -Inf, w = 1)
scmap.results <- data.frame(scmap_cluster = scmapCluster.results$combined_labs, scmap_cell = scmapCell_clusters$combined_labs, true.label = colData(sce.query)$celltype)
rownames(scmap.results) <- colnames(x = sce.query)
saveRDS(object = scmap.results, file = paste0("~/data/bipolar_projection/scmappred_norep", rep.to.remove, ".rds"))
results <- list(seurat.projection = seurat.projection, scmap.results = scmap.results)
saveRDS(object = results, file = paste0("~/data/bipolar_projection/all_prediction_results_norep", rep.to.remove, ".rds"))

seurat.accuracy <- table(seurat.projection$predicted.id == seurat.projection$true.label)[2] / 
  sum(table(seurat.projection$predicted.id == seurat.projection$true.label) )
seurat.results <- as.data.frame(table(seurat.projection$predicted.id, seurat.projection$true.label))
colnames(seurat.results) <- c("predicted", "actual", "count")
p1 <- ggplot(seurat.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
  geom_alluvium(aes(fill = actual), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
  ggtitle(paste0("Seurat: ",round(seurat.accuracy,2),"%")) +
  theme_cowplot() +
  NoAxes() +
  theme(legend.position = 'none')
ggsave(p1, filename = paste0("figures/bipolar_seurat_projection_norep", rep.to.remove , ".pdf"), width = 10, height = 8)


scmapcell.results <- as.data.frame(table(scmap.results$scmap_cell, scmap.results$true.label))
scmapcell.accuracy <- table(as.vector(scmap.results$scmap_cell) == scmap.results$true.label)[2] / 
  sum(table(as.vector(scmap.results$scmap_cell) == scmap.results$true.label))
colnames(scmapcell.results) <- c("predicted", "actual", "count")
p2 <- ggplot(scmapcell.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
  geom_alluvium(aes(fill = actual), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
  ggtitle(paste0("scmap_cell: ",round(scmapcell.accuracy,2),"%")) +
  theme_cowplot() +
  NoAxes() +
  theme(legend.position = 'none')
ggsave(p2, filename = paste0("figures/bipolar_scmapcell_projection_norep", rep.to.remove , ".pdf"), width = 10, height = 8)


scmapcluster.results <- as.data.frame(table(scmap.results$scmap_cluster, scmap.results$true.label))
scmapcluster.accuracy <- table(as.vector(scmap.results$scmap_cluster) == scmap.results$true.label)[2] / 
  sum(table(as.vector(scmap.results$scmap_cluster) == scmap.results$true.label))
colnames(scmapcluster.results) <- c("predicted", "actual", "count")
p3 <- ggplot(scmapcluster.results, aes(y = count, axis1 = actual, axis2 = predicted)) +
  geom_alluvium(aes(fill = actual), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = T,label.size = 0.1)+ 
  scale_x_continuous(breaks = 1:2, labels = c("Predicted", "Actual")) +
  ggtitle(paste0("scmap_cluster: ",round(scmapcluster.accuracy, 2),"%")) +
  theme_cowplot() +
  NoAxes() +
  theme(legend.position = 'none')
ggsave(p3, filename = paste0("figures/bipolar_scmapcluster_projection_norep", rep.to.remove , ".pdf"), width = 10, height = 8)
