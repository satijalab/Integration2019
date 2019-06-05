library(Seurat)
library(Matrix)
library(future)
library(future.apply)

plan(strategy = 'multicore', workers = 6)
options(future.globals.maxSize= 10^12)

# Load data
mmx <- readMM(file = "raw_data/atacmca/activity_scores.quantitative.mtx")
cells <- readLines("raw_data/atacmca/activity_scores.quantitative.cells.txt")
genes <- readLines("raw_data/atacmca/activity_scores.quantitative.genes.txt")
meta <- read.table("raw_data/atacmca/knn_results.txt", sep = "\t")
up56 <- readRDS("seurat_objects/allen_p56_uppercase.rds")
rownames(mmx) <- genes
colnames(mmx) <- cells

# build atac object
cell_md <- read.table("raw_data/atacmca/cell_metadata.txt",
                      header = TRUE, row.names = 1, sep="\t")
atac <- CreateSeuratObject(counts = mmx)
atac <- AddMetaData(atac, metadata = cell_md)
Idents(atac) <- 'cell_label'
pfc <- subset(atac, tissue == "PreFrontalCortex")
pfc$dataset <- "ATAC"
pfc <- AddMetaData(pfc, meta)
pfc <- NormalizeData(pfc, scale.factor = median(pfc$nCount_RNA))

# Make binary object for DA
atac_pks <- readMM("raw_data/atacmca/atac_matrix.binary.qc_filtered.mtx")
atac_genes <- readLines("raw_data/atacmca/atac_matrix.binary.qc_filtered.peaks.txt")
atac_cells <- readLines("raw_data/atacmca/atac_matrix.binary.qc_filtered.cells.txt")
rownames(atac_pks) <- atac_genes
colnames(atac_pks) <- atac_cells

pfc_pks <- atac_pks[,Cells(pfc)]
apfc <- CreateSeuratObject(counts = pfc_pks)

# filter low quality nuclei
pfc$nATAC <- apfc$nCount_RNA[Cells(pfc)]
VlnPlot(pfc, 'nATAC', log = TRUE)
pfc_full <- pfc
pfc <- subset(pfc_full, nATAC > 3000)

# Preprocessing
pfc <- FindVariableFeatures(pfc, nfeatures = 3000)
up56 <- FindVariableFeatures(up56, nfeatures = 3000)

ob.list <- list(up56, pfc)
genes.use <- SelectIntegrationFeatures(ob.list, nfeatures = 3000)
ob.list <- future_sapply(ob.list, ScaleData, features = genes.use)

# just for visualization
ob1 <- RunCCA(
  object1 = ob.list[[1]],
  object2 = ob.list[[2]],
  features = genes.use,
  rescale = FALSE,
  num.cc = 30
)

ob1 <- CosineCCA(ob1)
ob1 <- RunUMAP(
  object = ob1,
  dims = 1:30,
  metric = 'correlation',
  reduction = "cca.cosine"
)

DimPlot(ob1, group.by = 'dataset', reduction = 'umap')
FeaturePlot(ob1, c('SST', 'PVALB'), split.by = 'dataset', reduction = 'umap', max.cutoff = 'q99')

# Projection
predictions <- ProjectCells(
  reference = ob.list[[1]],
  query = ob.list[[2]],
  label.name = "celltype1",
  dims = 1:30,
  reduction = 'cca',
  k.neighbors = 100,
  k.weights = 50,
  k.filter = NA,
  k.mnn = 5,
  score = FALSE,
  score2 = TRUE,
  sd = 1,
  use.cosine = TRUE,
  features = genes.use
)

# # Visualize predictions
# x <- ob1[, Cells(ob.list[[2]])]
# x <- AddMetaData(x, metadata = predictions)
# Idents(x) <- 'predicted.id'
# x1 <- subset(x, prediction.score.max > 0.75)
# DimPlot(x1, reduction = 'umap', label = TRUE)
# FeaturePlot(x1, 'SST')
# VlnPlot(x1, "SST")
# table(x1$predicted.id)

# Find DE Peaks
apfc <- AddMetaData(apfc, metadata = predictions)
apfc <- subset(apfc, prediction.score.max > 0.75)

# save cell classifications
for(i in unique(apfc$predicted.id)) {
  write.table(colnames(apfc[, apfc$predicted.id == i]), paste0("analysis_data/atac/clusters/cluster_", i, ".txt"),
              row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# CGE vs MGE
ident.1 <- c(4, 7)
ident.2 <- c(1, 2)

Idents(apfc) <- 'predicted.id'

Idents(apfc, cells = WhichCells(apfc, idents = ident.1)) <- "A"
Idents(apfc, cells = WhichCells(apfc, idents = ident.2)) <- "B"

ob <- subset(apfc, select = c("A","B"))
cts <- ob[["RNA"]]@data
cts_sum <- colSums(cts)
ob_a <- rowMeans(cts[,WhichCells(ob,idents = "A")])
ob_b <- rowMeans(cts[,WhichCells(ob,idents = "B")])
ob_diff <- sort(ob_a-ob_b)

# ntest <- 25000
# ltest <- c(head(names(ob_diff),ntest),tail(names(ob_diff),ntest))
ltest <- names(ob_diff)
groups <- Idents(ob)
data <- data.frame(groups)
data$nread <- cts_sum

DAtest <- function(i) {
  data$GENE <- as.numeric(cts[i,])
  m1 <- glm(GENE ~ groups + nread, family = 'binomial',data = data)
  m2 <- glm(GENE ~ nread, family = 'binomial',data = data)
  return(lmtest::lrtest(m1,m2)$Pr[2])
}

pvals <- future_sapply(ltest, DAtest)
results <- data.frame(ob_diff[ltest], pvals)
results <- results[order(results$pvals),]
write.table(results, 'analysis_data/atac/mge_cge.tsv')

colnames(results)[1] <- 'avg_FC'
results <- results[results$avg_FC > 0, ]
# results <- results[order(results$avg_FC, decreasing = TRUE), ]
pks <- rownames(head(results, 1000))
pk1 <- sapply((pks), Seurat:::ExtractField, 1)
pk2 <- sapply((pks), Seurat:::ExtractField, 2)
pk3 <- sapply((pks), Seurat:::ExtractField, 3)
writeLines(sapply(1:length(pk1), function(x) paste(pk1[x],pk2[x],pk3[x],sep='\t')),
           con = "analysis_data/atacmca/interneuron_peaks_mge_cge.bed")

# Pvalb vs Sst
ident.1 <- 7
ident.2 <- 4

Idents(apfc) <- 'predicted.id'

Idents(apfc, cells = WhichCells(apfc, idents = ident.1)) <- "A"
Idents(apfc, cells = WhichCells(apfc, idents = ident.2)) <- "B"

ob <- subset(apfc, select = c("A","B"))
cts <- ob[["RNA"]]@data
cts_sum <- colSums(cts)
ob_a <- rowMeans(cts[,WhichCells(ob,idents = "A")])
ob_b <- rowMeans(cts[,WhichCells(ob,idents = "B")])
ob_diff <- sort(ob_a-ob_b)
# ntest <- 25000
# ltest <- c(head(names(ob_diff),ntest),tail(names(ob_diff),ntest))
ltest <- names(ob_diff)
groups <- Idents(ob)
data <- data.frame(groups)
data$nread <- cts_sum

pvals <- future_sapply(ltest, DAtest)
results <- data.frame(ob_diff[ltest], pvals)
results <- results[order(results$pvals),]
write.table(results, 'analysis_data/atac/pv_sst.tsv')

colnames(results)[1] <- 'avg_FC'
results <- results[results$avg_FC > 0, ]
# results <- results[order(results$avg_FC, decreasing = TRUE), ]
pks <- rownames(head(results, 1000))
pk1 <- sapply((pks), Seurat:::ExtractField, 1)
pk2 <- sapply((pks), Seurat:::ExtractField, 2)
pk3 <- sapply((pks), Seurat:::ExtractField, 3)
writeLines(sapply(1:length(pk1), function(x) paste(pk1[x],pk2[x],pk3[x],sep='\t')),
           con = "analysis_data/atacmca/interneuron_peaks_pv_sst.bed")

# Vip vs Id2
ident.1 <- 1
ident.2 <- 2

Idents(apfc) <- 'predicted.id'

Idents(apfc, cells = WhichCells(apfc, idents = ident.1)) <- "A"
Idents(apfc, cells = WhichCells(apfc, idents = ident.2)) <- "B"

ob <- subset(apfc, select = c("A","B"))
cts <- ob[["RNA"]]@data
cts_sum <- colSums(cts)
ob_a <- rowMeans(cts[,WhichCells(ob,idents = "A")])
ob_b <- rowMeans(cts[,WhichCells(ob,idents = "B")])
ob_diff <- sort(ob_a-ob_b)
# ntest <- 25000
# ltest <- c(head(names(ob_diff),ntest),tail(names(ob_diff),ntest))
ltest <- names(ob_diff)
groups <- Idents(ob)
data <- data.frame(groups)
data$nread <- cts_sum

pvals <- future_sapply(ltest, DAtest)
results <- data.frame(ob_diff[ltest], pvals)
results <- results[order(results$pvals),]
write.table(results, 'analysis_data/atac/vip_id2.tsv')

colnames(results)[1] <- 'avg_FC'
results <- results[results$avg_FC > 0, ]
# results <- results[order(results$avg_FC, decreasing = TRUE), ]
pks <- rownames(head(results, 1000))
pk1 <- sapply((pks), Seurat:::ExtractField, 1)
pk2 <- sapply((pks), Seurat:::ExtractField, 2)
pk3 <- sapply((pks), Seurat:::ExtractField, 3)
writeLines(sapply(1:length(pk1), function(x) paste(pk1[x],pk2[x],pk3[x],sep='\t')),
           con = "analysis_data/atacmca/interneuron_peaks_vip_id2.bed")

# save full set of peaks tested
pk1 <- sapply((atac_genes), Seurat:::ExtractField, 1)
pk2 <- sapply((atac_genes), Seurat:::ExtractField, 2)
pk3 <- sapply((atac_genes), Seurat:::ExtractField, 3)
writeLines(sapply(1:length(pk1), function(x) paste(pk1[x],pk2[x],pk3[x],sep='\t')),
           con = "analysis_data/atacmca/all_peaks.bed")