#library(Seurat)
library(Matrix)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
devtools::load_all(args[8])

starmap.1 <- readRDS("seurat_objects/20180505_BY3_1kgenes.rds")
starmap.2 <- readRDS("seurat_objects/20180410-BY3_1kgenes.rds")
starmap.2 <- RenameCells(starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)
allen <- readRDS("seurat_objects/allen_brain.rds")
dropseq <- readRDS("seurat_objects/dropseq_cortex.rds")

dir.create("figures/spatial/")

# emove HPC from starmap
class_labels.1 <- read.table(
  file = "raw_data/spatial/starmap/visual_1020/20180505_BY3_1kgenes/class_labels.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
class_labels.2 <- read.table(
  file = "raw_data/spatial/starmap/visual_1020/20180410-BY3_1kgenes/class_labels.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

class_labels.1$cellname <- paste0('starmap', rownames(class_labels.1))
class_labels.2$cellname <- paste0('r2_starmap', rownames(class_labels.2))

class_labels.1$ClusterName <- ifelse(is.na(class_labels.1$ClusterName), 'Other', class_labels.1$ClusterName)
class_labels.2$ClusterName <- ifelse(is.na(class_labels.2$ClusterName), 'Other', class_labels.2$ClusterName)

hpc.1 <- class_labels.1[class_labels.1$ClusterName == 'HPC', ]$cellname
hpc.2 <- class_labels.2[class_labels.2$ClusterName == 'HPC', ]$cellname

accept.cells.1 <- setdiff(colnames(starmap.1), hpc.1)
accept.cells.2 <- setdiff(colnames(starmap.2), hpc.2)

starmap.1 <- starmap.1[, accept.cells.1]
starmap.2 <- starmap.2[, accept.cells.2]

starmap.1@misc$spatial <- starmap.1@misc$spatial[starmap.1@misc$spatial$cell %in% accept.cells.1, ]
starmap.2@misc$spatial <- starmap.2@misc$spatial[starmap.2@misc$spatial$cell %in% accept.cells.2, ]

#Integrate the STARmap datasets together
i1 <- FindIntegrationAnchors(
  object.list = c(starmap.1, starmap.2),
  anchor.features = rownames(starmap.1),
  dims = 1:30
)

star <- IntegrateData(anchorset = i1)

DefaultAssay(star) <- 'integrated'
star <- ScaleData(star)
star <- RunPCA(star, npcs = 30)
star <- RunUMAP(star, dims = 1:30)

#Project on allen labels
i2 <- FindTransferAnchors(
  reference = allen,
  query = star,
  features = rownames(starmap.1),
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'integrated'
)

predictions <- TransferData(
  anchorset = i2,
  refdata = allen$subclass,
  weight.reduction = 'pca'
)

star <- AddMetaData(star, metadata = predictions)
starmap.1 <- AddMetaData(starmap.1, metadata = predictions)
starmap.2 <- AddMetaData(starmap.2, metadata = predictions)

refdata <- GetAssayData(
  object = allen,
  assay = 'RNA',
  slot = 'data'
)

imputation <- TransferData(
  anchorset = i2,
  refdata = refdata,
  weight.reduction = 'pca'
)

star[['ss2']] <- imputation
starmap.1[["ss2"]] <- star[, colnames(starmap.1)][['ss2']]
starmap.2[["ss2"]] <- star[, colnames(starmap.2)][['ss2']]

# Transfer from dropseq
i2 <- FindTransferAnchors(
  reference = dropseq,
  query = star,
  features = rownames(starmap.1),
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'integrated'
)

refdata <- GetAssayData(
  object = dropseq,
  assay = 'RNA',
  slot = 'data'
)

imputation <- TransferData(
  anchorset = i2,
  refdata = refdata,
  weight.reduction = 'pca'
)

star[['dropseq']] <- imputation
starmap.1[["dropseq"]] <- star[, colnames(starmap.1)][['dropseq']]
starmap.2[["dropseq"]] <- star[, colnames(starmap.2)][['dropseq']]

write.table(predictions, 'analysis_data/starmap_celltype_predictions.tsv', sep = '\t', quote = FALSE)
saveRDS(starmap.1, 'seurat_objects/20180505_BY3_1kgenes_imputed.rds')
saveRDS(starmap.2, 'seurat_objects/20180410-BY3_1kgenes_imputed.rds')
saveRDS(star, 'seurat_objects/integrated_starmap.rds')
write.table(GetAssayData(starmap.1, assay = 'ss2', slot = 'data'), file = 'analysis_data/starmap_r1_imputed_expression_ss2.tsv', sep = '\t', quote = FALSE)
write.table(GetAssayData(starmap.2, assay = 'ss2', slot = 'data'), file = 'analysis_data/starmap_r2_imputed_expression_ss2.tsv', sep = '\t', quote = FALSE)
write.table(GetAssayData(starmap.1, assay = 'dropseq', slot = 'data'), file = 'analysis_data/starmap_r1_imputed_expression_dropseq.tsv', sep = '\t', quote = FALSE)
write.table(GetAssayData(starmap.2, assay = 'dropseq', slot = 'data'), file = 'analysis_data/starmap_r2_imputed_expression_dropseq.tsv', sep = '\t', quote = FALSE)
prediction.max <- predictions[, c('predicted.id', 'prediction.score.max')]
colnames(prediction.max) <- c("predicted.id", 'prediction.score')
write.table(prediction.max, 'analysis_data/starmap_celltype_predictions_max.tsv', sep = '\t', quote = FALSE)

# coembed starmap
DefaultAssay(star) <- 'ss2'
genes.use <- VariableFeatures(allen)

star$technology <- 'STARmap'
allen$technology <- 'SMART-seq2'

genes.use <- VariableFeatures(dropseq)

coembed.assay <- merge(
  x = subset(x = allen[["RNA"]], features = genes.use), 
  y = subset(x = star[["ss2"]], features = genes.use)
)
coembed <- merge(
  x = allen[genes.use, ],
  y = star[genes.use, ]
)

new.assay <- new(
  Class = "Assay",
  counts = coembed.assay@counts[, colnames(coembed)],
  data = coembed.assay@data[, colnames(coembed)],
  scale.data = matrix()
)

coembed[['combined']] <- new.assay
DefaultAssay(object = coembed) <- 'combined'

coembed <- ScaleData(object = coembed, features = genes.use, scale = TRUE)
coembed <- RunPCA(object = coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(object = coembed, dims = 1:30)
coembed$celltype <- ifelse(is.na(coembed$subclass), coembed$predicted.id, coembed$subclass)

DimPlot(coembed, group.by = 'celltype', label = TRUE, repel = TRUE) +
   ggsave("figures/spatial/starmap_coembed_cluster.png", height = 6, width = 8)

DimPlot(coembed, group.by = 'technology') +
  ggsave("figures/spatial/starmap_coembed_dataset.png", height = 6, width = 8)
