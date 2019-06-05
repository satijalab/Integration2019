library(Seurat)
library(ggplot2)

set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)
dir.create("figures/atac/pbmc/")
dir.create("analysis_data/atac/10x/")
dir.create("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/")

# load 10x v3 pbmc data
# args <- c("seurat_objects/pbmc_10k_v3.rds", "seurat_objects/pbmc_10k_atac.rds")
rna <- readRDS(args[1])
atac <- readRDS(args[2])
save.prefix <- unlist(strsplit(unlist(strsplit(args[2], '/'))[[2]], '.', fixed = TRUE))[[1]]

DimPlot(rna, group.by = 'celltype', label = TRUE, repel = TRUE) +
  ggsave("figures/atac/pbmc/sn1_fig1.png", dpi = 500, height = 6, width = 7)

# Transfer
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  features = VariableFeatures(rna),
  reference.assay = 'RNA',
  query.assay = 'RNA',
  reduction = 'cca'
)

# transfer cluster ID
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = atac[['lsi']]
)

atac <- AddMetaData(atac, metadata = celltype.predictions)

DimPlot(atac, group.by = 'predicted.id', reduction = 'umap.lsi',
        label = TRUE, repel = TRUE) +
  theme_classic(base_size = 10) +
  ggsave(paste0("figures/atac/pbmc/", save.prefix, "_celltype_transfer.png"),
         height = 6, width = 6, dpi = 800)

FeaturePlot(atac, 'prediction.score.max', reduction = 'umap.lsi') +
  theme_classic(base_size = 10) +
  ggsave(paste0("figures/atac/pbmc/", save.prefix, "_celltype_scores.png"),
         height = 6, width = 6, dpi = 800)

write.table(celltype.predictions,
            paste0("analysis_data/atac/10x/", save.prefix, "_predictions.tsv"),
            sep = "\t", quote = FALSE)

# coembedding
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(
  object = rna,
  assay = 'RNA',
  slot = 'data'
)[genes.use, ]

imputation <- TransferData(
  anchorset = transfer.anchors,
  refdata = refdata,
  weight.reduction = atac[['lsi']]
)

DefaultAssay(atac) <- 'RNA'
atac[["imputed"]] <- imputation
atac[['ATAC']] <- NULL

coembed.assay <- merge(
  x = subset(x = rna[["RNA"]], features = genes.use), 
  y = subset(x = atac[["imputed"]], features = genes.use)
)
coembed <- merge(
  x = rna[genes.use, ],
  y = atac[genes.use, ]
)

new.assay <- new(
  Class = "Assay",
  counts = coembed.assay@counts[, colnames(coembed)],
  data = coembed.assay@data[, colnames(coembed)],
  scale.data = matrix()
)

coembed[['combined']] <- new.assay
DefaultAssay(object = coembed) <- 'combined'

coembed <- ScaleData(object = coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(object = coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(object = coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

saveRDS(coembed, "seurat_objects/pbmc_atac_coembed.rds")

DimPlot(coembed, reduction = 'umap', group.by = 'orig.ident', pt.size=0.1) +
  ggsave(paste0("figures/atac/pbmc/", save.prefix, "_coembed_tech.png"),
         height = 5, width = 6, dpi = 500)

DimPlot(coembed, reduction = 'umap', group.by = 'celltype', label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggsave(paste0("figures/atac/pbmc/", save.prefix, "_coembed_celltype.png"),
         height = 5, width = 6, dpi = 500)

# manually define groups and show vln plots (eg blacklist regions, etc)
clusterx_cells <- readLines("raw_data/10x_atac/clusterx.txt")
clustery_cells <- readLines("raw_data/10x_atac/clustery.txt")

atac$clusterx <- ifelse(colnames(atac) %in% clusterx_cells, "ClusterX", "Other")
atac$clustery <- ifelse(colnames(atac) %in% clustery_cells, "ClusterY", "Other")

DimPlot(coembed, cells.highlight = clusterx_cells, pt.size = 0.1) + 
  theme_classic(base_size = 20) +
  theme(legend.position = 'none') +
  ggsave("figures/atac/pbmc/clusterx_dimplot.png", height = 5, width = 5, dpi = 400)

DimPlot(coembed, cells.highlight = clustery_cells, pt.size = 0.1) + 
  theme_classic(base_size = 20) +
  theme(legend.position = 'none') +
  ggsave("figures/atac/pbmc/clustery_dimplot.png", height = 5, width = 5, dpi = 400)

VlnPlot(atac, "blacklist_region_fragments", group.by = 'clusterx') + 
  theme_classic(base_size = 20) +
  theme(legend.position = 'none') +
  ggsave("figures/atac/pbmc/clusterx_vlnplot.png", height = 5, width = 5, dpi = 400)

# volcano plot between clusterY and other monocyte cells in atac
Idents(atac) <- "predicted.id"
cd14.mono <- WhichCells(atac, idents = 'CD14+ Monocytes')
cd14.other <- setdiff(cd14.mono, clustery_cells)
cly <- rep("clustery", length(clustery_cells))
names(cly) <- clustery_cells
cd14 <- rep('cd14', length(cd14.other))
names(cd14) <- cd14.other

atac <- AddMetaData(atac, metadata = c(cly, cd14), col.name = 'clustery_mono')
Idents(atac) <- "clustery_mono"

# plot average expression correlation between monocyte populations
express <- GetAssayData(atac, assay = 'RNA', slot = 'data')
cd14.a <- express[, cd14.other]
cd14.b <- express[, clustery_cells]
cd14.a.av <- Matrix::rowMeans(cd14.a)
cd14.b.av <- Matrix::rowMeans(cd14.b)

cd14.df <- data.frame('cluster_y' = cd14.b.av, 'cd14' = cd14.a.av)

ggplot(cd14.df, aes(cluster_y, cd14)) +
  geom_point(size = 0.1) +
  xlab("Average activity Cluster Y cells") +
  ylab("Average activity CD14+ monocytes") +
  theme_classic() +
  ggsave("figures/atac/pbmc/cluster_y_cd14_correlation.png", height = 5, width = 5)

# save barcode lists for each cell type
Idents(atac) <- "predicted.id"
atac.subset <- subset(atac, subset = prediction.score.max > 0.5)

ct <- atac$predicted.id
for (i in 1:length(unique(ct))) {
  ct.use <- unique(ct)[[i]]
  cells.use <- names(which(ct == ct.use))
  valid.filename <- gsub(pattern = '[ /+]', replacement = '', x = ct.use)
  if (ct.use != "Unknown") {
    writeLines(cells.use, con = paste0("analysis_data/atac/10x/celltypes/single_cell/", save.prefix, "/", valid.filename, ".txt"))
  }
}
