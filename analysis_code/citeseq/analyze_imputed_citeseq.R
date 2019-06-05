library(Seurat)
library(ggplot2)
library(dplyr)
library(mixtools)
library(ape)
set.seed(1234)


`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

FindSpatialMarkers <- function(object, assay, features = NULL, cells = NULL,
                               dims = c(1, 2), sd = 10, binsize = NULL, expression.limit = 5,
                               reduction = 'spatial', verbose = TRUE, scale.factor = 1) {
  if (!(reduction %in% names(object@reductions))){
    stop('Object does not contain requested reduction')
  }
  assay <- assay %||% DefaultAssay(object)
  cells <- cells %||% colnames(object)
  features <- features %||% rownames(object)
  spatial_coords <- as.data.frame(Embeddings(object[[reduction]])[cells, dims])
  colnames(spatial_coords) <- paste0('spatial', dims)
  if(!is.null(binsize)){
    if (verbose)
      message("Binning spatial regions")
    spatial_coords <- as.data.frame(apply(spatial_coords, 2, function(x) round((x*scale.factor)/binsize)))
    if (verbose)
      message("Averaging gene expression in spatial bins")
    if (length(dims) == 2) {
      bins <- spatial_coords %>%
        unite(bin, spatial1, spatial2) %>%
        unique()
      exprs <- GetAssayData(object, assay = assay, slot = 'data')[features, cells]
      binned_expression <- matrix(0, nrow = length(features), ncol = nrow(bins))
      colnames(binned_expression) <- bins[, 1]
      rownames(binned_expression) <- features
      for(i in 0:max(spatial_coords[, 1])){
        for(j in 0:max(spatial_coords[, 2])){
          bin.cells <- rownames(spatial_coords[spatial_coords[, 1] == i & spatial_coords[, 2] == j,])
          if(length(bin.cells) == 1){
            binned_expression[,paste(i, j, sep = '_')] <- exprs[features, bin.cells]
          } else if(length(bin.cells) > 0){
            av_exp <- apply(exprs[features, bin.cells], 1, mean)
            binned_expression[,paste(i, j, sep = '_')] <- av_exp
          }
        }
      }
    } else if (length(dims) == 1) {
      bins <- unique(spatial_coords)
      exprs <- GetAssayData(object, assay = assay, slot = 'data')[features, cells]
      binned_expression <- matrix(0, nrow = length(features), ncol = nrow(bins))
      colnames(binned_expression) <- bins[, 1]
      rownames(binned_expression) <- features
      spatial_coords <- as.matrix(spatial_coords)
      for(i in 0:max(spatial_coords[, 1])){
        bin.cells <- names(spatial_coords[spatial_coords[, 1] == i, ])
        if(length(bin.cells) == 1){
          binned_expression[ , as.character(i)] <- exprs[features, bin.cells]
        } else if(length(bin.cells) > 0){
          av_exp <- apply(exprs[features, bin.cells], 1, mean)
          binned_expression[ , as.character(i)] <- av_exp
        }
      }
    } else {
      stop("Too many dimensions requested")
    }
  } else {
    binned_expression <- GetAssayData(object = object, assay = assay, slot = 'data')[features, cells]
  }
  if(verbose) message("Computing distance weights")
  if (!is.null(binsize)) {
    spatial_bins <- as.matrix(unique(spatial_coords[, dims]))
  } else {
    spatial_bins <- as.matrix(spatial_coords[, dims])
  }
  distances <- as.matrix(dist(spatial_bins))
  weights <- exp(-1*distances/(2*(1/sd))^2)
  diag(weights) <- 0
  weights[is.infinite(weights)] <- 0
  ivals <- data.frame()
  if(verbose) {
    message("Calculating spatial autocorrelation")
    pb <- txtProgressBar(min = 1, max = nrow(binned_expression), style = 3,  file = stderr())
  }
  for(i in 1:nrow(binned_expression)){
    if (sum(binned_expression[i, ]) > expression.limit) {
      mi <- Moran.I(binned_expression[i,], weight = weights)
      ivals[rownames(binned_expression)[i],'I'] <- mi$observed
      ivals[rownames(binned_expression)[i],'expected.I'] <- mi$expected
      ivals[rownames(binned_expression)[i],'sd'] <- mi$sd
      ivals[rownames(binned_expression)[i],'p.value'] <- mi$p.value
      ivals[rownames(binned_expression)[i],'gene'] <- rownames(binned_expression)[i]
    }
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  if(verbose) message("")
  return(ivals[order(ivals$I, decreasing = T),])
}

imputed <- readRDS("seurat_objects/imputed_citeseq.rds")
citeseq <- readRDS("seurat_objects/citeseq.rds")

# Featureplots for key imputed proteins
selected.proteins <- c('CD3', 'CD4', 'CD8a', 'CD69', 'CD14', 'CD16', 'CD34', 'CD57')

DefaultAssay(imputed) <- 'ADT'
plots <- FeaturePlot(
  object = imputed,
  features = selected.proteins,
  reduction = 'umap',
  pt.size = 0.1,
  label.size = 1,
  min.cutoff = 'q40',
  max.cutoff = 'q99',
  cols = c("grey", "darkgreen"),
  combine = FALSE
)

for(x in 1:length(plots)) {
  plots[[x]] <- plots[[x]] + theme_classic(base_size = 20) + theme(legend.position = 'none')
}

CombinePlots(plots, ncol = 4) +
  ggsave(filename = "figures/citeseq/key_featureplot_imputed_hca.png",
         height = 30, width = 55, units = 'cm')

# plot failed ADTs
plots <- FeaturePlot(
  object = citeseq,
  features = c('CD25', 'CD197-CCR7'),
  reduction = 'umap',
  pt.size = 0.1,
  label.size = 1,
  min.cutoff = 'q1', max.cutoff = 'q99',
  cols = c("grey", "darkgreen"),
  combine = FALSE
)

for(x in 1:length(plots)) {
  plots[[x]] <- plots[[x]] + theme_classic(base_size = 20)
}

CombinePlots(plots, ncol = 2) +
  ggsave(filename = "figures/citeseq/featureplot_failed_adt.png",
         height = 20, width = 40, units = 'cm')

# plot all ADTs
p <- FeaturePlot(
  object = imputed,
  features = rownames(imputed),
  reduction = 'umap',
  pt.size = 0.1,
  min.cutoff = 'q40',
  max.cutoff = 'q99',
  cols = c("grey", "darkgreen"),
  combine = FALSE
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + cowplot::theme_cowplot(font_size = 50) + theme(legend.position = 'none')
}
plots <- CombinePlots(p, ncol = 5)
ggsave(plot = plots, filename = "figures/citeseq/hca_adt_featureplots.png", height = 50, width = 50, dpi = 100, limitsize = FALSE)

# plot all marker genes
celltype.markers <- c(
  'AVP', 'LMO4', 'PF4', 'BLVRB', 'MME', 'DERL3', 'CLEC9A',
  'CD1C', 'MPO', 'AZU1', 'CD14', 'FCGR3A', 'VPREB3', 'MS4A1', 'CD79A',
  'IGKC', 'PF4', 'XCL1', 'CD8A', 'CD4', 'CCL5', 'CCR7', 'SH2D1A'
)
DefaultAssay(imputed) <- 'RNA'
p <- FeaturePlot(
  object = imputed,
  features = celltype.markers,
  reduction = 'umap',
  pt.size = 0.1,
  # min.cutoff = 'q40',
  max.cutoff = 'q99',
  combine = FALSE
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + cowplot::theme_cowplot(font_size = 50) + theme(legend.position = 'none')
}
plots <- CombinePlots(p, ncol = 5, legend = 'none')
ggsave(plot = plots, filename = "figures/citeseq/hca_markers_featureplots.png", height = 50, width = 50, dpi = 100, limitsize = FALSE)


# Cluster cells
DefaultAssay(imputed) <- "integrated"
imputed <- FindNeighbors(imputed, nn.eps = 0.5, dims = 1:30)
imputed <- FindClusters(
  object = imputed,
  nstart = 10,
  graph.name = "integrated_snn",
  resolution = 0.8
)

DimPlot(
  object = imputed,
  reduction = 'umap',
  group.by = 'integrated_snn_res.0.8',
  label = TRUE,
  pt.size = 0.1,
  label.size = 3
) + NoLegend() +
  ggsave("figures/citeseq/hca_clusters.png", height = 20, width = 20, units = 'cm')

# take cd69+ HCA cells and find markers
cd8.all <- subset(imputed, idents = 4)

# Fit mixture model to identify CD69+ population
nmix <- normalmixEM(
  x = FetchData(cd8.all, vars = "adt_CD69")[,1],
  k = 3
)
Idents(cd8.all, Cells(cd8.all)[which(nmix$posterior[,2]>0.1)]) <- 'CD69Low'
Idents(cd8.all, Cells(cd8.all)[which(nmix$posterior[,2]<=0.1)]) <- 'CD69Hi'

# Find markers that separate the population
cd69_markers <- FindMarkers(
  object = cd8.all,
  ident.1 = "CD69Hi",
  only.pos = FALSE,
  test.use = "LR",
  latent.vars = 'orig.ident'
)

cd69_markers <- cd69_markers[order(abs(cd69_markers$avg_logFC), decreasing = TRUE), ]
top.markers <- head(cd69_markers, 25)
top.markers <- top.markers[order(top.markers$avg_logFC, decreasing = TRUE), ]
write.table(x = top.markers, file = "analysis_data/cd69_markers.tsv", quote = FALSE, sep = '\t')

# average heatmap of CD69 markers split by replicate
cd8.all[['cd69.pos']] <- Idents(cd8.all)
cd8.all[['cd69_replicate']] <- paste0(cd8.all$cd69.pos, "_", cd8.all$orig.ident)
Idents(cd8.all) <- 'cd69_replicate'

DefaultAssay(cd8.all) <- "RNA"
object <- cd8.all
object@assays <- list('RNA' = object@assays$RNA)
object.list <- SplitObject(object, split.by = 'orig.ident')
for(i in 1:length(object.list)) {
  DefaultAssay(object.list[[i]]) <- 'RNA'
  object.list[[i]] <- ScaleData(object.list[[i]], features = rownames(top.markers))
}

markers.to.plot <- rownames(top.markers)
alldata <- matrix(nrow = length(x = markers.to.plot), ncol = 0)
for(i in 1:length(x = object.list)) {
  idents <- sort(x = names(x = which(x = table(Idents(object = object.list[[i]])) > 1)))
  newdata <- sapply(X = idents, function(x) {
    rowMeans(x = GetAssayData(object = object.list[[i]], slot = "scale.data", assay = "RNA")[, WhichCells(object = object.list[[i]], idents = x)])
  })[markers.to.plot, ]
  colnames(x = newdata) <- paste(idents, names(x = object.list)[i], sep = "_")
  alldata <- cbind(alldata, newdata)
}
rownames(alldata) <- make.unique(names = rownames(x = alldata))
alldata <- alldata[, sort(x = colnames(x = alldata))]
avg <- CreateSeuratObject(counts = alldata)
avg[["RNA"]]@scale.data <- alldata
hm_cd69 <- DoHeatmap(avg,features = markers.to.plot,disp.max = 0.75)+NoLegend()
ggsave(hm_cd69, filename = "figures/citeseq/avg_heatmap_cd69.png", height = 6, width = 5)


# order cells along CD69 gradient
# need to make this range from 1 to max expression, so binning works as expected
# this gives similar/same genes as above FindMarkers
cd69.data <- as.matrix(cd8.all[['ADT']]['CD69', ])
cd69.data <- cd69.data - min(cd69.data)
cd69.data <- cd69.data / max(cd69.data)
cd69.data <- cd69.data * 1000
cd69.data <- t(cd69.data)
colnames(cd69.data) <- "1"
dr.obj <- CreateDimReducObject(
  embeddings = cd69.data,
  assay = 'ADT',
  key = 'CD69_'
)
cd8.all[['CD69']] <- dr.obj

DefaultAssay(imputed) <- 'RNA'
cd69.markers <- FindSpatialMarkers(
  object = cd8.all,
  assay = 'RNA',
  features = rownames(imputed),
  reduction = 'CD69',
  dims = 1,
  binsize = 10,
  expression.limit = 0,
  sd = 1
)

# calculate for cite-seq cells
DefaultAssay(citeseq) <- 'RNA'
citeseq <- FindNeighbors(citeseq, dims = 1:30, nn.eps = 0.5)
citeseq <- FindClusters(citeseq, resolution = 1)
cd8.cite <- subset(citeseq, idents = 10)

cd69.data <- as.matrix(cd8.cite[['ADT']]['CD69', ])
cd69.data <- cd69.data - min(cd69.data)
cd69.data <- cd69.data / max(cd69.data)
cd69.data <- cd69.data * 1000
cd69.data <- t(cd69.data)
colnames(cd69.data) <- "1"
dr.obj <- CreateDimReducObject(
  embeddings = cd69.data,
  assay = 'ADT',
  key = 'CD69_'
)
cd8.cite[['CD69']] <- dr.obj
common.genes <- intersect(rownames(imputed), rownames(cd8.cite))
cd69.citeseq <- FindSpatialMarkers(
  object = cd8.cite,
  assay = 'RNA',
  features = common.genes,
  reduction = 'CD69',
  dims = 1,
  binsize = 10,
  expression.limit = 0,
  sd = 1
)

highlight.genes <- c('IFNG', "CCL4L2", "CCL3", "CCL4", "XCL2")

joint.moran <- left_join(cd69.markers, cd69.citeseq, by = 'gene')
joint.moran$color <- ifelse(joint.moran$gene %in% highlight.genes, 'red', 'black')

# save high moran genes for each experiment
high_hca <- filter(joint.moran, I.x > 0.5, I.y < 0.1)
high_cite <- filter(joint.moran, I.y > 0.25, I.x < 0.25)

write.table(x = high_hca, file = 'analysis_data/cd69_high_moran_hca.tsv', sep = "\t", quote = FALSE)
write.table(x = high_cite, file = 'analysis_data/cd69_high_moran_citeseq.tsv', sep = "\t", quote = FALSE)

ggplot(joint.moran, aes(I.x, I.y, color = color)) +
  geom_point(size = 0.1, alpha=0.5) +
  scale_color_manual(values = c('black', 'red')) +
  theme_classic(base_size = 10) +
  geom_text(aes(label=ifelse(gene %in% highlight.genes, gene, '')), hjust=0, vjust=0) +
  theme(legend.position = 'none') +
  xlab("HCA") +
  ylab("CITE-seq") +
  xlim(c(min(cd69.markers$I), 1)) + ylim(c(min(cd69.markers$I), 1)) +
  ggtitle("Moran's I: CD69") +
  ggsave("figures/citeseq/cd69_moran.png", height = 10, width = 10, units = 'cm', dpi = 600)

