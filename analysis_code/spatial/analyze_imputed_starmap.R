library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(ggrepel)
library(ape)
library(cowplot)
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

starmap.1 <- readRDS("seurat_objects/20180505_BY3_1kgenes_imputed.rds")
starmap.2 <- readRDS("seurat_objects/20180410-BY3_1kgenes_imputed.rds")
allstar <- readRDS("seurat_objects/integrated_starmap.rds")
allen <- readRDS("seurat_objects/allen_brain.rds")
osm <- readRDS("seurat_objects/osm_fish.rds")
drops <- readRDS("seurat_objects/dropseq_cortex.rds")

dir.create("figures/spatial/spatial_featureplots/")

# add starmap metadata
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

starmap1.type <- class_labels.1$ClusterName
names(starmap1.type) <- class_labels.1$cellname

starmap2.type <- class_labels.2$ClusterName
names(starmap2.type) <- class_labels.2$cellname

allstar <- AddMetaData(allstar, metadata = c(starmap1.type, starmap2.type), col.name = 'starmap.celltype')

# show predicted labels on starmap UMAP
allstar$predicted.celltype <- ifelse(allstar$prediction.score.max > 0.5, allstar$predicted.id, NA)
# cells.use <- rownames(allstar[[]][!is.na(allstar$predicted.celltype), ])
p1 <- DimPlot(allstar, reduction = 'umap', group.by = 'starmap.celltype', label = TRUE, repel = TRUE) + ggtitle('STARmap classifications')
p2 <- DimPlot(allstar, reduction = 'umap', group.by = 'predicted.celltype', label = TRUE, repel = TRUE) + ggtitle('Transferred celltype labels')
cowplot::plot_grid(p1, p2) +
  ggsave("figures/spatial/starmap_label_comparision.pdf", height = 5, width = 15)

# Show key imputed genes
DefaultAssay(starmap.1) <- "RNA"
DefaultAssay(starmap.2) <- "RNA"
starmap.genes <- rownames(starmap.1)

show.starmap <- c('Cux2', 'Lamp5', 'Rorb', 'Rab3c', 'Syt6', 'Sox2ot', 'Bsg', 'Sst')
show.osm.long <- c("Vtn", "Serpinf1", "Sox10", "Hexb", "Apln")
show.osm.short <- c("Vtn", "Sox10")
new.genes <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2', 'Tcrb')
osm.validation <- c('Rorb', 'Syt6', 'Lamp5', 'Sox10')

starmap.osm <- intersect(rownames(starmap.1), rownames(osm))
osm.unique <- setdiff(rownames(osm), rownames(starmap.1))

FormatPlots <- function(plots, ncol = NULL) {
  for (i in 1:length(plots)) {
    plots[[i]] <- plots[[i]] + 
      SpatialTheme() + theme(plot.title = element_blank())
  }
  plots <- CombinePlots(plots, ncol = ncol)
}

# show imputed using leave-out object for each gene
plots.ss2 <- list()
plots.ds <- list()
for(i in 1:length(show.starmap)) {
  dat <- readRDS(paste0('seurat_objects/spatial_leaveout/starmap_1_', show.starmap[[i]], '.rds'))
  if (show.starmap[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  DefaultAssay(dat) <- 'ss2'
  p <- PolyFeaturePlot(dat, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots.ss2[[i]] <- p
  DefaultAssay(dat) <- 'dropseq'
  p <- PolyFeaturePlot(dat, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots.ds[[i]] <- p
}

r1.imputed.ss2 <- plot_grid(plotlist = plots.ss2, ncol = length(plots.ss2))
r1.imputed.ds <- plot_grid(plotlist = plots.ds, ncol = length(plots.ds))

plots.ss2 <- list()
plots.ds <- list()
for(i in 1:length(show.starmap)) {
  dat <- readRDS(paste0('seurat_objects/spatial_leaveout/starmap_2_', show.starmap[[i]], '.rds'))
  if (show.starmap[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  DefaultAssay(dat) <- 'ss2'
  p <- PolyFeaturePlot(dat, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots.ss2[[i]] <- p
  DefaultAssay(dat) <- 'dropseq'
  p <- PolyFeaturePlot(dat, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots.ds[[i]] <- p
}

r2.imputed.ss2 <- plot_grid(plotlist = plots.ss2, ncol = length(plots.ss2))
r2.imputed.ds <- plot_grid(plotlist = plots.ds, ncol = length(plots.ds))

# plot measured values the same way (combining many plots) for consistent look
plots <- list()
for(i in 1:length(show.starmap)) {
  if (show.starmap[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap.1, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots[[i]] <- p
}

r1.measured <- plot_grid(plotlist = plots, ncol = length(plots))

plots <- list()
for(i in 1:length(show.starmap)) {
  if (show.starmap[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap.2, show.starmap[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  plots[[i]] <- p
}

r2.measured <- plot_grid(plotlist = plots, ncol = length(plots))

# show genes not in the starmap dataset
DefaultAssay(starmap.1) <- "ss2"
DefaultAssay(starmap.2) <- "ss2"

r1.imputed.osm.long <- PolyFeaturePlot(
  object = starmap.1,
  features = show.osm.long,
  ncol= length(show.osm.long),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r2.imputed.osm.long <- PolyFeaturePlot(
  object = starmap.2,
  features = show.osm.long,
  ncol= length(show.osm.long),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r1.imputed.osm.short <- PolyFeaturePlot(
  object = starmap.1,
  features = show.osm.short,
  ncol= length(show.osm.short),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r2.imputed.osm.short <- PolyFeaturePlot(
  object = starmap.2,
  features = show.osm.short,
  ncol= length(show.osm.short),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r1.new <- PolyFeaturePlot(
  object = starmap.1,
  features = new.genes,
  ncol = 2,
  flip.coords = TRUE,
  min.cutoff = 0,
  max.cutoff = 'q99' 
) + SpatialTheme()

r2.new <- PolyFeaturePlot(
  object = starmap.2,
  features = new.genes,
  ncol = 3,
  flip.coords = TRUE,
  min.cutoff = 0,
  max.cutoff = 'q99'
) + SpatialTheme()

# dropseq imputation
DefaultAssay(starmap.1) <- "dropseq"
DefaultAssay(starmap.2) <- "dropseq"
new.genes.ds <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2')

r1.imputed.osm.long.ds <- PolyFeaturePlot(
  object = starmap.1,
  features = show.osm.long,
  ncol= length(show.osm.long),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r2.imputed.osm.long.ds <- PolyFeaturePlot(
  object = starmap.2,
  features = show.osm.long,
  ncol= length(show.osm.long),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r1.imputed.osm.short.ds <- PolyFeaturePlot(
  object = starmap.1,
  features = show.osm.short,
  ncol= length(show.osm.short),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r2.imputed.osm.short.ds <- PolyFeaturePlot(
  object = starmap.2,
  features = show.osm.short,
  ncol= length(show.osm.short),
  min.cutoff = 0,
  max.cutoff = 'q99', 
  flip.coords = TRUE
) + SpatialTheme()

r1.new.ds <- PolyFeaturePlot(
  object = starmap.1,
  features = new.genes.ds,
  ncol = 2,
  flip.coords = TRUE,
  min.cutoff = 0,
  max.cutoff = 'q99' 
) + SpatialTheme()

r2.new.ds <- PolyFeaturePlot(
  object = starmap.2,
  features = new.genes.ds,
  ncol = 2,
  flip.coords = TRUE,
  min.cutoff = 0,
  max.cutoff = 'q99'
) + SpatialTheme()

osm.gene <- FeaturePlot(
  object = osm,
  features = osm.validation,
  reduction = 'spatial',
  coord.fixed = TRUE,
  pt.size = 0.5,
  cols = viridis::viridis(10),
  combine = FALSE,
  min.cutoff = 0,
  max.cutoff = 'q99'
)
osm.gene <- FormatPlots(osm.gene, ncol = 2)

# ss2
ggsave(plot = r1.imputed.ss2, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed.pdf")
ggsave(plot = r1.measured, filename = "figures/spatial/spatial_featureplots/starmap_r1_measured.pdf")
ggsave(plot = r1.imputed.osm.long, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed_osm_genes_long.pdf")
ggsave(plot = r1.imputed.osm.short, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed_osm_genes_short.pdf")
ggsave(plot = r1.new, filename = "figures/spatial/spatial_featureplots/starmap_r1_new.pdf")

ggsave(plot = r2.imputed.ss2, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed.pdf")
ggsave(plot = r2.measured, filename = "figures/spatial/spatial_featureplots/starmap_r2_measured.pdf")
ggsave(plot = r2.imputed.osm.long, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed_osm_genes_long.pdf")
ggsave(plot = r2.imputed.osm.short, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed_osm_genes_short.pdf")
ggsave(plot = r2.new, filename = "figures/spatial/spatial_featureplots/starmap_r2_new.pdf")

# dropseq
ggsave(plot = r1.imputed.ds, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed_dropseq.pdf")
ggsave(plot = r1.imputed.osm.long.ds, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed_osm_genes_long_dropseq.pdf")
ggsave(plot = r1.imputed.osm.short.ds, filename = "figures/spatial/spatial_featureplots/starmap_r1_imputed_osm_genes_short_dropseq.pdf")
ggsave(plot = r1.new.ds, filename = "figures/spatial/spatial_featureplots/starmap_r1_new_dropseq.pdf")

ggsave(plot = r2.imputed.ds, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed_dropseq.pdf")
ggsave(plot = r2.imputed.osm.long.ds, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed_osm_genes_long_dropseq.pdf")
ggsave(plot = r2.imputed.osm.short.ds, filename = "figures/spatial/spatial_featureplots/starmap_r2_imputed_osm_genes_short_dropseq.pdf")
ggsave(plot = r2.new.ds, filename = "figures/spatial/spatial_featureplots/starmap_r2_new_dropseq.pdf")

ggsave(plot = osm.gene, filename = "figures/spatial/spatial_featureplots/osm_genes.png",
       height = 11.5, width = 6, units = 'in', dpi = 500, bg = 'transparent')

# Find new spatial DE genes
DefaultAssay(starmap.1) <- 'ss2'
DefaultAssay(starmap.2) <- 'ss2'
starmap.1 <- FindVariableFeatures(starmap.1, selection.method = 'dispersion', nfeatures = 3000)
starmap.2 <- FindVariableFeatures(starmap.2, selection.method = 'dispersion', nfeatures = 3000)
genes.use <- unique(c(VariableFeatures(starmap.1), VariableFeatures(starmap.2)))

spatial.markers.1 <- FindSpatialMarkers(
  object = starmap.1,
  assay = 'ss2',
  features = genes.use,
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)

spatial.markers.2 <- FindSpatialMarkers(
  object = starmap.2,
  assay = 'ss2',
  features = genes.use,
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)

ds.genes.use <- intersect(genes.use, rownames(starmap.1[['dropseq']]))
spatial.markers.1.ds <- FindSpatialMarkers(
  object = starmap.1,
  assay = 'dropseq',
  features = ds.genes.use,
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)

spatial.markers.2.ds <- FindSpatialMarkers(
  object = starmap.2,
  assay = 'dropseq',
  features = ds.genes.use,
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)

# cross-tech correlation
sp.1 <- spatial.markers.1[rownames(spatial.markers.1.ds), ]
joint.spatial <- left_join(spatial.markers.1.ds, sp.1, by = 'gene')
r <- cor(joint.spatial$I.x, joint.spatial$I.y, use = 'complete.obs')
ggplot(joint.spatial, aes(I.x, I.y)) +
  geom_point(size=0.1, alpha=0.5) +
  theme_classic(base_size = 6) +
  geom_smooth(method = 'lm', se = FALSE) +
  xlab("Moran's I Drop-seq") +
  ylab("Moran's I SMART-seq2") +
  ggtitle(paste0('Pearson correlation: ', as.character(round(r, 3)))) +
  ggsave("figures/spatial/starmap_moran_rep1_ss2_ds.pdf", 
         height = 6, width = 4, units = 'cm')

sp.2 <- spatial.markers.2[rownames(spatial.markers.2.ds), ]
joint.spatial <- left_join(spatial.markers.2.ds, sp.2, by = 'gene')
r <- cor(joint.spatial$I.x, joint.spatial$I.y, use = 'complete.obs')
ggplot(joint.spatial, aes(I.x, I.y)) +
  geom_point(size=0.1, alpha=0.5) +
  theme_classic(base_size = 6) +
  geom_smooth(method = 'lm', se = FALSE) +
  xlab("Moran's I Drop-seq") +
  ylab("Moran's I SMART-seq2") +
  ggtitle(paste0('Pearson correlation: ', as.character(round(r, 3)))) +
  ggsave("figures/spatial/starmap_moran_rep2_ss2_ds.pdf", 
         height = 6, width = 4, units = 'cm')


# correlation between gene spatial dependence in both replicates
joint.spatial <- left_join(spatial.markers.1, spatial.markers.2, by = 'gene')
r <- cor(joint.spatial$I.x, joint.spatial$I.y, use = 'complete.obs')

divergent.spatial <- filter(joint.spatial, I.x > 0.4, I.y < 0.2)
p1 <- PolyFeaturePlot(starmap.1, head(divergent.spatial$gene, 8), flip.coords = TRUE, ncol = 4) + 
  theme_dark(base_size = 8) + DarkTheme() + NoGrid()
p2 <- PolyFeaturePlot(starmap.2, head(divergent.spatial$gene, 8), flip.coords = TRUE, ncol = 4) + 
  theme_dark(base_size = 8) + DarkTheme() + NoGrid()
cowplot::plot_grid(p1, p2) +
  ggsave("figures/spatial/spatial_outlier_genes.pdf", height = 20, width = 30, units = 'cm')

write.table(divergent.spatial, 'analysis_data/divergent_spatial_genes.tsv', sep = '\t', quote = FALSE)
write.table(spatial.markers.1, 'analysis_data/spatial_markers_starmap1.tsv', sep = '\t', quote = FALSE)
write.table(spatial.markers.2, 'analysis_data/spatial_markers_starmap2.tsv', sep = '\t', quote = FALSE)

# find VLMC markers and endo markers in allen data
Idents(allen) <- 'subclass'
markers.vlmc <- FindMarkers(allen, ident.1 = 'VLMC', min.pct = 0.5, min.diff.pct = 0.5)
markers.endo <- FindMarkers(allen, ident.1 = 'Endo', min.pct = 0.5, min.diff.pct = 0.5)
markers.peri <- FindMarkers(allen, ident.1 = 'Peri', min.pct = 0.5, min.diff.pct = 0.5)

top.vlmc <- rownames(markers.vlmc[markers.vlmc$avg_logFC > 1,])
top.endo <- rownames(markers.endo[markers.endo$avg_logFC > 1,])
top.peri <- rownames(markers.peri[markers.peri$avg_logFC > 1,])

joint.spatial$celltype <- ifelse(joint.spatial$gene %in% top.vlmc, "VLMC",
                                 ifelse(joint.spatial$gene %in% top.endo, "Endo",
                                        ifelse(joint.spatial$gene %in% top.peri, "Peri", "Other")))

ggplot(joint.spatial) +
  geom_point(aes(I.x, I.y, color = celltype), size=0.1, alpha=0.5) +
  theme_classic(base_size = 6) +
  geom_smooth(aes(I.x, I.y), method = 'lm', se = FALSE) +
  xlab("Moran's I replicate 1") +
  ylab("Moran's I replicate 2") +
  scale_color_manual(values = c("red", "black", "blue", "orange")) +
  theme(legend.position = 'none') +
  ggtitle(paste0('Pearson correlation: ', as.character(round(r, 3)))) +
  ggsave("figures/spatial/starmap_rep1_vs_rep2_moran.png", 
         height = 6, width = 5, units = 'cm', dpi = 500)

# plot cell types on collapsed spatial positions
nts1 <- starmap.1

# assign celltypes by taking highest 50% cells within each celltype
classifications <- nts1[[]] %>% 
  group_by(predicted.id) %>% 
  mutate(celltype = ifelse(prediction.score.max > quantile(prediction.score.max, 1/2), predicted.id, 'Unassigned'))

nts1$celltype <- classifications$celltype

Idents(nts1) <- "none"

celltypes <- c(
  "L2/3 IT", "L4", "L5 IT", "L5 PT","L6 IT", "L6 CT", "L6b", "NP",
  "Lamp5", "Vip", "Pvalb", "Sst",
  "Oligo", "VLMC", "Astro", "Endo", "Macrophage", "SMC"
)

colors <- c(
  rep('firebrick2', 8),
  rep('darkorchid4', 4),
  rep('orange1', 6)
)

plots <- list()
for(i in 1:length(celltypes)) {
  cell1 <- WhichCells(nts1, expression = (celltype == celltypes[i]))
  cell.remain <- setdiff(colnames(nts1), cell1)
  plots[[i]] <- DimPlot(
    object = nts1,
    reduction = 'spatial',
    cells = c(cell.remain, cell1),
    dims = c(2,1),
    cells.highlight = WhichCells(nts1, cells = cell1),
    cols = 'grey',
    pt.size = 0.5,
    cols.highlight = colors[i],
    sizes.highlight = 1.5
  ) + NoLegend() + NoAxes() + ggtitle(celltypes[i])
}

CombinePlots(plots, ncol = length(plots)) +
  ggsave("figures/spatial/celltype_distibution_1.png",
         height = 6, width = 10, dpi = 800)


star2 <- starmap.2

# assign celltypes by taking highest 50% cells within each celltype
classifications <- star2[[]] %>% 
  group_by(predicted.id) %>% 
  mutate(celltype = ifelse(prediction.score.max > quantile(prediction.score.max, 1/2), predicted.id, 'Unassigned'))

star2$celltype <- classifications$celltype

Idents(star2) <- "none"

plots <- list()
for(i in 1:length(celltypes)) {
  cell1 <- WhichCells(star2, expression = (celltype == celltypes[i]))
  cell.remain <- setdiff(colnames(star2), cell1)
  plots[[i]] <- DimPlot(
    object = star2,
    reduction = 'spatial',
    cells = c(cell.remain, cell1),
    dims = c(2,1),
    cells.highlight = WhichCells(star2, cells = cell1),
    cols = 'grey',
    pt.size = 0.5,
    cols.highlight = colors[i],
    sizes.highlight = 1.5
  ) + NoLegend() + NoAxes() + ggtitle(celltypes[i])
}

CombinePlots(plots, ncol = length(plots)) +
  ggsave("figures/spatial/celltype_distibution_2.png",
         height = 4, width = 10, dpi = 800)


# heatmap with projected labels
classifications <- allstar[[]] %>% 
  group_by(predicted.id) %>% 
  mutate(celltype = ifelse(prediction.score.max > quantile(prediction.score.max, 1/2), predicted.id, 'Unassigned'))

allstar$celltype <- classifications$celltype

star <- allstar
Idents(star) <- 'celltype'
star <- subset(star, idents = 'Unassigned', invert = TRUE)
star <- subset(star, idents = names(which(table(Idents(star))>5)))

ordering <- c(
  "L2/3 IT", "L4", "L5 IT", "L5 PT","L6 IT", "L6 CT", "L6b", "NP",
  "Lamp5", "Vip", "Pvalb", "Sst",
  "Oligo", "VLMC", "Astro", "Endo", "Macrophage"
)

Idents(star) <- factor(Idents(star), levels = ordering)

markers_all <- FindAllMarkers(
  object = star,
  assay = "RNA",
  test.use = 'LR',
  latent.vars = 'orig.ident',
  only.pos = TRUE,
  min.diff.pct = 0.1,
  logfc.threshold = 0,
  min.pct = 0.2
)

markers_hm <- as.character(data.frame(markers_all %>% group_by(cluster) %>% top_n(-10, p_val))$gene)
DefaultAssay(star) <- 'RNA'
star <- ScaleData(star, features = markers_hm)
e2 <- subset(star, downsample = 50)
DoHeatmap(
  object = e2,
  angle = 90,
  features = markers_hm,
  disp.max = 2.5
) + theme(axis.text = element_text(size=5)) + 
  ggsave("figures/spatial/starmap_marker_heatmap.png",
         height = 10, width = 12, dpi = 800)

write.table(markers_all, 'analysis_data/starmap_markers.tsv', sep = '\t', quote = FALSE)

# save colorscale
png("figures/colorscale.png")
image(matrix(seq(100)), col = PurpleAndYellow())
dev.off()

# histogram of correlation between ss2 and dropseq for marker genes
compare.genes <- intersect(markers_hm, intersect(rownames(star[['ss2']]), rownames(star[['dropseq']])))
ss2 <- GetAssayData(star, assay = 'ss2', slot = 'data')[compare.genes, ]
ds <- GetAssayData(star, assay = 'dropseq', slot = 'data')[compare.genes, ]

matrix.cor <- function(m1, m2) {
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  cA <- m1 - rowMeans(m1)
  cB <- m2 - rowMeans(m2)
  sA <- sqrt(rowMeans(cA^2))
  sB <- sqrt(rowMeans(cB^2))
  return(rowMeans(cA * cB) / (sA * sB))
}

r <- matrix.cor(ss2, ds)

# plot expression of genes that don't agree across tech
FeaturePlot(allen, names(head(sort(r))), pt.size = 0.1) +
  ggsave("figures/spatial/ss2_low_cor_genes.png", dpi = 500, height = 15, width = 10)

FeaturePlot(drops, names(head(sort(r))), pt.size = 0.1) +
  ggsave("figures/spatial/dropseq_low_cor_genes.png", dpi = 500, height = 15, width = 10)

cowplot::plot_grid(p1, p2, labels = c("SMART-seq2", 'Drop-seq'))

# scatterplots for some genes
ss2.df <- as.data.frame(ss2)
ss2.df$gene <- rownames(ss2.df)
ds.df <- as.data.frame(ds)
ds.df$gene <- rownames(ds.df)

ss2.df <- tidyr::gather(ss2.df, cell, expression.ss2, 1:(ncol(ss2.df)-1))
ds.df <- tidyr::gather(ds.df, cell, expression.dropseq, 1:(ncol(ds.df)-1))

all.tech <- dplyr::left_join(ss2.df, ds.df, by = c('cell', 'gene'))

genes.plot <- c('Lamp5', 'Sst', 'Pvalb', 'Bsg', 'Cux2')
plots <- list()
for (i in 1:length(genes.plot)) {
  data.plot <- all.tech[all.tech$gene == genes.plot[[i]], ]
  rval <- cor(data.plot$expression.ss2, data.plot$expression.dropseq)
  plots[[i]] <- ggplot(data.plot, aes(expression.ss2, expression.dropseq)) +
    geom_point(size = 0.5) +
    theme_classic(base_size = 10) +
    xlab("SMART-seq2 imputation") +
    ylab("Drop-seq imputation") +
    ggtitle(label = genes.plot[[i]], subtitle = paste0("Pearson correlation: ", as.character(round(rval, 3))))
}

rhist <- data.frame("pearson.cor" = r, "gene" = names(r))
plots[[i+1]] <- ggplot(rhist, aes(pearson.cor)) +
  geom_histogram() +
  theme_classic(base_size = 10) +
  ggtitle(label = "Pearson correlation", subtitle = paste0("Median: ", as.character(round(median(r), 3))))

cowplot::plot_grid(plotlist = plots, ncol = 3) +
  ggsave("figures/spatial/dropseq_vs_ss2_imputation_hist_scatter.pdf", height = 5, width = 8)

# Moran's I calculate on original STARmap data
DefaultAssay(starmap.1) <- 'RNA'
DefaultAssay(starmap.2) <- 'RNA'


moran.1 <- FindSpatialMarkers(
  object = starmap.1,
  assay = 'RNA',
  features = rownames(starmap.1),
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)
moran.2 <- FindSpatialMarkers(
  object = starmap.2,
  assay = 'RNA',
  features = rownames(starmap.2),
  verbose = TRUE,
  binsize = 10,
  expression.limit = 20,
  sd = 5
)

combined.moran <- left_join(moran.1, moran.2, by = 'gene')
r <- cor(combined.moran$I.x, combined.moran$I.y, use = 'complete.obs')

ggplot(combined.moran, aes(I.x, I.y)) +
  geom_point(size=0.1, alpha=0.5) +
  theme_classic(base_size = 6) +
  geom_smooth(method = 'lm', se = FALSE) +
  xlab("Moran's I replicate 1") +
  ylab("Moran's I replicate 2") +
  ggtitle(paste0('Pearson correlation: ', as.character(round(r, 3)))) +
  ggsave("figures/spatial/starmap_rep1_vs_rep2_moran_original_starmap.pdf", 
         height = 5, width = 5, units = 'cm')

combined.r1 <- left_join(moran.1, spatial.markers.1, by = 'gene')
combined.r1 <- combined.r1[!(is.na(combined.r1$I.x)) & !(is.na(combined.r1$I.y)), ]
r <- cor(combined.r1$I.x, combined.r1$I.y, use = 'complete.obs')

ggplot(combined.r1, aes(x=I.y, y=I.x)) +
  geom_point(size=0.1, alpha=0.5) +
  theme_classic(base_size = 6) +
  ylim(c(0, 0.8)) + xlim(c(0, 0.8)) +
  xlab("Moran's I imputed") +
  ylab("Moran's I measured") +
  ggtitle("Spatial autocorrelation") +
  ggsave("figures/spatial/imputed_vs_measured_moran.pdf",
         height = 3, width = 3)

combined.r1 <- combined.r1[, c('I.x', 'I.y', 'gene')]
colnames(combined.r1) <- c('Measured', "Imputed", "gene")

tidyr::gather(combined.r1, dataset, moran, Measured:Imputed) %>%
  ggplot(., aes(moran, fill = dataset)) +
    geom_density(alpha = 0.5, size = 0.1) +
    theme_classic() +
    xlab("Moran's I") +
    ggsave("figures/spatial/moran_density_imputed_measured.pdf",
         height = 4, width = 5)

# transfer allen labels to dropseq to check cux2 expression
anchors <- FindTransferAnchors(
  reference = allen,
  query = drops
)

predicted.labels <- TransferData(
  anchorset = anchors,
  refdata = allen$subclass
)

drops <- AddMetaData(drops, predicted.labels)
Idents(drops) <- "predicted.id"

DimPlot(drops, label = T)

Idents(drops) <- factor(Idents(drops), levels = levels(Idents(allen)))

p1 <- VlnPlot(drops, 'Cux2', pt.size = 0.1) + NoLegend() + ggtitle("Drop-seq")
p2 <- VlnPlot(allen, "Cux2", pt.size = 0.1) + NoLegend() + ggtitle("SMART-seq2")

cowplot::plot_grid(p1, p2, ncol = 1) +
  ggsave("figures/spatial/cux2_ss2_dropseq.png", dpi = 500, height = 8, width = 5)


# subset to obly l4 and l2/3
drops.subset <- subset(drops, subset = predicted.id %in% c("L2/3 IT", "L4"))
allen.subset <- subset(allen, subset = subclass %in% c("L2/3 IT", "L4"))

p1 <- VlnPlot(drops.subset, 'Cux2', pt.size = 0.1) + NoLegend() + ggtitle("Drop-seq")
p2 <- VlnPlot(allen.subset, "Cux2", pt.size = 0.1) + NoLegend() + ggtitle("SMART-seq2")

cowplot::plot_grid(p1, p2, ncol = 2) +
  ggsave("figures/spatial/cux2_ss2_dropseq_l4_23_only.png", dpi = 500, height = 5, width = 5)
