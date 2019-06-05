args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/20180505_BY3_1kgenes.rds", 
#           "~/Projects/muir/seurat_objects/20180410-BY3_1kgenes.rds",
#           "/home/stuartt/github/muir/seurat_objects/integrated_starmap.rds",
#           "/home/stuartt/github/muir/seurat_objects/allen_brain.rds",
#           "/home/stuartt/github/muir/seurat_objects/dropseq_cortex.rds",
#           "/home/stuartt/github/muir/analysis_data/starmap_markers.tsv")

devtools::load_all(args[10])
library(dplyr)
library(ggplot2)
library(cowplot)

# Downsample features with a guided strategy
# 
# Given a Seurat object and a marker list (from FindMarkers/FindAllMarkers), return a vector of 
# nfeatures. First select the N top markers for each cluster. If nfeatures is not a multiple of the 
# number of clusters, select the last set of cluster markers such that the most significant are 
# picked first up to nfeatures. If nfeatures is greater than the number of unique markers, then 
# start selecting from the possible features ranked by variance of the standardized feature values.
# 
# @param object Seurat object
# @param marker.df Dataframe containing the results of FindMarkers or FindAllMarkers
# @param possible.features Vector of features to restrict the possible values to
# @param ident Name of the meta.data column specifying the identity (cluster) of cells
# @param nfeatures How many features to return
# @return Vector of feature names of length nfeatures
#
# @export
# 
GuidedFeatureDownsample <- function(object, marker.df, possible.features = NULL, ident = NULL, nfeatures = 50) {
  if (!is.null(x = ident)) {
    Idents(object = object) <- ident
  }
  possible.features <- possible.features %||% rownames(x = object)
  marker.df <- marker.df[marker.df$gene %in% possible.features, ]
  # select marker genes first 
  ngroups <- length(x = levels(x = object))
  features.per.group <- floor(x = nfeatures / ngroups)
  features <- c()
  max.marker.features <- length(x = unique(x = marker.df$gene))
  while (length(x = features) < nfeatures & length(x = features) < max.marker.features) {
    features <- unique(x = unlist(x = lapply(
      X = unique(x = marker.df$cluster), 
      FUN = function(x) {
        cluster.markers <- marker.df[marker.df$cluster == x, ]
        cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE), ]
        return(head(x = cluster.markers, n = features.per.group)[, "gene"])
      }
    )))
    features.per.group <- features.per.group + 1
  }
  if (length(x = features) > nfeatures) {
    features.per.group <- features.per.group - 2
    features <- unique(x = unlist(x = lapply(
      X = unique(x = marker.df$cluster),
      FUN = function(x) {
        cluster.markers <- marker.df[marker.df$cluster == x, ]
        cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE), ]
        return(head(x = cluster.markers, n = features.per.group)[, "gene"])
      }
    )))
    features.per.group <- features.per.group + 1
    extra.features <- nfeatures - length(x = features)
    while (extra.features > 0 & length(x = features) <= max.marker.features) {
      features.to.add <- do.call(rbind, lapply(
        X = unique(x = marker.df$cluster),
        FUN = function(x) {
          cluster.markers <- marker.df[marker.df$cluster == x, ]
          cluster.markers <- cluster.markers[order(cluster.markers$avg_logFC, decreasing = TRUE), ]
          return(cluster.markers[features.per.group, , drop = FALSE])
        }
      ))
      features <- unique(x = c(features, head(features.to.add[order(features.to.add[, "avg_logFC"], decreasing = TRUE), "gene"], n = extra.features)))
      extra.features <- nfeatures - length(x = features)
      features.per.group <- features.per.group + 1
    }
  }
  # if exhausted marker list, pull top features ranked by standardized variance
  if (length(x = features) < nfeatures) {
    hvf <- HVFInfo(object = object)
    if (!"variance.standardized" %in% colnames(x = hvf)) {
      object <- FindVariableFeatures(object = object, selection.method = "vst", verbose = FALSE)
    }
    hvf <- HVFInfo(object = object)
    hvf <- hvf[order(hvf$variance.standardized), ]
    hvf <- hvf[(!rownames(x = hvf) %in% features), ]
    hvf <- hvf[rownames(x = hvf) %in% possible.features, ]
    extra.features <- nfeatures - length(x = features)
    features <- c(features, rownames(x = hvf)[1:extra.features]) 
  }
  return(features)
}

starmap.1 <- readRDS(file = args[1])
starmap.2 <- readRDS(file = args[2])
starmap.2 <- RenameCells(object = starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)

star <- readRDS(file = args[3])
allen <- readRDS(file = args[4])
dropseq <- readRDS(file = args[5])
markers <- read.table(file = args[6], sep = "\t", stringsAsFactors = FALSE)$gene

markers.ss2 <- intersect(markers, rownames(star[['ss2']]))
markers.ds <- intersect(markers, rownames(star[['dropseq']]))
keep.markers <- intersect(markers.ds, markers.ss2)
full.integration.ss2 <- GetAssayData(star, assay = 'ss2', slot = 'data')[keep.markers, ]
full.integration.ds <- GetAssayData(star, assay = 'dropseq', slot = 'data')[keep.markers, ]

DefaultAssay(star) <- 'integrated'

matrix.cor <- function(m1, m2) {
  # row-wise correlation between two matrices
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  cA <- m1 - Matrix::rowMeans(m1)
  cB <- m2 - Matrix::rowMeans(m2)
  sA <- sqrt(Matrix::rowMeans(cA^2))
  sB <- sqrt(Matrix::rowMeans(cB^2))
  return(Matrix::rowMeans(cA * cB) / (sA * sB))
}

run_imputation <- function(ref.obj, query.obj, features, full.integration) {
  message(paste0('using ', length(features), ' features'))
  anchors <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj,
    reference.assay = 'RNA',
    query.assay = 'integrated',
    features = features,
    dims = 1:30,
    reduction = 'cca'
  )
  refdata <- GetAssayData(
    object = ref.obj,
    assay = 'RNA',
    slot = 'data'
  )
  imputation <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    dims = 1:30,
    weight.reduction = 'pca',
    k.weight = 50,
    l2.norm = FALSE
  )
  marker.cor <- matrix.cor(full.integration, GetAssayData(imputation, slot = 'data')[rownames(full.integration), colnames(full.integration)])
  return(marker.cor)
}

# ss2
common.genes <- intersect(rownames(allen), rownames(star))

results_to_df <- function(results) {
  results.mat <- t(results)
  results.df <- as.data.frame(results.mat)
  results.df$ngene <- as.numeric(rownames(results.df))
  results.df <- tidyr::gather(results.df, gene, correlation, 1:(ncol(results.df)-1))
  
  results.df <- results.df %>% 
    group_by(gene) %>% 
    mutate(max_val = max(correlation, na.rm = TRUE)) %>% 
    ungroup()
  
  results.df$gene <- with(results.df, reorder(gene, max_val))
  return(results.df)
}

VariableFeatures(allen) <- unique(c(VariableFeatures(allen), common.genes))
allen <- ScaleData(allen)
allen <- RunPCA(allen, npcs = 30)


Idents(object = allen) <- "subclass"
#markers <- FindAllMarkers(object = allen, assay = "RNA", max.cells.per.ident = 200)
markers <- read.table(file = args[6], header = T, sep = "\t", stringsAsFactors = FALSE)

results.guided.ss2 <- list()
library(future.apply)
plan(strategy = "multicore", workers = 8)
options(future.globals.maxSize = 10^12)

#markers %>% group_by(cluster) %>% top_n(n = x, avg_logFC) %>% pull(gene)

sample.sizes <- seq(50, length(common.genes), 10)
results.guided.ss2 <- future_sapply(X = sample.sizes, FUN = function(x) {
  features.use <- GuidedFeatureDownsample(object = allen, marker.df = markers, possible.features = common.genes, nfeatures = x)
  run_imputation(ref.obj = allen,
                 query.obj = star,
                 features = features.use,
                 full.integration = full.integration.ss2
  )}
)

# create results matrix
colnames(x = results.guided.ss2) <- sample.sizes
results.guided.ss2.df <- results_to_df(results.guided.ss2)

ggplot(results.guided.ss2.df, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlab("Number of features")

ggplot(results.guided.ss2.df, aes(x = ngene, y = correlation)) +
  geom_line()

ggplot(results.guided.ss2.df, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlab("Number of features")
dir.create(paste0(getwd(), "/figures/sn2/"), showWarnings = FALSE)

ggsave(filename = args[7], #"~/xfer/starmap_downsample_guided_cluster_markers_ss2.png",
       height = 5, width = 5, units = 'in', dpi = 300, limitsize = FALSE)

sample.sizes <- seq(50, length(x = common.genes), 10)
set.seed(42)
random.order <- sample(x = common.genes)
results.random.ss2 <- future_sapply(X = sample.sizes, FUN = function(x) {
  features.use <- head(x = random.order, n = x)
  run_imputation(ref.obj = allen,
                 query.obj = star,
                 features = features.use,
                 full.integration = full.integration.ss2
  )}
)
colnames(x = results.random.ss2) <- sample.sizes
results.random.ss2.df <- results_to_df(results.random.ss2)

results.random.ss2.df$strat <- "random"
results.guided.ss2.df$strat <- "guided"


results.all <- rbind(results.random.ss2.df, results.guided.ss2.df)
results.all %>% filter(ngene %% 50 == 0 | ngene < 200 ) -> x

ggplot(results.all, aes(x = factor(ngene), y = correlation, color = strat))  + geom_boxplot(outlier.shape = NA, fatten = 12) +
  xlab("Number of Genes") + scale_color_discrete(name = "Downsampling Strategy") +
  theme(axis.text = element_text(size = 50), axis.title = element_text(size = 50), axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.5, "cm")) +
  #theme(legend.title = element_text(size = 40), legend.text = element_text(size = 30)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))) + 
  NoLegend()+
  guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_x_discrete(breaks = seq(50, 1000, 50)) +
  ggsave(filename = args[8], #"~/xfer/starmap_downsample_ss2_random_vs_guided_boxplots.png",
         height = 15, width = 35, units = "in", dpi = 300, limitsize = FALSE)

ggplot(results.random.ss2.df, aes(x = factor(ngene), y = correlation))  + geom_boxplot(outlier.shape = NA, fatten = 5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Number of Genes") +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40), axis.ticks = element_line(size = 1.5), axis.ticks.length = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20))  +  
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0))) + 
  theme(axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))) + 
  scale_x_discrete(breaks = seq(50, 1000, 50)) + 
  ggsave(filename = args[9], #"~/xfer/starmap_downsample_ss2_random_boxplots.png",
         height = 10, width = 15, units = "in", dpi = 300, limitsize = FALSE)




