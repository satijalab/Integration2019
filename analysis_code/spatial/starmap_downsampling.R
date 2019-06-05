library(Seurat)
library(ggplot2)
library(dplyr)


starmap.1 <- readRDS("seurat_objects/20180505_BY3_1kgenes.rds")
starmap.2 <- readRDS("seurat_objects/20180410-BY3_1kgenes.rds")
starmap.2 <- RenameCells(starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)
allen <- readRDS("seurat_objects/allen_brain.rds")
star <- readRDS("seurat_objects/integrated_starmap.rds")
dropseq <- readRDS("seurat_objects/dropseq_frontal_cortex.rds")
markers <- read.table("analysis_data/starmap_markers.tsv", sep = "\t", stringsAsFactors = FALSE)$gene
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
# randomize
random.order <- sample(common.genes)
results.ss2 <- list()
for(i in seq(50, length(common.genes), 10)) {
  features.use <- head(random.order, i)
  results.ss2[[as.character(length(common.genes)-i)]] <- run_imputation(ref.obj = allen, query.obj = star,
                                                                        features = features.use,
                                                                        full.integration = full.integration.ss2)
}

# dropseq
common.genes <- intersect(rownames(dropseq), rownames(star))
# randomize
random.order <- sample(common.genes)
results.ds <- list()
for(i in seq(50, length(common.genes), 10)) {
  features.use <- head(random.order, i)
  results.ds[[as.character(i)]] <- run_imputation(ref.obj = dropseq,
                                                  query.obj = star,
                                                  features = features.use,
                                                  full.integration = full.integration.ds)
}

results_to_df <- function(results) {
  results.mat <- do.call(what = rbind, args = results)
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

# create results matrix
results.df.ss2 <- results_to_df(results.ss2)
results.df.ds <- results_to_df(results.ds)
results.df.ds$gene <- with(results.df.ds, reorder(gene, results.df.ss2$max_val))

ggplot(results.df.ss2, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggsave("figures/spatial/starmap_downsampling_random_ss2.png",
         height = 5, width = 5, units = 'in', dpi = 300, limitsize = FALSE)

ggplot(results.df.ds, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggsave("figures/spatial/starmap_downsampling_random_dropseq.png",
         height = 5, width = 5, units = 'in', dpi = 300, limitsize = FALSE)
  
# guided downsampling
# should guide based on information in the RNA-seq
# need to re-run PCA including all the starmap genes
VariableFeatures(allen) <- unique(c(VariableFeatures(allen), common.genes))
allen <- ScaleData(allen)
allen <- RunPCA(allen, npcs = 30)

VariableFeatures(dropseq) <- unique(c(VariableFeatures(dropseq), common.genes))
dropseq <- ScaleData(dropseq)
dropseq <- RunPCA(dropseq, npcs = 30)

get_ranks <- function(obj, common.genes) {
  ranking <- abs(x = Loadings(obj[['pca']]) * obj[['pca']]@stdev)
  pca.genes <- intersect(rownames(ranking), common.genes)
  ranking <- ranking[pca.genes, ]
  ranking <- sort(x = rowSums(x = ranking), decreasing = TRUE)
  excluded.genes <- setdiff(common.genes, pca.genes)
  excl <- rep(0, length(excluded.genes))
  names(excl) <- excluded.genes
  ranking <- c(ranking, excl)
  return(ranking)
}

ranks.ss2 <- get_ranks(allen, common.genes)
ranks.ds <- get_ranks(dropseq, common.genes)

results.guided.ss2 <- list()
for(i in seq(50, length(common.genes), 10)) {
  features.use <- names(head(ranks.ss2, i))
  results.guided.ss2[[as.character(i)]] <- run_imputation(ref.obj = allen,
                                                      query.obj = star,
                                                      features = features.use,
                                                      full.integration = full.integration.ss2)
}

results.guided.ds <- list()
for(i in seq(50, length(common.genes), 10)) {
  features.use <- names(head(ranks.ds, i))
  results.guided.ds[[as.character(i)]] <- run_imputation(ref.obj = dropseq,
                                                         query.obj = star,
                                                         features = features.use,
                                                         full.integration = full.integration.ds)
}

# create results matrix
results.guided.ss2.df <- results_to_df(results.guided.ss2)
results.guided.ds.df <- results_to_df(results.guided.ds)
results.guided.ds.df$gene <- with(results.guided.ds.df, reorder(gene, results.guided.ss2.df$max_val))

ggplot(results.guided.ss2.df, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
ggsave("figures/spatial/starmap_downsample_guided_ss2.png",
       height = 5, width = 5, units = 'in', dpi = 300, limitsize = FALSE)

ggplot(results.guided.ds.df, aes(ngene, gene, fill = correlation)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_classic(base_size = 8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  ggsave("figures/spatial/starmap_downsample_guided_ds.png",
         height = 5, width = 5, units = 'in', dpi = 300, limitsize = FALSE)
