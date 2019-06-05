args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/20180505_BY3_1kgenes.rds", 
#           "~/Projects/muir/seurat_objects/20180410-BY3_1kgenes.rds",
#           "/home/stuartt/github/muir/seurat_objects/integrated_starmap.rds",
#           "/home/stuartt/github/muir/seurat_objects/allen_brain.rds",
#           "/home/stuartt/github/muir/seurat_objects/dropseq_cortex.rds",

devtools::load_all(args[7])
library(dplyr)
library(ggplot2)
library(cowplot)

starmap.1 <- readRDS(file = args[1])
starmap.2 <- readRDS(file = args[2])
starmap.2 <- RenameCells(starmap.2, add.cell.id = 'r2')
starmap.2@misc$spatial$cell <- paste0('r2_', starmap.2@misc$spatial$cell)
star <- readRDS(file = args[3])
allen <- readRDS(file = args[4])
dropseq <- readRDS(file = args[5])
common.genes.ss2 <- intersect(rownames(allen), rownames(star))
common.genes.ds <- intersect(rownames(dropseq), rownames(star))

DefaultAssay(star) <- 'integrated'

matrix.cor <- function(m1, m2) {
  # rowwise correaltion between matrices
  m1 <- as.matrix(m1)
  m2 <- as.matrix(m2)
  cA <- m1 - Matrix::rowMeans(m1)
  cB <- m2 - Matrix::rowMeans(m2)
  sA <- sqrt(Matrix::rowMeans(cA^2))
  sB <- sqrt(Matrix::rowMeans(cB^2))
  return(Matrix::rowMeans(cA * cB) / (sA * sB))
}

genes.use.ss2 <- unique(c(VariableFeatures(allen), common.genes.ss2))
full.integration.ss2 <- GetAssayData(star, assay = 'ss2', slot = 'data')[genes.use.ss2, ]
allen.expression <- GetAssayData(allen, assay = 'RNA', slot = 'data')[genes.use.ss2, ]
allen.expression <- t(as.matrix(allen.expression))
gene.gene.correlation <- cor(allen.expression)

redundancy <- function(gene.gene.correlation, selected.genes, cutoff=0.1) {
  assess.genes <- setdiff(rownames(gene.gene.correlation), selected.genes)
  #return(Matrix::rowSums(abs(gene.gene.correlation[assess.genes, selected.genes]) > cutoff))
  return(apply(
    X = abs(gene.gene.correlation[assess.genes, selected.genes]),
    FUN = function(x) {
      sort(x, decreasing = T)[3]
    },
    MARGIN = 1))
  #return(apply(X = abs(gene.gene.correlation[assess.genes, selected.genes]), FUN = max, MARGIN = 1))

  #return(Matrix::rowSums(abs(gene.gene.correlation[assess.genes, selected.genes]))/length(selected.genes))
}

find_mean_expression <- function(gene.gene.correlation, selected.genes, seq) {
  assess.genes <- setdiff(rownames(gene.gene.correlation), selected.genes)
  assess.cor <- abs(gene.gene.correlation[assess.genes, selected.genes])
  maxcor <- apply(assess.cor, 1, which.max)
  maxcor <- colnames(assess.cor)[maxcor]
  expression.data <- GetAssayData(seq, asssay = 'RNA', slot = 'data')[maxcor, ]
  binary <- expression.data
  binary[binary > 0] <- 1
  nonzero <- Matrix::rowSums(binary)
  mean.expressed <- Matrix::rowSums(expression.data) / nonzero
  return(mean.expressed)
}

run_imputation <- function(ref.obj, query.obj, features, full.integration, gene.gene.correlation) {
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
  return(imputation)
}

results <- list()
n <- 100
set.seed(42)
random.genes <- sapply(1:n, function(x) sample(x = common.genes.ss2, size = 300, replace = FALSE))
for(i in 1:1) {
  message(i)
  selected.genes <- random.genes[, i]
  imputed <- run_imputation(ref.obj = allen,
                            query.obj = star,
                            full.integration = full.integration.ss2,
                            features = selected.genes,
                            gene.gene.correlation = gene.gene.correlation)
  results[[i]] <- imputed
  #df <- data.frame('Accuracy' = imputed$accuracy, "Redundancy" = imputed$redundancy,
  #                 'Replicate' = i, 'Gene' = names(imputed$accuracy),
  #                 'Expression' = imputed$expression)
  #results <- rbind(results, df)
}

results2 <- data.frame()
for(i in 1:1){
  imputed <- results[[i]]
  marker.cor <- matrix.cor(full.integration.ss2, GetAssayData(imputed, slot = 'data')[rownames(full.integration.ss2), colnames(full.integration.ss2)])
  gene.redundancy <- redundancy(gene.gene.correlation = gene.gene.correlation, selected.genes = random.genes[, i])
  mean.expression <- find_mean_expression(gene.gene.correlation = gene.gene.correlation, selected.genes = random.genes[, i], seq = allen)
  
  df <- data.frame('Accuracy' = marker.cor[names(gene.redundancy)], 
                   "Redundancy" = gene.redundancy,
                  'Replicate' = i, 
                  'Gene' = names(marker.cor[names(gene.redundancy)]),
                  'Expression' = mean.expression)
  results2 <- rbind(results2, df)
}


# results.shash <- results

rsq <- cor(results2$Accuracy, results2$Redundancy, use='complete')

p <- ggplot(results2, aes(Redundancy, Accuracy)) +
  geom_point(size=0.1) +
  scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("3rd highest,  rsq = ", rsq)) + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25)) + 
  ggsave(filename = args[6], #"~/xfer/redundancy_vs_accuracy.png",
    height = 9, width = 9, units = "in", dpi = 300, limitsize = FALSE)


# p1 <- ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("Max,  r = ", rsq))
# 
# p2 <- ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("5th highest,  rsq = ", rsq))
# 
# p3 <- ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("10th highest,  rsq = ", rsq))


# p5 <-  ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("Average,  rsq = ", rsq))
# 
# p6 <-  ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("15th highest,  rsq = ", rsq))
# 
# p7 <-  ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("2nd highest,  rsq = ", rsq))
# 
# p8 <-  ggplot(results2, aes(Redundancy, Accuracy)) +
#   geom_point(size=0.1) +
#   scale_colour_distiller(palette = 'RdYlBu') + ggtitle(paste0("Mean Top 3,  rsq = ", rsq))
# 

# p1 <- results2 %>% 
#   filter(Gene == 'Creb3l4') %>% 
#   ggplot(., aes(Redundancy, Accuracy)) +
#   geom_point() + ggtitle("Creb3l4")
# 
# p2 <- results2 %>% 
#   filter(Gene == 'Fmr1nb') %>% 
#   ggplot(., aes(Redundancy, Accuracy)) +
#   geom_point() + ggtitle("Fmr1nb")
# 
# p3 <- results2 %>% 
#   filter(Gene == 'Lrr1') %>% 
#   ggplot(., aes(Redundancy, Accuracy)) +
#   geom_point() + ggtitle("Lrr1")
# 
# p4 <- results2 %>% 
#   filter(Gene == 'Mettl3') %>% 
#   ggplot(., aes(Redundancy, Accuracy)) +
#   geom_point() + ggtitle("Mettl3")
# 
# plot_grid(p1,p2,p3, p4)
