library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("seurat_objects/citeseq.rds", 
#           "seurat_objects/citeseq_downsampling/1000_feature_downsampling.rds", 
#           "figures/citeseq/nfeature_correlation.pdf",
#           "figures/citeseq/nfeature_select_correlation.pdf",
#           "figures/citeseq/nfeature_heatmap_correlation.pdf",
#           "figures/citeseq/nfeature_heatmap_correlation.rds")

files <- list.files("seurat_objects/citeseq_downsampling")
files <- paste0("seurat_objects/citeseq_downsampling/", files)
citeseq <- readRDS(file = args[1])
imputed.adt <- list()

for(i in files) {
  ds.level <- unlist(strsplit(i, "seurat_objects/citeseq_downsampling/", perl = FALSE))[[2]]
  ds.level <- unlist(strsplit(ds.level, "_", perl = FALSE))[[1]]
  obj <- readRDS(i)
  imputed.adt[[ds.level]] <- obj
}

ground.truth <- SubsetData(citeseq, cells = colnames(imputed.adt[[1]]))
cell.order <- colnames(ground.truth)

imputation.accuracy <- list()
x <- 1
for(i in names(imputed.adt)) {
  for(gene in 1:nrow(imputed.adt[[i]])) {
    correlation <- cor(ground.truth[['ADT']]@data[gene, cell.order], imputed.adt[[i]][gene, cell.order])
    gene.name <- rownames(imputed.adt[[i]])[gene]
    imputation.accuracy[[x]] <- c('ADT' = gene.name, 'pearson.cor' = correlation, 'nfeatures' = as.numeric(i))
    x <- x + 1
  }
}

results <- as.data.frame(do.call(what = rbind, args = imputation.accuracy), stringsAsFactors = FALSE)
results$pearson.cor <- as.numeric(results$pearson.cor)
results$nfeatures <- as.numeric(results$nfeatures)

results <- results %>% 
  group_by(ADT) %>% 
  mutate(max_val = max(pearson.cor)) %>% 
  ungroup()

results$ADT <- with(results, reorder(ADT, max_val))

select.genes <- c('CD69', 'CD56', 'CD16', 'CD4', 'CD8a', 'CD34')

ggplot(results, aes(nfeatures, pearson.cor)) +
  geom_line() +
  ylim(c(0, 1)) +
  theme_minimal(base_size = 8) +
  facet_wrap(~ADT) +
  xlab("Number of features") +
  ylab("Pearson correlation") + 
  ggsave(filename = args[3], height = 10, width = 18, units = 'cm')

ggplot(results[results$ADT %in% select.genes, ], aes(nfeatures, pearson.cor, color = ADT)) +
  geom_line() +
  ylim(c(0, 1)) +
  theme_minimal(base_size = 8) +
  xlab("Number of features") +
  ylab("Pearson correlation") +
  ggsave(filename = args[4], height = 5, width = 7, units = 'cm')

p <- ggplot(results, aes(nfeatures, ADT, fill = pearson.cor)) +
  geom_tile() +
  scale_fill_gradientn(colours = viridis::viridis(10), limits = c(0, 1)) +
  theme_classic(base_size = 8) 

ggsave(filename = args[5], plot = p, height = 6, width = 11, units = 'cm')
saveRDS(object = p, file = args[6])

