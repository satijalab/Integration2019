library(ggplot2)
library(tidyr)
library(dplyr)

corces <- read.table("raw_data/corces_atac/GSE74912_ATACseq_All_Counts.txt",
                     stringsAsFactors = FALSE, header = TRUE, sep = "\t")
cols.use <- colnames(corces)[grepl("CD4|CD8|Nk|Bcell", colnames(corces))]
corces <- corces[,  cols.use]
colnames(corces)[grepl('Bcell', colnames(corces))] <- paste0('B_', seq(4))
colnames(corces)[grepl('Nk', colnames(corces))] <- paste0('Nk_', seq(4))
colnames(corces)[grepl('CD8', colnames(corces))] <- paste0('CD8_', seq(5))
colnames(corces)[grepl('CD4', colnames(corces))] <- paste0('CD4_', seq(5))
sc.atac <- read.table("analysis_data/atac/10x/celltypes/single_cell/pbmc_10k_atac/celltype_coverage.bed",
                      stringsAsFactors = FALSE, header = FALSE, sep = "\t",
                      col.names = c('chr', 'start', 'stop', 'B.pro', 'B.pre', 'CD4.memory',
                                    'CD4.naive', 'CD8.effector', 'Double.negative', 'CD8.naive', 'NK',
                                    'CD14.monocytes', 'CD16.monocytes', 'Dendritic', 'pDC'))
sc.atac <- sc.atac[, 4:ncol(sc.atac)]
corces <- as.matrix(corces)
sc.atac <- as.matrix(sc.atac)

# need to filter peaks in corces to retain only those that are different between cell types
# otherwise most peaks are the same and boosts correlation between all cell types
find_da <- function(counts, group.1, group.2=NULL, limit=5, cutoff=50) {
  da <- c()
  if (is.null(group.2)) {
    group.2 <- setdiff(colnames(counts), group.1)
  }
  for (i in 1:nrow(counts)) {
    if (sum(counts[i, ]) > cutoff) {
      mn1 <- sum(counts[i, group.1]) / length(group.1)
      mn2 <- sum(counts[i, group.2]) / length(group.2)
      if (((mn1+1) / (mn2+1)) > limit) {
        da <- c(da, i)
      }
    }
  }
  return(da)
}

# Comparing to mean accessibility in all other cell types
cd8 <- find_da(corces, group.1 = paste0('CD8_', seq(5)))
cd4 <- find_da(corces, group.1 = paste0('CD4_', seq(5)))
nk <- find_da(corces, group.1 = paste0('Nk_', seq(4)))
b <- find_da(corces, group.1 = paste0('B_', seq(4)))

# remove peaks that are in multiple cell types
all.markers <- c(cd8, cd4, nk, b)
marker.counts <- table(all.markers)
markers <- as.numeric(names(which(marker.counts == 1)))

# subset matrix
corces.markers <- corces[markers, ]
sc.markers <- sc.atac[markers, ]

rvals <- cor(sc.markers, corces.markers)
rvals <- as.data.frame(rvals)
rvals$single.cell <- rownames(rvals)
r.tidy <- gather(rvals, bulk, correlation, 1:(ncol(rvals)-1))
r.tidy$single.cell <- factor(r.tidy$single.cell, levels = c('B.pre', 'B.pro', 'CD4.naive', 'CD4.memory',
                                                            'CD8.naive', 'Double.negative', 'CD8.effector',
                                                            'NK',
                                                            'Dendritic', 'pDC',
                                                            'CD14.monocytes', 'CD16.monocytes'))

ggplot(r.tidy, aes(single.cell, bulk, fill = correlation)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggsave("figures/atac/pbmc/bulk_correlation_all.pdf", height = 5, width = 6)


ct.use <- c('B.pre', 'B.pro', 'CD4.naive',
            'CD8.naive', 'CD8.effector', 'NK')

r.tidy %>% 
  filter(single.cell %in% ct.use) %>% 
  ggplot(., aes(single.cell, bulk, fill = correlation)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdYlBu") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggsave("figures/atac/pbmc/bulk_correlation.pdf", height = 5, width = 6)