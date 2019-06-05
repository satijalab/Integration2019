library(Seurat)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/citeseq.rds",
#           "~/Projects/muir/seurat_objects/citeseq_crossvalidation.rds",
#           "~/Projects/muir/figures/citeseq/cv_correlation.png",
#           "~/Projects/muir/figures/citeseq/cv_correlation_select.png",
#           "~/Projects/muir/figures/citeseq/cv_correlation_select.rds")

citeseq <- readRDS(file = args[1])
query.obj <- readRDS(file = args[2])

# all proteins
feature_scatter <- function(obj, feature1, feature2, title, base_size=10) {
  dat <- FetchData(object = obj, vars = c(feature1, feature2))
  r <- cor(dat[, 1], dat[, 2])
  df <- data_frame('Measured' = dat[, 1], 'Imputed' = dat[, 2])
  p <- ggplot(df, aes(Measured, Imputed)) +
    geom_point(size = 0.1, alpha = 0.2) +
    theme_classic(base_size = base_size) +
    theme(legend.position = 'none') +
    ggtitle(label = title,
            subtitle = paste("Pearson correlation:",
                             as.character(round(r, 3))))
  return(list(p, r))
}

p.list <- list()
r.vals <- c()
DefaultAssay(citeseq) <- "ADT"
for(i in 1:nrow(citeseq)) {
  gene <- rownames(citeseq)[[i]]
  p.data <- feature_scatter(query.obj, paste0('adt_', gene),
                            paste0('imputed_', gene),
                            gene)
  p.list[[i]] <- p.data[[1]]
  r.vals <- c(r.vals, p.data[[2]])
}
p1 <- CombinePlots(p.list, ncol = 5) 

ggsave(filename = args[3], plot = p1, height = 30, width = 30, units = 'cm')

# focus on set of 5
selected.proteins <- c('CD11a', 'CD8a', 'CD69', 'CD197-CCR7', 'CD123')

p.list <- list()
for(i in 1:length(selected.proteins)) {
  gene <- selected.proteins[[i]]
  p.data <- feature_scatter(query.obj, paste0('adt_', gene),
                            paste0('imputed_', gene),
                            gene, base_size = 18)
  p.list[[i]] <- p.data[[1]]
}

df <- data.frame('Pearson correlation' = r.vals)
p3 <- ggplot(df, aes(Pearson.correlation)) +
  geom_histogram(fill = "#31688EFF", bins = 20, size = 1/2) +
  xlim(c(0,1)) +
  ggtitle(label = "Accuracy",
          subtitle = paste0("Median: ", as.character(round(median(df$Pearson.correlation), 3)))) +
  theme_classic(base_size = 18)
p.list[[6]] <- p3

p2 <- CombinePlots(p.list, ncol = 3)
ggsave(filename = args[4], plot = p2, height = 8, width = 12, units = 'in', dpi = 500)
saveRDS(object = p2, file = args[5])
