library(Seurat, quietly = TRUE)
suppressMessages(library(cowplot, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/Projects/muir/seurat_objects/pbmc3k.rds",
#           "~/Projects/muir/seurat_objects/pbmc33k.rds",
#           "~/Projects/muir/analysis_data/pbmc3k_predictions.rds")

pbmc <- readRDS(file = args[1])
pbmc33k <- readRDS(file = args[2])
predictions <- readRDS(file = args[3])

pbmc <- AddMetaData(object = pbmc, metadata = predictions)
Idents(object = pbmc) <- "predicted.id"

markers33k <- FindAllMarkers(object = pbmc33k, test.use = "LR", only.pos = TRUE, max.cells.per.ident = 500)
markers33k %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) %>% pull(gene) -> top.markers

pbmc <- ScaleData(object = pbmc, features = top.markers)
pa <- DoHeatmap(object = pbmc, features = top.markers, size = 4, group.by = "cluster")
pb <- DoHeatmap(object = pbmc, features = top.markers, size = 4, group.by = "predicted.id")

saveRDS(object = list(pa, pb), file = args[4])