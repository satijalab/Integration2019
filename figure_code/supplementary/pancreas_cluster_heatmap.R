# average heatmap per dataset

library(Seurat)
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
args <- c("~/Projects/muir/seurat_objects/integrated_pancreas.rds", "~/Projects/muir/figures/supplementary/pancreas_heatmap.png")

pi <- readRDS(file = args[1])

DefaultAssay(object = pi) <- "RNA"
clusterIDs <- names(x = sort(x = table(Idents(object = pi)), decreasing = TRUE))
pi.sub <- pi[, WhichCells(object = pi, downsample = 200)]
Idents(object = pi.sub) <- factor(x = Idents(object = pi.sub), levels = clusterIDs)
Idents(object = pi) <- factor(x = Idents(object = pi), levels = clusterIDs)

all.markers <- FindAllMarkers(
  object = pi.sub, 
  logfc.threshold = 1,
  assay = "RNA",
  test.use = "LR",
  only.pos = TRUE,
  latent.vars = "tech",
  min.diff.pct = 0.2
)

all.markers %>% group_by(cluster) %>% top_n(-10, p_val) %>% pull(gene) -> hm.markers

object.list <- SplitObject(object = pi, split.by = "replicate")

sub.list <- list()
hm.list <- list()
for (i in 1:8) {
  sub.list[[i]] <- object.list[[i]][, WhichCells(object = object.list[[i]], downsample = 25)]
  sub.list[[i]] <- ScaleData(object = sub.list[[i]], features = hm.markers)
  hm.list[[i]] <- DoHeatmap(object = sub.list[[i]], features = hm.markers, disp.max = 3, size = 4, angle=90) 
}

names(x = hm.list) <- names(object.list)
names(x = hm.list[5:8]) <- paste0("indrop_", names(x = object.list)[5:8])
hm.ordered <- hm.list[c(5,6,7,8,2,4,1,3)]

for(i in 1:8) {
  hm.ordered[[i]] <- hm.ordered[[i]] + NoLegend() + ggtitle(names(x = hm.ordered)[[i]]) + theme(axis.text = element_text(size=3))
}

p.final <- plot_grid(plotlist = hm.ordered, ncol = 4)

save_plot(plot = p.final, filename = args[2], base_width = 30, base_height = 20)
