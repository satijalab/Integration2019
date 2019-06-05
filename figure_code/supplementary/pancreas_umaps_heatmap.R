suppressMessages(library(Seurat))
suppressMessages(library(cowplot))
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

args <- c("~/Projects/muir/seurat_objects/all_pancreas_no_integration.rds",
          "~/Projects/muir/seurat_objects/integrated_pancreas.rds",
          "~/Projects/muir/figures/supplementary/pancreas_full_umaps.pdf",
          "~/Projects/muir/figures/supplementary/pancreas_heatmap.png",
          "~/Projects/muir/figures/supplementary/pancreas_supp.pdf")

pi.noint <- readRDS(file = args[1])
pi <- readRDS(file = args[2])

p1 <- DimPlot(object = pi.noint, reduction = "umap", group.by = "replicate") + 
  theme(legend.position = "top", legend.justification = "center", legend.text = element_text(size = 18)) + 
  ggtitle("Unintegrated Datasets") + guides(colour = guide_legend(override.aes = list(size=8)))
p2 <- DimPlot(object = pi, reduction = "umap", group.by = "replicate") + 
  theme(legend.position = "top", legend.justification = "center", legend.text = element_text(size = 18)) +
  ggtitle("Integrated Datasets") + guides(colour = guide_legend(override.aes = list(size=8)))
p3 <- DimPlot(object = pi, reduction = "umap", group.by = "celltype") + 
  theme(legend.position = "top", legend.justification = "center", legend.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size=8)))
pfinal <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "lt", labels = LETTERS[1:3])

save_plot(plot = pfinal, filename = args[3], base_aspect_ratio = 2/1, base_height = 10)


# average heatmap per dataset


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

pb <- plot_grid(plotlist = hm.ordered, ncol = 4)

save_plot(plot = pb, filename = args[4], base_width = 30, base_height = 20)

ps <- plot_grid(pfinal, pb, nrow = 2, rel_heights = c(0.75, 1))

save_plot(plot = ps, filename = args[5], base_width = 30, base_height = 25)
