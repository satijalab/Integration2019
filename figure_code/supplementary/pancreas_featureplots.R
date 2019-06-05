suppressMessages(library(Seurat))
suppressMessages(library(cowplot))

args <- commandArgs(trailingOnly = TRUE)

#args <- c("~/Projects/muir/seurat_objects/integrated_pancreas.rds", "~/Projects/muir/figures/supplementary/pancreas_featureplots.pdf")
pi <- readRDS(file = args[1])

# activated vs quiscent stellate cells
# TFF1 ductal cells
# beta ER stress

DefaultAssay(object = pi) <- "RNA"
p1 <- DimPlot(object = pi, reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
p2 <- FeaturePlot(object = pi, features = "PDGFRA", reduction = "umap", min.cutoff = "q5", max.cutoff = "q95")
p3 <- FeaturePlot(object = pi, features = "RGS5", reduction = "umap", min.cutoff = "q5", max.cutoff = "q95")
p4 <- FeaturePlot(object = pi, features = "KRT19", reduction = "umap")
p5 <- FeaturePlot(object = pi, features = "TFF1", reduction = "umap")
p6 <- FeaturePlot(object = pi, features = "INS", reduction = "umap")
p7 <- FeaturePlot(object = pi, features = "HERPUD1", reduction = "umap", min.cutoff = "q3", max.cutoff = "q97")

r1 <- plot_grid(p1, p2, p3, labels = c("A", "B", "C"), nrow = 1, rel_widths = c(1, 0.5, 0.5))
r2 <- plot_grid(p4, p5, p6, p7, labels = c("D", "E", "F", "G"), nrow = 1)

pfinal <- plot_grid(
  r1, r2,
  nrow = 2
)

save_plot(plot = pfinal, filename = args[2], base_aspect_ratio = 1.5, base_height = 15)
